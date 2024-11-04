"""Module run of viperleed.calc.

Defines the main functionality for running a
viperleed calculation from a set of input files.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-11-12'  # Was originally tleedm.py
__license__ = 'GPLv3+'

import logging
import os
from pathlib import Path
import shutil

from viperleed import __version__
from viperleed.calc import LOGGER as logger
from viperleed.calc.classes import rparams
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.files import parameters, poscar
from viperleed.calc.files.tenserleed import get_tensorleed_path
from viperleed.calc.lib.log_utils import close_all_handlers
from viperleed.calc.lib.log_utils import prepare_calc_logger
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.sections.cleanup import cleanup
from viperleed.calc.sections.cleanup import prerun_clean
from viperleed.calc.sections.initialization import (
    warn_if_slab_has_atoms_in_multiple_c_cells
    )
from viperleed.calc.sections.run_sections import section_loop



def run_calc(system_name=None,
             console_output=True,
             slab=None,
             preset_params=None,
             source=None):
    """Run a ViPErLEED calculation.

    By default, a PARAMETERS and a POSCAR file are expected, but can be
    replaced by passing the `slab` and/or `present_params` kwargs.

    Parameters
    ----------
    system_name : str, optional
        Used as a comment in some output-file headers. If not
        specified, the name of the parent directory is used
        instead. Default is None.
    console_output : bool, optional
        If False, will not add a logging.StreamHandler. Output
        will only be printed to the log file. Default is True.
    slab : Slab, optional
        Start from a pre-existing slab, instead of reading from
        POSCAR. Default is None.
    preset_params : dict, optional
        Parameters to add to the Rparams object after a PARAMETERS file
        has been read. Keys should be attributes of Rparams. Values in
        `preset_params` will overwrite values read from the PARAMETERS
        file, if present in both. If no PARAMETERS file is read,
        parameters will be read exclusively from `preset_params`.
        Default is None.
    source : pathlike, optional
        Path where the tensor-LEED directory can be found, which
        contains all the TensErLEED source code. If not given
        or None, try taking it from the environment variable
        VIPERLEED_TENSORLEED. As a last resort, use the current
        directory. Default is None.

    Returns
    -------
    exit_code : int
        0: exit without errors.
        1: clean exit through KeyboardInterrupt
        2: exit due to Exception before entering main loop
        3: exit due to Exception during main loop
    state_recorder : CalcStateRecorder or None
        A collection of the intermediate states of the Slab and Rparams
        objects at the end of each executed section. None if the run
        terminates before any section is executed.
    """
    os.umask(0)
    # start logger, write to file:
    timestamp = DateTimeFormat.FILE_SUFFIX.now()

    log_name = f'{LOG_PREFIX}-{timestamp}.log'
    prepare_calc_logger(logger,
                        file_name=log_name,
                        with_console=console_output)
    logger.info(f"Starting new log: {log_name}\nTime of execution: "
                + DateTimeFormat.LOG_CONTENTS.now())
    logger.info(f"This is ViPErLEED version {__version__}\n")

    tmp_manifest = [DEFAULT_SUPP, DEFAULT_OUT, log_name]
    try:
        rp = parameters.read()
    except FileNotFoundError:
        if not preset_params:
            logger.error("No PARAMETERS file found, and no preset parameters "
                         "passed. Execution will stop.")
            cleanup(tmp_manifest)
            return 2, None
        rp = rparams.Rparams()
    except Exception:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmp_manifest)
        return 2, None

    # check if this is going to be a domain search
    domains = False
    if "DOMAIN" in rp.readParams:
        domains = True

    if domains:  # no POSCAR in main folder for domain searches
        slab = None
    elif slab is None:
        poscar_file = Path("POSCAR")
        if not poscar_file.is_file():
            logger.error("POSCAR not found. Stopping execution...")
            cleanup(tmp_manifest)
            return 2, None

        logger.info("Reading structure from file POSCAR")
        try:
            slab = poscar.read(filename=poscar_file)
        except Exception:
            logger.error("Exception while reading POSCAR", exc_info=True)
            cleanup(tmp_manifest)
            return 2, None

        if not slab.preprocessed:
            logger.info("The POSCAR file will be processed and overwritten. "
                        "Copying the original POSCAR to POSCAR_user...")
            try:
                shutil.copy2(poscar_file, "POSCAR_user")
            except OSError:
                logger.error("Failed to copy POSCAR to POSCAR_user. Stopping "
                             "execution...")
                cleanup(tmp_manifest)
                return 2, None
            tmp_manifest.append("POSCAR_user")
    try:
        # interpret the PARAMETERS file
        parameters.interpret(rp, slab=slab, silent=False)
    except (parameters.errors.ParameterNeedsSlabError,
            parameters.errors.SuperfluousParameterError):
        # Domains calculation is the only case in which slab is None
        logger.error('Main PARAMETERS file contains an invalid parameter '
                     'for a multi-domain calculation', exc_info=True)
        cleanup(tmp_manifest)
        return 2, None
    except parameters.errors.ParameterError:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmp_manifest)
        return 2, None

    rp.timestamp = timestamp
    rp.manifest = tmp_manifest

    # Load parameter presets, overriding those in PARAMETERS
    preset_params = preset_params or {}
    try:
        rp.update(preset_params)
    except (ValueError, TypeError):
        logger.warning(f"Error applying preset parameters: ", exc_info=True)

    # set logging level
    if 'LOG_LEVEL' in preset_params:
        logger.info(f'Overriding log level to {rp.LOG_LEVEL}.')
    logger.setLevel(rp.LOG_LEVEL)
    logger.debug("PARAMETERS file was read successfully")

    if not domains:
        warn_if_slab_has_atoms_in_multiple_c_cells(slab, rp)
        slab.full_update(rp)   # gets PARAMETERS data into slab
        rp.fileLoaded["POSCAR"] = True

    # set source directory
    try:
        rp.source_dir = get_tensorleed_path(source)
    except (ValueError, FileNotFoundError) as exc:
        logger.warning(f'{exc} This may cause errors.')
        rp.source_dir = Path(source or '').resolve()

    if system_name is None:
        system_name = _get_parent_directory_name()
    rp.systemName = system_name

    # check if halting condition is already in effect:
    if rp.halt >= rp.HALTING:
        logger.info("Halting execution...")
        cleanup(rp.manifest, rp)
        return 0, None

    rp.updateDerivedParams()
    logger.info(f"ViPErLEED is using TensErLEED version {str(rp.TL_VERSION)}.")

    prerun_clean(rp, log_name)
    exit_code, state_recorder = section_loop(rp, slab)

    # Finalize logging - if not done, will break unit testing
    close_all_handlers(logger)
    logging.shutdown()

    return exit_code, state_recorder



def _get_parent_directory_name():
    """Return the name of the directory above the current one."""
    _system_name = Path.cwd().resolve().parent.name
    logger.info('No system name specified. Using name of parent '
                f'directory: {_system_name}')
    return _system_name
