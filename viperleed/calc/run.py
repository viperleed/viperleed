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
import time

from viperleed import __version__
from viperleed.calc import LOGGER as logger
from viperleed.calc import LOG_PREFIX
from viperleed.calc.classes import rparams
from viperleed.calc.files import parameters, poscar
from viperleed.calc.lib.base import CustomLogFormatter
from viperleed.calc.sections.cleanup import prerun_clean, cleanup
from viperleed.calc.sections.initialization import (
    warn_if_slab_has_atoms_in_multiple_c_cells
    )
from viperleed.calc.sections.run_sections import section_loop


_DEFAULT_SYSTEM_NAME = 'unknown'


def run_calc(system_name=None,
             console_output=True,
             slab=None,
             preset_params=None,
             source=''):
    """Run a ViPErLEED calculation.

    By default, a PARAMETERS and a POSCAR file are expected, but can be
    replaced by passing the `slab` and/or `present_params` kwargs.

    Parameters
    ----------
    system_name : str, optional
        Used as a comment in some output-file headers. If not
        specified, 'unknown' is used instead. Default is None.
    console_output : bool, optional
        If False, will not add a logging.StreamHandler. Output will
        only be printed to the log file. Default is True.
    slab : Slab, optional
        Start from a pre-existing slab, instead of reading from POSCAR.
        Default is None.
    preset_params : dict, optional
        Parameters to add to the Rparams object after a PARAMETERS file
        has been read. Keys should be attributes of Rparams. Values in
        `preset_params` will overwrite values read from the PARAMETERS
        file, if present in both. If no PARAMETERS file is read,
        parameters will be read exclusively from `preset_params`.
        Default is None.
    source : pathlike, optional
        Path where the 'tensorleed' directory can be found, which
        contains all the TensErLEED source code. Default is an empty
        string, i.e., the current directory.

    Returns
    -------
    int
        0: exit without errors.
        1: clean exit through KeyboardInterrupt
        2: exit due to Exception before entering main loop
        3: exit due to Exception during main loop
    """
    os.umask(0)
    # start logger, write to file:
    timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())

    log_name = f'{LOG_PREFIX}-{timestamp}.log'
    logger.setLevel(logging.INFO)
    logFormatter = CustomLogFormatter()
    fileHandler = logging.FileHandler(log_name, mode="w")
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    if console_output:
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)

    logger.info(f"Starting new log: {log_name}\nTime of execution (UTC): "
                + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    logger.info(f"This is ViPErLEED version {__version__}\n")
    logger.info("! THIS VERSION IS A PRE-RELEASE NOT MEANT FOR PUBLIC "         # TODO: remove for v1.0
                "DISTRIBUTION !\n")

    tmp_manifest = ["SUPP", "OUT", log_name]
    try:
        rp = parameters.read()
    except FileNotFoundError:
        if not preset_params:
            logger.error("No PARAMETERS file found, and no preset parameters "
                         "passed. Execution will stop.")
            cleanup(tmp_manifest)
            return 2
        rp = rparams.Rparams()
    except Exception:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmp_manifest)
        return 2

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
            return 2

        logger.info("Reading structure from file POSCAR")
        try:
            slab = poscar.read(filename=poscar_file)
        except Exception:
            logger.error("Exception while reading POSCAR", exc_info=True)
            cleanup(tmp_manifest)
            return 2

        if not slab.preprocessed:
            logger.info("The POSCAR file will be processed and overwritten. "
                        "Copying the original POSCAR to POSCAR_user...")
            try:
                shutil.copy2(poscar_file, "POSCAR_user")
            except OSError:
                logger.error("Failed to copy POSCAR to POSCAR_user. Stopping "
                             "execution...")
                cleanup(tmp_manifest)
                return 2
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
        return 2
    except parameters.errors.ParameterError:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmp_manifest)
        return 2

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
    _source = Path(source).resolve()
    if not _source.is_dir():
        logger.warning(f"tensorleed directory {source} not found.")
    if _source.name == "tensorleed":
        rp.source_dir = _source
    elif _source.parent.name == "tensorleed":
        logger.warning(f"tensorleed directory found in {_source.parent}, "
                       f"using that instead of {_source}.")
        rp.source_dir = _source.parent
    elif (_source / "tensorleed").is_dir():
        logger.warning(f"tensorleed directory found in {_source}, using that "
                       f"instead of {_source}.")
        rp.source_dir = _source / "tensorleed"
    else:
        logger.warning(f"Could not find a tensorleed directory at {_source}. "
                       "This may cause errors.")
        rp.source_dir = _source

    rp.systemName = system_name or _DEFAULT_SYSTEM_NAME
    if not system_name:
        logger.info(f'No system name specified. Using name {rp.systemName!r}.')

    # check if halting condition is already in effect:
    if rp.halt >= rp.HALTING:
        logger.info("Halting execution...")
        cleanup(rp.manifest, rp)
        return 0

    rp.updateDerivedParams()
    logger.info(f"ViPErLEED is using TensErLEED version {rp.TL_VERSION_STR}.")

    prerun_clean(rp, log_name)
    exit_code = section_loop(rp, slab)

    # Finalize logging - if not done, will break unit testing
    logger.handlers.clear()
    logging.shutdown()

    return exit_code


def get_tensorleed_path(tensorleed_path=None):
    """Return the path to the TensErLEED source code.

    Parameters
    ----------
    tensorleed_path : pathlike, optional
        Path to the viperleed-tensorleed source code, by default None.
        If not given, tries to resolve the $VIPERLEED_TENSORLEED environment
        variable.

    Returns
    -------
    Path
        Path to the TensErLEED source code.

    Raises
    ------
    ValueError
        If neither the tensorleed_path argument nor the $VIPERLEED_TENSORLEED
        environment variable are set.
    """
    # if tensorleed arg is given, use that
    if tensorleed_path:
        return Path(tensorleed_path)
    # else check environment variable $VIPERLEED_TENSORLEED
    try:
        return Path(os.environ["VIPERLEED_TENSORLEED"])
    except KeyError as exc:
        # environment variable not set
        raise ValueError(
            "TensErLEED path not specified.\n"
            "Please either pass a path to the TensErLEED source code with "
            "the --tensorleed argument, or set the environment variable "
            "$VIPERLEED_TENSORLEED."
            ) from exc
