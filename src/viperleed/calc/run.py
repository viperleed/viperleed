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
from viperleed.calc import LOGGER
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.files import parameters
from viperleed.calc.files import poscar
from viperleed.calc.files.manifest import ManifestFile
from viperleed.calc.files.tenserleed import get_tensorleed_path
from viperleed.calc.lib.log_utils import close_all_handlers
from viperleed.calc.lib.log_utils import prepare_calc_logger
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.sections.cleanup import cleanup
from viperleed.calc.sections.cleanup import prerun_clean
from viperleed.calc.sections.cleanup import preserve_original_inputs
from viperleed.calc.sections.initialization import (
    warn_if_slab_has_atoms_in_multiple_c_cells
    )
from viperleed.calc.sections.run_sections import section_loop


def run_calc(
    system_name=None,
    console_output=True,
    slab=None,
    preset_params=None,
    source=None,
    home=None,
    ):
    """Run a ViPErLEED calculation in the current directory.

    By default, a PARAMETERS and a POSCAR file are expected, but can be
    replaced by passing the `slab` and/or `preset_params` kwargs.

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
    home : pathlike, optional
        Path to the folder in which viperleed.calc was originally
        executed, i.e., before moving to the work directory. If
        not given or None, take the parent of the current directory.
        Default is None.

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

    # Start logger, write to file:
    timestamp = DateTimeFormat.FILE_SUFFIX.now()
    log_name = f'{LOG_PREFIX}-{timestamp}.log'
    prepare_calc_logger(LOGGER,
                        file_name=log_name,
                        with_console=console_output)
    LOGGER.info(f'Starting new log: {log_name}\nTime of execution: '
                + DateTimeFormat.LOG_CONTENTS.now())
    LOGGER.info(f'This is ViPErLEED version {__version__}\n')

    manifest = ManifestFile(DEFAULT_SUPP, DEFAULT_OUT, log_name)
    try:
        # Read input files and load user arguments
        rpars, slab = _make_rpars_and_slab(manifest, preset_params, slab, home)
    except Exception:
        _finalize_on_early_exit(manifest)
        return 2, None

    # Load runtime information in rpars
    rpars.timestamp = timestamp
    rpars.manifest = manifest
    _set_tensorleed_source(rpars, source)
    _set_system_name(rpars, system_name)

    # Check if halting condition is already in effect
    if rpars.halt >= rpars.HALTING:
        LOGGER.info('Halting execution...')
        _finalize_on_early_exit(rpars)
        return 0, None

    rpars.updateDerivedParams()
    LOGGER.info(f'ViPErLEED is using TensErLEED version {rpars.TL_VERSION}.')

    prerun_clean(rpars, log_name)
    preserve_original_inputs(rpars)  # Store inputs BEFORE any edit!
    exit_code, state_recorder = section_loop(rpars, slab)

    # Prevent other sub-loggers from producing more
    # messages for the main log file of viperleed calc.
    close_all_handlers(LOGGER)
    return exit_code, state_recorder


def _finalize_on_early_exit(rpars_or_manifest):
    """Finish a calc execution before entering the `section_loop`."""
    cleanup(rpars_or_manifest)
    # Prevent other sub-loggers from producing more
    # messages for the main log file of viperleed calc.
    close_all_handlers(LOGGER)


def _get_parent_directory_name():
    """Return the name of the directory above the current one."""
    _system_name = Path.cwd().parent.name
    LOGGER.info('No system name specified. Using name of parent '
                f'directory: {_system_name}')
    return _system_name


def _make_rpars_and_slab(manifest, preset_params, slab, home):
    """Return an Rparams and a SurfaceSlab.

    Parameters
    ----------
    manifest : ManifestFile
        The manifest file for this calculation.
    preset_params : dict or None
        Values of PARAMETERS to replace those read from file.
    slab : SurfaceSlab or None
        A user-given slab. If given, no POSCAR file is read from the
        current directory.
    home : pathlike or None
        Path to the folder in which viperleed.calc was originally
        executed, i.e., before moving to the work directory. If
        not given or None, take the parent of the current directory.

    Returns
    -------
    rpars : Rparams
        The run parameters for this calculation, loaded with
        `preset_params` and ready to be used.
    slab : SurfaceSlab or None
        The slab for this calculation. Read from POSCAR and
        updated with the contents of `rpars`, unless this
        is a multi-domain calculation. None in the latter
        case.

    Raises
    ------
    Exception
        If reading PARAMETERS or, for a single-domain calculation,
        POSCAR fails.
    FileNotFoundError
        If PARAMETERS is not found in the current directory.
    FileNotFoundError
        If POSCAR is not found in the current directory for a
        single-domain calculation
    OSError
        If duplicating POSCAR to POSCAR_user fails.
    ParameterError
        If interpreting PARAMETERS fails.
    TypeError
        If a non-None `slab` is given in a multi-domain calculation.
    """
    rpars = _read_parameters_file(preset_params)

    # Check if this is going to be a domain search
    domains = 'DOMAIN' in rpars.readParams
    if not domains and slab is None:
        slab = _read_poscar_file(manifest)
    elif domains and slab is not None:
        # no POSCAR in main folder for domain searches
        raise TypeError('Cannot give slab argument '
                        'for a multi-domain calculation')

    # Store the directory in which the calculation originally started.
    # Must be done before interpreting PARAMETERS, as it is needed
    # for DOMAIN.
    if home is None:
        home = Path.cwd().parent
    rpars.paths.home = Path(home)

    # Interpret the PARAMETERS file and load presets
    _interpret_parameters(rpars, slab, preset_params or {})

    if not domains:
        warn_if_slab_has_atoms_in_multiple_c_cells(slab, rpars)
        slab.full_update(rpars)   # gets PARAMETERS data into slab
        rpars.fileLoaded['POSCAR'] = True
    return rpars, slab


def _interpret_parameters(rpars, slab, preset_params):
    """Interpret PARAMETERS read from file and apply presets."""
    try:
        # interpret the PARAMETERS file
        parameters.interpret(rpars, slab=slab, silent=False)
    except (parameters.errors.ParameterNeedsSlabError,
            parameters.errors.SuperfluousParameterError):
        # Domains calculation is the only case in which slab is None
        LOGGER.error('Main PARAMETERS file contains an invalid parameter '
                     'for a multi-domain calculation', exc_info=True)
        raise
    except parameters.errors.ParameterError:
        LOGGER.error('Exception while reading PARAMETERS file', exc_info=True)
        raise

    # Load parameter presets, overriding those in PARAMETERS
    try:
        rpars.update(preset_params)
    except (ValueError, TypeError):
        LOGGER.warning('Error applying preset parameters: ', exc_info=True)

    _set_log_level(rpars, preset_params)
    LOGGER.debug('PARAMETERS file was read successfully')


def _read_parameters_file(preset_params):
    """Return an Rparams read from PARAMETERS in the current directory."""
    try:
        return parameters.read()
    except FileNotFoundError:
        if preset_params:
            return Rparams()
        LOGGER.error('No PARAMETERS file found, and no preset parameters '
                     'passed. Execution will stop.')
        raise
    except Exception:
        LOGGER.error('Exception while reading PARAMETERS file', exc_info=True)
        raise


def _read_poscar_file(manifest):
    """Return a SurfaceSlab read from POSCAR in the current directory.

    If the POSCAR file in the current directory has not been used
    for a viperleed.calc execution before, it is duplicated as
    POSCAR_user.

    Parameters
    ----------
    manifest : ManifestFile
        The manifest file used for this calculation. Used to retain
        information on whether POSCAR_user was created.

    Returns
    -------
    slab : SurfaceSlab
        Slab read from POSCAR.

    Raises
    ------
    FileNotFoundError
        If no POSCAR file is found in the current directory.
    OSError
        If duplicating POSCAR to POSCAR_user fails.
    Exception
        If reading POSCAR fails.
    """
    poscar_file = Path('POSCAR')
    LOGGER.info('Reading structure from file POSCAR')
    try:
        slab = poscar.read(filename=poscar_file)
    except FileNotFoundError:
        LOGGER.error('POSCAR not found. Stopping execution...')
        raise
    except Exception:
        LOGGER.error('Exception while reading POSCAR', exc_info=True)
        raise

    if not slab.preprocessed:
        LOGGER.info('The POSCAR file will be processed and overwritten. '
                    'Copying the original POSCAR to POSCAR_user...')
        try:
            shutil.copy2(poscar_file, 'POSCAR_user')
        except OSError:
            LOGGER.error('Failed to copy POSCAR to POSCAR_user. Stopping '
                         'execution...')
            raise
        manifest.add('POSCAR_user')
    return slab


def _set_log_level(rpars, preset_params):
    """Assign a (user-defined) log level to the current logger."""
    # pylint: disable-next=magic-value-comparison
    if 'LOG_LEVEL' in preset_params:
        LOGGER.info(f'Overriding log level to {rpars.LOG_LEVEL}.')
    LOGGER.setLevel(rpars.LOG_LEVEL)


def _set_tensorleed_source(rpars, source):
    """Set `source` as the directory for tensor-LEED files/executables."""
    try:
        tensorleed = get_tensorleed_path(source).resolve()
    except (ValueError, FileNotFoundError) as exc:
        LOGGER.warning(f'{exc} This may cause errors.')
        tensorleed = Path(source or '').resolve()
    rpars.paths.tensorleed = tensorleed


def _set_system_name(rpars, system_name):
    """Assign a system name to rpars from `system_name` or the CWD."""
    if system_name is None:
        system_name = _get_parent_directory_name()
    rpars.systemName = system_name
