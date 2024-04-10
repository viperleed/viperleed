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

from viperleed import VIPERLEED_TENSORLEED_ENV
from viperleed import __version__
from viperleed.calc import LOGGER as logger
from viperleed.calc import LOG_PREFIX
from viperleed.calc.classes import rparams
from viperleed.calc.files import parameters, poscar
from viperleed.calc.lib.base import CustomLogFormatter
from viperleed.calc.sections.cleanup import cleanup
from viperleed.calc.sections.cleanup import prerun_clean
from viperleed.calc.sections.initialization import (
    warn_if_slab_has_atoms_in_multiple_c_cells
    )
from viperleed.calc.sections.run_sections import section_loop


_KNOWN_TENSORLEED_TOP_FOLDERS = {
    'viperleed-tensorleed',
    'tensorleed',  # For backwards compatibility
    }


def run_calc(system_name=None,
             console_output=True,
             slab=None,
             preset_params=None,
             source=None,
             inputs_dir=None,):
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

    # set input directory
    rp.inputs_dir = inputs_dir or Path.cwd().resolve()

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
        If not given, the $VIPERLEED_TENSORLEED environment variable
        is used, if defined. As a last resort, the current directory
        is used. In all cases, the path given may be the top-level tree
        containing the tensor-LEED code, its parent, or a first-level
        subfolder.

    Returns
    -------
    Path
        Path to the TensErLEED source code, that is, the directory
        that contains, e.g., potentially multiple TensErLEED sources
        as well as the EEASiSSS source code and compiled binaries.

    Raises
    ------
    ValueError
        If neither the tensorleed_path argument nor the
        VIPERLEED_TENSORLEED environment variable are set.
    FileNotFoundError
        If neither tensorleed_path (or VIPERLEED_TENSORLEED), its
        parent, or any of its children, point to one of the known
        tensor-LEED directory structures.
    """
    if tensorleed_path is not None:
        path_ = Path(tensorleed_path)
        return _verify_tensorleed_path(path_)

    # Check environment variable $VIPERLEED_TENSORLEED
    try:
        path_ = Path(os.environ[VIPERLEED_TENSORLEED_ENV])
    except KeyError:  # Environment variable not set. Try CWD below.
        pass
    else:
        return _verify_tensorleed_path(path_)

    # Last resort: try with the current directory
    try:
        return _verify_tensorleed_path(Path.cwd())
    except ValueError:
        raise ValueError(
            'TensErLEED path not specified.\n'
            'Either pass a path to the TensErLEED source code with the '
            '--tensorleed command-line argument, or set the environment '
            f'variable {VIPERLEED_TENSORLEED_ENV}.'
            ) from None


def _get_parent_directory_name():
    """Return the name of the directory above the current one."""
    _system_name = Path.cwd().resolve().parent.name
    logger.info('No system name specified. Using name of parent '
                f'directory: {_system_name}')
    return _system_name


def _verify_tensorleed_path(path_):
    """Return the path to the tensor-LEED folder, starting from `path_`.

    Parameters
    ----------
    path_ : Path
        Path where the tensor-LEED code is looked up.

    Returns
    -------
    verified_path : Path
        Path to the folder that contains tensor-LEED code. The
        following places are checked: `path_`, its containing
        folder, and all of its direct subfolders. `verified_path`
        is identified as a tensor-LEED containing folder if its
        name is one of the known tensor-LEED repository names
        or if it contains any subfolder whose name starts with
        'TensErLEED'.

    Raises
    ------
    FileNotFoundError
        If `path_` is not the path to an existing directory,
        or if finding a verified tensor-LEED path failed.
    """
    source = path_.resolve()
    if not source.is_dir():
        raise FileNotFoundError(f'Tensor-LEED directory {path_} not found.')

    # Collect various potential source-code places: the given path_,
    # its containing folder, and all of the path_'s direct subfolders
    potential_sources = (
        source,
        source.parent,
        *(d for d in source.iterdir() if d.is_dir())
        )

    # First look for an exact name match with
    # the known top-level tensor-LEED folders
    try:
        source_by_name = next(d for d in potential_sources
                              if d.name in _KNOWN_TENSORLEED_TOP_FOLDERS)
    except StopIteration:  # Not a known folder name
        pass
    else:
        if source_by_name != source:
            logger.warning(
                f'{source_by_name.name!r} directory found in '
                f'{source_by_name.parent}. Using {source_by_name} '
                f'instead of {source}.'
                )
        return source_by_name

    # Not found by name. Try to find if one contains a
    # TensErLEED directory instead. In principle here we
    # could use Rparams.get_tenserleed_directory() but we
    # do not care about versions and globbing is enough.
    try:
        with_tenserleed = (
            d for d in potential_sources
            if any(d.glob(f'{rparams.TENSERLEED_FOLDER_NAME}*'))
            )
    except StopIteration:
        raise FileNotFoundError('Could not find a known tensor-LEED '
                                f'source directory at {source}.') from None
    if with_tenserleed != source:
        logger.warning(
            f'TensErLEED code found in {with_tenserleed}. '
            f'Using the latter instead of {source}.'
            )
    return with_tenserleed
