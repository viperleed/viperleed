"""Module cleanup of viperleed.calc.sections.

Defines clean-up functions, to be used between
sections or before/after execution.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-06-04'
__license__ = 'GPLv3+'

import logging
from pathlib import Path
import re
import shutil
from zipfile import ZIP_DEFLATED
from zipfile import ZipFile

from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.constants import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.base import copytree_exists_ok
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.log_utils import close_all_handlers
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.calc_section import EXPBEAMS_NAMES


# Files to go in SUPP
_SUPP_FILES = (
    'AUXBEAMS',
    'AUXEXPBEAMS',
    'AUXGEO',
    'AUXLATGEO',
    'AUXNONSTRUCT',
    'BEAMLIST',
    'delta-input',
    'EEASISSS-input.txt',
    'eeasisss-input',
    'EEASISSS-log.txt',
    'muftin.f',
    'Phaseshifts_plots.pdf',
    'POSCAR_bulk_appended',
    'POSCAR_bulk',
    'POSCAR_mincell',
    'POSCAR_oricell',
    'POSCAR_vacuum_corrected',
    'refcalc-FIN',
    'refcalc-PARAM',
    'restrict.f',
    'rfactor-PARAM',
    'rfactor-WEXPEL',
    'search-PARAM',
    'search-rf.info',
    'search.steu',
    'searchpars.info',
    'superpos-CONTRIN',
    'superpos-PARAM',
    'VIBROCC_generated',
    )

_SUPP_DIRS = (
    ORIGINAL_INPUTS_DIR_NAME,
    'compile_logs',
    )

# Files that may be generated automatically and do not need
# storage into original_inputs.
OPTIONAL_INPUT_FILES = (
    'BEAMLIST',
    )

# Files to go in OUT
_OUT_FILES = (
    'Complex_amplitudes_imag.csv',
    'Complex_amplitudes_real.csv',
    'control.chem',
    'Errors_summary.csv',
    'Errors.pdf',
    'Errors.zip',
    'FD_Optimization_beams.pdf',
    'FD_Optimization.csv',
    'FD_Optimization.pdf',
    'FITBEAMS_norm.csv',
    'FITBEAMS.csv',
    'PatternInfo.tlm',
    'refcalc-amp.out',
    'Rfactor_analysis_refcalc.pdf',
    'Rfactor_analysis_superpos.pdf',
    'Rfactor_plots_refcalc.pdf',
    'Rfactor_plots_superpos.pdf',
    'SD.TL',
    'refcalc-fd.out',
    'Search-progress.csv',
    'Search-progress.pdf',
    'Search-report.pdf',
    'superpos-spec.out',
    'THEOBEAMS_norm.csv',
    'THEOBEAMS.csv',
    'THEOBEAMS.pdf',
    )

# Label given to workhistory folders when cleaning up stray remains
# from previous viperleed.calc executions from the work directory
PREVIOUS_LABEL = 'previous'

# Output files that may be inputs in future runs - keep during prerun
_IOFILES = (
    'control.chem',
    'refcalc-fd.out',
    'superpos-spec.out',
    )

logger = logging.getLogger(__name__)


def prerun_clean(rpars, logname=''):
    """Clean up the current directory before viperleed.calc starts.

    Delete workhistory, old executables, and old logfiles.
    Call move_oldruns if required.

    Parameters
    ----------
    rpars : Rparams
        The run parameters, needed for move_oldruns.
    logname : str, optional
        Name of the current log file, to be excluded from cleanup.

    Returns
    -------
    None.
    """
    _delete_old_root_directories()  # workhistory, SUPP, OUT, ...
    _delete_out_suffixed_files()    # POSCAR_OUT, etc.
    _delete_old_executables()       # refcalc, etc.

    # If there are old log files, move inputs/outputs to workhistory
    old_logs = (f for f in Path().glob('*.log')
                if f.is_file() and f.name != logname)
    if any(old_logs):
        try:
            move_oldruns(rpars, prerun=True)
        except OSError:
            logger.warning('Exception while trying to clean up from previous '
                           'run. Program will proceed, but old files may be '
                           'lost.', exc_info=True)

    # Get rid of other files that may be modified in place
    other_logs = (
        'fortran-compile.log',  # We append to this
        )
    _silently_remove_files(*other_logs)


def organize_workdir(rpars,
                     path,
                     delete_unzipped=False,
                     tensors=True,
                     deltas=True):
    """Reorganize files in path into SUPP, OUT, Tensors and Deltas.

    Tensors and Deltas folders are zipped and moved over. All other
    files are copied to appropriate locations in SUPP and OUT.

    Parameters
    ----------
    rpars : Rparams
        The run parameters associated with the calculation that
        is running in the directory to be cleaned up. Attributes
        accessed (read-only):
            TENSOR_INDEX:
                Picks the name of the Deltas_xxx folder.
            ZIP_COMPRESSION_LEVEL:
                Which compression level to use to create ZIP
                archives for Deltas/Tensors.
    path : pathlike
        The path to the work folder that contains the files to
        be reorganized.
    delete_unzipped : bool, optional
        Whether the original Delta- and Tensor-files should be
        deleted after making the archives. The default is False.
    tensors, deltas : bool, optional
        Whether the Tensor/Delta files contain new information
        and should be saved. The default is True.

    Returns
    -------
    None.
    """
    with execute_in_dir(path):
        _collect_delta_files(rpars.TENSOR_INDEX or 0)

        # Create ZIP files for Tensors and Deltas subfolders
        zip_args = delete_unzipped, rpars.ZIP_COMPRESSION_LEVEL
        if tensors or delete_unzipped:
            _zip_subfolders(DEFAULT_TENSORS, tensors, *zip_args)
        if deltas or delete_unzipped:
            _zip_subfolders(DEFAULT_DELTAS, deltas, *zip_args)

        _collect_supp_contents(rpars)
        _collect_out_contents(rpars)


def _collect_supp_contents(rpars):
    """Store relevant files/folder from the current directory to SUPP."""
    files_to_copy = set(Path(f) for f in _SUPP_FILES
                        if f not in rpars.files_to_out)
    directories_to_copy = (Path(d) for d in _SUPP_DIRS)

    # Also add log files into SUPP: skip calc logs (they go to
    # main dir), and compile logs (they go to compile_logs dir)
    logs_to_supp = (
        f for f in Path().glob('*.log')
        # pylint: disable-next=magic-value-comparison  # 'compile'
        if not f.name.startswith(LOG_PREFIX) and 'compile' not in f.name
        )
    files_to_copy.update(logs_to_supp)

    _copy_files_and_directories(files_to_copy,
                                directories_to_copy,
                                Path(DEFAULT_SUPP))


def _collect_out_contents(rpars):
    """Store relevant files/folder from the current directory to OUT."""
    out_path = Path(DEFAULT_OUT)
    out_files = set(Path(f) for f in _OUT_FILES)
    # Add R-factor output files
    out_files.update(Path().glob('R_*R=*'))
    # And POSCAR, PARAMETERS, VIBROCC files that we generated/edited.
    # They may be the ones created at initialization, or those from
    # an optimization.
    out_files.update(Path(f) for f in rpars.files_to_out)
    _copy_files_and_directories(out_files, (), out_path)


def _copy_files_and_directories(files, directories, target):
    """Copy files and directories to target, creating it if not existing."""
    try:
        target.mkdir(parents=True)
    except FileExistsError:
        pass
    except OSError:
        logger.error(f'Error creating {target.name} folder: ', exc_info=True)
        return

    for item in (*files, *directories):
        _copy = shutil.copy2 if item.is_file() else copytree_exists_ok
        try:
            _copy(item, target/item.name)
        except FileNotFoundError:
            pass
        except OSError:
            which = 'file' if item.is_file() else 'directory'
            logger.error(f'Error moving {target.name} {which} {item.name}: ',
                         exc_info=True)


def _zip_subfolders(at_path, archive, delete_unzipped, compression_level):
    """Archive all numbered subfolders `at_path`.

    Parameters
    ----------
    at_path : Path
        The folder containing the subfolders to be archived into
        a ZIP file. Only subfolders whose names begin with
        '<at_path.name>_ddd' are packed.
    archive : bool
        Whether subfolders should be archived or only deleted.
    delete_unzipped : bool
        Whether subfolders should be deleted after they have been
        successfully archived.
    compression_level : int
        Compression level to be applied while archiving.

    Returns
    -------
    None.
    """
    at_path = Path(at_path)
    if not at_path.is_dir():
        return
    root_name = at_path.name
    rgx = re.compile(rf'{root_name}_[0-9]{{3,}}')
    subfolders = (p for p in at_path.glob('*') if p.is_dir())
    for subfolder in subfolders:
        if not rgx.fullmatch(subfolder.name):
            continue
        if archive:
            try:
                _zip_folder(subfolder, compression_level)
            except OSError:
                continue

        if delete_unzipped:
            try:
                shutil.rmtree(subfolder)
            except OSError:
                logger.warning(
                    f'Error deleting unzipped {root_name} directory '
                    f'{subfolder}. This will increase the size of the '
                    'work folder, but not cause any problems.'
                    )


def _zip_folder(folder, compression_level):
    """Create a ZIP with the same name as folder.

    Parameters
    ----------
    folder : Path
        Path to the folder to be compressed. The archive will
        be saved in the parent of `folder`, with the same name.
        If the ZIP file already exists, the contents of `folder`
        are added to it.
    compression_level : int
        The level of compression to use when creating the ZIP.

    Raises
    ------
    OSError
        If creating the archive fails.
    """
    kwargs = {'compression': ZIP_DEFLATED, 'compresslevel': compression_level}
    arch_name = folder.with_suffix('.zip')
    logger.info(f'Packing {arch_name}...')
    # Don't pack the archive into itself
    to_pack = (f for f in folder.glob('*') if f != arch_name)
    try:  # pylint: disable=too-many-try-statements
        with ZipFile(arch_name, 'a', **kwargs) as archive:
            for item in to_pack:
                archive.write(item, item.relative_to(folder))
    except OSError:
        logger.error(f'Error packing {arch_name} file: ', exc_info=True)
        raise


def _collect_delta_files(tensor_index):
    """Move all 'DEL_' files in the current directory into a Deltas folder.

    Parameters
    ----------
    tensor_index : int
        The index of the Tensors that were used to generate these
        Deltas. Used to label the Deltas/Deltas_<index> folder.

    Returns
    -------
    None.
    """
    deltas = tuple(Path().glob('DEL_*'))
    if not deltas:
        return
    destination = Path(f'{DEFAULT_DELTAS}/{DEFAULT_DELTAS}_{tensor_index:03d}')
    try:
        destination.mkdir(parents=True)
    except FileExistsError:
        pass
    except OSError:
        logger.error(f'Failed to create {destination} folder: ',
                     exc_info=True)
        return
    errors = []
    for delta_file in deltas:
        try:
            delta_file.replace(destination/delta_file)
        except OSError as exc:
            errors.append(exc)
    if errors:
        logger.error(f'Error moving Delta files: {errors}')


def move_oldruns(rpars, prerun=False):
    """Copy relevant files to a new 'workhistory' subfolder.

    Files are copied from SUPP, OUT, and those in rpars.manifest.
    The main log file is excluded.

    Parameters
    ----------
    rpars : Rparams
        The run parameters.
    prerun : bool, optional
        If True, instead of using the manifest, all potentially
        interesting files will be copied. The new subfolder is
        indexed as 0. Then, SUPP, OUT, and old SUPP/OUT files
        are cleared from the main directory.

    Raises
    ------
    OSError
        If creation of workhistory or its subfolder fails.
    """
    worhistory_subfolder = _make_new_workhistory_subfolder(rpars, prerun)
    if not prerun:
        rpars.manifest.add(DEFAULT_WORK_HISTORY)
        kwargs = {'delete_unzipped': False,
                  'tensors': False,
                  'deltas': False}
        organize_workdir(rpars, path='', **kwargs)
        for domain in rpars.domainParams:
            organize_workdir(domain.rp, path=domain.workdir, **kwargs)
    _collect_worhistory_contents(rpars, prerun, worhistory_subfolder)


def _collect_worhistory_contents(rpars, prerun, to_path):
    """Copy or move files/directories to a workhistory subfolder."""
    files, directories = _find_next_workistory_contents(rpars, prerun)
    _copy = shutil.move if prerun else shutil.copy2
    for file in files:
        _copyfile = shutil.copy2 if file in _IOFILES else _copy
        try:
            _copyfile(file, to_path)
        except OSError:
            logger.warning(f'Error copying {file} to {to_path}. '
                           'File may get overwritten.')
    _copy = shutil.move if prerun else shutil.copytree
    for directory in directories:
        try:
            _copy(directory, to_path / directory)
        except OSError:
            logger.warning(f'Error copying {directory} to {to_path}. '
                           'Files in directory may get overwritten.')


def _find_next_workistory_contents(rpars, prerun):
    """Return files/folders for a fresh workhistory directory."""
    all_dirs = (f.name for f in Path().iterdir() if f.is_dir())
    all_files = (f.name for f in Path().iterdir() if f.is_file())
    if prerun:
        # Skip manifest, generated, and IO files. Take all logs, as
        # well as SUPP and OUT directories. Also take root files that
        # may have already been copied to SUPP/OUT: the sole purpose
        # is **removing them** from the root directory (via shutil.move
        # in _collect_worhistory_contents).
        files = [
            f for f in all_files
            if f not in rpars.manifest
            and f not in rpars.files_to_out
            and f not in _IOFILES
            and (f.endswith('.log') or f in _OUT_FILES or f in _SUPP_FILES)
            ]
        directories = [d for d in all_dirs if d in (DEFAULT_SUPP, DEFAULT_OUT)]
    else:
        # Take only files from manifest, and all directories
        # that are not potentially used in subsequent runs
        _calc_log = re.compile(rf'{LOG_PREFIX}.*\.log')
        _skip_dirs = {DEFAULT_TENSORS, DEFAULT_DELTAS, DEFAULT_WORK_HISTORY}
        files = [f for f in all_files
                 if f in rpars.manifest and not _calc_log.fullmatch(f)]
        directories = [d for d in all_dirs
                       if d in rpars.manifest and d not in _skip_dirs]
    return files, directories


def _find_next_workistory_dir_name(rpars, prerun):
    """Return the name of a fresh workhistory subfolder."""
    run_number = _find_next_workistory_run_number(rpars, prerun)
    dirname_prefix = f't{rpars.TENSOR_INDEX:03d}.r{run_number:03d}'
    if prerun:
        try:
            most_recent_log = max(
                f.name
                for f in Path().glob(f'{LOG_PREFIX}*.log')
                if f.is_file() and f.name not in rpars.manifest
                )
        except ValueError:  # No relevant log files
            old_timestamp = f'moved-{rpars.timestamp}'
        else:
            old_timestamp = most_recent_log[-17:-4]
        return f'{dirname_prefix}_{PREVIOUS_LABEL}_{old_timestamp}'

    sectionabbrv = {1: 'R', 2: 'D', 3: 'S'}
    new_segments = rpars.runHistory[len(rpars.lastOldruns):]
    abbreviations = ''.join(sectionabbrv.get(index, '')
                            for index in new_segments)
    if abbreviations:
        abbreviations = '_' + abbreviations
    rpars.lastOldruns = rpars.runHistory[:]
    return f'{dirname_prefix}{abbreviations}_{rpars.timestamp}'


def _find_next_workistory_run_number(rpars, prerun):
    """Return a numeric identifier for a fresh workhistory subfolder."""
    workhistory = Path(DEFAULT_WORK_HISTORY)
    try:
        subfolders = tuple(d.name for d in workhistory.iterdir() if d.is_dir())
    except FileNotFoundError:
        subfolders = ()
    # Keep only the subfolders for the specific TENSOR_INDEX
    rgx = re.compile(rf't{rpars.TENSOR_INDEX:03d}.r(?P<run>[0-9]{{3,}})_')
    matches = (rgx.match(f) for f in subfolders)
    run_numbers = (int(m['run']) for m in matches if m)
    try:
        return max(run_numbers) + 1
    except ValueError:  # No matching subfolder
        return 0 if prerun else 1


def _make_new_workhistory_subfolder(rpars, prerun):
    """Create a fresh subfolder of workhistory for storing previous results."""
    workhistory = Path(DEFAULT_WORK_HISTORY)
    try:
        workhistory.mkdir(exist_ok=True)
    except OSError:
        logger.error(f'Error creating {workhistory} folder: ', exc_info=True)
        raise

    dirname = _find_next_workistory_dir_name(rpars, prerun)
    subfolder = workhistory / dirname
    try:
        subfolder.mkdir()
    except OSError:
        logger.error(f'Error creating {subfolder}: ', exc_info=True)
        raise
    return subfolder


def cleanup(manifest, rpars=None):
    """Finalize a viperleed.calc execution.

    After a call to this function:
    - Files in the current directory and all domain work directories
      are organized into SUPP, OUT, Tensors, and Deltas folders
    - The manifest file is written in the current directory. It
      contains information about which files should be moved back
      to the original folder where calc was started.
    - Final messages are written to the log, including information
      about duration of the overall calculation, the segments that
      were executed, and the final R factors. Warnings concerning
      the bookkeeper and a checklist of items for user are also
      logged.

    The logging module is fully shut down after this function
    is executed, and should not be used any longer.

    Parameters
    ----------
    manifest : set of str
        The files and directories that should be preserved from
        the work folder.
    rpars : Rparams, optional
        The run parameters. If None, it is assumed that the run
        crashed before an Rparams object existed.

    Returns
    -------
    None.
    """
    logger.info('\nStarting cleanup...')
    if rpars is None:  # Make a dummy, essentially empty one
        rpars = Rparams()
        rpars.manifest = manifest
        rpars.timer = None  # To print the correct final message

    _organize_all_work_directories(rpars)
    _write_manifest_file(manifest)
    _write_final_log_messages(rpars)

    # Shut down logger
    close_all_handlers(logger)
    logging.shutdown()


def preserve_original_inputs(rpars):
    """Create the original_inputs directory and copy input files there.

    The original_inputs directory is created in the current directory
    if it does not exist yet. Input files are also taken from the
    current directory. Notice that all potentially relevant input
    files are stored, irrespective of whether they are used in the
    calculation or not.

    Parameters
    ----------
    rpars : Rparams
        The current run parameters. Used only for error reporting.

    Raises
    ------
    OSError
        If creating the directory fails.
    """
    orig_inputs = Path(ORIGINAL_INPUTS_DIR_NAME).resolve()
    try:
        orig_inputs.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise OSError(f'Could not create directory {orig_inputs}. '
                      'Check disk permissions.') from exc

    # We will copy all files that have potentially been used as
    # inputs. Make sure the correct version of EXPBEAMS is stored
    files_to_preserve = ALL_INPUT_FILES - set(EXPBEAMS_NAMES)
    try:
        files_to_preserve.add(
            next(f for f in EXPBEAMS_NAMES if Path(f).is_file())
            )
    except StopIteration:  # No EXPBEAMS
        pass

    for filename in files_to_preserve:
        file = Path(filename)
        if not file.is_file() and filename in OPTIONAL_INPUT_FILES:
            continue
        if not file.is_file():
            logger.warning(f'Could not find file {file}. It will not '
                           f'be stored in {ORIGINAL_INPUTS_DIR_NAME}.')
            rpars.setHaltingLevel(1)
            continue
        try:
            shutil.copy2(file, orig_inputs)
        except OSError:
            logger.warning(f'Could not copy file {file} to '
                           f'{ORIGINAL_INPUTS_DIR_NAME}.')
            rpars.setHaltingLevel(1)


def _delete_old_executables():
    """Remove compiled executables from the current directory."""
    executables = (  # They have a timestamp.
        'refcalc',
        'rfactor',
        'search',
        'superpos',
        )
    for file_name in executables:
        pattern = re.compile(file_name + r'-\d{6}-\d{6}')
        for file in Path().glob(f'{file_name}-*'):
            if not file.is_file() or not pattern.fullmatch(file.stem):
                continue
            try:
                file.unlink()
            except OSError:
                logger.debug(f'Failed to delete file {file}')


def _delete_old_root_directories():
    """Remove calc-created directories from the current directory."""
    directories = (
        DEFAULT_WORK_HISTORY,
        DEFAULT_SUPP,
        DEFAULT_OUT,
        )
    for directory in directories:
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            pass
        except OSError:
            logger.warning(f'Failed to clear {directory} folder.')


def _delete_out_suffixed_files():
    """Remove all files containing '_OUT' from the current directory."""
    for file in Path().glob('*_OUT*'):
        try:
            file.unlink()
        except OSError:
            logger.warning(f'Failed to delete previous {file} file.')


def _organize_all_work_directories(rpars):
    """Collect files from the current directory and all domain ones.

    After calling this function, files in both the current directory
    and those in the work directories of all subdomains are collected
    into their respective SUPP, OUT, Tensors, and Deltas folders.

    Parameters
    ----------
    rpars : Rparams
        The run parameters of the main calculation.

    Returns
    -------
    None.
    """
    rpars.closePdfReportFigs()
    to_sort = [{'tensors': DEFAULT_TENSORS in rpars.manifest,
                'deltas': DEFAULT_DELTAS in rpars.manifest,
                'rpars': rpars,
                'path': ''}]
    to_sort.extend(
        {'tensors': DEFAULT_TENSORS in dp.rp.manifest,
         'deltas': DEFAULT_DELTAS in dp.rp.manifest,
         'rpars': dp.rp,
         'path': dp.workdir}
        for dp in rpars.domainParams
        )
    for kwargs in to_sort:
        organize_workdir(delete_unzipped=True, **kwargs)


def _silently_remove_files(*files):
    """Delete files from this directory without complaining for errors."""
    for file in files:
        file = Path(file)
        try:
            file.unlink()
        except OSError:
            pass


def _write_final_log_messages(rpars):
    """Emit the last logging messages concerning the calculation."""
    elapsed = ('unknown' if not rpars.timer
               else rpars.timer.how_long(as_string=True))
    logger.info(f'\nFinishing execution at {DateTimeFormat.LOG_CONTENTS.now()}'
                f'\nTotal elapsed time: {elapsed}\n')

    # Write information about executed sections
    if rpars.runHistory:
        segments = ' '.join(str(s) for s in rpars.runHistory)
        logger.info(f'Executed segments: {segments}')

    # Write the final R factors, if any, including integer/fractional
    for section, r_factors in rpars.stored_R.items():
        if r_factors is None:
            continue
        overall, integer, fractional = r_factors
        msg = f'Final R ({section}): {overall:.4f}'
        if integer > 0 and fractional > 0:
            msg += f' ({integer:.4f} / {fractional:.4f})'
        logger.info(msg)

    # Warn about manually running bookkeeper for domain calculations
    if rpars.domainParams:
        logger.info(
            'Domain calculations have been run. Note that the bookkeeper will '
            'only run automatically in the top level calculation directory. '
            'To preserve optimizations for individual domains, please run '
            'bookkeeper manually in the respective domain directories. '
            'The command is: viperleed bookkeeper --archive.\n'
            )

    if rpars.checklist:
        logger.info('')
        logger.info('# The following issues should be '
                    'checked before starting again:')
        for item in rpars.checklist:
            logger.info(f'- {item}')
    logger.info('')


def _write_manifest_file(manifest_contents):
    """Write items in manifest_contents to file 'manifest'."""
    manifest_contents = set(manifest_contents)
    manifest = Path('manifest')
    try:
        manifest.write_text('\n'.join(manifest_contents) + '\n',
                            encoding='utf-8')
    except OSError:
        logger.error(f'Failed to write {manifest} file.')
    else:
        logger.info(f'Wrote {manifest} file successfully.')
