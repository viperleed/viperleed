"""Section Delta Amplitudes."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import hashlib
import itertools
import logging
from pathlib import Path
import shutil
import subprocess

import numpy as np

from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files import iodeltas
from viperleed.calc.files import iotensors
from viperleed.calc.files.displacements import readDISPLACEMENTS_block
from viperleed.calc.lib import leedbase
from viperleed.calc.lib import parallelization
from viperleed.calc.lib.checksums import validate_multiple_files
from viperleed.calc.lib.context import execute_in_dir

logger = logging.getLogger(__name__)

# TODO: would be nice to replace all os.path with pathlib


class DeltasError(Exception):
    """Base exception for delta-amplitudes-related errors."""


# TODO: see note in refcalc. #43
class DeltaCompileTask:
    """Information to compile a delta-amplitudes executable.

    Attributes
    ----------
    exename : str
        File name of the executable that will be compiled.
    foldername : str
        Name of the folder in which `exename` can be found
        after successful compilation.
    fortran_comp : tuple
        Compiler and compilation flags.
    hash : str
        A unique identifier for this task. Usually calculated
        from `param`.
    param : str
        Contents of the PARAM file, defining array dimensions
        for compilation.
    source_dir : Path
        Path to the folder containing the static Fortran
        source files to be compiled.
    """

    def __init__(self, param, source_dir, index):
        """Initialize instance.

        Parameters
        ----------
        param : str
            Contents of the PARAM file, defining array dimensions
            for compilation.
        source_dir : Path
            Path to the folder containing the static Fortran
            source files to be compiled.
        index : int
            A progressive identifier for the executable to be
            compiled. This value affects both the .foldername
            and the .exename attributes.

        Returns
        -------
        None.
        """
        self.exename = f'delta-{index}'
        self.foldername = f'Delta_Compile_{index}'
        self.fortran_comp = ['', '']
        self.param = param
        self.source_dir = Path(source_dir).resolve()  # where the fortran files are

        if os.name == 'nt':
            self.exename += '.exe'

    def __str__(self):
        """Return a string representation of this DeltaCompileTask."""
        return f'{type(self).__name__} {self.foldername}'

    @property
    def logfile(self):
        """Return the (relative) path to the compilation log file."""
        return Path(self.foldername) / 'fortran-compile.log'

    @property
    def compile_log_name(self):
        """Name of the log file as it should appear in compile_logs."""
        return self.foldername

    def get_source_files(self):
        """Return files needed for a delta-amplitude compilation."""
        srcpath = self.source_dir / 'src'
        srcname = next(srcpath.glob('delta*'), None)
        libpath = self.source_dir / 'lib'
        lib_tleed = next(libpath.glob('lib.tleed*'), None)
        lib_delta = next(libpath.glob('lib.delta*'), None)
        globalname = srcpath / 'GLOBAL'
        if any(f is None for f in (srcname, lib_tleed, lib_delta)):
            raise RuntimeError(f'Source files missing in {self.source_dir}')    # TODO: use a better custom exception in CompileTask (e.g., MissingSourceFileError)
        return srcname, lib_tleed, lib_delta, globalname

    def copy_source_files_to_local(self):
        """Copy delta source files to current directory."""
        for filepath in self.get_source_files():
            if filepath:
                shutil.copy2(filepath, filepath.name)


# TODO: see note in refcalc.run_refcalc. #43
class DeltaRunTask:
    """Information for executing a delta-amplitudes calculation.

    Attributes
    ----------
    comptask : DeltaCompileTask
        The task that was used to create an executable for this
        calculation.
    deltalogname : str
        Name of the log file to which runtime information is
        collected.
    deltaname : str
        Name of the main output file of this calculation.
    din : str
        Contents of the input piped to the delta-amplitudes
        calculation executable.
    din_short : str
        Similar to `din`, but with repeated input replaced with
        appropriate tags.
    tensorname : str
        File name for the Tensors file that this calculation uses.
    """

    def __init__(self, comptask):
        """Initialize instance from a compilation task."""
        self.comptask = comptask
        self.deltalogname = ''
        self.deltaname = ''
        self.din = ''
        self.din_short = ''
        self.tensorname = ''

    @property
    def foldername(self):
        """Return a name for the subfolder in which this should run."""
        return f'calculating_{self.name}'

    @property
    def name(self):
        """Return a name for this task."""
        return self.deltaname

    def __str__(self):
        """Return a string representation of this DeltaRunTask."""
        return f'{type(self).__name__} {self.name}'


# TODO: see note in refcalc.run_refcalc. #43
def run_delta(runtask):
    """Execute a single delta-amplitude DeltaRunTask.

    The calculation is executed in a subfolder of the current working
    directory. The subfolder is removed upon successful execution.
    This function is meant to be executed by parallelized workers.

    Parameters
    ----------
    runtask : DeltaRunTask
        Information about the delta-amplitudes calculation
        to be executed.

    Returns
    -------
    error_info : str
        A message with information about errors occurred during
        execution of `runtask`.
    """
    base = Path.cwd()
    workfolder = base / runtask.foldername
    # Make folder and go there:
    try:
        workfolder.mkdir()
    except FileExistsError:
        logger.warning(f'Folder {workfolder.name} already exists. '
                       'Contents may get overwritten.')
    with execute_in_dir(workfolder):
        # Collect tensor file
        tensor = base / DEFAULT_TENSORS / runtask.tensorname
        try:
            shutil.copy2(tensor, 'AMP')
        except FileNotFoundError:
            logger.error(f'{DEFAULT_TENSORS} file not '
                         f'found: {runtask.tensorname}')
            return (f'Error encountered by {runtask}: '
                    f'{DEFAULT_TENSORS} not found.')
        except OSError:
            logger.error(f'Error copying {DEFAULT_TENSORS} file: ',
                         exc_info=True)
            return (f'Error encountered by {runtask}: Error '
                    f'copying {DEFAULT_TENSORS} file.')

        # Fetch compiled executable
        exename = runtask.comptask.exename
        try:
            shutil.copy2(base / runtask.comptask.foldername / exename, '.')
        except OSError:
            logger.error('Error getting delta executable: ', exc_info=True)
            return (f'Error encountered by {runtask}: '
                    'Failed to get delta executable.')

        # Run delta-amplitudes calculation
        log_file = Path('delta.log')
        with log_file.open('w', encoding='utf-8') as log:
            try:
                subprocess.run(str(workfolder / exename),
                               input=runtask.din,
                               encoding='ascii',
                               stdout=log,
                               stderr=log,
                               check=False)
            except Exception:
                logger.error('Error while executing delta-amplitudes '
                             f'calculation for {runtask.name}. Also '
                             'check delta log file.',
                             exc_info=True)
                return (f'Error encountered by {runtask}: '
                        'Error during delta execution.')

        # Copy delta file out
        try:
            shutil.copy2('DELWV', base / runtask.deltaname)
        except OSError:
            logger.error('Failed to copy delta output file DELWV to main '
                         f'folder as {runtask.deltaname}', exc_info=True)
            return (f'Error encountered by {runtask}: '
                    'Failed to copy result file out.')

        # Collate log results
        log = ''
        try:
            log = log_file.read_text(encoding='utf-8')
        except OSError:
            logger.warning('Could not read local delta '
                           f'log for {runtask.name}')
        if log:
            deltalog = base / runtask.deltalogname
            try:  # pylint: disable=too-many-try-statements
                with deltalog.open('a', encoding='utf-8') as collated_log:
                    collated_log.write(
                        f'\n\n### STARTING LOG FOR {runtask.name} ###\n\n{log}'
                        )
            except OSError:
                logger.warning(
                    f'Error writing delta log part {runtask.name}: ',
                    exc_info=True
                    )

    # Clean up
    try:
        shutil.rmtree(workfolder)
    except OSError:
        logger.warning(f'Error deleting folder {workfolder.name}')
    return ''


# TODO: see note in refcalc.compile_refcalc. #43
def compile_delta(comptask):
    """Compile a delta-amplitudes calculation executable.

    Compilation is performed in the comptask.foldername
    subfolder of the current directory. This function may
    be executed by parallelized workers.

    Parameters
    ----------
    comptask : DeltaCompileTask
        Information about the compilation to be performed.

    Returns
    -------
    error_info : str
        Description of any error that occurred while compiling.
    """
    workfolder = Path(comptask.foldername).resolve()
    try:
        workfolder.mkdir()
    except FileExistsError:
        logger.warning(f'Folder {workfolder.name} already exists. '
                       'Contents may get overwritten.')
    with execute_in_dir(workfolder):
        # write PARAM:
        try:
            Path('PARAM').write_text(comptask.param, encoding='utf-8')
        except OSError:
            logger.error('Error writing PARAM file: ', exc_info=True)
            return (f'Error encountered by {comptask} '
                    'while trying to write PARAM file.')
        try:
            comptask.copy_source_files_to_local()
        except OSError:
            logger.error('Error getting TensErLEED files for '
                         'delta-amplitudes: ', exc_info=True)
            return (f'Error encountered by {comptask} while '
                    'trying to fetch fortran source files')

        # TODO: we could skip this, if we implemented
        # a general CompileTask (Issue #43)
        (srcname, lib_tleed, lib_delta, _) = (
             str(fname.name) if fname is not None else None
             for fname in comptask.get_source_files()
             )

        # compile
        compiler = comptask.fortran_comp
        compile_list = [
            (srcname, 'main.o'),
            (lib_tleed, 'lib.tleed.o'),
            (lib_delta, 'lib.delta.o'),
            ]
        ctasks = [(f'{compiler[0]} -o {oname} -c', fname, compiler[1])
                  for fname, oname in compile_list]
        _, object_files = zip(*compile_list)
        ctasks.append(
            (f'{compiler[0]} -o {comptask.exename}',
             ' '.join(object_files),
             compiler[1])
            )
        try:
            leedbase.fortran_compile_batch(ctasks)
        except Exception as exc:
            logger.error(f'Error compiling fortran files: {exc}')
            return f'Fortran compile error in {comptask}'
    return ''


def deltas(slab, rpars, subdomain=False):
    """Runs the delta-amplitudes calculation."""
    if rpars.domainParams:
        deltas_domains(rpars)
        return

    delta_tasks = _prepare_deltas_for_one_domain(slab, rpars, subdomain)

    # if execution is suppressed, stop here
    if rpars.SUPPRESS_EXECUTION and not subdomain:
        rpars.setHaltingLevel(3)
        return

    _, run_tasks = delta_tasks
    if subdomain and run_tasks:
        rpars.manifest.add(DEFAULT_DELTAS)

    if subdomain:  # Actual calculations done in deltas_domains
        return delta_tasks

    _compile_and_run_deltas_in_parallel(rpars, *delta_tasks)
    rpars.manifest.add(DEFAULT_DELTAS)


def deltas_domains(rpars):
    """Define and run delta calculations for all domains."""
    compile_tasks = []
    run_tasks = []
    # get input for all domains
    for domain in rpars.domainParams:
        logger.info(f'Getting input for delta calculations: {domain}')
        with execute_in_dir(domain.workdir):
            try:
                result = deltas(domain.slab, domain.rpars, subdomain=True)
            except Exception:
                logger.error(f'Error while creating delta input for {domain}')
                raise
        if type(result) == tuple:
            # if no deltas need to be calculated returns None
            compile_tasks.extend(result[0])
            run_tasks.extend(result[1])
        elif result is not None:
            raise RuntimeError('Unknown error while creating '
                               f'delta input for {domain}')

    # if execution is suppressed, stop here
    if rpars.SUPPRESS_EXECUTION:
        rpars.setHaltingLevel(3)
        return
    _compile_and_run_deltas_in_parallel(rpars, compile_tasks, run_tasks)


def _assemble_tasks(slab, rpars, atom_element_pairs, deltalogname):
    """Return compilation/execution tasks for atoms that need deltas.

    Parameters
    ----------
    slab : Slab
        The slab for which delta amplitudes should be calculated.
    rpars : Rparams
        The run parameters corresponding to `slab`.
    atom_element_pairs : Sequence of (Atom, str)
        Atoms and the corresponding element names for which delta
        amplitudes should be calculated.
    deltalogname : str
        Name of the log file where runtime information from the
        execution tasks should be collected.

    Returns
    -------
    compile_tasks : list of DeltaCompileTask
        Information about which executables need to be compiled
        to calculate delta amplitudes for `atom_element_pairs`.
        May have fewer items than `atom_element_pairs`.
    run_tasks : list of DeltaRunTask
        Information about which executions of `compile_tasks` are
        needed to produce delta amplitudes for `atom_element_pairs`.
        As many items as there are `atom_element_pairs`.
    """
    # Collect input files that are common to all deltas
    static_files_contents = iodeltas.collect_static_input_files(slab, rpars)

    # Assemble tasks
    comp_tasks_by_hash = {}  # Which versions to compile
    run_tasks = []           # Which deltas to run
    tensordir = Path(f'{DEFAULT_TENSORS}_{rpars.TENSOR_INDEX:03d}')
    tl_path = rpars.get_tenserleed_directory().path
    for atom, element in atom_element_pairs:
        *din, param = iodeltas.generateDeltaInput(atom,
                                                  element,
                                                  slab,
                                                  rpars,
                                                  *static_files_contents)
        hash_ = hashlib.md5(param.encode()).digest()
        compile_task = comp_tasks_by_hash.setdefault(
            hash_,
            DeltaCompileTask(param, tl_path, len(comp_tasks_by_hash)),
            )
        run_tasks.append(
            _create_runtask(atom, element, compile_task, din, tensordir)
            )
        run_tasks[-1].deltalogname = deltalogname
    return list(comp_tasks_by_hash.values()), run_tasks


def _create_runtask(atom, element, compile_task, din, tensordir):
    """Return a new DeltaRunTask for an (atom, element) pair."""
    run_task = DeltaRunTask(compile_task)
    run_task.din, run_task.din_short = din
    run_task.tensorname = str(tensordir / f'T_{atom.num}')
    run_task.deltaname = _get_unique_delta_file_name(atom, element)
    atom.current_deltas.append(run_task.deltaname)
    return run_task


def _compile_and_run_deltas_in_parallel(rpars, compile_tasks, run_tasks):
    """Compile and run multiple delta-amplitude calculations in parallel."""
    rpars.updateCores()  # In case it was not known yet
    _compile_deltas_in_parallel(rpars, compile_tasks)
    if rpars.STOP:
        return
    ran_anything = _run_deltas_in_parallel(rpars, run_tasks)
    if rpars.STOP:
        return
    if ran_anything:
        logger.info('Delta calculations finished.')

    # Clean up compile folders
    for comp_task in compile_tasks:
        leedbase.copy_compile_log(rpars,
                                  comp_task.logfile,
                                  comp_task.compile_log_name)
        try:
            shutil.rmtree(comp_task.foldername)
        except Exception:
            logger.warning('Error deleting delta '
                           f'compile folder {comp_task.foldername}')


def _compile_deltas_in_parallel(rpars, compile_tasks):
    """Compile multiple delta-amplitudes in a parallelized manner."""
    if not compile_tasks:  # Nothing to compile
        return False

    # Make sure there's a compiler ready
    if rpars.FORTRAN_COMP[0] == "":                                             # TODO: this may mask problems of PARAMETERS settings when running a delta without a refcalc if the specified compiler does not exist.
        try:
            rpars.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("No Fortran compiler")

    for task in compile_tasks:
        task.fortran_comp = rpars.FORTRAN_COMP

    # Validate TensErLEED checksums
    if not rpars.TL_IGNORE_CHECKSUM:
        validate_multiple_files(compile_tasks[0].get_source_files(),
                                logger, "delta calculations",
                                rpars.TL_VERSION)

    # Compile files
    logger.info('Compiling fortran files...')
    poolsize = min(len(compile_tasks), rpars.N_CORES)
    try:
        parallelization.monitoredPool(rpars, poolsize,
                                      compile_delta,
                                      compile_tasks)
    except Exception:  # Save log files in case of error
        for task in compile_tasks:
            leedbase.copy_compile_log(rpars,
                                      task.logfile,
                                      task.compile_log_name)
        raise
    return True


def _ensure_tensors_loaded(slab, rpars):
    """Unpack the current Tensors and ensure `slab` is up to date with them."""
    tensors_dir = Path(DEFAULT_TENSORS)
    if not tensors_dir.is_dir():
        logger.error(f'No {tensors_dir} directory found.')
        raise FileNotFoundError(f'{tensors_dir} not found')
    iotensors.fetch_unpacked_tensor(rpars.TENSOR_INDEX)
    if 1 not in rpars.runHistory:
        load_from = tensors_dir / f'{DEFAULT_TENSORS}_{rpars.TENSOR_INDEX:03d}'
        logger.debug(
            'Running without reference calculation, checking input files '
            f'in {load_from.name} to determine original configuration.'
            )
        iotensors.getTensorOriStates(slab, load_from)
        slab.restoreOriState(keepDisp=True)


def _find_atoms_that_need_deltas(slab, rpars):
    """Return information about atoms that need delta calculations.

    Parameters
    ----------
    slab : Slab
        The slab for which delta-amplitudes are being calculated.
    rpars : Rparams
        The current PARAMETERS.

    Returns
    -------
    attodo : list of Atom
        All atoms in `slab` that have some sort of variation, and that
        thus require delta-amplitudes to be available.
    atElTodo : list of tuples
        (Atom, element) pairs of all atoms in `slab` that require a
        new delta-amplitude calculation.
    vaclist : list of Atom
        All atoms in `attodo` that have partial occupation.
    """
    # go through atoms, remove those that have no variation whatsoever:
    attodo = [at for at in slab if not at.is_bulk]
    j = 0
    while j < len(attodo):
        found = False
        at = attodo[j]
        for el in at.disp_occ.keys():
            at.mergeDisp(el)
        for d in [at.disp_occ, at.disp_geo, at.disp_vib]:
            for el in d:
                if len(d[el]) > 1:
                    found = True
                    break
        if not found:
            for el in at.disp_vib:
                if abs(at.disp_vib[el][0]) >= 1e-4:
                    found = True
                    break
        if not found:
            for el in at.disp_geo:
                if np.linalg.norm(at.disp_geo[el][0]) >= 1e-4:
                    found = True
                    break
        if not found:
            occlists = []
            for k in at.disp_occ:
                occlists.append(at.disp_occ[k])
            for i in range(0, len(occlists[0])):
                totalocc = 0.
                for ol in occlists:
                    if len(ol) <= i:
                        break  # error - will pop up again later...
                    totalocc += ol[i]
                if totalocc < 1 - 1e-4:
                    found = True
                    break
        if not found:
            attodo.pop(j)
        else:
            j += 1

    vaclist = []    # atoms for which a vacancy delta file is needed
    for at in attodo:
        occlists = []
        for k in at.disp_occ:
            occlists.append(at.disp_occ[k])
        for i in range(0, len(occlists[0])):
            totalocc = 0.
            for ol in occlists:
                if len(ol) <= i:
                    err_msg = f'Inconsistent occupancy lists for {at}'
                    logger.error(err_msg)
                    raise ValueError(err_msg)
                totalocc += ol[i]
            if totalocc < 1 - 1e-4:
                vaclist.append(at)
                break

    # check existing delta files
    countExisting = 0
    atElTodo = []
    for at in attodo:
        checkEls = list(at.disp_occ.keys())
        if at in vaclist:
            checkEls.append('vac')
        for el in checkEls:
            dfiles = (f.name for f in Path.cwd().glob(f'DEL_{at.num}_{el}*'))
            found = False
            for df in dfiles:
                if iodeltas.checkDelta(df, at, el, rpars):
                    found = True
                    at.current_deltas.append(df)
                    countExisting += 1
                    break
            if not found:
                atElTodo.append((at, el))

    if len(atElTodo) == 0:
        logger.info('All Delta files specified in DISPLACEMENTS are '
                    'already present in the Deltas.zip file. Skipping new '
                    'calculations.')
        return attodo, atElTodo, vaclist
    if countExisting > 0:
        logger.info(f'{countExisting} of {len(atElTodo) + countExisting} '
                    'required Delta-files are already present. '
                    'Generating remaining {len(atElTodo)} files...')
    return attodo, atElTodo, vaclist


def _get_unique_delta_file_name(atom, element):
    """Return a unique name for a delta file for an (atom, element) pair."""
    def _get_delta_index(file):
        try:
            return int(file.name.split('_')[-1])
        except ValueError:
            return 0

    delta_name_prefix = f'DEL_{atom.num}_{element}'
    delta_indices = (_get_delta_index(f)
                     for f in Path.cwd().glob(f'{delta_name_prefix}_*'))
    return f'{delta_name_prefix}_{max(delta_indices, default=0) + 1}'


def _prepare_deltas_for_one_domain(slab, rpars, subdomain):
    """Prepare input files and return compile/run tasks for a single slab.
    
    Parameters
    ----------
    slab : Slab
        The slab for which delta-amplitude calculations should be done.
    rpars : Rparams
        The parameters corresponding to `slab`.
    subdomain : bool
        Whether `slab` is a domain of a multi-domain calculation.
    
    Returns
    -------
    comp_tasks : list of DeltaCompileTask
        Information about which executables need to be compiled
        to calculate delta amplitudes for `slab`.
    run_tasks : list of DeltaRunTask
        Information about which executions of `compile_tasks`
        are needed to produce delta amplitudes for `slab`.
    """
    # read DISPLACEMENTS block
    if not rpars.disp_block_read:
        readDISPLACEMENTS_block(rpars,
                                slab,
                                rpars.disp_blocks[rpars.search_index])
        rpars.disp_block_read = True

    _ensure_tensors_loaded(slab, rpars)

    # If there are old deltas, pull them in here
    leedbase.getDeltas(rpars.TENSOR_INDEX, required=False)

    (atoms_with_displacements,
     atom_element_pairs_requiring_new_deltas,
     atoms_with_vacancies) = _find_atoms_that_need_deltas(slab, rpars)
    if not atom_element_pairs_requiring_new_deltas:
        # Nothing to calculate
        return [], []

    _remove_old_param_file()

    # Create compilation and running tasks, as well as a log
    # file for collating the runtime information of the latter.
    collated_logs_name = _prepare_log_file(rpars, subdomain)
    delta_tasks = _assemble_tasks(slab,
                                  rpars,
                                  atom_element_pairs_requiring_new_deltas,
                                  collated_logs_name)
    # Ensure stable sorting of the stored delta files. Do so
    # only after assembling the tasks, as this modifies the
    # current_deltas of atom_elements_todo.
    _sort_current_deltas_by_element(atoms_with_displacements,
                                    atoms_with_vacancies)
    iodeltas.write_delta_input_file(*delta_tasks)
    return delta_tasks


def _prepare_log_file(rpars, subdomain):
    """Create a log file for collating those of multiple calculations."""
    log_name = f'delta-{rpars.timestamp}.log'
    if not subdomain:
        logger.info('Generating delta files...\n'
                    'Delta log will be written to local subfolders, '
                    f'and collected in {log_name}')
    log_header = ('Logs from multiple delta calculations are collected '
                  'here. Their order may not be preserved.\n')
    try:
        Path(log_name).write_text(log_header, encoding='utf-8')
    except OSError:
        logger.warning('Error creating delta log file. This will not '
                       'affect execution, proceeding...')
    return log_name


def _remove_old_param_file():
    """Rename or delete a PARAM file in the current directory."""
    param = Path('PARAM')
    try:
        param.replace('PARAM-old')
    except FileNotFoundError:
        return  # Nothing to do
    except OSError:
        pass
    try:
        param.unlink()
    except OSError:
        logger.warning(
            'Section Delta-Amplitudes: Cannot rename/remove old PARAM file. '
            'This might cause the Delta generation to fail!'
            )


def _run_deltas_in_parallel(rpars, run_tasks):
    """Execute multiple delta-amplitude calculations in parallel."""
    if not run_tasks:  # Nothing to run
        return False
    logger.info('Running delta calculations...')
    poolsize = min(len(run_tasks), rpars.N_CORES)
    parallelization.monitoredPool(rpars, poolsize, run_delta, run_tasks)
    return True


def _sort_current_deltas_by_element(atoms, atoms_with_vacancies):
    """Re-sort the current_deltas for each `atoms` by element."""
    # See bf776a37e14ece851cc8c783a7c376a44e55145f for details on
    # why sorting is needed.
    # Each (atom, element) pair should have, by now, exactly one delta
    # file. In fact, .current_deltas is cleared (i) before each refcalc
    # (in refcalc._reinitialize_deltas) and (ii) after each search
    # segment (in run_sections, after each SUPERPOS R, via calls to
    # slab.restoreOriState(keepDisp=False). The sorting we do here
    # would not be needed if we would fix #447. However, we would
    # still probably want to raise if some elements have no delta.

    def _by_element(delta_file):
        """Extract the "element" from a `delta_file` name."""
        return delta_file.split('_')[-2].lower()

    for atom in atoms:
        current_deltas = sorted(atom.current_deltas, key=_by_element)
        deltas_by_element = {
            el: tuple(fnames)
            for el, fnames in itertools.groupby(current_deltas,
                                                key=_by_element)
            }
        duplicate_deltas = {el: fnames
                            for el, fnames in deltas_by_element.items()
                            if len(fnames) > 1}
        if duplicate_deltas:
            raise DeltasError(
                f'Found multiple delta files for {atom}, elements '
                + ', '.join(duplicate_deltas)
                )
        element_names = list(atom.disp_occ)
        if atom in atoms_with_vacancies:
            element_names.append('vac')
        sorted_by_element = []
        for element in element_names:
            try:
                sorted_by_element.append(deltas_by_element[element.lower()][0])
            except KeyError:
                raise DeltasError(f'Found no delta files for {atom}, '
                                  f'element {element}') from None
        atom.current_deltas = sorted_by_element
