"""Section Delta Amplitudes."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import hashlib
import logging
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files import iodeltas
from viperleed.calc.files import iotensors
from viperleed.calc.files.beams import writeAUXBEAMS
from viperleed.calc.files.displacements import readDISPLACEMENTS_block
from viperleed.calc.lib import leedbase
from viperleed.calc.lib import parallelization
from viperleed.calc.lib.checksums import validate_multiple_files
from viperleed.calc.lib.context import execute_in_dir

logger = logging.getLogger(__name__)

# TODO: would be nice to replace all os.path with pathlib


# TODO: see note in refcalc. #43
class DeltaCompileTask():
    """Stores information for a worker to compile a delta file, and keeps
    track of the folder that the compiled file is in afterwards."""

    def __init__(self, param, hash_, source_dir, index):
        self.exename = f'delta-{index}'
        self.foldername = f'Delta_Compile_{index}'
        self.fortran_comp = ["", ""]
        self.hash = hash_
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
        globalname = srcpath / "GLOBAL"
        if any(f is None for f in (srcname, lib_tleed, lib_delta)):
            raise RuntimeError("Source files missing in {self.source_dir}")     # TODO: use a better custom exception in CompileTask (e.g., MissingSourceFileError)
        return srcname, lib_tleed, lib_delta, globalname

    def copy_source_files_to_local(self):
        """Copy delta source files to current directory."""
        for filepath in self.get_source_files():
            if filepath:
                shutil.copy2(filepath, filepath.name)


# TODO: see note in refcalc.run_refcalc. #43
class DeltaRunTask():
    """Stores information needed to copy the correct delta executable and
    tensor file to a subfolder, execute the delta-calculation there, and copy
    results back."""

    def __init__(self, comptask):
        self.comptask = comptask
        self.deltalogname = ""
        self.deltaname = ""
        self.din = ""
        self.din_short = ""
        self.tensorname = ""

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


def deltas(sl, rp, subdomain=False):
    """Runs the delta-amplitudes calculation."""

    if rp.domainParams:
        deltas_domains(rp)
        return

    # read DISPLACEMENTS block
    if not rp.disp_block_read:
        readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
        rp.disp_block_read = True

    # get Tensors
    if not Path(DEFAULT_TENSORS).is_dir():
        logger.error(f'No {DEFAULT_TENSORS} directory found.')
        raise RuntimeError(f'{DEFAULT_TENSORS} not found')                      # TODO: FileNotFoundError?
    iotensors.fetch_unpacked_tensor(rp.TENSOR_INDEX)
    if 1 not in rp.runHistory:
        load_from = Path(DEFAULT_TENSORS)
        load_from /= f'{DEFAULT_TENSORS}_{rp.TENSOR_INDEX:03d}'
        logger.debug(
            'Running without reference calculation, checking input files '
            f'in {load_from.name} to determine original configuration.'
            )
        iotensors.getTensorOriStates(sl, load_from)
        sl.restoreOriState(keepDisp=True)

    # if there are old deltas, fetch them
    leedbase.getDeltas(rp.TENSOR_INDEX, required=False)
    dbasic = iodeltas.generateDeltaBasic(sl, rp)

    # get AUXBEAMS; if AUXBEAMS is not in work folder, check SUPP folder
    auxbeams_file = Path('AUXBEAMS')
    if not auxbeams_file.is_file() and (DEFAULT_SUPP/auxbeams_file).is_file():
        try:
            shutil.copy2(DEFAULT_SUPP/auxbeams_file, auxbeams_file.name)
        except OSError:
            logger.warning(f'Failed to copy {auxbeams_file.name} from '
                           f'{DEFAULT_SUPP} folder. Generating new file...')
    if not auxbeams_file.is_file():
        try:
            writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
        except Exception:                                                       # TODO: better exception
            logger.error('Exception during writeAUXBEAMS: ')
            raise
    try:
        auxbeams = auxbeams_file.read_text(encoding='utf-8')
    except OSError:
        logger.error('Could not read {auxbeams_file.name} for delta input')
        raise
    if not auxbeams.endswith('\n'):
        auxbeams += '\n'

    # get PHASESHIFTS
    try:
        phaseshifts = Path('PHASESHIFTS').read_text(encoding='utf-8')
    except OSError:
        logger.error('Could not read PHASESHIFTS for delta-input')
        raise
    if not phaseshifts.endswith('\n'):
        phaseshifts += '\n'

    # go through atoms, remove those that have no variation whatsoever:
    attodo = [at for at in sl if not at.is_bulk]
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
                    else:
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
                    logger.error("Inconsistent occupancy lists for {} "
                                 .format(at))
                    raise ValueError("Inconsistent occupancy lists for {}"
                                     .format(at))
                else:
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
            checkEls.append("vac")
        for el in checkEls:
            dfiles = [f for f in os.listdir(".")
                      if f.startswith(f'DEL_{at.num}_{el}')]
            found = False
            for df in dfiles:
                if iodeltas.checkDelta(df, at, el, rp):
                    found = True
                    at.known_deltas.append(df)
                    countExisting += 1
                    break
            if not found:
                atElTodo.append((at, el))

    if len(atElTodo) == 0:
        logger.info("All Delta files specified in DISPLACEMENTS are "
                    "already present in the Deltas.zip file. Skipping new "
                    "calculations.")
        return
    if countExisting > 0:
        logger.info("{} of {} required Delta-files are already present. "
                    "Generating remaining {} files..."
                    .format(countExisting, len(atElTodo) + countExisting,
                            len(atElTodo)))
    # create log file:
    deltaname = "delta-"+rp.timestamp
    deltalogname = deltaname+".log"
    if not subdomain:
        logger.info(
            "Generating delta files...\n"
            "Delta log will be written to local subfolders, and collected in "
            + deltalogname)
    try:
        with open(deltalogname, "w") as wf:
            wf.write("Logs from multiple delta calculations are collected "
                     "here. Their order may not be preserved.\n")
    except Exception:
        logger.warning("Error creating delta log file. This will not "
                       "affect execution, proceeding...")

    # move PARAM file
    if os.path.isfile("PARAM"):
        try:
            os.rename("PARAM", "PARAM-old")
        except Exception:
            try:
                os.remove(os.path.join(".", "PARAM"))
            except Exception:
                logger.warning(
                    "Section Delta-Amplitudes: Cannot rename/remove old PARAM "
                    "file. This might cause the Delta generation to fail!")
    # assemble tasks
    deltaCompTasks = []  # keep track of what versions to compile
    deltaRunTasks = []   # which deltas to run
    tensordir = f"{DEFAULT_TENSORS}_{rp.TENSOR_INDEX:03d}"
    tl_source = rp.get_tenserleed_directory()
    tl_path = tl_source.path
    for (at, el) in atElTodo:
        din, din_short, param = iodeltas.generateDeltaInput(
            at, el, sl, rp, dbasic, auxbeams, phaseshifts)
        h = hashlib.md5(param.encode()).digest()
        found = False
        for ct in deltaCompTasks:
            if ct.hash == h:
                found = True
                rt = DeltaRunTask(ct)
                break
        if not found:
            index = len(deltaCompTasks)
            ct = DeltaCompileTask(param, h, tl_path, index)
            deltaCompTasks.append(ct)
            rt = DeltaRunTask(ct)
        deltaRunTasks.append(rt)
        rt.din = din
        rt.din_short = din_short
        rt.tensorname = os.path.join(tensordir, f'T_{at.num}')
        nameBase = f'DEL_{at.num}_{el}'
        n = 1
        nums = []
        for fn in [f for f in os.listdir(".") if f.startswith(nameBase)]:
            try:
                nums.append(int(fn.split("_")[-1]))
            except Exception:
                pass
        if nums:
            n = max(nums) + 1
        rt.deltaname = nameBase + "_{}".format(n)
        rt.deltalogname = deltalogname
        at.known_deltas.append(rt.deltaname)

    # sort known_deltas
    for at in attodo:
        checkEls = list(at.disp_occ.keys())
        if at in vaclist:
            checkEls.append("vac")
        copydel = at.known_deltas[:]
        at.known_deltas = []
        for el in checkEls:
            at.known_deltas.append(
                [df for df in copydel
                 if df.split("_")[-2].lower() == el.lower()][0])
        if len(at.known_deltas) != len(copydel):
            logger.error("Failed to sort delta files for {}".format(at))
            raise RuntimeError("Inconsistent delta files")

    # write delta-input file
    dinput = ("""# ABOUT THIS FILE:
# Input for the delta-calculations is collected here. The blocks of data are
# new 'PARAM' files, which are used to recompile the fortran code, and input
# for generation of specific DELTA files. Lines starting with '#' are comments
# on the function of the next block of data.
# In the DELTA file blocks, [AUXBEAMS] and [PHASESHIFTS] denote where the
# entire contents of the AUXBEAMS and PHASESHIFTS files should be inserted.
""")
    for ct in deltaCompTasks:
        dinput += ("\n#### NEW 'PARAM' FILE: ####\n\n" + ct.param + "\n")
        for rt in [t for t in deltaRunTasks if t.comptask == ct]:
            dinput += ("\n#### INPUT for new DELTA file {}: ####\n\n"
                       .format(rt.deltaname) + rt.din_short + "\n")
    try:
        with open("delta-input", "w") as wf:
            wf.write(dinput)
    except Exception:
        logger.warning("Failed to write file 'delta-input'. This will "
                       "not affect TensErLEED execution, proceeding...")

    # if execution is suppressed, stop here
    if rp.SUPPRESS_EXECUTION and not subdomain:
        rp.setHaltingLevel(3)
        return

    # make sure there's a compiler ready:
    if rp.FORTRAN_COMP[0] == "" and not subdomain:                              # TODO: this may mask problems of PARAMETERS settings when running a delta without a refcalc if the specified compiler does not exist.
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("No Fortran compiler")

    for ct in deltaCompTasks:
        ct.fortran_comp = rp.FORTRAN_COMP

    if subdomain and deltaRunTasks:
        rp.manifest.add(DEFAULT_DELTAS)

    if subdomain:  # Actual calculations done in deltas_domains
        return deltaCompTasks, deltaRunTasks

    rp.updateCores()

    # Validate TensErLEED checksums
    if not rp.TL_IGNORE_CHECKSUM:
        validate_multiple_files(deltaCompTasks[0].get_source_files(),
                                logger, "delta calculations",
                                rp.TL_VERSION)

    # compile files
    logger.info("Compiling fortran files...")
    poolsize = min(len(deltaCompTasks), rp.N_CORES)
    try:
        parallelization.monitoredPool(rp, poolsize,
                                      compile_delta,
                                      deltaCompTasks)
    except Exception:
        # save log files in case of error:
        for ct in deltaCompTasks:
            leedbase.copy_compile_log(rp, ct.logfile, ct.compile_log_name)
        raise
    if rp.STOP:
        return

    # run executions
    logger.info("Running delta calculations...")
    poolsize = min(len(deltaRunTasks), rp.N_CORES)
    parallelization.monitoredPool(rp, poolsize, run_delta, deltaRunTasks)
    if rp.STOP:
        return
    logger.info("Delta calculations finished.")

    # clean up compile folders
    for ct in deltaCompTasks:
        leedbase.copy_compile_log(rp, ct.logfile, ct.compile_log_name)
        try:
            shutil.rmtree(ct.foldername)
        except Exception:
            logger.warning("Error deleting delta compile folder "
                           + ct.foldername)
    rp.manifest.add(DEFAULT_DELTAS)


def deltas_domains(rp):
    """Define and run delta calculations for all domains."""
    deltaCompTasks = []
    deltaRunTasks = []
    # get input for all domains
    for dp in rp.domainParams:
        logger.info(f'Getting input for delta calculations: {dp}')
        with execute_in_dir(dp.workdir):
            try:
                r = deltas(dp.slab, dp.rpars, subdomain=True)
            except Exception:
                logger.error(f'Error while creating delta input for {dp}')
                raise
        if type(r) == tuple:  # if no deltas need to be calculated returns None
            deltaCompTasks.extend(r[0])
            deltaRunTasks.extend(r[1])
        elif r is not None:
            raise RuntimeError('Unknown error while creating '
                               f'delta input for {dp}')

    # if execution is suppressed, stop here
    if rp.SUPPRESS_EXECUTION:
        rp.setHaltingLevel(3)
        return

    # make sure there's a compiler ready, and we know the number of cores:
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("No Fortran compiler")
    for ct in deltaCompTasks:
        ct.fortran_comp = rp.FORTRAN_COMP
    rp.updateCores()  # if number of cores is not defined, try to find it

    # compile files
    if len(deltaCompTasks) > 0:
        logger.info("Compiling fortran files...")
        poolsize = min(len(deltaCompTasks), rp.N_CORES)
        parallelization.monitoredPool(rp, poolsize,
                                      compile_delta,
                                      deltaCompTasks)
        if rp.STOP:
            return

    # run executions
    if len(deltaRunTasks) > 0:
        logger.info("Running delta calculations...")
        poolsize = min(len(deltaRunTasks), rp.N_CORES)
        parallelization.monitoredPool(rp, poolsize, run_delta, deltaRunTasks)
        if rp.STOP:
            return
        logger.info("Delta calculations finished.")

    # clean up
    for ct in deltaCompTasks:
        leedbase.copy_compile_log(rp, ct.logfile, ct.compile_log_name) # copy compile folder
        try:
            shutil.rmtree(ct.foldername)
        except Exception:
            logger.warning('Error deleting delta '
                           f'compile folder {ct.foldername}')
