"""Section Reference Calculation."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path
import shutil
import subprocess

import fortranformat as ff
import numpy as np

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.files import beams
from viperleed.calc.files import iorefcalc
from viperleed.calc.files.ivplot import plot_iv
from viperleed.calc.lib import fs_utils
from viperleed.calc.lib import leedbase
from viperleed.calc.lib import parallelization
from viperleed.calc.lib.checksums import validate_multiple_files
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version

logger = logging.getLogger(__name__)
_TENSOR_INPUT_FILES = (
    # These are input/output files of the calculation that are stored
    # inside Tensors_XXX.zip files for purposes of reverting to the
    # same state without re-running a full reference calculation
    'IVBEAMS',
    'PARAMETERS',
    'PHASESHIFTS',
    'POSCAR',
    'VIBROCC',
    )

# TODO: we should have a parent class for compile tasks (issue #43)
# CompileTask subclasses would need a class-level list of
# glob patterns for retrieving source files. Probably allowing
# a special "{__sourcedir__}" format specifier to be formatted
# with .format(__sourcedir__=self.source_dir) before globbing.
# Similar considerations regarding the base names for foldername
# and exename.
# TODO: when implementing the compile task class, we should also
# refactor the tenserleed source class such that we can directly use
# zipped source directories
class RefcalcCompileTask():
    """Information to compile a reference calculation executable.

    Attributes
    ----------
    exename : str
        File name of the executable that will be compiled.
    foldername : str
        Name of the folder in which `exename` can be found
        after successful compilation.
    fortran_comp : tuple
        Compiler and compilation flags.
    lmax : int
        The maximum angular momentum quantum number used
        for the reference calculation.
    param : str
        Contents of the PARAM file, defining array dimensions
        for compilation.
    source_dir : Path
        Path to the folder containing the static Fortran
        source files to be compiled.

    Notes
    -----
    It is important to create instances of this class while
    the current directory is the **main** work directory for
    the reference calculation to be compiled.
    """

    def __init__(self, param, lmax, fortran_comp, sourcedir):
        """Initialize instance.

        Parameters
        ----------
        param : str
            Contents of the PARAM file, defining array dimensions
            for compilation.
        lmax : int
            The maximum angular momentum quantum number used
            for the reference calculation. This value affects
            both the .foldername and the .exename attributes.
        fortran_comp : tuple
            Compiler and compilation flags.
        sourcedir : Path
            Path to the folder containing the static Fortran
            source files to be compiled.

        Returns
        -------
        None.
        """
        self.exename = f'refcalc-{lmax}'
        self.foldername = f'refcalc-compile_LMAX{lmax}'
        self.fortran_comp = fortran_comp
        self.lmax = lmax
        self.param = param
        self.source_dir = Path(sourcedir).resolve()  # where the fortran files are

        if os.name == 'nt':
            self.exename += '.exe'

        self._basedir = Path.cwd()  # The main work folder

    def __str__(self):
        """Return a string representation of this RefcalcCompileTask."""
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
        """Return a tuple of source files needed for running a refcalc."""
        sourcedir = Path(self.source_dir).resolve()
        libpath = sourcedir / 'lib'
        srcpath = sourcedir / 'src'
        lib_tleed = next(libpath.glob('lib.tleed*'), None)
        srcname = next(srcpath.glob('ref-calc*'), None)
        globalname = srcpath / "GLOBAL"
        _muftin = self._basedir / "muftin.f"
        muftinname =_muftin if _muftin.is_file() else None
        if any(f is None for f in (srcname, lib_tleed)):
            raise RuntimeError(f"Source files missing in {sourcedir}")          # TODO: use a more appropriate custom exception in CompileTask (e.g., MissingSourceFileError)
        return lib_tleed, srcname, globalname, muftinname

    def copy_source_files_to_local(self):
        """Copy ref-calc files to current directory."""
        for filepath in self.get_source_files():
            if filepath:
                shutil.copy2(filepath, filepath.name)


# TODO: similar to RefcalcCompileTask, also RefcalcRunTask could
# profit from some refactoring to collect portions of code similar
# to those in deltas.DeltaRunTask (to be done in #43).
class RefcalcRunTask():
    """Stores information for a worker to create a subfolder, copy input there,
    compile and run a reference calculation, and copy results back."""

    def __init__(self, fin, energy, comptask, logname,
                 collect_at="", single_threaded=False, tl_version=0.):
        self.collect_at = collect_at
        self.comptask = comptask
        self.energy = energy
        self.fin = fin
        self.foldername = f'refcalc-part_{energy:.2f}eV'
        self.logname = logname
        self.single_threaded = single_threaded
        self.tl_version = tl_version

    @property
    def name(self):
        """Return a name for this task."""
        return '(single-threaded)' if self.single_threaded else self.foldername

    def __str__(self):
        """Return a string representation of this RefcalcRunTask."""
        return f'{type(self).__name__} {self.name}'


# TODO: this is largely overlapping with deltas.compile_delta. It could
# become a method of RefcalcCompileTask (or, better, its ABC once we
# do #43.
def compile_refcalc(comptask):
    """Compile a reference calculation executable.

    Compilation is performed in the comptask.foldername
    subfolder of the current directory. This function may
    be executed by parallelized workers.

    Parameters
    ----------
    comptask : RefcalcCompileTask
        Information about the compilation to be performed.

    Returns
    -------
    error_info : str
        Description of any error that occurred while compiling.
    """
    workfolder = Path(comptask.foldername).resolve()
    # Make compilation subfolder and go there
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
            logger.error('Error getting TensErLEED files for refcalc: ',
                         exc_info=True)
            return (f'Error encountered by {comptask} while '
                    'trying to fetch fortran source files')

        # TODO: we could skip this, if we implemented
        # a general CompileTask (Issue #43)
        (libname, srcname, _, muftinname) = (
             str(fname.name) if fname is not None else None
             for fname in comptask.get_source_files()
             )

        # compile
        compiler = comptask.fortran_comp
        compile_list = [
            (srcname, 'main.o'),
            (libname, 'lib.tleed.o'),
            ]
        if muftinname:
            compile_list.append((muftinname, 'muftin.o'))
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


# TODO: there is some repeated code here and in deltas.run_delta.
# This function could become a concrete method of RefcalcRunTask
# when we solve #43
def run_refcalc(runtask):
    """Run (a part of) the reference calculation.

    The calculation is executed in the current working directory
    (if in single-threaded mode) or a temporary subfolder (if only
    a subset of energies should be calculated). The subfolder is
    created when needed, and removed upon successful execution.

    Parameters
    ----------
    runtask : RefcalcRunTask
        The worker that contains information about the (part of
        the) reference calculation to execute.

    Returns
    -------
    error_info : str
        A message with information about errors occurred during
        execution of `runtask`.
    """
    base = Path.cwd()
    workfolder = base if runtask.single_threaded else base/runtask.foldername
    if not runtask.single_threaded:
        try:
            workfolder.mkdir()
        except FileExistsError:
            logger.warning(f'Folder {workfolder.name} already '
                           'exists. Contents may get overwritten.')

    with execute_in_dir(workfolder):
        if runtask.single_threaded:
            log_file = Path(runtask.logname)
            fin = runtask.fin
        else:
            log_file = Path('refcalc.log')
            fin = edit_fin_energy_lmax(runtask)
            try:
                Path('refcalc-FIN').write_text(fin, encoding='utf-8')
            except OSError:
                pass  # local FIN is just for information...

        # get executable
        exename = runtask.comptask.exename
        try:
            shutil.copy2(base/runtask.comptask.foldername/exename, workfolder)
        except OSError:
            logger.error('Error getting refcalc executable: ', exc_info=True)
            return (f'Error encountered by {runtask}: '
                    'Failed to get refcalc executable.')
        # run execution
        with log_file.open('w', encoding='utf-8') as log:
            try:
                subprocess.run(str(workfolder/exename),
                               input=fin,
                               encoding='ascii',
                               stdout=log,
                               stderr=log,
                               check=False)
            except Exception:
                logger.error('Error while executing reference calculation '
                             f'{runtask.name}. Also check refcalc log file.',
                             exc_info=True)
                return (f'Error encountered by {runtask}: '
                        'Error during refcalc execution.')

        if runtask.single_threaded:
            return ''

        # move/copy files out
        targetpath = (Path(runtask.collect_at).resolve() if runtask.collect_at
                      else base)
        energy = f'_{runtask.energy:.2f}eV'
        tensorfiles = (f for f in Path().glob('T_*') if f.is_file())
        for tensor in tensorfiles:
            try:   # move instead of copy to not duplicate the large files
                fs_utils.move(tensor, targetpath/f'{tensor.name}{energy}')
            except OSError:
                logger.error('Failed to copy refcalc output file '
                             f'{tensor} to main folder.', exc_info=True)
                return (f'Error encountered by {runtask}: '
                        'Failed to copy Tensor file out.')
        try:
            shutil.copy2('fd.out', targetpath/f'fd{energy}.out')
        except OSError:
            logger.error('Failed to copy refcalc output file fd.out '
                         'to main folder.', exc_info=True)
            return (f'Error encountered by {runtask}: '
                    'Failed to copy fd.out file out.')
        try:
            shutil.copy2('amp.out', targetpath/f'amp{energy}.out')
        except FileNotFoundError:
            if runtask.tl_version >= Version('1.7.3'):
                logger.warning('Refcalc output file amp.out not found.')
        except OSError as exc:      # warn but continue
            logger.warning('Failed to copy refcalc output file '
                           f'amp.out to main folder: {exc}')
        # append log
        log = ''
        try:
            log = log_file.read_text(encoding='utf-8')
        except OSError:
            logger.warning(f'Could not read local refcalc log {log_file}')
        if log:
            global_log_path = base / runtask.logname
            try:  # pylint: disable=too-many-try-statements
                with global_log_path.open('a', encoding='utf-8') as global_log:
                    global_log.write(
                        f'\n\n### STARTING LOG FOR {runtask.name} ###\n\n{log}'
                        )
            except OSError:
                logger.warning('Error writing refcalc log part '
                               f'{runtask.name}: ', exc_info=True)
    # clean up
    try:
        shutil.rmtree(workfolder)
    except OSError:
        logger.warning(f'Error deleting folder {runtask.foldername}')
    return ''


def edit_fin_energy_lmax(runtask):
    """modify FIN: replace the energy range (second line)"""
    comment, _, rest = runtask.fin.split("\n", maxsplit=2)
    if runtask.tl_version < Version('1.7.0'):
        eformatter = ff.FortranRecordWriter('3F7.2')
        lj = 24
    else:
        eformatter = ff.FortranRecordWriter('3F9.2')
        lj = 30
    energy = (eformatter.write([runtask.energy, runtask.energy + 0.01, 1.0])
              .ljust(lj) + 'EI,EF,DE')
    # now replace LMAX - now makes sure LMAX can be part of directory name
    # this works because even if the directory were to be named LMAX, there is
    #   a timestap after it rather than a \n
    before_LMAX, after_LMAX = rest.split("   LMAX", maxsplit=1)
    before_LMAX, *_ = before_LMAX.rsplit('\n', maxsplit=1)
    after_LMAX = (str(runtask.comptask.lmax).rjust(3).ljust(45)
                  + "LMAX" + after_LMAX)
    # fin = finparts[0] + "\n" + nl + finparts[1]
    fin = "\n".join((comment, energy, before_LMAX, after_LMAX))
    return fin


def refcalc(sl, rp, subdomain=False, parent_dir=Path()):
    """Main function to execute the reference calculation segment."""
    if rp.domainParams:
        refcalc_domains(rp)
        return
    sl.update_cartesian_from_fractional(update_origin=True)
    sl.update_layer_coordinates()
    # delete old refcalc-fd.out if present - not earlier because can be input
    #   for r-factor calculation if no refcalc is executed
    try:
        os.remove("refcalc-fd.out")
    except FileNotFoundError:
        pass
    try:
        iorefcalc.writeAUXLATGEO(sl, rp)
    except Exception:
        logger.error("Exception during writeAUXLATGEO: ")
        raise
    try:
        iorefcalc.writeAUXNONSTRUCT(sl, rp)
    except Exception:
        logger.error("Exception during writeAUXNONSTRUCT: ")
        raise
    if rp.TL_VERSION < Version('1.7.3'):
        try:
            beams.writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
        except Exception:
            logger.error("Exception during writeAUXBEAMS: ")
            raise
    try:
        iorefcalc.writeAUXGEO(sl, rp)
    except Exception:
        logger.error("Exception during writeAUXGEO: ")
        raise
    try:
        fin = iorefcalc.collectFIN(version=rp.TL_VERSION)
    except Exception:
        logger.error("Exception while trying to collect input for "
                     "refcalc FIN: ")
        raise
    try:
        with open("refcalc-FIN", "w") as wf:
            wf.write(fin)
        logger.debug("Wrote input for refcalc as file refcalc-FIN.")
    except Exception:
        logger.error(
            "Exception while trying to write refcalc-FIN file. "
            "Execution will proceed. The exception was: ",
            exc_info=True
            )
    if rp.TL_VERSION < Version('1.7.0'):   # muftin.f deprecated in version 1.7
        try:
            iorefcalc.writeMuftin(rp)
        except Exception:
            logger.error("Exception during writeMuftin: ")
            raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Reference "
                       "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return

    energies = np.arange(rp.THEO_ENERGIES.start, rp.THEO_ENERGIES.stop+0.01,         # TODO: use better arange
                         rp.THEO_ENERGIES.step)
    tl_source = rp.get_tenserleed_directory()
    tl_path = tl_source.path
    rp.updateCores()
    single_threaded = (rp.N_CORES <= 1)
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("Fortran compile error")

    # first, figure out for which LMAX to compile:
    uses_one_lmax = (
        single_threaded
        or rp.LMAX.has_single_value
        or rp.TL_VERSION <= Version('1.6')
        )
    if uses_one_lmax:
        which_lmax = {rp.LMAX.max,}
    else:    # find appropriate LMAX per energy
        ps_en = [(i, ps[0]*leedbase.HARTREE_TO_EV)
                 for (i, ps) in enumerate(rp.phaseshifts)]
        lmax = {}  # lmax as a function of energy
        warn_small = True
        warn_large = True
        _, _max = rp.get_limits("LMAX")
        for en in energies:
            try:
                ps_ind = [pe[0] for pe in ps_en if pe[1] >= en][0]
            except IndexError:
                if ps_en[-1][1] > en and ps_en[-1][1] - en > 1.:
                    logger.warning(
                        "No approriate phaseshifts found for {:.2f} eV. Will "
                        "try highest available phaseshift energy {:.2f} eV "
                        "instead.".format(en, ps_en[-1][1]))
                ps_ind = len(rp.phaseshifts) - 1
            lmax_cands = [1]
            for site_ps in rp.phaseshifts[ps_ind][1]:
                try:
                    lmax_cands.append(max(
                        [site_ps.index(v) for v in site_ps
                         if abs(v) > rp.PHASESHIFT_EPS]) + 1)
                except (IndexError, ValueError):
                    pass
            lmax[en] = min(max((rp.LMAX.min, *lmax_cands)), rp.LMAX.max)
            if lmax[en] < 6 and warn_small:
                warn_small = False
                logger.debug("Found small LMAX value based on PHASESHIFT_EPS "
                             f"parameter (LMAX = {lmax[en]}, E = {en:.2f} eV)")
            if lmax[en] > _max:
                lmax[en] = _max
                if warn_large:
                    warn_large = False
                    logger.info("The LMAX found based on the PHASESHIFT_EPS "
                                f"parameter is greater than {_max}, which is "
                                "currently not supported. LMAX was set to "
                                f"{_max}.")
        which_lmax = set(lmax.values())

    # collect compile tasks
    comp_tasks = []
    collect_param = ""
    for lm in which_lmax:
        try:
            param = iorefcalc.writePARAM(sl, rp, lmax=lm)
        except Exception:
            logger.error("Exception during writePARAM: ",
                         exc_info=rp.is_debug_mode)
            raise
        comp_tasks.append(
            RefcalcCompileTask(param, lm, rp.FORTRAN_COMP, tl_path)
            )
        collect_param += f"### PARAM file for LMAX = {lm} ###\n\n{param}\n\n"
    try:
        with open("refcalc-PARAM", "w") as wf:
            wf.write(collect_param)
    except Exception:
        # just for information, can continue
        logger.warning("Error writing refcalc-PARAM file: ", exc_info=True)

    # set up log
    logname = "refcalc-"+rp.timestamp+".log"
    if not single_threaded:
        try:
            with open(logname, "w") as wf:
                wf.write(
                    "Logs from multiple reference calculations are "
                    "collected  here. Their order may not be preserved.\n")
        except Exception:
            logger.warning("Error creating refcalc log file. This will not "
                           "affect execution, proceeding...")
    # set up collection directory
    if not single_threaded:
        collection_dir = Path("refcalc-out").resolve()
        if collection_dir.is_dir():
            try:
                shutil.rmtree(collection_dir)
            except Exception:
                logger.warning(
                    "Failed to delete existing folder "
                    + os.path.basename(collection_dir) + ". This may cause "
                    f"old data to end up in the final {DEFAULT_TENSORS}, "
                    "check results!"
                    )
        collection_dir.mkdir(parents=True, exist_ok=True)
    # collect run tasks
    ref_tasks = []
    if not single_threaded:
        for en in sorted(energies, reverse=True):
            if len(which_lmax) == 1:
                ct = comp_tasks[0]
            else:
                ct = next(ct for ct in comp_tasks if ct.lmax == lmax[en])
            ref_tasks.append(RefcalcRunTask(fin, en, ct, logname,
                                            collect_at=collection_dir,
                                            single_threaded=False,
                                            tl_version=rp.TL_VERSION))
    else:
        ct = comp_tasks[0]
        ref_tasks.append(RefcalcRunTask(fin, -1, ct, logname,
                                        single_threaded=True,
                                        tl_version=rp.TL_VERSION))

    # Validate TensErLEED checksums
    if not rp.TL_IGNORE_CHECKSUM:
        # @issue #43: this could be a class method
        files_to_check = (file for file in comp_tasks[0].get_source_files()
                          if file and 'muftin' not in file.name.lower())
        validate_multiple_files(files_to_check,
                                logger, "reference calculation",
                                rp.TL_VERSION)

    if single_threaded:
        logger.info("Compiling fortran files...")
        try:
            r = compile_refcalc(comp_tasks[0])
        except Exception:
            # if something goes wrong copy log file to compile logs
            leedbase.copy_compile_log(rp, comp_tasks[0].logfile,
                                      comp_tasks[0].compile_log_name)
            raise
        if r:
            logger.error(r)
            raise RuntimeError("Error compiling fortran files.")
        logger.info("Starting reference calculation...\n"
                    "Refcalc log will be written to file "+logname)
        logger.info("Reference calculation running without parallelization. "
                    "Set the N_CORES parameter to speed it up.")
        r = run_refcalc(ref_tasks[0])
        if r:
            logger.error(r)
            raise RuntimeError("Error in reference calculation.")
        logger.info("Reference calculation finished. Processing files...")
    else:
        # compile files
        logger.info("Compiling fortran files...")
        poolsize = min(len(comp_tasks), rp.N_CORES)
        try:
            parallelization.monitoredPool(rp, poolsize,
                                          compile_refcalc,
                                          comp_tasks,
                                          update_from=parent_dir)
        except Exception:
            # save log files in case of error:
            for ct in comp_tasks:
                leedbase.copy_compile_log(rp, ct.logfile, ct.compile_log_name)
            raise
        if rp.STOP:
            return
        # run executions
        logger.info("Running reference calculations...")
        poolsize = min(len(ref_tasks), rp.N_CORES)
        parallelization.monitoredPool(rp, poolsize, run_refcalc, ref_tasks,
                                      update_from=parent_dir)
        if rp.STOP:
            return
        logger.info("Reference calculations finished. Processing files...")

    # clean up compile folders - AMI: move logs first to compile_logs !
    for ct in comp_tasks:
        leedbase.copy_compile_log(rp, ct.logfile, ct.compile_log_name)
        try:
            shutil.rmtree(ct.foldername)
        except Exception:
            logger.warning("Error deleting refcalc compile folder "
                           + ct.foldername)

    if not single_threaded:
        iorefcalc.combine_fdout(oripath=collection_dir)
        if 1 in rp.TENSOR_OUTPUT:
            iorefcalc.combine_tensors(oripath=collection_dir)
        try:
            shutil.rmtree(collection_dir)
        except Exception:
            logger.warning("Failed to delete empty directory "
                           + os.path.basename(collection_dir))

    try:
        rp.theobeams["refcalc"], rp.refcalc_fdout = iorefcalc.readFdOut()
    except FileNotFoundError:
        logger.error("fd.out not found after reference calculation. "
                     "Check settings and refcalc log.")
        raise
    except Exception:
        logger.error("Error reading fd.out after reference calculation. "
                     "Check settings and refcalc log.")
        raise
    if rp.theobeams["refcalc"] is None:
        logger.error("No data found in fd.out. Check if file is empty.")
        raise RuntimeError                                                      # TODO: better exception
    # clear oriState for atoms and sites, current state will be new origin
    for at in sl:
        at.oriState = None
    for site in sl.sitelist:
        site.oriState = None
    # compare beam sets:
    eq = True
    eps = 1e-3
    if len(rp.ivbeams) != len(rp.theobeams["refcalc"]):
        eq = False
        message = "Number of beams is inconsitent."
    else:
        eq = all([rp.ivbeams[i].isEqual(rp.theobeams["refcalc"][i], eps=eps)
                  for i in range(0, len(rp.ivbeams))])
        message = "Beam labels are inconsistent."
    if not eq:
        logger.error("The list of beams read from IVBEAMS is not "
                     "equivalent to the list of beams in the fd.out file "
                     "produced by the reference calculation! " + message)
        rp.setHaltingLevel(2)
    # check for beams with very low values
    if not subdomain:
        small_intensities = []
        for b in [b for b in rp.theobeams["refcalc"]
                  if max(b.intens.values()) < 1e-10]:
            small_intensities.append(b.label)
        if small_intensities:
            logger.warning(
                "Some calculated beams only contain very small intensities. "
                "This may indicate that the beams do not exist for this "
                "structure. Consider removing them from IVBEAMS: "
                + ", ".join(small_intensities))
    try:
        beams.writeOUTBEAMS(rp.theobeams["refcalc"], filename="THEOBEAMS.csv")
        theobeams_norm = copy.deepcopy(rp.theobeams["refcalc"])
        for b in theobeams_norm:
            b.normMax()
        beams.writeOUTBEAMS(theobeams_norm, filename="THEOBEAMS_norm.csv")
    except Exception:
        logger.error("Error writing THEOBEAMS after reference "
                     "calculation: ", exc_info=True)
        rp.setHaltingLevel(2)
    if len(rp.theobeams["refcalc"][0].complex_amplitude) != 0:
        try:
            beams.writeOUTBEAMS(rp.theobeams["refcalc"],
                                filename="Complex_amplitudes_real.csv",
                                which="amp_real")
            beams.writeOUTBEAMS(rp.theobeams["refcalc"],
                                filename="Complex_amplitudes_imag.csv",
                                which="amp_imag")
        except Exception:
            logger.error("Error writing complex amplitudes after reference "
                         "calculation.", exc_info=rp.is_debug_mode)
    try:
        plot_iv(theobeams_norm, "THEOBEAMS.pdf", formatting=rp.PLOT_IV)
    except Exception:
        logger.warning("Error writing THEOBEAMS.pdf after reference "
                       "calculation.")

    try:
        os.rename('fd.out', 'refcalc-fd.out')
    except Exception:
        logger.warning("Failed to rename refcalc output file fd.out to "
                       "refcalc-fd.out")
    try:
        os.rename('amp.out', 'refcalc-amp.out')
    except FileNotFoundError:
        pass
    except Exception:
        logger.warning("Failed to rename refcalc output file amp.out to "
                       "refcalc-amp.out")
    if 1 not in rp.TENSOR_OUTPUT:
        return

    # Move and zip tensor files
    rp.TENSOR_INDEX = leedbase.getMaxTensorIndex() + 1
    rp.manifest.add(DEFAULT_TENSORS)
    tensor_folder = Path(DEFAULT_TENSORS)
    tensor_folder /= f'{DEFAULT_TENSORS}_{rp.TENSOR_INDEX:03d}'
    tensor_folder.mkdir(parents=True, exist_ok=True)
    for tensor_file in Path().glob('T_*'):
        try:
            tensor_file.replace(tensor_folder / tensor_file.name)
        except OSError:
            logger.error('Error moving Tensor files: ')
            raise
    for input_file in _TENSOR_INPUT_FILES:
        try:
            shutil.copy2(input_file, tensor_folder)
        except FileNotFoundError:
            continue
        except OSError:
            logger.warning(f'Failed to add input file {input_file} to '
                           f'{DEFAULT_TENSORS} folder {tensor_folder.name}.')
    try:
        shutil.copy2('refcalc-fd.out', tensor_folder / 'refcalc-fd.out')
    except Exception:
        logger.warning('Failed to copy refcalc-fd.out '
                       f'to {DEFAULT_TENSORS} folder.')

    # remove references to Deltas from old tensors
    _reinitialize_deltas(sl)


def _reinitialize_deltas(slab):
    """Remove references to deltas from previous tensors.

    Delete old delta files in main work folder, if necessary
    (there should not be any, unless there was an error). Also,
    empty all atom.known_deltas because they would refer to
    previous tensors.

    Parameters
    ----------
    slab : Slab
    """
    # delete old delta files in main work folder, if necessary
    deltas_to_remove = (f for f in Path().glob('DEL_*') if f.is_file())
    for delta_file in deltas_to_remove:
        try:
            delta_file.unlink()
        except Exception:
            logger.warning(
                "Error deleting old Delta file in work directory. This may "
                "cause the delta file to incorrectly be labelled as belonging "
                "with the new set of tensors.")

    # empty atom.known_deltas
    for at in slab:
        at.known_deltas = []


def run_refcalc_for_one_domain(domain):
    """Run the reference calculation for a single domain.

    Parameters
    ----------
    domain : DomainParameters
        Information about the domain for which the reference
        calculation should be executed.

    Raises
    ------
    Exception
        Should the reference calculation for `domain` fail.
    """
    with execute_in_dir(domain.workdir):
        try:
            refcalc(domain.slab, domain.rpars, subdomain=True)
        except Exception:
            logger.error('Exception during reference calculation '
                         f'for {domain}: ', exc_info=True)
            raise


def refcalc_domains(rp):
    """Runs reference calculations for the domains that require them."""
    rr = [dp for dp in rp.domainParams if dp.refcalc_required]
    if not rr:
        logger.info("Found no domain which requires a reference calculation.")
        return
    # make sure there's a compiler ready, and we know the number of cores:
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise
    for dp in rp.domainParams:
        dp.rpars.FORTRAN_COMP = rp.FORTRAN_COMP
    rp.updateCores()  # if number of cores is not defined, try to find it
    rr = [dp for dp in rp.domainParams if dp.refcalc_required]
    logger.info("Running reference calculations in subfolders for domains: "
                + ", ".join([d.name for d in rr]))

    for dp in rr:
        logger.info(f'Starting reference calculation for {dp}')
        run_refcalc_for_one_domain(dp)
    logger.info("Domain reference calculations finished.")

    if len(rr) < len(rp.domainParams):
        return
    # if refcalcs were done for all domains, get averaged beams
    if 3 in rp.runHistory and rp.searchResultConfig:
        weights = [rp.searchResultConfig[0][i][0]
                   for i in range(0, len(rp.domainParams))]
        logger.info(
            "Reference calculations were done for all domains. Getting "
            "weighted average over domain beams, using weights from last "
            "search result...")
    else:
        weights = None
        logger.info(
            "Reference calculations were done for all domains, but no "
            "area weights for the different domains are available yet. "
            "Getting unweighted average over domain beams...")
    rp.theobeams["refcalc"] = beams.averageBeams([dp.rpars.theobeams["refcalc"]
                                                  for dp in rp.domainParams],
                                                 weights=weights)
    try:
        beams.writeOUTBEAMS(rp.theobeams["refcalc"], filename="THEOBEAMS.csv")
        theobeams_norm = copy.deepcopy(rp.theobeams["refcalc"])
        for b in theobeams_norm:
            b.normMax()
        beams.writeOUTBEAMS(theobeams_norm, filename="THEOBEAMS_norm.csv")
    except Exception:
        logger.error("Error writing THEOBEAMS after reference calculation.",
                     exc_info=rp.is_debug_mode)
    try:
        rp.superpos_specout = beams.writeFdOut(rp.theobeams["refcalc"],
                                               rp.beamlist,
                                               filename="refcalc-fd.out",
                                               header=rp.systemName)
    except Exception:
        logger.error("Error writing averaged refcalc-fd.out for R-factor "
                     "calculation.", exc_info=rp.is_debug_mode)
