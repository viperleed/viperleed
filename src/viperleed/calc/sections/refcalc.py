"""Section Reference Calculation."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
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
from viperleed.calc.files import parameters
from viperleed.calc.files.ivplot import plot_iv
from viperleed.calc.lib import leedbase
from viperleed.calc.lib import parallelization
from viperleed.calc.lib.checksums import validate_multiple_files
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
    """Stores information for a worker to compile a refcalc file, and keeps
    track of the folder that the compiled file is in afterwards."""

    def __init__(self, param, lmax, fortran_comp, sourcedir,
                 basedir=Path()):
        self.param = param
        self.lmax = lmax
        self.fortran_comp = fortran_comp
        self.source_dir = Path(sourcedir).resolve()  # where the fortran files are
        self.basedir = Path(basedir)  # where the calculation is based
        self.foldername = f'refcalc-compile_LMAX{lmax}'
        self.exename = f'refcalc-{lmax}'

        if os.name == 'nt':
            self.exename += '.exe'

    def __str__(self):
        """Return a string representation of this RefcalcCompileTask."""
        return f'{type(self).__name__} {self.foldername}'

    @property
    def logfile(self):
        return self.basedir / self.foldername / "fortran-compile.log"

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
        _muftin = self.basedir / "muftin.f"
        muftinname =_muftin if _muftin.is_file() else None
        if any(f is None for f in (srcname, lib_tleed)):
            raise RuntimeError(f"Source files missing in {sourcedir}")          # TODO: use a more appropriate custom exception in CompileTask (e.g., MissingSourceFileError)
        return lib_tleed, srcname, globalname, muftinname

    def copy_source_files_to_local(self):
        """Copy ref-calc files to current directory."""
        for filepath in self.get_source_files():
            if filepath:
                shutil.copy2(filepath, filepath.name)


class RefcalcRunTask():
    """Stores information for a worker to create a subfolder, copy input there,
    compile and run a reference calculation, and copy results back."""

    def __init__(self, fin, energy, comptask, logname,
                 collect_at="", single_threaded=False, tl_version=0.):
        self.fin = fin
        self.energy = energy
        self.comptask = comptask
        self.logname = logname
        self.foldername = "refcalc-part_{:.2f}eV".format(energy)
        self.collect_at = collect_at
        self.single_threaded = single_threaded
        self.tl_version = tl_version


# TODO: this is largely overlapping with deltas.compile_delta. It could
# become a method of RefcalcCompileTask (or, better, its ABC once we
# do #43.
def compile_refcalc(comptask):
    """Compile a reference calculation executable.

    This function may be executed by parallelized workers.

    Parameters
    ----------
    comptask : RefcalcCompileTask
        Information about the compilation to be performed.

    Returns
    -------
    error_info : str
        Description of any error that occurred while compiling.
    """
    workfolder = Path(comptask.basedir) / comptask.foldername
    # Make compilation subfolder and go there
    try:
        workfolder.mkdir()
    except FileExistsError:
        logger.warning(f'Folder {workfolder} already exists. '
                       'Contents may get overwritten.')
    os.chdir(workfolder)
    # write PARAM:
    try:
        with open('PARAM', 'w', encoding='utf-8') as param_file:
            param_file.write(comptask.param)
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
        (libname, 'lib.tleed.o'),
        (srcname, 'main.o'),
        ]
    if muftinname:
        compile_list.append((muftinname, 'muftin.o'))
    ctasks = [(f'{compiler[0]} -o {oname} -c', fname, compiler[1])
              for (fname, oname) in compile_list]
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
    os.chdir(comptask.basedir)
    return ''


def run_refcalc(runtask):
    """Runs a part of a reference calculation in a subfolder, or the whole
    refcalc here if in single-threaded mode."""
    base = runtask.comptask.basedir
    workfolder = base
    task_name = "(single-threaded)"
    if not runtask.single_threaded:
        # make folder and go there:
        workfolder = os.path.join(base, runtask.foldername)
        task_name = runtask.foldername
        if os.path.isdir(workfolder):
            logger.warning("Folder "+runtask.foldername+" already exists. "
                           "Contents may get overwritten.")
        else:
            os.mkdir(workfolder)
        os.chdir(workfolder)

    if runtask.single_threaded:
        fin = runtask.fin
        logname = runtask.logname
    else:
        logname = "refcalc.log"
        fin = edit_fin_energy_lmax(runtask)
        try:
            with open("refcalc-FIN", "w") as wf:
                wf.write(fin)
        except Exception:
            pass  # local FIN is just for information...

    # get executable
    exename = runtask.comptask.exename
    try:
        shutil.copy2(os.path.join(base, runtask.comptask.foldername, exename),
                     os.path.join(workfolder, exename))
    except Exception:
        logger.error("Error getting refcalc executable: ", exc_info=True)
        return ("Error encountered by RefcalcRunTask " + task_name
                + ": Failed to get refcalc executable.")
    # run execution
    try:
        with open(logname, "w") as log:
            subprocess.run(os.path.join(workfolder, exename),
                           input=fin, encoding="ascii",
                           stdout=log, stderr=log)
    except Exception:
        logger.error("Error while executing reference calculation "
                     + task_name + ". Also check refcalc log file.",
                     exc_info=True)
        return ("Error encountered by RefcalcRunTask " + task_name
                + ": Error during refcalc execution.")

    if runtask.single_threaded:
        return ""

    # move/copy files out
    if runtask.collect_at:
        targetpath = os.path.abspath(runtask.collect_at)
    else:
        targetpath = base
    en_str = "_{:.2f}eV".format(runtask.energy)
    tensorfiles = [f for f in os.listdir() if f.startswith("T_")
                   and os.path.isfile(f)]
    for tf in tensorfiles:
        try:   # move instead of copy to not duplicate the large files
            shutil.move(os.path.join(workfolder, tf),
                        os.path.join(targetpath, tf + en_str))
        except Exception:
            logger.error("Failed to copy refcalc output file " + tf +
                         " to main folder.", exc_info=True)
            return ("Error encountered by RefcalcRunTask " + task_name
                    + ": Failed to copy Tensor file out.")
    try:
        shutil.copy2(os.path.join(workfolder, "fd.out"),
                     os.path.join(targetpath, "fd" + en_str + ".out"))
    except Exception:
        logger.error("Failed to copy refcalc output file fd.out "
                     " to main folder.", exc_info=True)
        return ("Error encountered by RefcalcRunTask " + task_name
                + ": Failed to copy fd.out file out.")
    try:
        shutil.copy2(os.path.join(workfolder, "amp.out"),
                     os.path.join(targetpath, "amp" + en_str + ".out"))
    except FileNotFoundError:
        if runtask.tl_version >= Version('1.7.3'):
            logger.warning("Refcalc output file amp.out not found.")
    except Exception as e:      # warn but continue
        logger.warning("Failed to copy refcalc output file amp.out "
                       "to main folder: " + str(e))
    # append log
    log = ""
    try:
        with open(logname, "r") as rf:
            log = rf.read()
    except Exception:
        logger.warning("Could not read local refcalc log " + task_name)
    if log != "":
        globallog = os.path.join(base, runtask.logname)
        try:
            with open(globallog, "a") as wf:
                wf.write("\n\n### STARTING LOG FOR " + task_name
                         + " ###\n\n" + log)
        except Exception:
            logger.warning("Error writing refcalc log part "
                           + task_name + ": ", exc_info=True)
    # clean up
    os.chdir(base)
    try:
        shutil.rmtree(workfolder)
    except Exception:
        logger.warning("Error deleting folder " + runtask.foldername)
    return ""


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
    if single_threaded or rp.LMAX.has_single_value or rp.TL_VERSION <= Version('1.6'):
        which_lmax = {rp.LMAX.max,}
    else:    # find appropriate LMAX per energy
        ps_en = [(i, ps[0]*leedbase.HARTREE_TO_EV) for (i, ps) in enumerate(rp.phaseshifts)]
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
        # TODO: see if we can remove paths.work here. Was introduced in         # TODO
        # 7c6018c455d95440ae2f0b5b5c18ad1379600d85 due to a problem of
        # not-reverting-back-to-home in case of failures. The better
        # solution there seems to be use a try...finally block, or the
        # execute_in_dir context manager.
        comp_tasks.append(RefcalcCompileTask(param, lm, rp.FORTRAN_COMP,
                                             tl_path, basedir=rp.paths.work))
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
        collection_dir = os.path.join(os.getcwd(), "refcalc-out")
        if os.path.isdir(collection_dir):
            try:
                shutil.rmtree(collection_dir)
            except Exception:
                logger.warning(
                    "Failed to delete existing folder "
                    + os.path.basename(collection_dir) + ". This may cause "
                    f"old data to end up in the final {DEFAULT_TENSORS}, "
                    "check results!"
                    )
        os.makedirs(collection_dir, exist_ok=True)
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
        home = os.getcwd()
        try:
            r = compile_refcalc(comp_tasks[0])
        except Exception:
            # if something goes wrong copy log file to compile logs
            leedbase.copy_compile_log(rp, comp_tasks[0].logfile,
                                      comp_tasks[0].compile_log_name)
            raise
        finally:
            os.chdir(home)
        if r:
            logger.error(r)
            raise RuntimeError("Error compiling fortran files.")
        logger.info("Compiling fortran files...")
        logger.info("Starting reference calculation...\n"
                    "Refcalc log will be written to file "+logname)
        logger.info("Reference calculation running without parallelization. "
                    "Set the N_CORES parameter to speed it up.")
        try:
            r = run_refcalc(ref_tasks[0])
        finally:
            os.chdir(home)
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
            shutil.rmtree(os.path.join(ct.basedir, ct.foldername))
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
    if DEFAULT_TENSORS not in rp.manifest:
        rp.manifest.append(DEFAULT_TENSORS)
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
        input_file = Path(input_file)
        if not input_file.is_file():
            # If there was no input, there is also no output
            continue
        output_file = Path(input_file.name + '_OUT')
        should_take_out_suffixed = (
            input_file.name in {'POSCAR', 'VIBROCC', 'PARAMETERS'}
            and 3 in rp.runHistory  # Search
            and output_file.is_file()
            )
        if not should_take_out_suffixed:
                output_file = input_file
        try:
            shutil.copy2(output_file, tensor_folder / input_file.name)
        except OSError:
            logger.warning(f'Failed to add input file {input_file.name} to '
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


def runDomainRefcalc(dp):
    """Runs the reference calculation for one domain, based on the
    DomainParameters object."""
    home = os.getcwd()
    try:
        os.chdir(dp.workdir)
        refcalc(dp.sl, dp.rp, subdomain=True)
    except Exception:
        logger.error("Exception during reference calculation for domain {}: "
                     .format(dp.name), exc_info=True)
        raise
    finally:
        os.chdir(home)
    return


def refcalc_domains(rp):
    """Runs reference calculations for the domains that require them."""
    rr = [dp for dp in rp.domainParams if dp.refcalcRequired]
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
        dp.rp.FORTRAN_COMP = rp.FORTRAN_COMP
    rp.updateCores()  # if number of cores is not defined, try to find it
    rr = [dp for dp in rp.domainParams if dp.refcalcRequired]
    logger.info("Running reference calculations in subfolders for domains: "
                + ", ".join([d.name for d in rr]))

    for dp in rr:
        logger.info("Starting reference calculation for domain {}"
                    .format(dp.name))
        runDomainRefcalc(dp)
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
    rp.theobeams["refcalc"] = beams.averageBeams([
        dp.rp.theobeams["refcalc"] for dp in rp.domainParams], weights=weights)
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
