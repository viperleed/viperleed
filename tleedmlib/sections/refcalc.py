# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Reference Calculation
"""

import os
import logging
import copy
import shutil
import subprocess
import numpy as np

from viperleed import fortranformat as ff
from viperleed.tleedmlib.leedbase import (
    fortranCompile, getTLEEDdir, getMaxTensorIndex, monitoredPool)
from viperleed.tleedmlib.base import splitMaxRight
from viperleed.tleedmlib.files.parameters import modifyPARAMETERS
import viperleed.tleedmlib.files.beams as beams
import viperleed.tleedmlib.files.iorefcalc as io

logger = logging.getLogger("tleedm.refcalc")


class RefcalcCompileTask():
    """Stores information for a worker to compile a refcalc file, and keeps
    track of the folder that the compiled file is in afterwards."""

    def __init__(self, param, lmax, fortran_comp, sourcedir,
                 basedir=os.getcwd()):
        self.param = param
        self.lmax = lmax
        self.fortran_comp = fortran_comp
        self.sourcedir = sourcedir  # where the fortran files are
        self.basedir = basedir    # where the calculation is based
        self.foldername = "refcalc-compile_LMAX{}".format(lmax)
        self.exename = "refcalc-{}".format(lmax)


class RefcalcRunTask():
    """Stores information for a worker to create a subfolder, copy input there,
    compile and run a reference calculation, and copy results back."""

    def __init__(self, fin, energy, comptask, logname,
                 collect_at="", single_threaded=False):
        self.fin = fin
        self.energy = energy
        self.comptask = comptask
        self.logname = logname
        self.foldername = "refcalc-part_{:.2f}eV".format(energy)
        self.collect_at = collect_at
        self.single_threaded = single_threaded


def compile_refcalc(comptask):
    """Function meant to be executed by parallelized workers. Executes a
    RefcalcCompileTask."""
    logger = logging.getLogger("tleedm.refcalc")
    home = os.getcwd()
    workfolder = os.path.join(comptask.basedir, comptask.foldername)
    # make folder and go there:
    if os.path.isdir(workfolder):
        logger.warning("Folder "+comptask.foldername+" already exists. "
                       "Contents may get overwritten.")
    else:
        os.mkdir(workfolder)
    os.chdir(workfolder)
    # write PARAM:
    try:
        with open("PARAM", "w") as wf:
            wf.write(comptask.param)
    except Exception:
        logger.error("Error writing PARAM file: ", exc_info=True)
        return ("Error encountered by RefcalcCompileTask "
                + comptask.foldername + "while trying to write PARAM file.")
    # get Fortran source files
    try:
        libpath = os.path.join(comptask.sourcedir, 'lib')
        libname = [f for f in os.listdir(libpath)
                   if f.startswith('lib.tleed')][0]
        shutil.copy2(os.path.join(libpath, libname), libname)
        srcpath = os.path.join(comptask.sourcedir, 'src')
        srcname = [f for f in os.listdir(srcpath)
                   if f.startswith('ref-calc')][0]
        shutil.copy2(os.path.join(srcpath, srcname), srcname)
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath, globalname), globalname)
        muftinname = "muftin.f"
        shutil.copy2(os.path.join(comptask.basedir, muftinname), muftinname)
    except Exception:
        logger.error("Error getting TensErLEED files for "
                     "refcalc-amplitudes: ", exc_info=True)
        return ("Error encountered by RefcalcCompileTask "
                + comptask.foldername + " while trying to fetch fortran "
                "source files")
    # compile
    try:
        for (fname, oname) in [(muftinname, "muftin.o"),
                               (libname, "lib.tleed.o"),
                               (srcname, "main.o")]:
            fortranCompile(
                comptask.fortran_comp[0] + " -o " + oname + " -c",
                fname, comptask.fortran_comp[1])
        fortranCompile(
            comptask.fortran_comp[0] + " -o " + comptask.exename,
            "muftin.o lib.tleed.o main.o", comptask.fortran_comp[1])
    except Exception:
        logger.error("Error compiling fortran files: ")
        return ("Fortran compile error in RefcalcCompileTask "
                + comptask.foldername)
    os.chdir(home)
    return ""


def run_refcalc(runtask):
    """Runs a part of a reference calculation in a subfolder, or the whole
    refcalc here if in single-threaded mode."""
    logger = logging.getLogger("tleedm.refcalc")
    home = os.getcwd()
    base = runtask.comptask.basedir
    workfolder = home
    task_name = "(single-threaded)"
    if not runtask.single_threaded:
        # make folder and go there:
        workfolder = os.path.join(home, runtask.foldername)
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
        # modify FIN: replace the energy range (second line)
        finparts = runtask.fin.split("\n", maxsplit=2)
        f72x3 = ff.FortranRecordWriter('3F7.2')
        nl = (f72x3.write([runtask.energy, runtask.energy + 0.01, 1.0])
              .ljust(24) + 'EI,EF,DE')
        fin = "\n".join((finparts[0], nl, finparts[2]))
        # now replace LMAX
        finparts = fin.split("LMAX", maxsplit=1)
        finparts = [splitMaxRight(finparts[0], "\n")[0], finparts[1]]
        nl = str(runtask.comptask.lmax).rjust(3).ljust(45) + "LMAX"
        fin = finparts[0] + "\n" + nl + finparts[1]
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
    os.chdir(home)
    try:
        shutil.rmtree(workfolder)
    except Exception:
        logger.warning("Error deleting folder " + runtask.foldername)
    return ""


def refcalc(sl, rp, subdomain=False):
    """Main function to execute the reference calculation segment."""
    if rp.domainParams:
        refcalc_domains(rp)
        return
    sl.getCartesianCoordinates(updateOrigin=True)
    sl.updateLayerCoordinates()
    try:
        io.writeAUXLATGEO(sl, rp)
    except Exception:
        logger.error("Exception during writeAUXLATGEO: ")
        raise
    try:
        io.writeAUXNONSTRUCT(sl, rp)
    except Exception:
        logger.error("Exception during writeAUXNONSTRUCT: ")
        raise
    try:
        beams.writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
    except Exception:
        logger.error("Exception during writeAUXBEAMS: ")
        raise
    try:
        io.writeAUXGEO(sl, rp)
    except Exception:
        logger.error("Exception during writeAUXGEO: ")
        raise
    try:
        fin = io.collectFIN()
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
            "Execution will proceed. The exception was: ", exc_info=True)
    try:
        io.writeMuftin(sl, rp)
    except Exception:
        logger.error("Exception during writeMuftin: ")
        raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Reference "
                       "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return

    energies = np.arange(rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+0.01,
                         rp.THEO_ENERGIES[2])
    tldir = os.path.abspath(getTLEEDdir(home=rp.workdir,
                                        version=rp.TL_VERSION))
    if not tldir:
        raise RuntimeError("TensErLEED code not found.")
    rp.updateCores()
    single_threaded = (rp.N_CORES <= 1)
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("Fortran compile error")

    # first, figure out for which LMAX to compile:
    if single_threaded or rp.LMAX[0] == rp.LMAX[1] or rp.TL_VERSION <= 1.6:
        which_lmax = set([rp.LMAX[1]])
    else:    # find appropriate LMAX per energy
        if rp.PHASESHIFT_EPS == 0:
            rp.PHASESHIFT_EPS = 0.01
        ps_en = [(i, ps[0]*27.2116) for (i, ps) in enumerate(rp.phaseshifts)]
        lmax = {}  # lmax as a function of energy
        warn_small = True
        warn_large = True
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
            lmax[en] = min(max(max(lmax_cands), rp.LMAX[0]), rp.LMAX[1])
            if lmax[en] < 8 and warn_small:
                warn_small = False
                logger.debug(
                    "Found small LMAX value based on PHASESHIFT_EPS parameter "
                    "(LMAX = {}, E = {:.2f} eV)".format(lmax[en], en))
            if lmax[en] > 18:
                lmax[en] = 18
                if warn_large:
                    warn_large = False
                    logger.info(
                        "The LMAX found based on the PHASESHIFT_EPS "
                        "parameter is greater than 18, which is currently "
                        "not supported. LMAX was set to 18.")
        which_lmax = set(lmax.values())

    # collect compile tasks
    comp_tasks = []
    collect_param = ""
    for lm in which_lmax:
        try:
            param = io.writePARAM(sl, rp, lmax=lm)
        except Exception:
            logger.error("Exception during writePARAM: ",
                         exc_info=rp.LOG_DEBUG)
            raise
        comp_tasks.append(RefcalcCompileTask(param, lm, rp.FORTRAN_COMP,
                                             tldir, basedir=os.getcwd()))
        collect_param += ("### PARAM file for LMAX = {} ###\n\n".format(lm)
                          + param + "\n\n")
    try:
        with open("refcalc-PARAM", "w") as wf:
            wf.write(collect_param)
    except Exception:
        # just for information, can continue
        logger.warning("Error writing refcalc-PARAM file: ", exc_info=True)

    # set up log
    logname = "refcalc-"+rp.timestamp+".log"
    rp.manifest.append(logname)
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
                    "old data to end up in the final Tensors, check results!")
        os.makedirs(collection_dir, exist_ok=True)
    # collect run tasks
    ref_tasks = []
    if not single_threaded:
        for en in energies:
            if len(which_lmax) == 1:
                ct = comp_tasks[0]
            else:
                ct = [ct for ct in comp_tasks if ct.lmax == lmax[en]][0]
            ref_tasks.append(RefcalcRunTask(fin, en, ct, logname,
                                            collect_at=collection_dir,
                                            single_threaded=False))
    else:
        ct = comp_tasks[0]
        ref_tasks.append(RefcalcRunTask(fin, -1, ct, logname,
                                        single_threaded=True))

    if single_threaded:
        home = os.getcwd()
        try:
            r = compile_refcalc(comp_tasks[0])
        except Exception:
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
        except Exception:
            raise
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
        monitoredPool(rp, poolsize, compile_refcalc, comp_tasks)
        if rp.STOP:
            return
        # run executions
        logger.info("Running reference calculations...")
        poolsize = min(len(ref_tasks), rp.N_CORES)
        monitoredPool(rp, poolsize, run_refcalc, ref_tasks)
        if rp.STOP:
            return
        logger.info("Reference calculations finished. Processing files...")

    # clean up compile folders
    for ct in comp_tasks:
        try:
            shutil.rmtree(os.path.join(ct.basedir, ct.foldername))
        except Exception:
            logger.warning("Error deleting refcalc compile folder "
                           + ct.foldername)

    if not single_threaded:
        io.combine_tensors(oripath=collection_dir)
        io.combine_fdout(oripath=collection_dir)
        try:
            shutil.rmtree(collection_dir)
        except Exception:
            logger.warning("Failed to delete empty directory "
                           + os.path.basename(collection_dir))

    try:
        rp.theobeams["refcalc"], rp.refcalc_fdout = io.readFdOut()
    except FileNotFoundError:
        logger.error("fd.out not found after reference calculation. "
                     "Check settings and refcalc log.")
        raise
    except Exception:
        logger.error("Error reading fd.out after reference calculation. "
                     "Check settings and refcalc log.")
        raise
    # clear oriState for atoms and sites, current state will be new origin
    for at in sl.atlist:
        at.oriState = None
    for site in sl.sitelist:
        site.oriState = None
    # compare beam sets:
    eq = True
    eps = 1e-3
    if len(rp.ivbeams) != len(rp.theobeams["refcalc"]):
        eq = False
    else:
        eq = all([rp.ivbeams[i].isEqual(rp.theobeams["refcalc"][i], eps=eps)
                  for i in range(0, len(rp.ivbeams))])
    if not eq:
        logger.error("The list of beams read from IVBEAMS is not "
                     "equivalent to the list of beams in the fd.out file "
                     "produced by the reference calculation!")
        rp.setHaltingLevel(2)
    # check for beams with very low values
    if not subdomain:
        for b in [b for b in rp.theobeams["refcalc"]
                  if max(b.intens.values()) < 1e-10]:
            logger.warning(
                "Beam {} only contains very small intensities. This may "
                "indicate that the beam does not exist for this structure. "
                "Consider removing it from IVBEAMS.".format(b.label))
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
    # rename and move files
    try:
        os.rename('fd.out', 'refcalc-fd.out')
    except Exception:
        logger.warning("Failed to rename refcalc output file fd.out to "
                       "refcalc-fd.out")
    # move and zip tensor files
    rp.TENSOR_INDEX = getMaxTensorIndex() + 1
    rp.manifest.append("Tensors")
    dn = "Tensors_"+str(rp.TENSOR_INDEX).zfill(3)
    os.makedirs(os.path.join(".", "Tensors", dn), exist_ok=True)
    try:
        for tf in [f for f in os.listdir('.') if f.startswith("T_")]:
            shutil.move(tf, os.path.join(".", "Tensors", dn, tf))
    except Exception:
        logger.error("Error moving Tensor files: ")
        raise
    tInputFiles = ["POSCAR", "PARAMETERS", "VIBROCC", "IVBEAMS",
                   "PHASESHIFTS"]
    for f in [f for f in tInputFiles if f in os.listdir('.')]:
        of = f
        for fn in ["POSCAR", "VIBROCC"]:
            if (f == fn and 3 in rp.runHistory
                    and os.path.isfile(fn + "_OUT_" + rp.timestamp)):
                of = fn + "_OUT_" + rp.timestamp
        try:
            shutil.copy2(of, os.path.join("Tensors", dn, f))
        except Exception:
            logger.warning("Failed to add input file " + f
                           + " to Tensors folder.")
    try:
        shutil.copy2("refcalc-fd.out",
                     os.path.join("Tensors", dn, "refcalc-fd.out"))
    except Exception:
        logger.warning("Failed to copy refcalc-fd.out to Tensors folder.")
    # modify PARAMETERS to contain the energies and LMAX that were really used
    if os.path.isfile(os.path.join("Tensors", dn, "PARAMETERS")):
        modifyPARAMETERS(rp, "THEO_ENERGIES",
                         new=" ".join(["{:.4g}".format(v)
                                       for v in rp.THEO_ENERGIES]),
                         path=os.path.join("Tensors", dn),
                         suppress_ori=True)
        modifyPARAMETERS(rp, "LMAX", new=str(rp.LMAX),
                         path=os.path.join("Tensors", dn), suppress_ori=True)
    # delete old delta files in main work folder, if necessary
    #   (there should not be any, unless there was an error)
    for df in [f for f in os.listdir(".") if f.startswith("DEL_") and
               os.path.isfile(os.path.join(".", f))]:
        try:
            os.remove(df)
        except Exception:
            logger.warning(
                "Error deleting old Delta file in work directory. This may "
                "cause the delta file to incorrectly be labelled as belonging "
                "with the new set of tensors.")
    return


def runDomainRefcalc(dp):
    """Runs the reference calculation for one domain, based on the
    DomainParameters object."""
    logger = logging.getLogger("tleedm.refcalc")
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
        weights = []
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
                     exc_info=rp.LOG_DEBUG)
    try:
        rp.superpos_specout = beams.writeFdOut(rp.theobeams["refcalc"],
                                               rp.beamlist,
                                               filename="refcalc-fd.out",
                                               header=rp.systemName)
    except Exception:
        logger.error("Error writing averaged refcalc-fd.out for R-factor "
                     "calculation.", exc_info=rp.LOG_DEBUG)
    return
