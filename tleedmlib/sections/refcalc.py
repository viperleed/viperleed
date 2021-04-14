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

from viperleed.tleedmlib.leedbase import (
    fortranCompile, getTLEEDdir, getMaxTensorIndex)
from viperleed.tleedmlib.files.parameters import modifyPARAMETERS
import viperleed.tleedmlib.files.beams as beams
import viperleed.tleedmlib.files.iorefcalc as io

logger = logging.getLogger("tleedm.refcalc")


def refcalc(sl, rp, subdomain=False):
    """Runs the reference calculation."""
    if rp.domainParams:
        refcalc_domains(rp)
        return
    sl.getCartesianCoordinates(updateOrigin=True)
    sl.updateLayerCoordinates()
    try:
        io.writePARAM(sl, rp)
    except Exception:
        logger.error("Exception during writePARAM: ")
        raise
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
    try:
        tldir = getTLEEDdir(home=rp.workdir, version=rp.TL_VERSION)
        if not tldir:
            raise RuntimeError("TensErLEED code not found.")
        libpath = os.path.join(tldir, 'lib')
        libname = [f for f in os.listdir(libpath)
                   if f.startswith('lib.tleed')][0]
        shutil.copy2(os.path.join(libpath, libname), libname)
        srcpath = os.path.join(tldir, 'src')
        srcname = [f for f in os.listdir(srcpath)
                   if f.startswith('ref-calc')][0]
        shutil.copy2(os.path.join(srcpath, srcname), srcname)
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath, globalname), globalname)
    except Exception:
        logger.error("Error getting TensErLEED files for refcalc: ")
        raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Reference "
                       "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return
    logger.info("Compiling fortran input files...")
    rcname = "refcalc-"+rp.timestamp
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("Fortran compile error")
    try:
        fortranCompile(rp.FORTRAN_COMP[0]+" -o muftin.o -c",
                       "muftin.f", rp.FORTRAN_COMP[1])
        fortranCompile(rp.FORTRAN_COMP[0]+" -o lib.tleed.o -c",
                       libname, rp.FORTRAN_COMP[1])
        fortranCompile(rp.FORTRAN_COMP[0]+" -o main.o -c", srcname,
                       rp.FORTRAN_COMP[1])
        fortranCompile(rp.FORTRAN_COMP[0]+" -o "+rcname, "muftin.o "
                       "lib.tleed.o main.o", rp.FORTRAN_COMP[1])
        logger.debug("Compiled fortran files successfully")
    except Exception:
        logger.error("Error compiling fortran files: ", exc_info=True)
        raise
    rclogname = rcname+".log"
    logger.info("Starting reference calculation...\n"
                "Refcalc log will be written to file "+rclogname)
    rp.manifest.append(rclogname)
    try:
        with open(rclogname, "w") as log:
            subprocess.run(os.path.join('.', rcname), input=fin,
                           encoding="ascii", stdout=log, stderr=log)
    except Exception:
        logger.error("Error during reference calculation. Also check "
                     "refcalc log file.")
        raise
    logger.info("Finished reference calculation. Processing files...")
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
    try:
        os.rename('PARAM', 'refcalc-PARAM')
    except Exception:
        logger.warning("Failed to rename refcalc input file PARAM to "
                       "refcalc-PARAM")
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
    #    PARALELLIZED VERSION - doesn't work like this, because the dp objects
    #       don't get written to (probably blocked by being attached to rp?).
    #       Could be made to work, but serial execution works just as well.
    # poolsize = min(len(rr), rp.N_CORES)
    # loglevel = logger.level
    # try:
    #     logger.setLevel(logging.WARNING)
    #     logging.getLogger("tleedm.files.iorefcalc").setLevel(logging.WARNING)
    #     logging.getLogger("tleedm.files.beams").setLevel(logging.WARNING)
    #     with multiprocessing.Pool(poolsize) as pool:
    #         r = pool.map(runDomainRefcalc, rr)
    # except:
    #     raise
    # finally:
    #     logger.setLevel(loglevel)
    #     logging.getLogger("tleedm.files.iorefcalc").setLevel(loglevel)
    #     logging.getLogger("tleedm.files.beams").setLevel(loglevel)
    # for v in [v for v in r if v != 0]:
    #     logger.error(v)
    #     raise RuntimeError("Error during domain reference calculations")

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
