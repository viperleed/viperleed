# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Rfactor
"""

import os
import logging
import shutil
import subprocess

from viperleed.tleedmlib.files.iorefcalc import readFdOut
from viperleed.tleedmlib.leedbase import (fortranCompile, getTLEEDdir,
                                          getTensors)
import viperleed.tleedmlib.files.iorfactor as io

logger = logging.getLogger("tleedm.rfactor")


def rfactor(sl, rp, index, for_error=False):
    """Runs the r-factor calculation for either the reference calculation
    (index 11) or the superpos (index 12)."""
    if int((rp.THEO_ENERGIES[1]-rp.THEO_ENERGIES[0])
           / rp.THEO_ENERGIES[2]) + 1 < 2:
        logger.info("Only one theoretical energy found: Cannot calculate "
                    "a meaningful R-Factor. Stopping...")
        return
    if for_error and rp.best_v0r is None:
        logger.error("Cannot calculate R-factor for error without a stored"
                     "V0r value. Execute a normal R-factor calculation first.")
    if index == 11:
        name = "refcalc"
    else:
        name = "superpos"
    if ((index == 11 and len(rp.refcalc_fdout) == 0)
            or (index == 12 and len(rp.superpos_specout) == 0)):
        if index == 11:
            fn = "refcalc-fd.out"
        else:
            fn = "superpos-spec.out"
        if os.path.isfile(os.path.join(".", fn)):
            logger.warning(
                "R-factor calculation was called without "
                "stored spectrum data. Reading from file " + fn
                + " in work folder...")
            path = os.path.join(".", fn)
        elif os.path.isfile(os.path.join(".", "OUT", fn)):
            logger.warning(
                "R-factor calculation was called without "
                "stored spectrum data. Reading from file " + fn
                + "in OUT folder...")
            path = os.path.join(".", "OUT", fn)
        else:
            path = ""
            if index == 11:
                # try getting from Tensors
                getTensors(rp.TENSOR_INDEX, required=False)
                dn = "Tensors_"+str(rp.TENSOR_INDEX).zfill(3)
                if os.path.isfile(os.path.join(".", "Tensors", dn, fn)):
                    logger.warning(
                        "R-factor calculation was called without "
                        "stored spectrum data. Reading from file " + fn
                        + "in " + dn + "folder...")
                    path = os.path.join(".", "Tensors", dn, fn)
            if not path:
                logger.error("Cannot execute R-factor calculation: no stored "
                             "spectrum data and no "+fn+" file was "
                             "found.")
                raise RuntimeError("No spectrum data found")
        try:
            theobeams, theospec = readFdOut(readfile=path)
            if index == 11:
                rp.refcalc_fdout = theospec
            else:
                rp.superpos_specout = theospec
            rp.theobeams[name] = theobeams
        except Exception:
            logger.error("Failed to read "+path)
            raise
        # if we haven't returned, then it was read. check data vs ivbeams:
        eq = True
        eps = 1e-3
        if len(rp.ivbeams) != len(theobeams):
            eq = False
        else:
            eq = all([rp.ivbeams[i].isEqual(theobeams[i], eps=eps)
                      for i in range(0, len(rp.ivbeams))])
        if not eq:
            logger.error(
                "The list of beams read from IVBEAMS is not "
                "equivalent to the list of beams in "+path+". R-Factor "
                "calculation cannot proceed.")
            raise ValueError("Contradiction in beam sets")
    if index == 11:
        theospec = rp.refcalc_fdout
    elif index == 12:
        theospec = rp.superpos_specout
    theobeams = rp.theobeams[name]
    # WEXPEL before PARAM, to make sure number of exp. beams is correct
    try:
        io.writeWEXPEL(sl, rp, theobeams, for_error=for_error)
    except Exception:
        logger.error("Exception during writeWEXPEL: ")
        raise
    try:
        io.writeRfactPARAM(rp, theobeams, for_error=for_error)
    except Exception:
        logger.error("Exception during writeRfactPARAM: ")
        raise
    # get fortran files and compile
    try:
        tldir = getTLEEDdir(home=rp.workdir, version=rp.TL_VERSION)
        if not tldir:
            raise RuntimeError("TensErLEED code not found.")
        libpath = os.path.join(tldir, 'lib')
        libname = [f for f in os.listdir(libpath)
                   if f.startswith('rfacsb')][0]
        shutil.copy2(os.path.join(libpath, libname), libname)
        srcpath = os.path.join(tldir, 'src')
        srcname = [f for f in os.listdir(srcpath)
                   if f.startswith('rfactor.')][0]
        shutil.copy2(os.path.join(srcpath, srcname), srcname)
    except Exception:
        logger.error("Error getting TensErLEED files for r-factor "
                     "calculation: ")
        raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. R-factor "
                       "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return
    logger.info("Compiling fortran input files...")
    rfacname = "rfactor-"+rp.timestamp
    if rp.FORTRAN_COMP[0] == "":
        rp.getFortranComp()
    try:
        fortranCompile(rp.FORTRAN_COMP[0]+" -o rfacsb.o -c",
                       libname, rp.FORTRAN_COMP[1])
        fortranCompile(rp.FORTRAN_COMP[0]+" -o main.o -c", srcname,
                       rp.FORTRAN_COMP[1])
        fortranCompile(rp.FORTRAN_COMP[0]+" -o "+rfacname, "rfacsb.o "
                       "main.o", rp.FORTRAN_COMP[1])
        logger.debug("Compiled fortran files successfully")
    except Exception:
        logger.error("Error compiling fortran files: ", exc_info=True)
        raise
    rfaclogname = rfacname+".log"
    logger.info("Starting R-factor calculation...\n"
                "R-factor log will be written to file "+rfaclogname)
    rp.manifest.append(rfaclogname)
    try:
        with open(rfaclogname, "w") as log:
            subprocess.run(os.path.join('.', rfacname),
                           input=theospec, encoding="ascii",
                           stdout=log, stderr=log)
    except Exception:
        logger.error("Error during R-factor calculation. Also check "
                     "R-factor log file.")
        raise
    logger.info("Finished R-factor calculation. Processing files...")
    # rename and move files
    try:
        os.rename('WEXPEL', 'rfactor-WEXPEL')
    except Exception:
        logger.warning("Failed to rename R-factor input file WEXPEL to "
                       "rfactor-WEXPEL")
    try:
        os.rename('PARAM', 'rfactor-PARAM')
    except Exception:
        logger.warning("Failed to rename R-factor input file PARAM to "
                       "rfactor-PARAM")
    if not os.path.isfile(os.path.join(".", "ROUT")):
        logger.error("No ROUT file was found after R-Factor calculation!")
        rp.setHaltingLevel(2)
        return

    # read output
    if for_error:
        try:
            rfaclist = io.readROUTSHORT()
        except Exception:
            logger.error("Error reading ROUTSHORT file", exc_info=rp.LOG_DEBUG)
            rp.setHaltingLevel(2)
        return rfaclist

    try:
        (rfac, r_int, r_frac), v0rshift, rfaclist = io.readROUT()
    except Exception:
        logger.error("Error reading ROUT file", exc_info=rp.LOG_DEBUG)
        rp.setHaltingLevel(2)
    else:
        logger.info("With inner potential shift of {:.2f} eV: "
                    "R = {:.4f}\n".format(v0rshift, rfac))
        rp.best_v0r = v0rshift
        dl = ["."]
        if os.path.isdir("OUT"):
            dl.append("OUT")
        for dn in dl:
            for fn in [f for f in os.listdir(dn)
                       if f.startswith("R_OUT_"+rp.timestamp)
                       and os.path.isfile(os.path.join(dn, f))]:
                try:  # delete old R_OUT files
                    os.remove(os.path.join(dn, fn))
                except Exception:
                    pass
        if rfac == 0:
            logger.error(
                "ROUT reports R-Factor as zero. This means "
                "something went wrong in the reference "
                "calculation or in the R-factor calculation.")
            rp.setHaltingLevel(2)
            fn = "R_OUT_"+name+"_"+rp.timestamp
        else:
            fn = "R_OUT_"+name+"_"+rp.timestamp+"_R={:.4f}".format(rfac)
            rp.last_R = rfac
            rp.stored_R[name] = (rfac, r_int, r_frac)
        try:
            os.rename("ROUT", fn)
        except Exception:
            logger.warning("Failed to rename R-factor output file "
                           "ROUT to "+fn)
        if len(rfaclist) != len(rp.expbeams):
            logger.warning("Failed to read R-Factors per beam from "
                           "R-factor output file ROUT.")
            rfaclist = [-1]*len(rp.expbeams)
        outname = "Rfactor_plots_{}.pdf".format(name)
        aname = "Rfactor_analysis_{}.pdf".format(name)
        try:
            io.writeRfactorPdf([(b.label, rfaclist[i]) for (i, b)
                                in enumerate(rp.expbeams)],
                               plotcolors=rp.PLOT_COLORS_RFACTOR,
                               outName=outname, analysisFile=aname,
                               v0i=rp.V0_IMAG)
        except Exception:
            logger.warning("Error plotting R-factors.", exc_info=True)
    return
