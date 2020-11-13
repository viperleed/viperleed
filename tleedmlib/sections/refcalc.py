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

from tleedmlib.leedbase import fortranCompile, getTLEEDdir, getMaxTensorIndex
import tleedmlib.files.beams as beams
import tleedmlib.files.iorefcalc as io

logger = logging.getLogger("tleedm.refcalc")

def refcalc(sl, rp):
    """Runs the reference calculation. Returns 0 when finishing without 
    errors, or an error message otherwise."""
    sl.getCartesianCoordinates(updateOrigin=True)
    sl.updateLayerCoordinates()
    try:
        io.writePARAM(sl, rp)
    except:
        logger.error("Exception during writePARAM: ")
        raise
    try:
        io.writeAUXLATGEO(sl, rp)
    except:
        logger.error("Exception during writeAUXLATGEO: ")
        raise
    try:
        io.writeAUXNONSTRUCT(sl, rp)
    except:
        logger.error("Exception during writeAUXNONSTRUCT: ")
        raise
    try:
        beams.writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
    except:
        logger.error("Exception during writeAUXBEAMS: ")
        raise
    try:
        io.writeAUXGEO(sl, rp)
    except:
        logger.error("Exception during writeAUXGEO: ")
        raise
    try:
        fin = io.collectFIN()
    except:
        logger.error("Exception while trying to collect input for "
                      "refcalc FIN: ")
        raise
    try:
        with open("refcalc-FIN", "w") as wf:
            wf.write(fin)
        logger.debug("Wrote input for refcalc as file refcalc-FIN.")
    except:
        logger.error("Exception while trying to write refcalc-FIN file. "
            "Execution will proceed. The exception was: ", exc_info=True)
    try:
        io.writeMuftin(sl, rp)
    except:
        logger.error("Exception during writeMuftin: ")
        raise
    try:
        tldir = getTLEEDdir(home=rp.workdir)
        libpath = os.path.join(tldir,'lib')
        libname = [f for f in os.listdir(libpath) 
                      if f.startswith('lib.tleed')][0]
        shutil.copy2(os.path.join(libpath,libname), libname)
        srcpath = os.path.join(tldir,'src')
        srcname = [f for f in os.listdir(srcpath) 
                      if f.startswith('ref-calc')][0]
        shutil.copy2(os.path.join(srcpath,srcname), srcname)
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath,globalname), globalname)
    except:
        logger.error("Error getting TensErLEED files for refcalc: ")
        raise
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Reference "
            "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return 0
    logger.info("Compiling fortran input files...")
    rcname = "refcalc-"+rp.timestamp
    if rp.FORTRAN_COMP[0] == "":
        if rp.getFortranComp() != 0:    #returns 0 on success
            logger.error("No fortran compiler found, cancelling...")
            return ("Fortran compile error")
    try:
        r=fortranCompile(rp.FORTRAN_COMP[0]+" -o muftin.o -c", 
                         "muftin.f", rp.FORTRAN_COMP[1])
        if r:
            logger.error("Error compiling muftin.f, cancelling...")
            return ("Fortran compile error")
        r=fortranCompile(rp.FORTRAN_COMP[0]+" -o lib.tleed.o -c", 
                                     libname, rp.FORTRAN_COMP[1])
        if r:
            logger.error("Error compiling "+libname+", cancelling...")
            return ("Fortran compile error")
        r=fortranCompile(rp.FORTRAN_COMP[0]+" -o main.o -c", srcname,
                            rp.FORTRAN_COMP[1])
        if r:
            logger.error("Error compiling "+srcname+", cancelling...")
            return ("Fortran compile error")
        r=fortranCompile(rp.FORTRAN_COMP[0]+" -o "+rcname, "muftin.o "
                          "lib.tleed.o main.o", rp.FORTRAN_COMP[1])
        if r:
            logger.error("Error compiling fortran files, cancelling...")
            return ("Fortran compile error")
        logger.debug("Compiled fortran files successfully")
    except:
        logger.error("Error compiling fortran files: ")
        raise
    rclogname = rcname+".log"
    logger.info("Starting reference calculation...\n"
        "Refcalc log will be written to file "+rclogname)
    rp.manifest.append(rclogname)
    try:
        with open(rclogname, "w") as log:
            subprocess.run(os.path.join('.',rcname), input=fin, 
                           encoding="ascii", stdout=log, stderr=log)
    except:
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
    except:
        logger.error("Error reading fd.out after reference calculation. "
                      "Check settings and refcalc log.")
        raise
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
    try:
        beams.writeOUTBEAMS(rp.theobeams["refcalc"], filename="THEOBEAMS.csv")
        theobeams_norm = copy.deepcopy(rp.theobeams["refcalc"])
        for b in theobeams_norm:
            b.normMax()
        beams.writeOUTBEAMS(theobeams_norm,filename="THEOBEAMS_norm.csv")
    except:
        logger.error("Error writing THEOBEAMS after reference "
                      "calculation: ", exc_info = True)
        rp.setHaltingLevel(2)
    # rename and move files
    try:
        os.rename('fd.out','refcalc-fd.out')
    except:
        logger.warning("Failed to rename refcalc output file fd.out to "
                        "refcalc-fd.out")
    try:
        os.rename('PARAM','refcalc-PARAM')
    except:
        logger.warning("Failed to rename refcalc input file PARAM to "
                        "refcalc-PARAM")
    # move and zip tensor files
    if not os.path.isdir(os.path.join(".","Tensors")):
        os.mkdir(os.path.join(".","Tensors"))
    rp.TENSOR_INDEX = getMaxTensorIndex() + 1
    rp.manifest.append("Tensors")
    dn = "Tensors_"+str(rp.TENSOR_INDEX).zfill(3)
    if not os.path.isdir(os.path.join(".","Tensors",dn)):
        os.mkdir(os.path.join(".","Tensors",dn))
    try:
        for tf in [f for f in os.listdir('.') if f.startswith("T_")]:
            shutil.move(tf, os.path.join(".","Tensors",dn,tf))
    except:
        logger.error("Error moving Tensor files: ")
        raise
    tInputFiles = ["POSCAR", "PARAMETERS", "VIBROCC", "IVBEAMS", 
                   "_PHASESHIFTS"]
    for f in [f for f in tInputFiles if f in os.listdir('.')]:
        of = f
        for fn in ["POSCAR","VIBROCC"]:
            if (f == fn and 3 in rp.runHistory 
                        and os.path.isfile(fn + "_OUT_" + rp.timestamp)):
                of = fn + "_OUT_" + rp.timestamp
        try:
            shutil.copy2(of, os.path.join(".","Tensors",dn,f))
        except:
            logger.warning("Failed to add input file " + f
                            + " to Tensors folder.")
    # delete old delta files in main work folder, if necessary
    #   (there should not be any, unless there was an error)
    for df in [f for f in os.listdir(".") if f.startswith("DEL_") and 
               os.path.isfile(os.path.join(".", f))]:
        try:
            os.remove(df)
        except:
            logger.warning("Error deleting old Delta file in work "
                    "directory. This may cause the delta file to "
                    "incorrectly be labelled as belonging with the new "
                    "set of tensors.")
    return 0