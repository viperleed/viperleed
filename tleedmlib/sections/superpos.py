# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Superpos
"""

import os
import logging
import copy
import shutil
import subprocess

import tleedmlib as tl

logger = logging.getLogger("tleedm.superpos")

def superpos(sl, rp):
    """Runs the superpos calculation. Returns 0 when finishing without 
    errors, or an error message otherwise."""
    # read DISPLACEMENTS block
    if not rp.disp_block_read:
        tl.readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
        rp.disp_block_read = True
    if not (2 in rp.runHistory or 3 in rp.runHistory):
        try:
            r = tl.getDeltas(rp.TENSOR_INDEX, required=True)
        except:
            raise
        if r != 0:
            return r
    # check whether there is anything to evaluate
    if rp.searchResultConfig != None:
        config = rp.searchResultConfig[0]
    else:
        config = None
        #check for an SD.TL file
        sdtl = None
        if os.path.isfile("SD.TL"):
            try:
                sdtl = tl.readSDTL_end(filename="SD.TL")
            except:
                logger.error("Superpos: Error reading SD.TL")
                rp.setHaltingLevel(2)
                return 0
        elif os.path.isfile(os.path.join("OUT","SD.TL")):
            try:
                sdtl = tl.readSDTL_end(filename = os.path.join("OUT",
                                                               "SD.TL"))
            except:
                logger.error("Superpos: Error reading SD.TL")
                rp.setHaltingLevel(2)
                return 0
        else:
            logger.error("Superpos: Found no stores results from recent "
                          "search and no SD.TL file. Cancelling...")
            rp.setHaltingLevel(2)
            return 0
        if sdtl == None:
            logger.error("Superpos: No data found in SD.TL")
            rp.setHaltingLevel(2)
            return 0
        else:
            config = None
            for line in sdtl:
                if "|" in line and line[:3] != "IND":
                    # this is a line with data in it
                    try:
                        valstring = line.split("|")[-1].rstrip()
                        pars = []
                        while len(valstring) > 0:
                            pars.append(int(valstring[:4]))
                            valstring = valstring[4:]
                        config = pars
                    except:
                        logger.error("Could not read line in SD.TL:\n"
                                      +line)
                    if config:
                        break # we're only interested in the best config
        if not config:
            logger.error("Superpos: Failed to read best "
                          "configuration from SD.TL")
            rp.setHaltingLevel(2)
            return 0
    # make sure search parameters are initialized
    if not 3 in rp.runHistory:
        logger.debug("Superpos calculation executed without search. "
                     "Search parameters will be inferred from input files.")
        r = rp.generateSearchPars(sl, rp)
        if r != 0:
            logger.error("Error getting search parameters. Superpos will "
                          "stop.")
            return 0
    # now we have configuration and parameters, create input:
    contrin = ""
    try:
        contrin = tl.writeCONTRIN(sl, rp, config)
        logger.debug("Wrote to Superpos-CONTRIN successfully")
    except:
        logger.error("Error getting input data for Superpos: ", 
                      exc_info = True)
        rp.setHaltingLevel(2)
        return 0
    if contrin == "":
        logger.error("Error getting input data for Superpos: "
                     "writeCONTRIN returned empty. Cancelling Superpos...")
        rp.setHaltingLevel(2)
        return 0
    try:
        tl.writeSuperposPARAM(rp)
    except:
        logger.error("Error writing PARAM file for Superpos: ",
                      exc_info = True)
        rp.setHaltingLevel(2)
        return 0
    # if execution is suppressed, stop here:
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Superpos "
            "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return 0
    if rp.FORTRAN_COMP[0] == "":
        if rp.getFortranComp() != 0:    #returns 0 on success
            logger.error("No fortran compiler found, cancelling...")
            return ("Fortran compile error")
    # get fortran files
    try:
        tldir = tl.getTLEEDdir()
        srcpath = os.path.join(tldir,'src')
        srcname = [f for f in os.listdir(srcpath) 
                      if f.startswith('superpos')][0]
        shutil.copy2(os.path.join(srcpath,srcname), srcname)
        libpath = os.path.join(tldir,'lib')
        libname = [f for f in os.listdir(libpath) 
                      if f.startswith('lib.superpos')][0]
        shutil.copy2(os.path.join(libpath,libname), libname)           
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath,globalname), globalname)
    except:
        logger.error("Error getting TensErLEED files for superpos: ")
        raise
    # compile fortran files
    sposname = "superpos-"+rp.timestamp
    logger.info("Compiling fortran input files...")
    try:
        r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o", sposname+" "
                          + srcname + " " + libname, rp.FORTRAN_COMP[1])
        if r:
            logger.error("Error compiling fortran files, cancelling...")
            return ("Fortran compile error")
        logger.debug("Compiled fortran files successfully")
    except:
        logger.error("Error compiling fortran files: ")
        raise
    logger.info("Starting Superpos calculation...")
    outname = "superpos-spec.out"
    try:
        with open(outname, "w") as out:
            subprocess.run(os.path.join('.',sposname), 
                           input=contrin, encoding="ascii",
                           stdout=out, stderr=subprocess.STDOUT)
    except:
        logger.error("Error during Superpos calculation.")
        raise
    logger.info("Finished Superpos calculation. Processing files...")
    try:
        theobeams, rp.superpos_specout = tl.readFdOut(outname)
    except FileNotFoundError:
        logger.error(outname + " not found after superpos calculation.")
        raise
    except:
        logger.error("Error reading "+ outname +" after superpos "
                      " calculation.")
        raise
    try:
        tl.writeOUTBEAMS(theobeams, filename="FITBEAMS.csv")
        theobeams_norm = copy.deepcopy(theobeams)
        for b in theobeams_norm:
            b.normMax()
        tl.writeOUTBEAMS(theobeams_norm,filename="FITBEAMS_norm.csv")
    except:
        logger.error("Error writing FITBEAMS after superpos calculation.")
    # rename and move files
    try:
        os.rename('PARAM','superpos-PARAM')
    except:
        logger.warning("Failed to rename superpos input file PARAM to "
                        "superpos-PARAM")
    try:
        os.rename('DOC','superpos-DOC')
    except:
        pass
    return 0