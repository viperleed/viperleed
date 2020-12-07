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

import tleedmlib.files.iosuperpos as io
from tleedmlib.leedbase import getDeltas, getTLEEDdir, fortranCompile
from tleedmlib.files.beams import writeOUTBEAMS, averageBeams, writeFdOut
from tleedmlib.files.displacements import readDISPLACEMENTS_block
from tleedmlib.files.iosearch import readSDTL_end, readSDTL_blocks
from tleedmlib.files.iorefcalc import readFdOut

logger = logging.getLogger("tleedm.superpos")

def superpos(sl, rp, subdomain=False):
    """Runs the superpos calculation. Returns 0 when finishing without 
    errors, or an error message otherwise."""
    # check whether there is anything to evaluate
    if rp.searchResultConfig != None:
        config = rp.searchResultConfig[0]
    else:
        config = None
        #check for an SD.TL file
        sdtl = None
        if os.path.isfile("SD.TL"):
            try:
                sdtl = readSDTL_end(filename="SD.TL")
            except:
                logger.error("Superpos: Error reading SD.TL")
                rp.setHaltingLevel(2)
                return 0
        elif os.path.isfile(os.path.join("OUT","SD.TL")):
            try:
                sdtl = readSDTL_end(filename = os.path.join("OUT", "SD.TL"))
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
        sdtlContent = readSDTL_blocks("\n".join(sdtl), 
                                      whichR = rp.SEARCH_BEAMS)
        if not sdtlContent:
            logger.error("Superpos: No data found in SD.TL")
            rp.setHaltingLevel(2)
            return 0
        try:
            config = sdtlContent[0][2][0]  # first block, config, best only
        except:
            logger.error("Superpos: Failed to read best "
                          "configuration from SD.TL")
            rp.setHaltingLevel(2)
            return 0
    
    if rp.domainParams:
        try:
            r = superpos_domains(rp, config)
        except:
            raise
        if r != 0:
            return r
        return 0
    
    # read DISPLACEMENTS block and fetch deltas
    if not rp.disp_block_read:
        readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
        rp.disp_block_read = True
    if not any([ind in rp.runHistory for ind in [2, 3]]):
        try:
            r = getDeltas(rp.TENSOR_INDEX, required=True)
        except:
            raise
        if r != 0:
            return r
    # make sure search parameters are initialized
    if not 3 in rp.runHistory and not subdomain:
        logger.debug("Superpos calculation executed without search. "
                     "Search parameters will be inferred from input files.")
        r = rp.generateSearchPars(sl)
        if r != 0:
            logger.error("Error getting search parameters. Superpos will "
                          "stop.")
            return 0
    # now we have configuration and parameters, create input:
    contrin = ""
    try:
        contrin = io.writeSuperposInput(sl, rp, config[0][1])
        logger.debug("Wrote Superpos input successfully")
    except:
        logger.error("Error getting input data for Superpos: ", 
                      exc_info = True)
        rp.setHaltingLevel(2)
        return 0
    if contrin == "":
        logger.error("Error getting input data for Superpos: "
                     "writeSuperposInput returned empty. Cancelling "
                     "Superpos...")
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
        tldir = getTLEEDdir(home=rp.workdir, version=rp.TL_VERSION)
        if not tldir:
            return("TensErLEED code not found.")
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
        r=fortranCompile(rp.FORTRAN_COMP[0]+" -o", sposname+" "
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
        rp.theobeams["superpos"], rp.superpos_specout = readFdOut(outname)
    except FileNotFoundError:
        logger.error(outname + " not found after superpos calculation.")
        raise
    except:
        logger.error("Error reading "+ outname +" after superpos "
                      " calculation.")
        raise
    try:
        writeOUTBEAMS(rp.theobeams["superpos"], filename="FITBEAMS.csv")
        theobeams_norm = copy.deepcopy(rp.theobeams["superpos"])
        for b in theobeams_norm:
            b.normMax()
        writeOUTBEAMS(theobeams_norm,filename="FITBEAMS_norm.csv")
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

def superpos_domains(rp, config):
    """Runs the normal superpos function for each subdomain, collects the 
    results and averages over the beams to get an overall result."""
    # make sure search parameters are initialized
    if not 3 in rp.runHistory:
        logger.debug("Superpos calculation executed without search. "
                     "Search parameters will be inferred from input files.")
        r = rp.generateSearchPars(None)
        if r != 0:
            logger.error("Error getting search parameters. Superpos will "
                          "stop.")
            return 0
    home = os.getcwd()
    percentages = []
    for (i, dp) in enumerate(rp.domainParams):
        (percent, params) = config[i]
        percentages.append(percent)
        dp.rp.searchResultConfig = [[(100, params)]]
        logger.info("Running superpos calculation for domain {}"
                    .format(dp.name))
        try:
            os.chdir(dp.workdir)
            r = superpos(dp.sl, dp.rp, subdomain=True)
        except:
            logger.error("Error while running superpos calculation for domain "
                         "{}".format(dp.name))
            raise
        finally:
            os.chdir(home)
        if r != 0:
            return r
    logger.info("Getting weighted average over domain beams...")
    rp.theobeams["superpos"] = averageBeams([dp.rp.theobeams["superpos"] 
                            for dp in rp.domainParams], weights=percentages)
    try:
        writeOUTBEAMS(rp.theobeams["superpos"], filename="FITBEAMS.csv")
        theobeams_norm = copy.deepcopy(rp.theobeams["superpos"])
        for b in theobeams_norm:
            b.normMax()
        writeOUTBEAMS(theobeams_norm,filename="FITBEAMS_norm.csv")
    except:
        logger.error("Error writing FITBEAMS after superpos calculation.",
                     exc_info = rp.LOG_DEBUG)
    try:
        rp.superpos_specout = writeFdOut(rp.theobeams["superpos"], rp.beamlist,
                                         filename="superpos-spec.out", 
                                         header=rp.systemName)
    except:
        logger.error("Error writing averaged superpos-spec.out for R-factor "
                     "calculation.", exc_info = rp.LOG_DEBUG)
    return 0