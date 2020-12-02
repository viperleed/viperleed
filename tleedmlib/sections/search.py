# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Search
"""

import os
import logging
import shutil
import subprocess
import time
from timeit import default_timer as timer
import numpy as np
import signal
import re

import tleedmlib.files.iosearch as io
import tleedmlib as tl
from tleedmlib.leedbase import fortranCompile
from tleedmlib.files.parameters import updatePARAMETERS
from tleedmlib.files.displacements import readDISPLACEMENTS_block
from tleedmlib.files.searchpdf import (writeSearchProgressPdf, 
                                       writeSearchReportPdf)

logger = logging.getLogger("tleedm.search")

def processSearchResults(sl, rp, final=True):
    """Looks at the SD.TL file and gets the data from the best structure in 
    the last generation into the slab. The "final" tag defines whether this is 
    the final interpretation that *must* work and should go into the log, or 
    if this will be repeated anyway. Then calls writeSearchOutput to write 
    POSCAR_OUT, VIBROCC_OUT, and modify data in sl and rp."""
    # get the last block from SD.TL:
    try:
        lines = io.readSDTL_end()
    except FileNotFoundError:
        logger.error("Could not process Search results: SD.TL file not "
                      "found.")
        raise
    except:
        logger.error("Failed to get last block from SD.TL file.")
        raise
    # read the block
    sdtlContent = io.readSDTL_blocks("\n".join(lines), 
                                     whichR = rp.SEARCH_BEAMS, logInfo = final)
    if not sdtlContent:
        if final:
            logger.error("No data found in SD.TL file!")
            rp.setHaltingLevel(2)
        return("No data in SD.TL")
    (generation, rfacs, configs) = sdtlContent[0]
    # collect equal entries
    pops = [] # tuples (r, pars) with 
              #   r..rfactor, pars..list of parameter indices
    popcount = []   # how often a given population is there
    dparlist = []    # for control.chem
    for i in range(0, len(rfacs)):
        if (rfacs[i], configs[i]) in pops:
            popcount[pops.index((rfacs[i], configs[i]))] += 1
        else:
            pops.append((rfacs[i], configs[i]))
            popcount.append(1)
        dparlist.append(configs[i])
    writeControlChem = False
    if not os.path.isfile("control.chem"):
        writeControlChem = True
    else:
        # check file, which generation
        try:
            with open("control.chem", "r") as rf:
                rf.readline() # first line is empty.
                s = rf.readline()
                s = s.split("No.")[-1].split(":")[0]
                if int(s) != generation:
                    writeControlChem = True
                if len(rf.readlines()) < rp.SEARCH_POPULATION:
                    writeControlChem = True
        except:
            # try overwriting, just to be safe
            writeControlChem = True
    # write backup file
    output = ("\nParameters of generation No." + str(generation).rjust(6)
              + ":\n")
    if rp.domainParams:
        astep = rp.DOMAIN_STEP
    else:
        astep = 100
    for dpars in dparlist:
        ol = ""
        areapars = []
        for (percent, l) in dpars:
            for ind in l:
                ol += str(ind).rjust(3)
            areapars.append(int(percent/astep)+1)
        for ap in areapars:
            ol += str(ap).rjust(3)
        output += ol + "\n"
    rp.controlChemBackup = output
    if writeControlChem:
        try:
            with open("control.chem", "w") as wf:
                wf.write(output)
        except:
            # if final:
            logger.error("Failed to write control.chem")
            rp.setHaltingLevel(1)
        
    # info for log:
    maxpop = max(popcount)
    if final:
        if len(pops) == 1:
            logger.info("All trial structures converged to the same "
                    "configuration (R = {:.4f})".format(pops[0][0]))
        elif any([v > rp.SEARCH_POPULATION/2 for v in popcount]):
            info = ("The search outcome is dominated by one configuration "
                    "(population {} of {}, R = {:.4f})\n".format(maxpop, 
                    rp.SEARCH_POPULATION, pops[popcount.index(maxpop)][0]))
        else:
            info = ("The search outcome is not dominated by any "
                    "configuration. The best result has population {} (of {})"
                    ", R = {:.4f}\n".format(popcount[0], rp.SEARCH_POPULATION, 
                                            pops[0][0]))
        if len(pops) > 1:
            info += ("The best configurations are:\nPOP       R |")
            if rp.domainParams:
                info += " area |"
            info += " PARAMETERS\n"
            for i in range(0,min(5,len(pops))):
                for (j, (percent, pars)) in enumerate(pops[i][1]):
                    if j == 0:
                        info += "{:>3}  {:.4f} |".format(popcount[i], 
                                                         pops[i][0])
                    else:
                        info += "            |"
                    if rp.domainParams:
                        info += " {:>3}% |".format(percent)
                    for v in pars:
                        info += "{:>3}".format(v)
                    info += "\n"
            logger.info(info)
    # now writeSearchOutput:
    if not rp.domainParams:
        result = io.writeSearchOutput(sl, rp, pops[0][1][0][1], 
                                      silent=(not final))
    else:
        home = os.getcwd()
        result = 0
        for (i, dp) in enumerate(rp.domainParams):
            try:
                os.chdir(dp.workdir)
                r = io.writeSearchOutput(dp.sl, dp.rp, pops[0][1][i][1], 
                                         silent=(not final))
            except:
                logger.error("Error while writing search output for domain {}"
                             .format(dp.name), exc_info = rp.LOG_DEBUG)
                rp.setHaltingLevel(2)
            finally:
                os.chdir(home)
            if r != 0:
                result = r
    return result

def search(sl, rp):
    """Runs the search. Returns 0 when finishing without errors, or an error 
    message otherwise."""
    rp.searchResultConfig = None
    if rp.domainParams:
        initToDo = [(dp.rp, dp.sl, dp.workdir) for dp in rp.domainParams]
    else:
        initToDo = [(rp, sl, ".")]
    for (rpt, slt, path) in initToDo:
        # read DISPLACEMENTS block
        if not rpt.disp_block_read:
            readDISPLACEMENTS_block(rpt, slt, 
                                    rpt.disp_blocks[rpt.search_index])
            rpt.disp_block_read = True
        # get Deltas
        if not 2 in rpt.runHistory:
            if "Tensors" in rpt.manifest:
                logger.error("New tensors were calculated, but no new delta "
                              "files were generated. Cannot execute search.")
                return ("Delta calculations was not run for current tensors.")
            try:
                r = tl.leedbase.getDeltas(rpt.TENSOR_INDEX, basedir=path, 
                                          targetdir=path, required=True)
            except:
                raise
            if r != 0:
                return r
    r = rp.updateCores()
    if r != 0:
        return r
    # generate rf.info
    try:
        rfinfo = io.writeRfInfo(sl, rp, filename="rf.info")
    except:
        logger.error("Error generating search input file rf.info")
        raise
    # generate PARAM and search.steu
    #   needs to go AFTER rf.info, as writeRfInfo may remove expbeams!
    try:
        r = io.generateSearchInput(sl, rp)
        if r != 0:
            logger.error("Error generating search input")
            return ("generateSearchInput failed")
    except:
        logger.error("Error generating search input")
        raise
    if rp.indyPars == 0:  # never for calculations with domains # !!! CHECK
        logger.info("Found nothing to vary in search. Will proceed "
                "directly to writing output and starting SUPERPOS.")
        rp.searchResultConfig = [(100, [1] * (len(rp.searchpars))-1)]
        for (i, sp) in enumerate(rp.searchpars):
            if type(sp.restrictTo) == int:
                rp.searchResultsConfig[i] = sp.restrictTo
            elif type(sp.restrictTo) == tl.SearchPar:
                rp.searchResultsConfig[i] = (rp.searchpars.index(
                                                    sp.restrictTo) + 1)
            elif type(sp.linkedTo) == tl.SearchPar:
                rp.searchResultsConfig[i] = (rp.searchpars.index(
                                                    sp.linkedTo) + 1)
        io.writeSearchOutput(sl, rp)
        return 0
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Search "
            " will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return 0
    # check for mpirun, decide whether to use parallelization
    usempi = True
    
    if (shutil.which("mpirun", os.X_OK) == None 
            or shutil.which("mpiifort", os.X_OK) == None):
        usempi = False
        logger.warning("mpirun / mpiifort are not present. Search "
            "will be compiled and executed without parallelization. "
            "This will be much slower!")
        if rp.FORTRAN_COMP[0] == "":
            if rp.getFortranComp() != 0:    #returns 0 on success
                logger.error("No fortran compiler found, cancelling...")
                return ("Fortran compile error")
    else:
        if rp.FORTRAN_COMP_MPI[0] == "":
            rp.FORTRAN_COMP_MPI[0] = "mpiifort -Ofast"
            
            
    if usempi and shutil.which("mpirun", os.X_OK) == None:
        usempi = False
        logger.warning("mpirun is not present. Search will be compiled "
            "and executed without parallelization. This will be much "
            "slower!")
        if rp.FORTRAN_COMP[0] == "":
            if rp.getFortranComp() != 0:    #returns 0 on success
                logger.error("No fortran compiler found, cancelling...")
                return ("Fortran compile error")
    else:
        if rp.FORTRAN_COMP_MPI[0] == "":
            if rp.getFortranMpiComp() != 0:    #returns 0 on success
                logger.error("No fortran mpi compiler found, "
                              "cancelling...")
                return ("Fortran compile error")
    # get fortran files
    try:
        tldir = tl.leedbase.getTLEEDdir(home=rp.workdir)
        srcpath = os.path.join(tldir,'src')
        srcname = [f for f in os.listdir(srcpath) 
                      if f.startswith('search.mpi')][0]
        shutil.copy2(os.path.join(srcpath,srcname), srcname)
        libpath = os.path.join(tldir,'lib')
        libname = [f for f in os.listdir(libpath) 
                      if f.startswith('lib.search.mpi')][0]
        shutil.copy2(os.path.join(libpath,libname), libname)           
        if usempi: # these are short C scripts - use pre-compiled versions
            randnamefrom = "MPIrandom_.o"
        else:
            randnamefrom = "random_.o"
        randname = "random_.o"
        shutil.copy2(os.path.join(libpath,randnamefrom), randname)
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath,globalname), globalname)
    except:
        logger.error("Error getting TensErLEED files for search: ")
        raise
    # compile fortran files
    searchname = "search-"+rp.timestamp
    if usempi:
        fcomp = rp.FORTRAN_COMP_MPI
    else:
        fcomp = rp.FORTRAN_COMP
    logger.info("Compiling fortran input files...")
    try:
        r=fortranCompile(fcomp[0]+" -o lib.search.o -c", 
                            libname, fcomp[1])
        if r:
            logger.error("Error compiling "+libname+", cancelling...")
            return ("Fortran compile error")
        r=fortranCompile(fcomp[0]+" -o restrict.o -c", 
                            "restrict.f", fcomp[1])
        if r:
            logger.error("Error compiling restrict.f, cancelling...")
            return ("Fortran compile error")
        r=fortranCompile(fcomp[0]+" -o search.o -c", srcname,
                            fcomp[1])
        if r:
            logger.error("Error compiling "+srcname+", cancelling...")
            return ("Fortran compile error")
        # combine
        r=fortranCompile(fcomp[0]+" -o "+ searchname, "search.o "
                            "random_.o lib.search.o restrict.o", fcomp[1])
        if r:
            logger.error("Error compiling fortran files, cancelling...")
            return ("Fortran compile error")
        logger.debug("Compiled fortran files successfully")
    except:
        logger.error("Error compiling fortran files: ")
        raise
    if rp.LOG_SEARCH:
        searchlogname = searchname+".log"
        logger.info("Search log will be written to file "+searchlogname)
        rp.manifest.append(searchlogname)
    if rp.N_CORES == 1:
        logger.warning("The N_CORES parameter is set to 1. The search "
                "will be run without multiprocessing. This will be much "
                "slower!")
                # TODO: This shouldn't need MPICOMPILE, but I'm not sure 
                #   if the non-parallelized version of the search is 
                #   equivalent, so currently, both are *compiled* with mpi,
                #   and this only switches off mpirun
        usempi = False
    # if there is an old SD.TL file, it needs to be removed
    if os.path.isfile("SD.TL"):
        try:
            os.remove("SD.TL")
        except:
            logger.warning("Failed to delete old SD.TL file. This may "
                            "cause errors in the interpretation of search "
                            "progress.")
    # start search process
    repeat = True
    first = True
    genOffset = 0
    gens = []           # generation numbers in SD.TL, but continuous if 
                        #   search restarts
    markers = []
    rfaclist = []
    realLastConfig = {"all": [], "best": [], "dec": []}
    realLastConfigGen = {"all": 0, "best": 0, "dec": 0}
    convergedConfig = {"all": None, "best": None, "dec": None}
    lastconfig = None
    rp.searchMaxGenInit = rp.SEARCH_MAX_GEN
    while repeat:
        if first:
            logger.info("Starting search. See files Search-progress.pdf "
                         "and SD.TL for progress information.")
            first = False
        repeat = False
        interrupted = False
        proc = None
        if usempi:
            command = ["mpirun", "-n", str(rp.N_CORES), 
                       os.path.join(".",searchname)]
        else:
            command = os.path.join('.',searchname)
        try:
            if not rp.LOG_SEARCH:
                proc = subprocess.Popen(command, 
                            encoding="ascii", 
                            stdout=subprocess.DEVNULL, 
                            stderr=subprocess.STDOUT,
                            preexec_fn=os.setsid)
            else:
                logExists = os.path.isfile(searchlogname)
                with open(searchlogname, "a") as log:
                    if logExists:
                        log.write("\n\n-------\nRESTARTING\n-------\n\n")
                    proc = subprocess.Popen(command, 
                            encoding="ascii", stdout=log, stderr=log,
                            preexec_fn=os.setsid)
        except:
            logger.error("Error starting search. Check SD.TL file.")
            raise
        if proc == None:
            logger.error("Error starting search subprocess... Stopping.")
            return("Error running search")
        # FEED INPUT
        try:
            proc.communicate(input = rfinfo, timeout = 0.1)
        except subprocess.TimeoutExpired:
            pass # started successfully; monitoring below
        except:
            logger.error("Error starting search. Check SD.TL file.")
        # MONITOR SEARCH
        searchStartTime = timer()
        filepos = 0
        timestep = 1 # time step to check files
        evaluationTime = 30 # how often should SD.TL be evaluated
        lastEval = 0 # last evaluation time, counting in seconds from 
                     #   searchStartTime
        comment = ""
        sdtlGenNum = 0
        gaussianWidthOri = rp.GAUSSIAN_WIDTH
        
        try:
            while proc.poll() == None:
                time.sleep(timestep)
                updatePARAMETERS(rp)
                # check convergence criteria
                stop = False
                checkrepeat = True
                if rp.STOP == True:
                    stop = True
                    checkrepeat = False
                    logger.info("Search stopped by STOP command.")
                    if not os.path.isfile("SD.TL"):
                        # try saving by waiting for SD.TL to be created...
                        logger.warning("SD.TL file not found. Trying to "
                                        "wait, maximum 5 minutes...")
                        i = 0
                        while not os.path.isfile("SD.TL") and not i >= 300:
                            time.sleep(1)
                            i += 1
                else:
                    for k in ["dec", "best", "all"]:
                        if (rp.SEARCH_MAX_DGEN[k] > 0 and len(gens) > 1 
                              and rp.GAUSSIAN_WIDTH_SCALING != 1
                              and gens[-1] - realLastConfigGen[k] >= 
                                                  rp.SEARCH_MAX_DGEN[k]):
                            stop = True
                            o = {"all": "all structures", 
                                 "best": "best structure",
                                 "dec": "best 10% of structures"}
                            logger.info("Search convergence criterion "
                                "reached: max. generations without change "
                                "({}): {}/{}."
                                .format(o[k], 
                                        gens[-1] - realLastConfigGen[k], 
                                        int(rp.SEARCH_MAX_DGEN[k])))
                            break
                if rp.GAUSSIAN_WIDTH != gaussianWidthOri:
                    stop = True
                    repeat = True
                    comment = ("GAUSSIAN_WIDTH = {}"
                               .format(rp.GAUSSIAN_WIDTH))
                    logger.info("GAUSSIAN_WIDTH parameter changed. "
                            "Search will restart.")
                t = timer() - searchStartTime
                if (t - lastEval > evaluationTime) or stop:
                    # evaluate
                    lastEval = t
                    newData = []
                    if os.path.isfile("SD.TL"):
                        filepos, content = io.readSDTL_next(
                                                         offset = filepos)
                        if content:
                            newData = io.readSDTL_blocks(content, 
                                                  whichR = rp.SEARCH_BEAMS)
                    for (gen, rfacs, configs) in newData:
                        gens.append(gen + genOffset)
                        sdtlGenNum = gen
                        rfaclist.append(np.array(rfacs))
                        if gen % 1000 == 0:
                            logger.debug("R = {:.4f} (Generation {})"
                                  .format(min(rfacs), gens[-1]))
                        if configs != realLastConfig["all"]:
                            realLastConfig["all"] = configs
                            realLastConfigGen["all"] = gens[-1]
                        if (configs[:int(np.ceil(
                                            rp.SEARCH_POPULATION * 0.1))] 
                                                != realLastConfig["dec"]):
                            realLastConfig["dec"] = configs[:int(np.ceil(
                                            rp.SEARCH_POPULATION * 0.1))]
                            realLastConfigGen["dec"] = gens[-1]
                        if configs[0] != realLastConfig["best"]:
                            realLastConfig["best"] = configs[0]
                            realLastConfigGen["best"] = gens[-1]
                    if len(newData) > 0:
                        lastconfig = newData[-1][2]
                    if len(gens) > 1:
                        try:
                            writeSearchProgressPdf(rp, gens, rfaclist, 
                                                lastconfig, markers=markers)
                        except KeyboardInterrupt:
                            raise
                        except:
                            logger.warning("Error writing Search-progress.pdf", 
                                            exc_info = rp.LOG_DEBUG)
                        try:
                            writeSearchReportPdf(rp)
                        except KeyboardInterrupt:
                            raise
                        except:
                            logger.warning("Error writing Search-report.pdf",
                                           exc_info = rp.LOG_DEBUG)
                    if (len(gens) > 1 and os.path.isfile("SD.TL") and 
                                                (repeat or not stop)):
                        try:
                            r = processSearchResults(sl, rp, final=False)
                        except Exception as e:
                            r = repr(e)
                        if r != 0:
                            logger.warning("Failed to update POSCAR_OUT "
                                            "and VIBROCC_OUT: " + str(r))
                if stop:
                    logger.info("Stopping search...")
                    pgid = os.getpgid(proc.pid)
                    proc.kill()
                    proc.wait()
                    try:
                        os.killpg(pgid, signal.SIGTERM)  
                                        # needed to kill mpirun children
                    except ProcessLookupError:
                        pass # already dead
                    if (not repeat and not rp.GAUSSIAN_WIDTH_SCALING == 1
                                   and checkrepeat):
                        repeat = True
                        block = False
                        for k in [k for k in ["dec", "best", "all"] 
                                  if rp.SEARCH_MAX_DGEN[k] > 0]:
                            if convergedConfig[k] != realLastConfig[k]:
                                convergedConfig[k] = realLastConfig[k][:]
                                block = True
                            else:
                                repeat = False
                                o = {"all": "any structure", 
                                     "best": "best structure",
                                     "dec": "best 10% of structures"}
                                logger.info("Convergence reached: No "
                                    "improvement to " + o[k] + " since "
                                    "changing GAUSSIAN_WIDTH.")
                        if not repeat and block:
                            repeat = True
                            logger.info("Multiple convergence "
                                "criteria are defined, not all are "
                                "met. Search continues.")
                        if repeat:
                            rp.GAUSSIAN_WIDTH *= rp.GAUSSIAN_WIDTH_SCALING
                            if rp.GAUSSIAN_WIDTH < 0.0001:
                                rp.GAUSSIAN_WIDTH = 0.0001
                            for k in ["dec", "best", "all"]:
                                rp.SEARCH_MAX_DGEN[k] *= (
                                            rp.SEARCH_MAX_DGEN_SCALING[k])
                                realLastConfigGen[k] = gens[-1]
                            logger.info("Reducing GAUSSIAN_WIDTH "
                                "parameter to {} and restarting search..."
                                .format(round(rp.GAUSSIAN_WIDTH,4)))
                            comment = ("GAUSSIAN_WIDTH = {}"
                                   .format(round(rp.GAUSSIAN_WIDTH,4)))
                
        except KeyboardInterrupt:
            if not os.path.isfile("SD.TL"):
                # try saving by waiting for SD.TL to be created...
                logger.warning("SD.TL file not found. Trying to wait, "
                                "interrupt again to override...")
                try:
                    i = 0
                    while not os.path.isfile("SD.TL") and not i >= 60:
                        time.sleep(1)
                        i += 1
                except KeyboardInterrupt:
                    pass # user insisted, give up
            interrupted = True
            rp.STOP = True
            try:
                pgid = os.getpgid(proc.pid)
                proc.kill()
                proc.wait()
                os.killpg(pgid, signal.SIGTERM)
                                    # needed to kill mpirun children
            except ProcessLookupError:
                pass # already dead
            logger.warning("Search interrupted by user. Attempting "
                            "analysis of results...")
        except:
            logger.error("Error during search. Check SD.TL file.")
            raise
        
        if repeat:
            rp.SEARCH_START = "control"
            if gens:
                genOffset = gens[-1]
                rp.SEARCH_MAX_GEN -= sdtlGenNum
                markers.append((genOffset, comment))
            try:
                r = io.generateSearchInput(sl, rp, steuOnly=True, 
                                           cull=True)
                if r != 0:
                    logger.error("Error re-generating search input")
                    return ("generateSearchInput failed")
            except:
                logger.error("Error re-generating search input")
                raise
            if os.path.isfile("SD.TL"):
                try:
                    os.remove("SD.TL")
                except:
                    logger.warning("Failed to delete old SD.TL file. "
                                "This may cause errors in the "
                                "interpretation of search progress.")
    if not interrupted:
        logger.info("Finished search. Processing files...")
    else:
        logger.info("Processing files...")
    # write pdf one more time
    if len(gens) > 1:
        try:
            writeSearchProgressPdf(rp, gens, rfaclist, lastconfig,
                                      markers=markers)
        except:
            logger.warning("Error writing Search-progress.pdf", 
                            exc_info = True)
        try:
            writeSearchReportPdf(rp)
        except:
            logger.warning("Error writing Search-report.pdf",
                            exc_info = True)
    # process SD.TL to get POSCAR_OUT, VIBROCC_OUT
    r = processSearchResults(sl, rp)
    if r != 0:
        logger.error("Error processing search results: "+r)
        return r
    # if deltas were copied from domain folders, clean them up
    if rp.domainParams:
        rgx = re.compile(r'D\d+_DEL_')
        for file in [f for f in os.listdir() if rgx.match(f)]:
            try:
                os.remove(file)
            except:
                logger.warning('Failed to deleted redundant domain delta file '
                               +file)
    # process files
    try:
        os.rename('PARAM','search-PARAM')
    except:
        logger.warning("Failed to rename search input file PARAM to "
                        "search-PARAM")
    try:
        os.rename('rf.info','search-rf.info')
    except:
        logger.warning("Failed to rename search input file rf.info to "
                        "search-rf.info")
    if lastconfig != None:
        rp.searchResultConfig = lastconfig
    return 0