# -*- coding: utf-8 -*-
"""
Created on Nov 12 2019

@author: Florian Kraushofer

Master script running the TensErLEED Manager.
"""

import sys
import time
from timeit import default_timer as timer
import logging
import os
import shutil
import copy
import subprocess
#import fortranformat as ff
import numpy as np
import re
import multiprocessing
import hashlib
import signal

cd = os.path.realpath(os.path.dirname(__file__))  # !!! check this.
tleedmap_path = os.path.realpath(os.path.join(cd, '..'))
if tleedmap_path not in sys.path:
    sys.path.append(tleedmap_path)

import tleedmlib as tl

logger = logging.getLogger("tleedm")


starttime = 0.0

class CustomLogFormatter(logging.Formatter):
    """Logging Formatter for level-dependent message formatting"""

    FORMATS = {
#        logging.DEBUG: "DBG: %(module)s: %(lineno)d: %(msg)s",
        logging.DEBUG: "dbg: %(msg)s",
        logging.INFO: "%(msg)s",
        logging.WARNING: "# WARNING: %(msg)s",
        logging.ERROR: "### ERROR ### in %(module)s:%(funcName)s:%(lineno)s\n"
                       "# %(msg)s \n#############",
        logging.CRITICAL: "### CRITICAL ### in %(module)s:%(funcName)s:"
                          "%(lineno)s\n# %(msg)s \n################",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
    
class DeltaCompileTask():
    """Stores information for a worker to compile a delta file, and keeps 
    track of the folder that the compiled file is in afterwards."""
    def __init__(self, param, hash_, index):
        self.param = param
        self.hash = hash_
        self.foldername = "Delta_Compile_{}".format(index)
        self.exename = "delta-{}".format(index)
        self.fortran_comp = ["",""]
        
class DeltaRunTask():
    """Stores information needed to copy the correct delta executable and 
    tensor file to a subfolder, execute the delta-calculation there, and copy 
    results back."""
    def __init__(self, comptask):
        self.tensorname = ""
        self.din = ""
        self.din_short = ""
        self.deltaname = ""
        self.comptask = comptask
        self.deltalogname = ""

def runDelta(runtask):
    """Function meant to be executed by parallelized workers. Executes a 
    DeltaRunTask."""
    logger = logging.getLogger("tleedm")
    home = os.getcwd()
    workname = "calculating_"+runtask.deltaname
    workfolder = os.path.join(home, workname)
    # make folder and go there:
    if os.path.isdir(workfolder):
        logger.warning("Folder "+workname+" already exists. "
                        "Contents may get overwritten.")
    else:
        os.mkdir(workfolder)
    os.chdir(workfolder)
    # get tensor file
    if os.path.isfile(os.path.join(home,"Tensors",runtask.tensorname)):
        try:
            shutil.copy2(os.path.join(home,"Tensors",runtask.tensorname),"AMP")
        except:
            logger.error("Error copying Tensor file: ", exc_info = True)
            return ("Error encountered by DeltaRunTask " + runtask.deltaname
                    + ": Error copying Tensor file.")
    else:
        logger.error("Tensor file not found: " + runtask.tensorname)
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Tensor not found.")
    # get executable
    exename = runtask.comptask.exename
    try:
        shutil.copy2(os.path.join(home, runtask.comptask.foldername, exename),
                     os.path.join(workfolder, exename))
    except:
        logger.error("Error getting delta executable: ", exc_info = True)
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Failed to get delta executable.")
    # run execution
    try:
        with open("delta.log", "w") as log:
            subprocess.run(os.path.join(workfolder, exename), 
                           input=runtask.din, encoding="ascii", 
                           stdout=log, stderr=log)
    except:
        logger.error("Error while executing delta-amplitudes "
            "calculation for " + runtask.deltaname + ". Also check delta "
            "log file.")
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Error during delta execution.")
    # copy delta file out
    try:
        shutil.copy2(os.path.join(workfolder, "DELWV"), 
                     os.path.join(home, runtask.deltaname))
    except:
        logger.error("Failed to copy delta output file DELWV"
                        " to main folder as" + runtask.deltaname)
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Failed to copy result file out.")
    # append log
    log = ""
    try:
        with open("delta.log", "r") as rf:
            log = rf.read()
    except:
        logger.warning("Could not read local delta log for "
                        + runtask.deltaname)
    if log != "":
        deltalog = os.path.join(home, runtask.deltalogname)
        try:
            with open(deltalog, "a") as wf:
                wf.write("\n\n### STARTING LOG FOR " + runtask.deltaname 
                         + " ###\n\n" + log)
        except:
            logger.warning("Error writing delta log part "
                            + runtask.deltaname + ". The exception was: ",
                            exc_info = True)
    # clean up
    os.chdir(home)
    try:
        shutil.rmtree(workfolder)
    except:
        logger.warning("Error deleting folder " + workname)
    return 0

def compileDelta(comptask):
    """Function meant to be executed by parallelized workers. Executes a 
    DeltaCompileTask."""
    home = os.getcwd()
    workfolder = os.path.join(home, comptask.foldername)
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
    except:
        logger.error("Error writing PARAM file: ", exc_info = True)
        return ("Error encountered by DeltaCompileTask " + comptask.foldername
                + "while trying to write PARAM file.")
    # get Fortran source files
    try:
        tldir = tl.getTLEEDdir(home=home)
        libpath = os.path.join(tldir,'lib')
        libname1 = [f for f in os.listdir(libpath) 
                      if f.startswith('lib.tleed')][0]
        shutil.copy2(os.path.join(libpath,libname1), libname1)
        libname2 = [f for f in os.listdir(libpath) 
                      if f.startswith('lib.delta')][0]
        shutil.copy2(os.path.join(libpath,libname2), libname2)
        srcpath = os.path.join(tldir,'src')
        srcname = [f for f in os.listdir(srcpath) 
                      if f.startswith('delta')][0]
        shutil.copy2(os.path.join(srcpath,srcname), srcname)
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath,globalname), globalname)
    except:
        logger.error("Error getting TensErLEED files for "
                      "delta-amplitudes: ", exc_info = True)
        return ("Error encountered by DeltaCompileTask " + comptask.foldername
                + "while trying to fetch fortran source files")
    # compile
    try:
        for (fname, oname) in [(srcname, "main.o"), 
                               (libname1, "lib.tleed.o"), 
                               (libname2, "lib.delta.o")]:
            r = tl.fortranCompile(comptask.fortran_comp[0]+" -o "
                    +oname+" -c", fname, comptask.fortran_comp[1])
            if r:
                logger.error("Error compiling "+srcname)
                return ("Fortran compile error in DeltaCompileTask "
                        + comptask.foldername)
        r=tl.fortranCompile(comptask.fortran_comp[0]+" -o " + comptask.exename
                            +" main.o lib.tleed.o lib.delta.o",
                            comptask.fortran_comp[1])
        if r:
            logger.error("Error compiling fortran files")
            return ("Fortran compile error in DeltaCompileTask "
                        + comptask.foldername)
    except:
        logger.error("Error compiling fortran files: ")
        return ("Fortran compile error in DeltaCompileTask "
                    + comptask.foldername)
    os.chdir(home)
    return 0

def getTensors(index, required=True):
    """Fetches Tensor files from Tensors or archive with specified tensor 
    index. If required is set True, an error will be printed if no Delta files 
    are found."""
    dn = "Tensors_"+str(index).zfill(3)
    if not os.path.isdir(os.path.join(".","Tensors",dn)):
        if os.path.isfile(os.path.join(".","Tensors",dn+".zip")):
            try:
                logger.info("Unpacking {}.zip...".format(dn))
                os.mkdir(os.path.join(".","Tensors",dn))
                shutil.unpack_archive(os.path.join(".","Tensors",
                                                   dn+".zip"),
                                      os.path.join(".","Tensors",dn))
            except:
                logger.error("Failed to unpack {}.zip".format(dn))
                raise
        else:
            logger.error("Tensors not found")
            return ("Tensors not found")
    return 0

def getDeltas(index, required=True):
    """Fetches Delta files from Deltas or archive with specified tensor index. 
    If required is set True, an error will be printed if no Delta files are 
    found."""
    dn = "Deltas_"+str(index).zfill(3)
    if os.path.isdir(os.path.join(".","Deltas",dn)):
        for f in [f for f in os.listdir(os.path.join(".","Deltas",dn))
                  if (os.path.isfile(os.path.join(".","Deltas",dn,f)) 
                      and f.startswith("DEL_"))]:
            try:
                shutil.copy2(os.path.join(".","Deltas",dn,f), ".")
            except:
                logger.error("Could not copy existing delta files to "
                              "work directory")
                raise
    elif os.path.isfile(os.path.join(".","Deltas",dn+".zip")):
        try:
            logger.info("Unpacking {}.zip...".format(dn))
            shutil.unpack_archive(os.path.join(".","Deltas",dn+".zip"),
                                  ".")
        except:
            logger.error("Failed to unpack {}.zip".dn)
            raise
    elif required:
        logger.error("Deltas not found")
        return ("Deltas not found")
    return 0


def runSection(index, sl, rp):
    """Runs a specific part of the program. Returns 0 when finishing without
    errors, or an error message otherwise."""
    sectionNames = {0: "INITIALIZATION", 
                    1: "REFERENCE CALCULATION", 
                    2: "DELTA-AMPLITUDES", 
                    3: "SEARCH",
                    11: "R-FACTOR CALCULATION", 
                    12: "R-FACTOR CALCULATION", 
                    31: "SUPERPOS"}
    requiredFiles = {0: ["POSCAR", "PARAMETERS", "VIBROCC", "IVBEAMS"],
                     1: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS", 
                         "IVBEAMS", "VIBROCC"],
                     2: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS"],
                     3: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS","EXPBEAMS"],
                     11: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "EXPBEAMS"],
                     12: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "EXPBEAMS"],
                     31: ["BEAMLIST", "POSCAR", "PARAMETERS", "IVBEAMS", 
                          "VIBROCC", "DISPLACEMENTS"]}
                # files that need to be there for the different parts to run
    o = "\nSTARTING SECTION: "+sectionNames[index]
    if index == 3 and rp.disp_blocks[rp.search_index][1]:
        o += " "+rp.disp_blocks[rp.search_index][1]  # displacement block name
    logger.info(o)
    sectionStartTime = timer()
    rp.runHistory.append(index)
    i = 0
    while i < len(requiredFiles[index]):
        filename = requiredFiles[index][i]
        ignoreError = False
        if not rp.fileLoaded[filename]:
            # try loading files
            if filename == "EXPBEAMS":
                if len(rp.THEO_ENERGIES) == 0:
                    er = []
                else:
                    er = rp.THEO_ENERGIES[:2]
                try:
                    rp.expbeams = tl.readOUTBEAMS("EXPBEAMS.csv", enrange=er)
                    rp.fileLoaded["EXPBEAMS"] = True
                except FileNotFoundError:
                    try:
                        # try without the .csv extension
                        rp.expbeams = tl.readOUTBEAMS("EXPBEAMS", enrange=er)
                        rp.fileLoaded["EXPBEAMS"] = True
                    except:
                        logger.error("Error while reading required file "
                                      "EXPBEAMS.csv")
                        raise
                except:
                    logger.error("Error while reading required file EXPBEAMS")
                    raise
                if index != 0:
                    tl.checkEXPBEAMS(sl, rp)
            elif filename == "IVBEAMS":
                try:
                    rp.ivbeams = tl.readIVBEAMS()
                    rp.ivbeams_sorted = False
                    rp.fileLoaded["IVBEAMS"] = True
                except FileNotFoundError:
                    if (os.path.isfile("EXPBEAMS") or 
                             os.path.isfile("EXPBEAMS.csv")):
                        requiredFiles[index].insert(i+1, "EXPBEAMS")
                        logger.warning("IVBEAMS file not found. Will attempt "
                                "generating IVBEAMS from EXBEAMS.")
                        ignoreError = True
                    else:
                        logger.error("Neither IVBEAMS not EXPBEAMS file "
                                      "found.")
                        raise
                except:
                    logger.error("Error while reading required file IVBEAMS")
                    raise
            elif filename == "BEAMLIST":
                try:
                    rp.beamlist = tl.readBEAMLIST()
                    rp.fileLoaded["BEAMLIST"] = True
                except:
                    logger.error("Error while reading required file "
                                  "_BEAMLIST")
                    raise
            elif filename == "VIBROCC":
                try:
                    changeVIBROCC = tl.readVIBROCC(rp,sl)
                    rp.fileLoaded["VIBROCC"] = True
                except:
                    logger.error("Error while reading required file VIBROCC")
                    raise
                sl.fullUpdate(rp)
                if changeVIBROCC:
                    if os.path.isfile("VIBROCC"):
                        os.rename("VIBROCC", "VIBROCC_user")
                        rp.manifest.append("VIBROCC_user")
                        logger.info("VIBROCC file was modified with "
                            "automatically generated vibrational amplitudes.")
                    tl.writeVIBROCC_OUT(sl, rp, "VIBROCC")
                    rp.manifest.append("VIBROCC")
                if rp.T_EXPERIMENT is not None:
                    tl.modifyPARAMETERS(rp, "T_EXPERIMENT", new="")
                if rp.T_DEBYE is not None:
                    tl.modifyPARAMETERS(rp, "T_DEBYE", new="")
                if len(rp.VIBR_AMP_SCALE) > 0:
                    tl.modifyPARAMETERS(rp, "VIBR_AMP_SCALE", new="")
            elif filename == "PHASESHIFTS":
                try:
                    (rp.phaseshifts_firstline, rp.phaseshifts,
                     newpsGen, newpsWrite) = tl.readPHASESHIFTS(sl, rp)
                    if newpsGen:
                        logging.critical("_PHASESHIFT file generation is only "
                            "supported during initialization. Stopping "
                            "execution...")
                        return ("Inconsistent _PHASESHIFT file")
                    elif newpsWrite:
                        logger.warning("Writing a new _PHASESHIFT file is "
                            "only supported during initialization. The "
                            "data in the provided file will be used, but "
                            "running the initialization is recommended.")
                        rp.fileLoaded["PHASESHIFTS"] = True
                    else:
                        rp.fileLoaded["PHASESHIFTS"] = True
                except:
                    logger.error("Error while reading required file "
                                  "_PHASESHIFTS")
                    raise
            elif filename == "DISPLACEMENTS":
                try:
                    r = tl.readDISPLACEMENTS(rp)
                    rp.fileLoaded["DISPLACEMENTS"] = True
                except:
                    logger.error("Error while reading required file "
                                  "DISPLACEMENTS")
                    raise
                if r != 0:
                    return ("Error while reading DISPLACEMENTS file")
            if not rp.fileLoaded[filename] and not ignoreError:
                # and if that didn't work, stop:
                logger.error("Step '"+sectionNames[index]+"' requires file "
                              +filename+". Stopping execution...")
                return ("File not found or file loading error: "+filename)
        i += 1
    
    ####    INITIALIZATION    ####
    
    if index == 0:
        # check for experimental beams:
        expbeamsname = ""
        for fn in ["EXPBEAMS.csv", "EXPBEAMS"]:
            if os.path.isfile(fn):
                expbeamsname = fn
                break
        if expbeamsname:
            if len(rp.THEO_ENERGIES) == 0:
                er = []
            else:
                er = rp.THEO_ENERGIES[:2]
            if not rp.fileLoaded["EXPBEAMS"]:
                try:
                    rp.expbeams = tl.readOUTBEAMS(fn, enrange=er)
                    rp.fileLoaded["EXPBEAMS"] = True
                except:
                    logger.error("Error while reading file "+fn, exc_info=True)
        rp.initTheoEnergies()  # may be initialized based on exp. beams
        # check whether _PHASESHIFTS are present & consistent:
        newpsGen, newpsWrite = True, True
                              # new phaseshifts need to be generated/written
        if os.path.isfile(os.path.join(".","_PHASESHIFTS")):
            try:
                (rp.phaseshifts_firstline, rp.phaseshifts,
                     newpsGen, newpsWrite) = tl.readPHASESHIFTS(sl, rp)
            except:
                logger.warning("Found a _PHASESHIFTS file but could not "
                    "read it. A new _PHASESHIFTS file will be generated."
                    "The exception during read was: ", exc_info=True)
                rp.setHaltingLevel(1)
        if newpsGen:
            try:
                rundgrenpath = os.path.join('.', 'source', 'EEASiSSS.x')
                serneliuspath = os.path.join('.', 'source', 'seSernelius')
                logger.info("Generating phaseshifts data... ")
                (rp.phaseshifts_firstline, 
                            rp.phaseshifts) = tl.runPhaseshiftGen(sl, rp,
                                                   psgensource = rundgrenpath, 
                                                   excosource = serneliuspath)
                logger.debug("Finished generating phaseshift data")
            except:
                logger.error("Exception while calling phaseshiftgen: ")
                raise
        if newpsWrite:
            try:
                tl.writePHASESHIFTS(rp.phaseshifts_firstline, rp.phaseshifts)
            except:
                logger.error("Exception during writePHASESHIFTS: ")
                raise
        rp.fileLoaded["PHASESHIFTS"] = True
        rp.updateDerivedParams()
        rp.manifest.append("_PHASESHIFTS")

        # if necessary, run findSymmetry:
        if sl.planegroup == "unknown":
            sl.findSymmetry(rp)
            sl.enforceSymmetry(rp)
        
        # generate new POSCAR
        tmpslab = copy.deepcopy(sl)
        tmpslab.sortOriginal()
        try:
            tl.writeCONTCAR(tmpslab, filename='POSCAR', comments='all')
        except:
            logger.error("Exception occurred while writing new POSCAR")
            raise
        rp.manifest.append('POSCAR')
        # generate POSCAR_oricell
        tmpslab.revertUnitCell()
        try:
            tl.writeCONTCAR(tmpslab, filename='POSCAR_oricell', 
                            comments='nodir')
        except:
            logger.error("Exception occurred while writing POSCAR_oricell, "
                          "execution will continue...")

        # create bulk slab:
        if sl.bulkslab == tl.DEFAULT:
            sl.bulkslab = sl.makeBulkSlab(rp)
        bsl = sl.bulkslab
        # find minimum in-plane unit cell for bulk:
        logger.info("Checking bulk unit cell...")
        changecell, mincell = bsl.getMinUnitCell(rp)
        if changecell:
            sl.changeBulkCell(rp, mincell)
            bsl = sl.bulkslab
        if not rp.superlattice_defined:
            ws = tl.writeWoodsNotation(rp.SUPERLATTICE) 
                    # !!! replace the writeWoodsNotation from baselib with 
                    #   the one from guilib
            si = rp.SUPERLATTICE.astype(int)
            if ws:
                logger.info("Found SUPERLATTICE = "+ws)
            else:
                logger.info("Found SUPERLATTICE M = {} {}, {} {}".format(
                                        si[0,0], si[0,1], si[1,0], si[1,1]))
        
        # bulk plane group detection:
        logger.info("Initializing bulk symmetry search...")
        bsl.findSymmetry(rp, bulk=True, output=False)
        bsl.revertUnitCell() # keep origin matched with main slab
        logger.info("Found bulk plane group: "+bsl.foundplanegroup)
        bsl.findBulkSymmetry(rp)
        
        # write POSCAR_bulk
        bsl = copy.deepcopy(sl.bulkslab)
        bsl.sortOriginal()
        try:
            tl.writeCONTCAR(bsl, filename='POSCAR_bulk', comments='bulk')
        except:
            logger.error("Exception occurred while writing POSCAR_bulk")
            raise
        
        # generate beamlist
        logger.info("Generating _BEAMLIST...")
        try:
            bgenpath = os.path.join('.', 'source', 'beamgen3.out')
            tl.runBeamGen(sl,rp,beamgensource = bgenpath)
            # this does NOT read the resulting file!
        except:
            logger.error("Exception occurred while calling beamgen.")
            raise
        rp.manifest.append("_BEAMLIST")
        try:
            rp.beamlist = tl.readBEAMLIST()
            rp.fileLoaded["BEAMLIST"] = True
        except:
            logger.error("Error while reading required file _BEAMLIST")
            raise
        
        tl.writePatternInfo(sl, rp)
        
        # if EXPBEAMS was loaded, it hasn't been check yet - check now
        if rp.fileLoaded["EXPBEAMS"]:
            tl.checkEXPBEAMS(sl, rp)
        # write and sort IVBEAMS
        if not rp.fileLoaded["IVBEAMS"]:
            try:
                rp.ivbeams = tl.writeIVBEAMS(sl, rp)
                rp.ivbeams_sorted = False
                rp.fileLoaded["IVBEAMS"] = True
                rp.manifest.append("IVBEAMS")
            except:
                logger.error("Error while writing IVBEAMS file based on "
                              "EXPBEAMS data.")
                raise
        if not rp.ivbeams_sorted:
            rp.ivbeams = tl.sortIVBEAMS(sl, rp)
            rp.ivbeams_sorted = True

    ####    REFERENCE CALCULATION    ####
   
    elif index == 1:
        sl.getCartesianCoordinates(updateOrigin=True)
        sl.updateLayerCoordinates()
        try:
            tl.writePARAM(sl, rp)
        except:
            logger.error("Exception during writePARAM: ")
            raise
        try:
            tl.writeAUXLATGEO(sl, rp)
        except:
            logger.error("Exception during writeAUXLATGEO: ")
            raise
        try:
            tl.writeAUXNONSTRUCT(sl, rp)
        except:
            logger.error("Exception during writeAUXNONSTRUCT: ")
            raise
        try:
            tl.writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
        except:
            logger.error("Exception during writeAUXBEAMS: ")
            raise
        try:
            tl.writeAUXGEO(sl, rp)
        except:
            logger.error("Exception during writeAUXGEO: ")
            raise
        try:
            fin = tl.collectFIN()
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
            tl.writeMuftin(sl, rp)
        except:
            logger.error("Exception during writeMuftin: ")
            raise
        try:
            tldir = tl.getTLEEDdir()
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
            r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o muftin.o -c", 
                                "muftin.f", rp.FORTRAN_COMP[1])
            if r:
                logger.error("Error compiling muftin.f, cancelling...")
                return ("Fortran compile error")
            r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o lib.tleed.o -c", 
                                libname, rp.FORTRAN_COMP[1])
            if r:
                logger.error("Error compiling "+libname+", cancelling...")
                return ("Fortran compile error")
            r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o main.o -c", srcname,
                                rp.FORTRAN_COMP[1])
            if r:
                logger.error("Error compiling "+srcname+", cancelling...")
                return ("Fortran compile error")
            r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o "+rcname, "muftin.o "
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
            theobeams, rp.refcalc_fdout = tl.readFdOut()
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
        if len(rp.ivbeams) != len(theobeams):
            eq = False
        else:
            eq = all([rp.ivbeams[i].isEqual(theobeams[i], eps=eps) for i in 
                                                range(0, len(rp.ivbeams))])
        if not eq:
            logger.error("The list of beams read from IVBEAMS is not "
                "equivalent to the list of beams in the fd.out file "
                "produced by the reference calculation!")
            rp.setHaltingLevel(2)
        try:
            tl.writeOUTBEAMS(theobeams, filename="THEOBEAMS.csv")
            theobeams_norm = copy.deepcopy(theobeams)
            for b in theobeams_norm:
                b.normMax()
            tl.writeOUTBEAMS(theobeams_norm,filename="THEOBEAMS_norm.csv")
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
        rp.TENSOR_INDEX = tl.getMaxTensorIndex() + 1
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

    ####    R-FACTOR    ####
    
    elif index == 11 or index == 12:
        if int((rp.THEO_ENERGIES[1]-rp.THEO_ENERGIES[0]) 
                       / rp.THEO_ENERGIES[2]) + 1 < 2:
            logger.info("Only one theoretical energy found: Cannot calculate "
                         "a meaningful R-Factor. Stopping...")
            return 0
        if ((index == 11 and len(rp.refcalc_fdout) == 0) 
                    or (index == 12 and len(rp.superpos_specout) == 0)):
            if index == 11:
                fn = "refcalc-fd.out"
            else:
                fn = "superpos-spec.out"
            if os.path.isfile(os.path.join(".",fn)):
                logger.warning("R-factor calculation was called without "
                    "stored spectrum data. Reading from file " + fn
                    + " in work folder...")
                path = os.path.join(".",fn)
                try:
                    theobeams, theospec = tl.readFdOut(readfile = path)
                    if index == 11: 
                        rp.refcalc_fdout = theospec
                    else:
                        rp.superpos_specout = theospec
                except:
                    logger.error("Failed to read "+path)
                    raise
            elif os.path.isfile(os.path.join(".","OUT",fn)):
                logger.warning("R-factor calculation was called without "
                    "stored spectrum data. Reading from file " + fn
                    + "in OUT folder...")
                path = os.path.join(".","OUT",fn)
                try:
                    theobeams, theospec = tl.readFdOut(readfile = path)
                    if index == 11: 
                        rp.refcalc_fdout = theospec
                    else:
                        rp.superpos_specout = theospec
                except:
                    logger.error("Failed to read "+path)
                    raise
            else:
                logger.error("Cannot execute R-factor calculation: no stored "
                              "spectrum data and no "+fn+" file was "
                              "found.")
                return("No spectrum data found")
            # if we haven't returned, then it was read. check data vs ivbeams:
            eq = True
            eps = 1e-3
            if len(rp.ivbeams) != len(theobeams):
                eq = False
            else:
                eq = all([rp.ivbeams[i].isEqual(theobeams[i], eps=eps) for i 
                                      in range (0, len(rp.ivbeams))])
            if not eq:
                logger.error("The list of beams read from IVBEAMS is not "
                    "equivalent to the list of beams in "+path+". R-Factor "
                    "calculation cannot proceed.")
                return("Contradiction in beam sets")
        if index == 11:
            theospec = rp.refcalc_fdout
        elif index == 12:
            theospec = rp.superpos_specout
        # WEXPEL before PARAM, to make sure number of exp. beams is correct
        try:
            tl.writeWEXPEL(sl, rp)
        except:
            logger.error("Exception during writeWEXPEL: ")
            raise
        try:
            tl.writeRfactPARAM(rp)
        except:
            logger.error("Exception during writeRfactPARAM: ")
            raise
        # get fortran files and compile
        try:
            tldir = tl.getTLEEDdir()
            libpath = os.path.join(tldir,'lib')
            libname = [f for f in os.listdir(libpath) 
                          if f.startswith('rfacsb')][0]
            shutil.copy2(os.path.join(libpath,libname), libname)
            srcpath = os.path.join(tldir,'src')
            srcname = [f for f in os.listdir(srcpath) 
                          if f.startswith('rfactor.')][0]
            shutil.copy2(os.path.join(srcpath,srcname), srcname)
        except:
            logger.error("Error getting TensErLEED files for r-factor "
                          "calculation: ")
            raise
        if rp.SUPPRESS_EXECUTION:
            logger.warning("SUPPRESS_EXECUTION parameter is on. R-factor "
                "calculation will not proceed. Stopping...")
            rp.setHaltingLevel(3)
            return 0
        logger.info("Compiling fortran input files...")
        rfacname = "rfactor-"+rp.timestamp
        if rp.FORTRAN_COMP[0] == "":
            if rp.getFortranComp() != 0:    #returns 0 on success
                logger.error("No fortran compiler found, cancelling...")
                return ("Fortran compile error")
        try:
            r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o rfacsb.o -c", 
                                libname, rp.FORTRAN_COMP[1])
            if r:
                logger.error("Error compiling "+libname+", cancelling...")
                return ("Fortran compile error")
            r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o main.o -c", srcname,
                                rp.FORTRAN_COMP[1])
            if r:
                logger.error("Error compiling "+srcname+", cancelling...")
                return ("Fortran compile error")
            r=tl.fortranCompile(rp.FORTRAN_COMP[0]+" -o "+rfacname, "rfacsb.o "
                              "main.o", rp.FORTRAN_COMP[1])
            if r:
                logger.error("Error compiling fortran files, cancelling...")
                return ("Fortran compile error")
            logger.debug("Compiled fortran files successfully")
        except:
            logger.error("Error compiling fortran files: ")
            raise
        rfaclogname = rfacname+".log"
        logger.info("Starting R-factor calculation...\n"
            "R-factor log will be written to file "+rfaclogname)
        rp.manifest.append(rfaclogname)
        try:
            with open(rfaclogname, "w") as log:
                subprocess.run(os.path.join('.',rfacname), 
                               input=theospec, encoding="ascii",
                               stdout=log, stderr=log)
        except:
            logger.error("Error during R-factor calculation. Also check "
                "R-factor log file.")
            raise
        logger.info("Finished R-factor calculation. Processing files...")
        if not os.path.isfile(os.path.join(".","ROUT")):
            logger.error("No ROUT file was found after R-Factor calculation!")
            rp.setHaltingLevel(2)
        else:
            try:
                rfac, v0rshift, rfaclist = tl.readROUT()
            except:
                logger.error("Error reading ROUT file")
                rp.setHaltingLevel(2)
            else:
                logger.info("With inner potential shift of {:.2f} eV: "
                             "R = {:.4f}\n"
                             .format(v0rshift, rfac))
                dl = ["."]
                if os.path.isdir("OUT"):
                    dl.append("OUT")
                for dn in dl:
                    for fn in [f for f in os.listdir(dn) 
                               if f.startswith("R_OUT_"+rp.timestamp)
                               and os.path.isfile(os.path.join(dn, f))]:
                        try:  # delete old R_OUT files
                            os.remove(os.path.join(dn, fn))
                        except:
                            pass
                if rfac == 0:
                    logger.error("ROUT reports R-Factor as zero. This means "
                            "something went wrong in the reference "
                            "calculation or in the R-factor calculation.")
                    rp.setHaltingLevel(2)
                    fn = "R_OUT_"+rp.timestamp
                else:
                    fn = "R_OUT_"+rp.timestamp+"_R={:.4f}".format(rfac)
                    rp.last_R = rfac
                try:
                    os.rename("ROUT",fn)
                except:
                    logger.warning("Failed to rename R-factor output file "
                                    "ROUT to "+fn)
                if len(rfaclist) != len(rp.expbeams):
                    logger.warning("Failed to read R-Factors per beam from "
                                    "R-factor output file ROUT.")
                    rfaclist = [-1]*len(rp.expbeams)
                if index == 11:
                    outname = "Rfactor_plots_refcalc.pdf"
                else:
                    outname = "Rfactor_plots_superpos.pdf"
                try:
                    r = tl.writeRfactorPdf([(b.label, rfaclist[i]) for (i, b) 
                                                    in enumerate(rp.expbeams)],
                                           plotcolors = rp.PLOT_COLORS_RFACTOR,
                                           outName = outname)
                except:
                    logger.warning("Error plotting R-factors.")
                else:
                    if r != 0:
                        logger.warning("Error plotting R-factors.")
        # rename and move files
        try:
            os.rename('WEXPEL','rfactor-WEXPEL')
        except:
            logger.warning("Failed to rename R-factor input file WEXPEL to "
                            "rfactor-WEXPEL")
        try:
            os.rename('PARAM','rfactor-PARAM')
        except:
            logger.warning("Failed to rename R-factor input file PARAM to "
                            "rfactor-PARAM")
            
    ####    DELTA-AMPLITUDES    ####
    
    elif index == 2:
        # read DISPLACEMENTS block
        if not rp.disp_block_read:
            tl.readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
            rp.disp_block_read = True
        # get Tensors
        if not os.path.isdir(os.path.join(".","Tensors")):
            logger.error("No Tensors directory found.")
            return("Tensors not found")
        try:
            getTensors(rp.TENSOR_INDEX)
        except:
            raise
        if not 1 in rp.runHistory:
            dn = "Tensors_"+str(rp.TENSOR_INDEX).zfill(3)
            logger.debug("Running without reference calculation, checking "
                "input files in "+dn+" to determine original configuration.")
            r = tl.getTensorOriStates(sl, os.path.join(".","Tensors",dn))
            if r != 0:
                return r
            sl.restoreOriState(keepDisp=True)
        # if there are old deltas, fetch them
        try:
            getDeltas(rp.TENSOR_INDEX, required=False)
        except:
            raise
        dbasic = tl.generateDeltaBasic(sl, rp)
        # get AUXBEAMS; if AUXBEAMS is not in work folder, check AUX folder 
        if not os.path.isfile(os.path.join(".","AUXBEAMS")):
            if os.path.isfile(os.path.join(".","AUX","AUXBEAMS")):
                try:
                    shutil.copy2(os.path.join(".","AUX","AUXBEAMS"), 
                                 "AUXBEAMS")
                except:
                    logger.warning("Failed to copy AUXBEAMS from AUX folder. "
                                    "Generating new file...")
                    try:
                        tl.writeAUXBEAMS(ivbeams=rp.ivbeams, 
                                         beamlist=rp.beamlist)
                    except:
                        logger.error("Exception during writeAUXBEAMS: ")
                        raise
            else:
                try:
                    tl.writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
                except:
                    logger.error("Exception during writeAUXBEAMS: ")
                    raise
        try:
            with open("AUXBEAMS", "r") as rf:
                auxbeams = rf.read()
            if auxbeams[-1] != "\n":
                auxbeams += "\n"
        except:
            logger.error("Could not read AUXBEAMS for delta-input")
            raise
        # get _PHASESHIFTS
        try:
            with open("_PHASESHIFTS", "r") as rf:
                phaseshifts = rf.read()
            if phaseshifts[-1] != "\n":
                phaseshifts += "\n"
        except:
            logger.error("Could not read _PHASESHIFTS for delta-input")
            raise
        
        # go through atoms, remove those that have no variation whatsoever:
        attodo = [at for at in sl.atlist if not at.layer.isBulk]
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
                    if at.disp_vib[el][0] != 0.:
                        found = True
                        break
            if not found:
                for el in at.disp_geo:
                    if np.linalg.norm(at.disp_geo[el][0]) > 0.:
                        found = True
                        break
            if not found:
                occlists = []
                for k in at.disp_occ:
                    occlists.append(at.disp_occ[k])
                for i in range(0,len(occlists[0])):
                    totalocc = 0.
                    for l in occlists:
                        if len(l) <= i:
                            break #error - will pop up again later...
                        else:
                            totalocc += l[i]
                    if totalocc < 1 - 1e-4:
                        found = True
                        break
            if not found:
                attodo.pop(j)
            else:
                j += 1
                
        vaclist = []    #atoms for which a vacancy delta file is needed
        for at in attodo:
            occlists = []
            for k in at.disp_occ:
                occlists.append(at.disp_occ[k])
            for i in range(0,len(occlists[0])):
                totalocc = 0.
                for l in occlists:
                    if len(l) <= i:
                        logger.error("Inconsistent occupancy lists for atom "
                                      +str(at.oriN))
                        return ("Inconsistent input")
                    else:
                        totalocc += l[i]
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
                          if f.startswith("DEL_{}_".format(at.oriN) + el)]
                found = False
                for df in dfiles:
                    if tl.checkDelta(df, at, el, rp):
                        found = True
                        at.deltasGenerated.append(df)
                        countExisting += 1
                        break
                if not found:
                    atElTodo.append((at, el))
        
        if len(atElTodo) == 0:
            logger.info("All Delta files specified in DISPLACEMENTS are "
                    "already present in the Deltas.zip file. Skipping new "
                    "calculations.")
            return 0
        elif countExisting > 0:
            logger.info("{} of {} required Delta-files are already present. "
                         "Generating remaining {} files..."
                         .format(countExisting, len(atElTodo)+countExisting,
                                 len(atElTodo)))
        # create log file:
        deltaname = "delta-"+rp.timestamp
        deltalogname = deltaname+".log"
        logger.info("Generating delta files...\n"
            "Delta log will be written to local subfolders, and collected in "
            +deltalogname)
        rp.manifest.append(deltalogname)
        try:
            with open(deltalogname, "w") as wf:
                wf.write("Logs from multiple delta calculations are collected "
                         "here. Their order may not be preserved.\n")
        except:
            logger.warning("Error creating delta log file. This will not "
                            "affect execution, proceeding...")
        
        # move PARAM file
        if os.path.isfile("PARAM"):
            try:
                os.rename("PARAM", "PARAM-old")
            except:
                try:
                    os.remove(os.path.join(".","PARAM"))
                except:
                    logger.warning("Section Delta-Amplitudes: Cannot rename/"
                        "remove old PARAM file. This might cause the Delta "
                        "generation to fail!")
        # assemble tasks
        deltaCompTasks = [] # keep track of what versions to compile
        deltaRunTasks = [] # which deltas to run
        tensordir = "Tensors_"+str(rp.TENSOR_INDEX).zfill(3)
        for (at, el) in atElTodo:
            din, din_short, param = tl.generateDeltaInput(at, el, sl, rp, 
                                            dbasic, auxbeams, phaseshifts)
            h = hashlib.md5(param.encode()).digest()
            found = False
            for ct in deltaCompTasks:
                if ct.hash == h:
                    found = True
                    rt = DeltaRunTask(ct)
                    break
            if not found:
                index = len(deltaCompTasks)
                ct = DeltaCompileTask(param, h, index)
                deltaCompTasks.append(ct)
                rt = DeltaRunTask(ct)
            deltaRunTasks.append(rt)
            rt.din = din
            rt.din_short = din_short
            rt.tensorname = os.path.join(tensordir, "T_{}".format(at.oriN))
            nameBase = "DEL_{}_".format(at.oriN) + el
            n = 1
            nums = []
            for fn in [f for f in os.listdir(".") if f.startswith(nameBase)]:
                try:
                    nums.append(int(fn.split("_")[-1]))
                except:
                    pass
            if nums:
                n = max(nums) + 1
            rt.deltaname = nameBase + "_{}".format(n)
            rt.deltalogname = deltalogname
            at.deltasGenerated.append(rt.deltaname)
            
        # write delta-input file
        dinput = ("""# ABOUT THIS FILE:
# Input for the delta-calculations is collected here. The blocks of data are
# new 'PARAM' files, which are used to recompile the fortran code, and input  
# for generation of specific DELTA files. Lines starting with '#' are comments 
# on the function of the next block of data.
# In the DELTA file blocks, [AUXBEAMS] and [_PHASESHIFTS] denote where the  
# entire contents of the AUXBEAMS and _PHASESHIFTS files should be inserted.
""")
        for ct in deltaCompTasks:
    #         wf.write("\n#### NEW 'PARAM' FILE: ####\n\n")
    #         wf.write(param+"\n")
            dinput += ("\n#### NEW 'PARAM' FILE: ####\n\n" + ct.param + "\n")
            for rt in [t for t in deltaRunTasks if t.comptask == ct]:
                dinput += ("\n#### INPUT for new DELTA file {}: ####\n\n"
                           .format(rt.deltaname) + rt.din_short + "\n")
        try:
            with open("delta-input", "w") as wf:
                wf.write(dinput)
        except:
            logger.warning("Failed to write file 'delta-input'. This will "
                            "not affect TensErLEED execution, proceeding...")
            
        # if execution is suppressed, stop here
        if rp.SUPPRESS_EXECUTION:
            rp.setHaltingLevel(3)
            return 0
        
        # make sure there's a compiler ready:
        if rp.FORTRAN_COMP[0] == "":
            if rp.getFortranComp() != 0:    #returns 0 on success
                logger.error("No fortran compiler found, "
                              "cancelling...")
                return ("No Fortran compiler")
        for ct in deltaCompTasks:
            ct.fortran_comp = rp.FORTRAN_COMP
            
        # if number of cores is not defined, try to find it
        if rp.N_CORES == 0:
            try:
                rp.N_CORES = tl.available_cpu_count()
            except:
                logger.error("Failed to detect number of cores.")
            logger.info("Automatically detected number of available CPUs: {}"
                         .format(rp.N_CORES))
        if rp.N_CORES == 0:
            logger.error("Failed to detect number of cores.")
            return("N_CORES undefined, automatic detection failed")
        
        # compile files
        logger.info("Compiling fortran files...")
        poolsize = min(len(deltaCompTasks), rp.N_CORES)
        with multiprocessing.Pool(poolsize) as pool:
            r = pool.map(compileDelta, deltaCompTasks)
        for v in [v for v in r if v != 0]:
            logger.error(v)
            return ("Fortran compile error")
        
        # run executions
        logger.info("Running delta calculations...")
        poolsize = min(len(deltaRunTasks), rp.N_CORES)
        with multiprocessing.Pool(poolsize) as pool:
            r = pool.map(runDelta, deltaRunTasks)
        for v in [v for v in r if v != 0]:
            logger.error(v)
            return ("Error during delta execution")
        logger.info("Delta calculations finished.")
        
        # clean up
        for ct in deltaCompTasks:
            try:
                shutil.rmtree(os.path.join(".", ct.foldername))
            except:
                logger.warning("Error deleting delta compile folder "
                                + ct.foldername)
        rp.manifest.append("Deltas")

    ####   SEARCH    ####
            
    elif index == 3:
        # read DISPLACEMENTS block
        if not rp.disp_block_read:
            tl.readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
            rp.disp_block_read = True
        rp.searchResultConfig = None
        # get Deltas
        if not 2 in rp.runHistory:
            if "Tensors" in rp.manifest:
                logger.error("New tensors were calculated, but no new delta "
                              "files were generated. Cannot execute search.")
                return ("Delta calculations was not run for current tensors.")
            try:
                r = getDeltas(rp.TENSOR_INDEX, required=True)
            except:
                raise
            if r != 0:
                return r
        # if number of cores is not defined, try to find it
        if rp.N_CORES == 0:
            try:
                rp.N_CORES = tl.available_cpu_count()
            except:
                logger.error("Failed to detect number of cores.")
            logger.info("Automatically detected number of available CPUs: {}"
                         .format(rp.N_CORES))
        if rp.N_CORES == 0:
            logger.error("Failed to detect number of cores.")
            return("N_CORES undefined, automatic detection failed")
        # generate rf.info
        try:
            rfinfo = tl.writeRfInfo(sl, rp, filename="rf.info")
        except:
            logger.error("Error generating search input file rf.info")
            raise
        # generate PARAM and search.steu
        #   needs to go AFTER rf.info, as writeRfInfo may remove expbeams!
        try:
            r = tl.generateSearchInput(sl, rp)
            if r != 0:
                logger.error("Error generating search input")
                return ("generateSearchInput failed")
        except:
            logger.error("Error generating search input")
            raise
        if rp.indyPars == 0:
            logger.info("Found nothing to vary in search. Will proceed "
                    "directly to writing output and starting SUPERPOS.")
            rp.searchResultConfig = [[1] * len(rp.searchpars)]
            for (i, sp) in enumerate(rp.searchpars):
                if type(sp.restrictTo) == int:
                    rp.searchResultsConfig[i] = sp.restrictTo
                elif type(sp.restrictTo) == tl.SearchPar:
                    rp.searchResultsConfig[i] = (rp.searchpars.index(
                                                        sp.restrictTo) + 1)
                elif type(sp.linkedTo) == tl.SearchPar:
                    rp.searchResultsConfig[i] = (rp.searchpars.index(
                                                        sp.linkedTo) + 1)
            tl.writeSearchOutput(sl, rp)
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
                
                
        if shutil.which("mpirun", os.X_OK) == None:
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
            tldir = tl.getTLEEDdir()
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
            r=tl.fortranCompile(fcomp[0]+" -o lib.search.o -c", 
                                libname, fcomp[1])
            if r:
                logger.error("Error compiling "+libname+", cancelling...")
                return ("Fortran compile error")
            r=tl.fortranCompile(fcomp[0]+" -o restrict.o -c", 
                                "restrict.f", fcomp[1])
            if r:
                logger.error("Error compiling restrict.f, cancelling...")
                return ("Fortran compile error")
            r=tl.fortranCompile(fcomp[0]+" -o search.o -c", srcname,
                                fcomp[1])
            if r:
                logger.error("Error compiling "+srcname+", cancelling...")
                return ("Fortran compile error")
            # combine
            r=tl.fortranCompile(fcomp[0]+" -o "+ searchname, "search.o "
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
                    # re-read PARAMETERS
                    tl.updatePARAMETERS_searchOnly(rp)
                    # check convergence criteria
                    stop = False
                    checkrepeat = True
                    if rp.SEARCH_KILL == True:
                        stop = True
                        checkrepeat = False
                        logger.info("Search stopped by SEARCH_KILL command.")
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
                            filepos, content = tl.readSDTL_next(
                                                             offset = filepos)
                            if content != "":
                                newData = tl.readSDTL_blocks(content, 
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
                                tl.writeSearchProgressPdf(rp, gens, rfaclist, 
                                                   lastconfig, markers=markers)
                            except:
                                logger.warning("Error writing "
                                                "Search-progress.pdf")
                            try:
                                tl.writeSearchReportPdf(rp)
                            except:
                                logger.warning("Error writing "
                                                "Search-report.pdf")
                        if (len(gens) > 1 and os.path.isfile("SD.TL") and 
                                                    (repeat or not stop)):
                            try:
                                r = tl.processSearchResults(sl, rp, 
                                                            final=False)
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
                rp.SEARCH_KILL = True
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
                    r = tl.generateSearchInput(sl, rp, steuOnly=True, 
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
                tl.writeSearchProgressPdf(rp, gens, rfaclist, lastconfig,
                                          markers=markers)
            except:
                logger.warning("Error writing Search-progress.pdf", 
                                exc_info = True)
            try:
                tl.writeSearchReportPdf(rp)
            except:
                logger.warning("Error writing Search-report.pdf",
                                exc_info = True)
        # process SD.TL to get POSCAR_OUT, VIBROCC_OUT
        r = tl.processSearchResults(sl, rp)
        if r != 0:
            logger.error("Error processing search results: "+r)
            return r
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



    ####   SUPERPOS    ####
            
    elif index == 31:
        # read DISPLACEMENTS block
        if not rp.disp_block_read:
            tl.readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
            rp.disp_block_read = True
        if not (2 in rp.runHistory or 3 in rp.runHistory):
            try:
                r = getDeltas(rp.TENSOR_INDEX, required=True)
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
    elapsedTimeStr = getElapsedTimeString(timer() - sectionStartTime)
    logger.info("Finishing section at " + time.strftime("%H:%M:%S", 
                                                         time.localtime())
                 + ". Section took " + elapsedTimeStr + ".")
    return 0

def sortfiles(tensorIndex, delete_unzipped = False, tensors = True, 
              deltas = True):
    """Makes Tensors and Deltas zip files. Copies files to AUX and OUT folders 
    as appropriate. If delete_unzipped is set to True, deletes unzipped Deltas 
    and Tensors directories."""
    # move files to AUX and OUT folders
    auxfiles = ["AUXBEAMS", "AUXGEO", "AUXLATGEO", "AUXNONSTRUCT", 
                "POSCAR_oricell", "POSCAR_bulk", "muftin.f", 
                "refcalc-PARAM", "refcalc-FIN", "rfactor-WEXPEL", 
                "rfactor-PARAM", "delta-input", "search.steu", 
                "search-rf.info", "seach-PARAM", "AUXEXPBEAMS", 
                "eeasisss-input", "searchpars.info", "superpos-PARAM",
                "superpos-CONTRIN"]
    outfiles = ["THEOBEAMS.csv", "THEOBEAMS_norm.csv", 
                "PatternInfo.tlm", "SD.TL", "refcalc-fd.out",
                "Rfactor_plots_refcalc.pdf", "control.chem", 
                "Search-progress.pdf", "Search-progress.csv", 
                "Search-report.pdf", "FITBEAMS.csv", "FITBEAMS_norm.csv", 
                "superpos-spec.out", "Rfactor_plots_superpos.pdf"]
    # outfiles with variable names:
    outfiles.extend([f for f in os.listdir(".") if 
                         (f.startswith("POSCAR_OUT") or
                          f.startswith("VIBROCC_OUT") or
                          f.startswith("R_OUT"))])
    # clean up deltas
    deltalist = [f for f in os.listdir('.') if f.startswith("DEL_")]
    if len(deltalist) > 0:
        fn = "Deltas_"+str(tensorIndex).zfill(3)
        if not os.path.isdir(os.path.join(".","Deltas")):
            os.mkdir(os.path.join(".","Deltas"))
        if not os.path.isdir(os.path.join(".","Deltas",fn)):
            os.mkdir(os.path.join(".","Deltas",fn))
        try:
            for df in deltalist:
                shutil.move(df, os.path.join(".","Deltas",fn,df))
        except:
            logger.error("Error moving Delta files: ", exc_info = True)
    # if there are unzipped Tensors or Deltas directories, zip them:
    if os.path.isdir(os.path.join(".","Tensors")) and (tensors or 
                                                       delete_unzipped):
        rgx = re.compile(r'Tensors_[0-9]{3}')
        for d in [d for d in os.listdir(os.path.join(".","Tensors")) 
                  if (os.path.isdir(os.path.join(".","Tensors",d))
                      and rgx.match(d))]:
            if not rgx.match(d).span()[1] == 11:
                continue
            delete = delete_unzipped
            if tensors:
                logger.info("Packing {}.zip...".format(d))
                try:
                    shutil.make_archive(os.path.join("Tensors",d),"zip",
                                        os.path.join("Tensors",d))
                except:
                    logger.error("Error packing {}.zip file: ".format(d))
                    delete = False
            if delete:
                try:
                    shutil.rmtree(os.path.join(".","Tensors",d))
                except:
                    logger.warning("Error deleting unzipped Tensors "
                        "directory. This will increase the size of the work "
                        "folder, but not cause any problems.")
    if os.path.isdir(os.path.join(".","Deltas")) and (deltas or 
                                                      delete_unzipped):
        rgx = re.compile(r'Deltas_[0-9]{3}')
        for d in [d for d in os.listdir(os.path.join(".","Deltas"))
                  if (os.path.isdir(os.path.join(".","Deltas",d))
                      and rgx.match(d))]:
            if not rgx.match(d).span()[1] == 10:
                continue
            delete = delete_unzipped
            if deltas:
                logger.info("Packing {}.zip...".format(d))
                try:
                    shutil.make_archive(os.path.join("Deltas",d),"zip",
                                        os.path.join("Deltas",d))
                except:
                    logger.error("Error packing {}.zip file.".format(d))
                    delete = False
            if delete:
                try:
                    shutil.rmtree(os.path.join(".","Deltas",d))
                except:
                    logger.warning("Error deleting unzipped Deltas "
                        "directory. This will increase the size of the work "
                        "folder, but not cause any problems.")
    # sort AUX and OUT files:
    try:
        if not os.path.isdir(os.path.join(".","AUX")):
            os.mkdir(os.path.join(".","AUX"))
        if not os.path.isdir(os.path.join(".","OUT")):
            os.mkdir(os.path.join(".","OUT"))
    except:
        logger.error("Error creating AUX and OUT folders: ", exc_info = True)
    for f in auxfiles:
        if os.path.isfile(os.path.join(".",f)):
            try:
                shutil.copy2(f, os.path.join(".","AUX",f))
            except:
                logger.error("Error moving AUX file "+f+": ", exc_info = True)
    for f in outfiles:
        if os.path.isfile(os.path.join(".",f)):
            try:
                shutil.copy2(f, os.path.join(".","OUT",f))
            except:
                logger.error("Error copying OUT file "+f+": ", 
                              exc_info = True)


###############################################
#            CLEANUP FUNCTIONS                #
###############################################
    
def moveoldruns(rp, prerun = False):
    """Makes a new folder in 'workhistory'. Copies AUX, OUT and files in 
    manifest (except main log) to that new folder. If prerun is set True, then 
    instead of using the manifest, all potentially interesting files will be 
    copied, and the new folder will get index 0."""
    sectionabbrv = {1: "R", 2: "D", 3: "S"}
    if not os.path.isdir(os.path.join(".","workhistory")):
        try:
            os.mkdir(os.path.join(".","workhistory"))
        except:
            logger.error("Error creating workhistory folder: ", 
                          exc_info = True)
            return 1
    if not prerun:
        rp.manifest.append("workhistory")
    dl = [n for n in os.listdir("workhistory") 
                      if os.path.isdir(os.path.join("workhistory",n))]
    maxnum = -1
    rgx = re.compile(r't'+'{:03d}'.format(rp.TENSOR_INDEX)+r'.r[0-9]{3}_')
    for d in dl:
        m = rgx.match(d)
        if m:
            try:
                i = int(d[6:9])
                if i > maxnum:
                    maxnum = i
            except:
                pass
    if maxnum == -1:
        if prerun:
            num = 0
        else:
            num = 1
    else:
        num = maxnum + 1
    if prerun:
        oldlogfiles = sorted([f for f in os.listdir() if os.path.isfile(f) and
                              f.endswith(".log") and f.startswith("tleedm-")
                              and not f in rp.manifest])
        if len(oldlogfiles) > 0:
            oldTimeStamp = oldlogfiles[-1][7:20]
        else:
            oldTimeStamp = "moved-" + rp.timestamp
        dirname = ("t{:03d}.r{:03d}_previous_".format(rp.TENSOR_INDEX, num) 
                   + oldTimeStamp)
    else:
        dirname = "t{:03d}.r{:03d}_".format(rp.TENSOR_INDEX, num)
        for ind in rp.runHistory[len(rp.lastOldruns):]:
            if ind in sectionabbrv:
                dirname += sectionabbrv[ind]
        rp.lastOldruns = rp.runHistory[:]
        dirname += "_" + rp.timestamp
    dirpath = os.path.join(".","workhistory",dirname)
    try:
        os.mkdir(dirpath)
    except:
        logger.error("Error creating workhistory subfolder: ", 
                      exc_info = True)
        return 1
    if not prerun:
        sortfiles(rp.TENSOR_INDEX, delete_unzipped=False, 
                  tensors = False, deltas = False)
    if prerun:
        filelist = [f for f in os.listdir() if os.path.isfile(f) and 
                    f.endswith(".log") and f not in rp.manifest]
        dirlist = ["AUX", "OUT"]
    else:
        filelist = [f for f in rp.manifest if os.path.isfile(f) and not
                    (f.startswith("tleedm-") and f.endswith(".log"))]
        dirlist = [d for d in rp.manifest if os.path.isdir(d) and not
                   d in ["Tensors", "Deltas", "workhistory"]]
    for f in filelist:
        try:
            if not prerun:
                shutil.copy2(f, os.path.join(dirpath, f))
            else:
                shutil.move(f, os.path.join(dirpath, f))
        except:
            logger.warning("Error copying "+f+" to " 
                            + os.path.join(dirpath, f)
                            + ". File may get overwritten.")
    for d in dirlist:
        try:
            shutil.copytree(d, os.path.join(dirpath, d))
        except:
            logger.warning("Error copying "+d+" to " 
                            + os.path.join(dirpath, d)
                            + ". Files in directory may get overwritten.")
    return 0

def getElapsedTimeString(t):
    """Takes a time in seconds, returns a string giving the elapsed times 
    either in minutes or hours, as appropriate."""
    if t >= 3600:
        elapsedTimeStr = (str(int(t/3600)) + ":"+(str(int(t/60)%60)).zfill(2)
                          +" hours")
    elif t >= 60:
        elapsedTimeStr = (str(int(t/60)) + ":"+(str(int(t)%60)).zfill(2) 
                          +" minutes")
    else:
        elapsedTimeStr = str(round(t, 2)) +" seconds"
    return elapsedTimeStr

def cleanup(manifest, rp = None):
    """Moves files to AUX and OUT folders, writes manifest, adds a final 
    message to the log, then shuts down everything."""
    global starttime
    
    logger.info("\nStarting cleanup...")
    if rp is None:
        history = []
        newTensors = False
        newDeltas = False
        tind = 0
    else:
        history = rp.runHistory
        newTensors = ("Tensors" in rp.manifest)
        newDeltas = ("Deltas" in rp.manifest)
        tl.closePdfReportFigs(rp)
        tind = rp.TENSOR_INDEX

    try:
        sortfiles(tind, delete_unzipped=True, tensors=newTensors, 
                                                         deltas=newDeltas)
    except:
        logger.warning("Error sorting files to AUX/OUT folders: ", 
                        exc_info = True)
    # write manifest
    written = []
    try:
        with open("manifest","w") as wf:
            for fn in manifest:
                if not fn in written:
                    wf.write(fn+"\n")
                    written.append(fn)
        logger.info("Wrote manifest file successfully.")
    except:
        logger.error("Failed to write manifest file.")
    
    # write final log message 
    elapsedTimeStr = getElapsedTimeString(timer() - starttime)
    logger.info("\nFinishing execution at "+time.strftime("%Y-%m-%d %H:%M:%S", 
                                                           time.localtime())
                 +"\nTotal elapsed time: "+elapsedTimeStr+"\n")
    if len(history) > 0:
        s = ""
        for ind in history:
            s += (str(ind)+" ")
        logger.info("Executed segments: "+s[:-1]+"\n")
    # shut down logger
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])
    logging.shutdown()


###############################################
#                  MAIN                       #
###############################################

def main():
    global starttime
    #start logger, write to file:
    timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
    logname = 'tleedm-'+timestamp+'.log'

    logger = logging.getLogger("tleedm")
    logger.setLevel(logging.INFO)
    logFormatter = CustomLogFormatter()

    fileHandler = logging.FileHandler(logname, mode="w")
    fileHandler.setFormatter(logFormatter)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)

    logger.addHandler(consoleHandler)
    logger.addHandler(fileHandler)

       
    logger.info("Starting new log: "+logname+"\nTime of execution (UTC): "
                 +time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+"\n")
    starttime = timer()
    
    tmpmanifest = ["AUX","OUT",logname]
    
    # read POSCAR and PARAMETERS first; they need to be there for everything
    #   anyway, and PARAMETERS contains information on what to run

    if os.path.isfile(os.path.join(".","POSCAR")):
        poscarfile = os.path.join(".","POSCAR")
        logger.info("Reading structure from file POSCAR")
        try:
            sl = tl.readPOSCAR(filename = poscarfile)
        except:
            logger.error("Exception while reading POSCAR", exc_info=True)
            cleanup(tmpmanifest)
            return 2
    else:
        logger.error("POSCAR not found. Stopping execution...")
        cleanup(tmpmanifest)
        return 1
    
    if not sl.preprocessed:
        logger.info("The POSCAR file will be processed and overwritten. "
                     "Copying the original POSCAR to POSCAR_user...")
        try:
            shutil.copy2(poscarfile, "POSCAR_user")
            tmpmanifest.append("POSCAR_user")
        except:
            logger.error("Failed to copy POSCAR to POSCAR_user. Stopping "
                          "execution...")
            cleanup(tmpmanifest)
            return 1

    try:
        rp = tl.readPARAMETERS(slab=sl)
        if rp.LOG_DEBUG:
            logger.setLevel(logging.DEBUG)
            logger.debug("PARAMETERS file was read successfully")
    except:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmpmanifest)
        return 2
    
    rp.timestamp = timestamp
    rp.manifest = tmpmanifest
    sl.fullUpdate(rp)   #gets PARAMETERS data into slab
    rp.fileLoaded["POSCAR"] = True
    if sl.preprocessed:
        rp.SYMMETRY_FIND_ORI = False
    
    # get name of parent folder, use as system name (for file headers):
    rp.systemName = os.path.basename(os.path.abspath(os.path.join(os.getcwd(),
                                                                  os.pardir)))
    # check if halting condition is already in effect:
    if rp.halt >= rp.HALTING:
        logger.info("Halting execution...")
        cleanup(rp.manifest, rp)
        return 0
   
    rp.updateDerivedParams()
    # clean out the workhistory folder, if there is one
    if os.path.isdir(os.path.join(".","workhistory")):
        try:
            shutil.rmtree(os.path.join(".","workhistory"))
        except:
            logger.warning("Failed to clear workhistory folder.")
    # see if there are old logfiles
    oldlogs = [f for f in os.listdir() if os.path.isfile(f) and 
               f.endswith(".log") and f != logname]
    if len(oldlogs) > 0:
        try:
            moveoldruns(rp, prerun = True)
        except:
            logger.warning("Exception while trying to clean up from previous "
                            "run. Program will proceed, but old files may be "
                            "lost.", exc_info=True)
    if os.path.isfile("fortran-compile.log"):
        try:
            os.remove("fortran-compile.log")
        except:
            pass
        
    sectionorder = [0, 1, 11, 2, 3, 31, 12]
    # searchLoopRfacs = []
    searchLoopR = None
    searchLoopLevel = 0
    initHalt = False
    while len(rp.RUN) > 0:
        try:
            sec = rp.RUN.pop(0)
            if rp.runHistory and (sectionorder.index(sec) 
                                  < sectionorder.index(rp.runHistory[-1])):
                logger.info("\nExecution repeats. Moving old output to "
                             "workhistory folder.")
                moveoldruns(rp)
            r = runSection(sec, sl, rp)
            if r != 0:
                logger.error("Error in tleedm execution: "+str(r))
                cleanup(rp.manifest, rp)
                return 1
            elif (sec == 0 and not sl.preprocessed and rp.HALTING <= 2 
                  and len(rp.RUN) > 0):
                logger.info("Initialization finished. Execution will stop. "
                    "Please check whether comments in POSCAR are correct, "
                    "then restart.")
                rp.setHaltingLevel(2)
                initHalt = True
            elif (sec == 1 and rp.fileLoaded["EXPBEAMS"]):
                if rp.RUN[:1] != [11]:   # r-factor after refcalc
                    rp.RUN.insert(0, 11)
            elif (sec == 3 and rp.fileLoaded["EXPBEAMS"]):
                if rp.RUN[:1] != [31]:  # superpos after search
                    rp.RUN.insert(0, 31)
            elif sec == 31 and rp.fileLoaded["EXPBEAMS"]:
                if rp.RUN[:1] != [12]:   # r-factor after superpos
                    rp.RUN.insert(0, 12)
            elif sec == 12 and not rp.SEARCH_KILL:
                loops = [t for t in rp.disp_loops if t[1] == rp.search_index]
                if loops:
                    if searchLoopLevel == 0 or searchLoopR > rp.last_R:
                        searchLoopR = rp.last_R
                        searchLoopLevel = len(loops)
                    elif searchLoopR <= rp.last_R:
                        searchLoopLevel -= 1
                    if searchLoopLevel != 0:
                        rp.search_index = sorted(loops)[searchLoopLevel-1][0]
                        logger.info("Search loop: repeating at block "
                                     +rp.disp_blocks[rp.search_index][1])
                    else:
                        rp.search_index += 1
                        o = "Search loop ends."
                        if len(rp.disp_blocks) > rp.search_index:
                            o += (" Continuing at block "
                                  + rp.disp_blocks[rp.search_index][1])
                        logger.info(o)
                else:
                    rp.search_index += 1
                if len(rp.disp_blocks) > rp.search_index:
                    sl.restoreOriState()
                    rp.resetSearchConv()
                    if rp.SEARCH_START == "control":
                        rp.SEARCH_START = "crandom"
                    if rp.RUN[:2] != [2,3]:
                        rp.RUN = [2,3] + rp.RUN
        except KeyboardInterrupt:
            logger.warning("Stopped by keyboard interrupt, attempting "
                            "clean exit...")
            cleanup(rp.manifest, rp)
            return 2
        except:
            logger.error("Exception during tleedm execution: ", exc_info=True)
            cleanup(rp.manifest, rp)
            return 2
        if rp.halt >= rp.HALTING:
            if not initHalt:
                logger.info("# An exception occured that meets the halting "
                    "criteria defined by the HALTING parameter. Execution "
                    "will stop, check log for warnings and errors.")
            break
    cleanup(rp.manifest, rp)
    return 0

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()