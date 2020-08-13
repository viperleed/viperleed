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
import re
import multiprocessing

cd = os.path.realpath(os.path.dirname(__file__))  # !!! check this.
tleedmap_path = os.path.realpath(os.path.join(cd, '..'))
if tleedmap_path not in sys.path:
    sys.path.append(tleedmap_path)

import tleedmlib as tl
import tleedmlib.sections as sections

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
    if index == 3 and rp.disp_blocks and rp.disp_blocks[rp.search_index][1]:
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

    r = 0
    if index == 0:
        sections.initialization(sl, rp)  
    elif index == 1:
        r = sections.refcalc(sl, rp)
    elif index == 11 or index == 12:
        r = sections.rfactor(sl, rp, index)
    elif index == 2:
        r = sections.deltas(sl, rp)
    elif index == 3:
        r = sections.search(sl, rp)
    elif index == 31:
        r = sections.superpos(sl, rp)
    if r != 0:
        return r

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
                "superpos-spec.out", "Rfactor_plots_superpos.pdf",
                "Rfactor_analysis_refcalc.pdf", 
                "Rfactor_analysis_superpos.pdf"]
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
        rp.closePdfReportFigs()
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
        logger.info("Executed segments: "+s[:-1])
    if rp is not None:
        for t in ["refcalc", "superpos"]:
            if rp.stored_R[t] is not None:
                o = "Final R ({}): {:.4f}".format(t, rp.stored_R[t][0])
                if rp.stored_R[t][1] > 0 and rp.stored_R[t][2] > 0:
                    o += " ({:.4f} / {:.4f})".format(rp.stored_R[t][1], 
                                                     rp.stored_R[t][2])
                logger.info(o)
    logger.info("")
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