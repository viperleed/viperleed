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


cd = os.path.realpath(os.path.dirname(__file__))
# NB: it's necessary to add vpr_path to sys.path so that viperleed
#     can be loaded correctly at the top-level package
vpr_path = os.path.realpath(os.path.join(cd, '..'))
for import_path in (cd, vpr_path):
    if import_path not in sys.path:
        sys.path.append(import_path)


import viperleed.tleedmlib.sections as sections
from viperleed.vprglobals import GLOBALS
from viperleed.tleedmlib.files.parameters import (
    readPARAMETERS, interpretPARAMETERS, modifyPARAMETERS, updatePARAMETERS)
from viperleed.tleedmlib.files.phaseshifts import readPHASESHIFTS
from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.files.beams import (
    readBEAMLIST, readIVBEAMS, readOUTBEAMS, checkEXPBEAMS)
from viperleed.tleedmlib.files.vibrocc import readVIBROCC, writeVIBROCC
from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS

logger = logging.getLogger("tleedm")
starttime = 0.0

auxfiles = ["AUXBEAMS", "AUXGEO", "AUXLATGEO", "AUXNONSTRUCT",
            "POSCAR_oricell", "POSCAR_bulk", "muftin.f",
            "refcalc-PARAM", "refcalc-FIN", "rfactor-WEXPEL",
            "rfactor-PARAM", "delta-input", "search.steu",
            "search-rf.info", "search-PARAM", "AUXEXPBEAMS",
            "eeasisss-input", "searchpars.info", "superpos-PARAM",
            "superpos-CONTRIN", "POSCAR_bulk_appended", "POSCAR_mincell",
            "restrict.f"]
outfiles = ["THEOBEAMS.csv", "THEOBEAMS_norm.csv",
            "PatternInfo.tlm", "SD.TL", "refcalc-fd.out",
            "Rfactor_plots_refcalc.pdf", "control.chem",
            "Search-progress.pdf", "Search-progress.csv",
            "Search-report.pdf", "FITBEAMS.csv", "FITBEAMS_norm.csv",
            "superpos-spec.out", "Rfactor_plots_superpos.pdf",
            "Rfactor_analysis_refcalc.pdf",
            "Rfactor_analysis_superpos.pdf", "Errors.csv", "Errors.pdf"]


class CustomLogFormatter(logging.Formatter):
    """Logging Formatter for level-dependent message formatting"""
    FORMATS = {
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
    """
    Runs a specific part of the program.

    Parameters
    ----------
    index : int
        Index of the section, see sectionNames variable below.
    sl : Slab
        Slab object containing atom information.
    rp : Rparams
        The run parameters.

    Returns
    -------
    None.

    """
    sectionNames = {0: "INITIALIZATION",
                    1: "REFERENCE CALCULATION",
                    2: "DELTA-AMPLITUDES",
                    3: "SEARCH",
                    11: "R-FACTOR CALCULATION",
                    12: "R-FACTOR CALCULATION",
                    31: "SUPERPOS",
                    5:  "ERROR CALCULATION"}
    # files that need to be there for the different parts to run
    requiredFiles = {0: ["POSCAR", "PARAMETERS", "VIBROCC", "IVBEAMS"],
                     1: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC"],
                     2: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS"],
                     3: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS", "EXPBEAMS"],
                     11: ["BEAMLIST", "PARAMETERS", "IVBEAMS", "EXPBEAMS"],
                     12: ["BEAMLIST", "PARAMETERS", "IVBEAMS", "EXPBEAMS"],
                     31: ["BEAMLIST", "POSCAR", "PARAMETERS", "IVBEAMS",
                          "VIBROCC", "DISPLACEMENTS"],
                     5: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS", "EXPBEAMS"]}

    checkfiles = requiredFiles[index][:]
    o = "\nSTARTING SECTION: "+sectionNames[index]
    if index == 3 and rp.disp_blocks and rp.disp_blocks[rp.search_index][1]:
        o += " "+rp.disp_blocks[rp.search_index][1]  # displacement block name
    if rp.domainParams or rp.DOMAINS:
        o += " (DOMAINS)"
        for fn in ["POSCAR", "VIBROCC", "PHASESHIFTS"]:
            try:
                checkfiles.remove(fn)
            except Exception:
                pass
    logger.info(o)
    sectionStartTime = timer()
    rp.runHistory.append(index)
    for dp in rp.domainParams:
        dp.rp.runHistory = rp.runHistory
    i = 0
    while i < len(checkfiles):
        filename = checkfiles[i]
        ignoreError = False
        if not rp.fileLoaded[filename]:
            # try loading files
            if filename == "EXPBEAMS":
                if len(rp.THEO_ENERGIES) == 0:
                    er = []
                else:
                    er = rp.THEO_ENERGIES[:2]
                try:
                    rp.expbeams = readOUTBEAMS("EXPBEAMS.csv", enrange=er)
                    rp.fileLoaded["EXPBEAMS"] = True
                except FileNotFoundError:
                    try:
                        # try without the .csv extension
                        rp.expbeams = readOUTBEAMS("EXPBEAMS", enrange=er)
                        rp.fileLoaded["EXPBEAMS"] = True
                    except Exception:
                        logger.error("Error while reading required file "
                                     "EXPBEAMS.csv")
                except Exception as e:
                    logger.error("Error while reading required file EXPBEAMS",
                                 exc_info=(type(e) != FileNotFoundError))
                if index != 0:
                    checkEXPBEAMS(sl, rp)
            elif filename == "IVBEAMS":
                try:
                    rp.ivbeams = readIVBEAMS()
                    rp.ivbeams_sorted = False
                    rp.fileLoaded["IVBEAMS"] = True
                except FileNotFoundError:
                    if (os.path.isfile("EXPBEAMS")
                            or os.path.isfile("EXPBEAMS.csv")):
                        checkfiles.insert(i+1, "EXPBEAMS")
                        logger.warning("IVBEAMS file not found. Will attempt "
                                       "generating IVBEAMS from EXBEAMS.")
                        ignoreError = True
                    else:
                        logger.error("Neither IVBEAMS not EXPBEAMS file "
                                     "found.")
                except Exception as e:
                    logger.error("Error while reading required file IVBEAMS",
                                 exc_info=(type(e) != FileNotFoundError))
            elif filename == "BEAMLIST":
                try:
                    rp.beamlist = readBEAMLIST()
                    rp.fileLoaded["BEAMLIST"] = True
                except Exception as e:
                    logger.error("Error while reading required file "
                                 "_BEAMLIST", exc_info=(type(e) !=
                                                        FileNotFoundError))
            elif filename == "VIBROCC":
                changeVIBROCC = False
                try:
                    changeVIBROCC = readVIBROCC(rp, sl)
                    rp.fileLoaded["VIBROCC"] = True
                except Exception as e:
                    logger.error("Error while reading required file VIBROCC",
                                 exc_info=(type(e) != FileNotFoundError))
                sl.fullUpdate(rp)
                if changeVIBROCC:
                    if os.path.isfile("VIBROCC"):
                        os.rename("VIBROCC", "VIBROCC_user")
                        rp.manifest.append("VIBROCC_user")
                        logger.info(
                            "VIBROCC file was modified with automatically "
                            "generated vibrational amplitudes.")
                    writeVIBROCC(sl, rp, "VIBROCC")
                    rp.manifest.append("VIBROCC")
                if rp.T_EXPERIMENT is not None:
                    modifyPARAMETERS(rp, "T_EXPERIMENT", new="")
                if rp.T_DEBYE is not None:
                    modifyPARAMETERS(rp, "T_DEBYE", new="")
                if len(rp.VIBR_AMP_SCALE) > 0:
                    modifyPARAMETERS(rp, "VIBR_AMP_SCALE", new="")
            elif filename == "PHASESHIFTS":
                try:
                    (rp.phaseshifts_firstline, rp.phaseshifts,
                     newpsGen, newpsWrite) = readPHASESHIFTS(sl, rp)
                    if newpsGen:
                        logger.error(
                            "_PHASESHIFT file generation is only supported "
                            "during initialization. Stopping execution...")
                        raise RuntimeError("Inconsistent _PHASESHIFT file")
                    elif newpsWrite:
                        logger.warning(
                            "Writing a new _PHASESHIFT file is "
                            "only supported during initialization. The "
                            "data in the provided file will be used, but "
                            "running the initialization is recommended.")
                        rp.fileLoaded["PHASESHIFTS"] = True
                    else:
                        rp.fileLoaded["PHASESHIFTS"] = True
                except Exception as e:
                    logger.error("Error while reading required file "
                                 "_PHASESHIFTS", exc_info=(type(e) !=
                                                           FileNotFoundError))
            elif filename == "DISPLACEMENTS":
                try:
                    readDISPLACEMENTS(rp)
                    rp.fileLoaded["DISPLACEMENTS"] = True
                except Exception as e:
                    logger.error("Error while reading required file "
                                 "DISPLACEMENTS", exc_info=(type(e) !=
                                                            FileNotFoundError))
            if not rp.fileLoaded[filename] and not ignoreError:
                # and if that didn't work, stop:
                logger.error("Step '" + sectionNames[index] + "' requires "
                             "file " + filename + ". Stopping execution...")
                rp.setHaltingLevel(3)
                return
        i += 1
    try:
        if index == 0:
            sections.initialization(sl, rp)
        elif index == 1:
            sections.refcalc(sl, rp)
        elif index in [11, 12]:
            sections.rfactor(sl, rp, index)
        elif index == 2:
            sections.deltas(sl, rp)
        elif index == 3:
            sections.search(sl, rp)
        elif index == 31:
            sections.superpos(sl, rp)
        elif index == 5:
            sections.errorcalc(sl, rp)
    except Exception:
        logger.error("Error in section {}".format(sectionNames[index]))
        raise
    elapsedTimeStr = getElapsedTimeString(timer() - sectionStartTime)
    logger.info("Finishing section at " + time.strftime("%H:%M:%S",
                                                        time.localtime())
                + ". Section took " + elapsedTimeStr + ".")
    return


def sortfiles(tensorIndex, delete_unzipped=False, tensors=True,
              deltas=True, path=""):
    """
    Makes Tensors and Deltas zip files. Copies files to AUX and OUT folders
    as appropriate. If delete_unzipped is set to True, deletes unzipped Deltas
    and Tensors directories.

    Parameters
    ----------
    tensorIndex : int
        Which Delta and Tensor files should be considered.
    delete_unzipped : bool, optional
        Whether the original Delta- and Tensor-files should be deleted
        after making the archives. The default is False.
    tensors, deltas : bool, optional
        Whether the Tensor/Delta files contain new information and should be
        saved. The default is True.
    path : str, optional
        The base path to check files in. The default is "".

    Returns
    -------
    None.

    """
    # move files to AUX and OUT folders

    # outfiles with variable names:
    if not path:
        path = "."
    outfiles.extend([f for f in os.listdir(path)
                     if (f.startswith("POSCAR_OUT") or
                         f.startswith("VIBROCC_OUT") or
                         f.startswith("R_OUT"))])
    # clean up deltas
    deltalist = [f for f in os.listdir(path) if f.startswith("DEL_")]
    if len(deltalist) > 0:
        fn = "Deltas_"+str(tensorIndex).zfill(3)
        os.makedirs(os.path.join(path, "Deltas", fn), exist_ok=True)
        try:
            for df in deltalist:
                shutil.move(os.path.join(path, df),
                            os.path.join(path, "Deltas", fn, df))
        except Exception:
            logger.error("Error moving Delta files: ", exc_info=True)

    # if there are unzipped Tensors or Deltas directories, zip them:
    for t in ["Tensors", "Deltas"]:
        if t == "Tensors":
            do = tensors
        else:
            do = deltas
        rgx = re.compile(t+r'_[0-9]{3}')
        if not os.path.isdir(os.path.join(path, t)):
            continue
        if not (do or delete_unzipped):
            continue
        for d in [d for d in os.listdir(os.path.join(path, t))
                  if (os.path.isdir(os.path.join(path, t, d))
                      and rgx.match(d))]:
            if not rgx.match(d).span()[1] == len(t)+4:
                continue
            delete = delete_unzipped
            if do:
                if not path:
                    o = d
                else:
                    o = os.path.relpath(os.path.join(path, t, d))
                logger.info("Packing {}.zip...".format(o))
                try:
                    shutil.make_archive(os.path.join(path, t, d), "zip",
                                        os.path.join(path, t, d))
                except Exception:
                    logger.error("Error packing {}.zip file: ".format(o))
                    delete = False
            if delete:
                try:
                    shutil.rmtree(os.path.join(path, t, d))
                except Exception:
                    logger.warning(
                        "Error deleting unzipped {} directory. "
                        "This will increase the size of the work folder, "
                        "but not cause any problems.".format(t))
    # sort AUX and OUT files:
    for t in ["AUX", "OUT"]:
        try:
            os.makedirs(os.path.join(path, t), exist_ok=True)
        except Exception:
            logger.error("Error creating {} folder: ".format(t), exc_info=True)
        if t == "AUX":
            filelist = auxfiles
        else:
            filelist = outfiles
        for f in [f for f in filelist
                  if os.path.isfile(os.path.join(path, f))]:
            try:
                shutil.copy2(os.path.join(path, f), os.path.join(path, t, f))
            except Exception:
                logger.error("Error moving {} file {}: ".format(t, f),
                             exc_info=True)


###############################################
#            CLEANUP FUNCTIONS                #
###############################################

def moveoldruns(rp, prerun=False):
    """
    Makes a new folder in 'workhistory'. Copies AUX, OUT and files in
    manifest (except main log) to that new folder.

    Parameters
    ----------
    rp : Rparams
        The run parameters.
    prerun : bool, optional
        If True, then instead of using the manifest, all potentially
        interesting files will be copied, and the new folder will get index 0.
        Then clears out the AUX and OUT folders and old AUX and OUT files from
        the main directory.

    Returns
    -------
    None

    """
    sectionabbrv = {1: "R", 2: "D", 3: "S"}
    try:
        os.makedirs(os.path.join(".", "workhistory"), exist_ok=True)
    except Exception:
        logger.error("Error creating workhistory folder: ", exc_info=True)
        raise
    if not prerun:
        rp.manifest.append("workhistory")
    dl = [n for n in os.listdir("workhistory")
          if os.path.isdir(os.path.join("workhistory", n))]
    maxnum = -1
    rgx = re.compile(r't'+'{:03d}'.format(rp.TENSOR_INDEX)+r'.r[0-9]{3}_')
    for d in dl:
        m = rgx.match(d)
        if m:
            try:
                i = int(d[6:9])
                if i > maxnum:
                    maxnum = i
            except Exception:
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
                              and f not in rp.manifest])
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
    dirpath = os.path.join(".", "workhistory", dirname)
    try:
        os.mkdir(dirpath)
    except Exception:
        logger.error("Error creating workhistory subfolder: ", exc_info=True)
        raise
    if not prerun:
        sortfiles(rp.TENSOR_INDEX, delete_unzipped=False,
                  tensors=False, deltas=False)
        for dp in rp.domainParams:
            sortfiles(dp.rp.TENSOR_INDEX, delete_unzipped=False,
                      tensors=False, deltas=False)
    if prerun:
        filelist = [f for f in os.listdir() if os.path.isfile(f) and
                    (f.endswith(".log") or f in outfiles or f in auxfiles)
                    and f not in rp.manifest]
        dirlist = ["AUX", "OUT"]
    else:
        filelist = [f for f in rp.manifest if os.path.isfile(f) and not
                    (f.startswith("tleedm-") and f.endswith(".log"))]
        dirlist = [d for d in rp.manifest if os.path.isdir(d) and
                   d not in ["Tensors", "Deltas", "workhistory"]]
    for f in filelist:
        try:
            if not prerun:
                shutil.copy2(f, os.path.join(dirpath, f))
            else:
                shutil.move(f, os.path.join(dirpath, f))
        except Exception:
            logger.warning("Error copying "+f+" to "
                           + os.path.join(dirpath, f)
                           + ". File may get overwritten.")
    for d in dirlist:
        try:
            if not prerun:
                shutil.copytree(d, os.path.join(dirpath, d))
            else:
                shutil.move(d, os.path.join(dirpath, d))
        except Exception:
            logger.warning("Error copying "+d+" to "
                           + os.path.join(dirpath, d)
                           + ". Files in directory may get overwritten.")
    return


def getElapsedTimeString(t):
    """
    Takes a time in seconds, returns a string giving the elapsed times
    either in minutes or hours, as appropriate.
    """
    if t >= 3600:
        elapsedTimeStr = (str(int(t/3600)) + ":"+(str(int(t/60) % 60)).zfill(2)
                          + " hours")
    elif t >= 60:
        elapsedTimeStr = (str(int(t/60)) + ":"+(str(int(t) % 60)).zfill(2)
                          + " minutes")
    else:
        elapsedTimeStr = str(round(t, 2)) + " seconds"
    return elapsedTimeStr


def cleanup(manifest, rp=None):
    """
    Moves files to AUX and OUT folders, writes manifest, adds a final
    message to the log, then shuts down everything.

    Parameters
    ----------
    manifest : list of str
        The files and directories that should be preserved from the work
        folder.
    rp : Rparams, optional
        The run parameters. If None, it is assumed that the run crashed before
        an Rparams object existed.

    Returns
    -------
    None.

    """
    global starttime

    logger.info("\nStarting cleanup...")
    if rp is None:
        history = []
        to_sort = [{"newTensors": False, "newDeltas": False, "tind": 0,
                    "path": ""}]
    else:
        history = rp.runHistory
        rp.closePdfReportFigs()
        if not rp.domainParams:
            to_sort = [{"newTensors": ("Tensors" in rp.manifest),
                        "newDeltas": ("Deltas" in rp.manifest),
                        "tind": rp.TENSOR_INDEX, "path": ""}]
        else:
            to_sort = [{"newTensors": False, "newDeltas": False, "tind": 0,
                        "path": ""}]
            for dp in rp.domainParams:
                to_sort.append({"newTensors": ("Tensors" in dp.rp.manifest),
                                "newDeltas": ("Deltas" in dp.rp.manifest),
                                "tind": dp.rp.TENSOR_INDEX,
                                "path": dp.workdir})
    for d in to_sort:
        try:
            sortfiles(d["tind"], delete_unzipped=True,
                      tensors=d["newTensors"],
                      deltas=d["newDeltas"], path=d["path"])
        except Exception:
            logger.warning("Error sorting files to AUX/OUT folders: ",
                           exc_info=True)
    # write manifest
    written = []
    try:
        with open("manifest", "w") as wf:
            for fn in manifest:
                if fn not in written:
                    wf.write(fn+"\n")
                    written.append(fn)
        logger.info("Wrote manifest file successfully.")
    except Exception:
        logger.error("Failed to write manifest file.")

    # write final log message
    elapsedTimeStr = getElapsedTimeString(timer() - starttime)
    logger.info("\nFinishing execution at "+time.strftime("%Y-%m-%d %H:%M:%S",
                                                          time.localtime())
                + "\nTotal elapsed time: "+elapsedTimeStr+"\n")
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
    if rp:
        if rp.checklist:
            logger.info("")
            logger.info("# The following issues should be checked before "
                        "starting again:")
            for s in rp.checklist:
                logger.info("- "+s)
    logger.info("")
    # shut down logger
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])
    logging.shutdown()
    return


###############################################
#                  MAIN                       #
###############################################

def main():
    """
    Main of TensErLEED Manager.

    Returns
    -------
    int
        0: exit without errors.
        1: clean exit through KeyboardInterrupt
        2: exit due to Exception before entering main loop
        3: exit due to Exception during main loop

    """
    global starttime
    os.umask(0)
    # start logger, write to file:
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

    logger.info("Starting new log: " + logname + "\nTime of execution (UTC): "
                + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    logger.info("This is ViPErLEED version " + GLOBALS["version"] + "\n")
    logger.info("! THIS VERSION IS A PRE-RELEASE NOT MEANT FOR PUBLIC "
                "DISTRIBUTION !\n")
    starttime = timer()

    tmpmanifest = ["AUX", "OUT", logname]
    try:
        rp = readPARAMETERS()
    except Exception:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmpmanifest)
        return 2

    # check if this is going to be a domain search
    domains = False
    if "DOMAIN" in rp.readParams:
        domains = True

    if domains:  # no POSCAR in main folder for domain searches
        sl = None
    else:
        if os.path.isfile(os.path.join(".", "POSCAR")):
            poscarfile = os.path.join(".", "POSCAR")
            logger.info("Reading structure from file POSCAR")
            try:
                sl = readPOSCAR(filename=poscarfile)
            except Exception:
                logger.error("Exception while reading POSCAR", exc_info=True)
                cleanup(tmpmanifest)
                return 2
        else:
            logger.error("POSCAR not found. Stopping execution...")
            cleanup(tmpmanifest)
            return 2

        if not sl.preprocessed:
            logger.info("The POSCAR file will be processed and overwritten. "
                        "Copying the original POSCAR to POSCAR_user...")
            try:
                shutil.copy2(poscarfile, "POSCAR_user")
                tmpmanifest.append("POSCAR_user")
            except Exception:
                logger.error("Failed to copy POSCAR to POSCAR_user. Stopping "
                             "execution...")
                cleanup(tmpmanifest)
                return 2
    try:
        interpretPARAMETERS(rp, slab=sl)
        if rp.LOG_DEBUG:
            logger.setLevel(logging.DEBUG)
            logger.debug("PARAMETERS file was read successfully")
    except Exception:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmpmanifest)
        return 2
    rp.timestamp = timestamp
    rp.manifest = tmpmanifest
    if not domains:
        sl.fullUpdate(rp)   # gets PARAMETERS data into slab
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
    if os.path.isdir(os.path.join(".", "workhistory")):
        try:
            shutil.rmtree(os.path.join(".", "workhistory"))
        except Exception:
            logger.warning("Failed to clear workhistory folder.")
    # get rid of old POSCAR_OUT, VIBROCC_OUT and R_OUT files:
    for d in [".", os.path.join(".", "OUT")]:
        if os.path.isdir(d):
            for s in ["POSCAR_OUT", "VIBROCC_OUT", "R_OUT"]:
                for f in [fn for fn in os.listdir(d) if fn.startswith(s)]:
                    try:
                        os.remove(os.path.join(d, f))
                    except Exception:
                        logger.debug("Failed to delete file {}"
                                     .format(os.path.join(d, f)))
    # clean up old executable files:
    for fn in ["refcalc", "rfactor", "search", "superpos"]:
        p = re.compile(fn+r'-\d{6}-\d{6}')
        for f in [f for f in os.listdir()
                  if len(f) == len(fn) + 14 and p.match(f)]:
            try:
                os.remove(f)
            except Exception:
                logger.debug("Failed to delete file {}".format(f))
    # see if there are old logfiles
    oldlogs = [f for f in os.listdir() if os.path.isfile(f) and
               f.endswith(".log") and f != logname]
    if len(oldlogs) > 0:
        try:
            moveoldruns(rp, prerun=True)
        except Exception:
            logger.warning("Exception while trying to clean up from previous "
                           "run. Program will proceed, but old files may be "
                           "lost.", exc_info=True)
    if os.path.isfile("fortran-compile.log"):
        try:
            os.remove("fortran-compile.log")
        except Exception:
            pass

    sectionorder = [0, 1, 11, 2, 3, 31, 12, 4, 5]
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
                try:
                    moveoldruns(rp)
                except Exception:
                    logger.warning(
                        "Exception while trying to clean up earlier segments. "
                        "Program will proceed, but old files may be lost.",
                        exc_info=True)
            runSection(sec, sl, rp)
            if domains and sl is None:
                sl = rp.pseudoSlab
            if rp.domainParams:
                rp.setHaltingLevel(max([dp.rp.halt for dp in rp.domainParams]))
            if (sec == 0 and not domains and not sl.preprocessed
                    and rp.HALTING <= 2 and len(rp.RUN) > 0):
                logger.info(
                    "Initialization finished. Execution will stop. Please "
                    "check whether comments in POSCAR are correct, then "
                    "restart.")
                rp.checklist.append(
                    "Check whether comments in POSCAR are correct")
                rp.setHaltingLevel(2)
                initHalt = True
            elif (sec == 1 and rp.fileLoaded["EXPBEAMS"]):
                if (rp.RUN[:1] != [11] and          # r-factor after refcalc
                        (not rp.domainParams or 3 in rp.runHistory)):
                    rp.RUN.insert(0, 11)
            elif (sec == 3 and rp.fileLoaded["EXPBEAMS"]):
                if rp.RUN[:1] != [31]:  # superpos after search
                    rp.RUN.insert(0, 31)
            elif sec == 31 and rp.fileLoaded["EXPBEAMS"]:
                if rp.RUN[:1] != [12]:   # r-factor after superpos
                    rp.RUN.insert(0, 12)
            elif sec == 12 and not rp.STOP:
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
                                    + rp.disp_blocks[rp.search_index][1])
                    else:
                        rp.search_index += 1
                        o = "Search loop ends."
                        if len(rp.disp_blocks) > rp.search_index:
                            o += (" Continuing at block "
                                  + rp.disp_blocks[rp.search_index][1])
                        logger.info(o)
                else:
                    rp.search_index += 1
                for dp in rp.domainParams:
                    dp.rp.search_index = rp.search_index
                if len(rp.disp_blocks) > rp.search_index:
                    if not rp.domainParams:
                        sl.restoreOriState()
                    rp.resetSearchConv()
                    for dp in rp.domainParams:
                        dp.sl.restoreOriState()
                        dp.rp.resetSearchConv()
                    if rp.SEARCH_START == "control":
                        rp.SEARCH_START = "crandom"
                    if rp.RUN[:2] != [2, 3]:
                        rp.RUN = [2, 3] + rp.RUN
        except KeyboardInterrupt:
            logger.warning("Stopped by keyboard interrupt, attempting "
                           "clean exit...")
            cleanup(rp.manifest, rp)
            return 1
        except Exception:
            logger.error("Exception during tleedm execution: ", exc_info=True)
            cleanup(rp.manifest, rp)
            return 3
        if rp.halt >= rp.HALTING:
            if not initHalt:
                logger.info(
                    "# An exception occured that meets the halting "
                    "criteria defined by the HALTING parameter. Execution "
                    "will stop, check log for warnings and errors.")
            break
        updatePARAMETERS(rp)
        if rp.RUN and rp.STOP and not rp.RUN[0] in [11, 12, 31]:
            logger.info("# Stopped by user STOP command.")
            break
    cleanup(rp.manifest, rp)
    return 0


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
