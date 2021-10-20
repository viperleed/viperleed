# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 16:03:36 2021

@author: Florian Kraushofer

Cleanup functions, to be used between sections or before/after tleedm
execution.
"""

import time
from timeit import default_timer as timer
import logging
import os
import shutil
import re

from viperleed.tleedmlib.base import get_elapsed_time_str

suppfiles = ["AUXBEAMS", "AUXGEO", "AUXLATGEO", "AUXNONSTRUCT", "BEAMLIST",
             "POSCAR_oricell", "POSCAR_bulk", "muftin.f",
             "refcalc-PARAM", "refcalc-FIN", "rfactor-WEXPEL",
             "rfactor-PARAM", "delta-input", "search.steu",
             "search-rf.info", "search-PARAM", "AUXEXPBEAMS",
             "eeasisss-input", "searchpars.info", "superpos-PARAM",
             "superpos-CONTRIN", "POSCAR_bulk_appended", "POSCAR_mincell",
             "restrict.f", "Phaseshifts_plots.pdf"]

supp_dirs = ["original_inputs", "compile_logs"]

outfiles = ["THEOBEAMS.csv", "THEOBEAMS_norm.csv",
            "PatternInfo.tlm", "SD.TL", "refcalc-fd.out",
            "Rfactor_plots_refcalc.pdf", "control.chem",
            "Search-progress.pdf", "Search-progress.csv",
            "Search-report.pdf", "FITBEAMS.csv", "FITBEAMS_norm.csv",
            "superpos-spec.out", "Rfactor_plots_superpos.pdf",
            "Rfactor_analysis_refcalc.pdf",
            "Rfactor_analysis_superpos.pdf", "Errors.csv", "Errors.pdf"]

logger = logging.getLogger("tleedm.sections.cleanup")


def prerun_clean(rp, logname=""):
    """
    Cleans up the work directory before tleedm starts. Deletes workhistory,
    old executables, and old logfiles. Calls move_oldruns if required.

    Parameters
    ----------
    rp : Rparams
        The run parameters, needed for move_oldruns.
    logname : str, optional
        Name of the current log file, to be excluded from cleanup.

    """
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
            move_oldruns(rp, prerun=True)
        except Exception:
            logger.warning("Exception while trying to clean up from previous "
                           "run. Program will proceed, but old files may be "
                           "lost.", exc_info=True)
    if os.path.isfile("fortran-compile.log"):
        try:
            os.remove("fortran-compile.log")
        except Exception:
            pass
    return


def sortfiles(tensorIndex, delete_unzipped=False, tensors=True,
              deltas=True, path=""):
    """
    Makes Tensors and Deltas zip files. Copies files to SUPP and OUT folders
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
    # move files to SUPP and OUT folders

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
    # sort SUPP and OUT files:
    for t in ["SUPP", "OUT"]:
        try:
            os.makedirs(os.path.join(path, t), exist_ok=True)
        except Exception:
            logger.error("Error creating {} folder: ".format(t), exc_info=True)
        if t == "SUPP":
            filelist = suppfiles
            directory_list = supp_dirs # move directories original_inputs and compile_logs to SUPP
        else:
            filelist = outfiles
            directory_list = []
        for f in [f for f in filelist
                  if os.path.isfile(os.path.join(path, f))]: # copies files into SUPP and OUT directories
            try:
                shutil.copy2(os.path.join(path, f), os.path.join(path, t, f))
            except Exception:
                logger.error("Error moving {} file {}: ".format(t, f),
                             exc_info=True)
        for d in directory_list:
            if os.path.isdir(os.path.join(path, d)):
                try:
                    shutil.copytree(os.path.join(path, d), os.path.join(path, t, d))
                except Exception:
                    logger.error("Error moving {} directory {}: ".format(t, d),
                                 exc_info=True)


def move_oldruns(rp, prerun=False):
    """
    Makes a new folder in 'workhistory'. Copies SUPP, OUT and files in
    manifest (except main log) to that new folder.

    Parameters
    ----------
    rp : Rparams
        The run parameters.
    prerun : bool, optional
        If True, then instead of using the manifest, all potentially
        interesting files will be copied, and the new folder will get index 0.
        Then clears out the SUPP and OUT folders and old SUPP and OUT files
        from the main directory.

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
                    (f.endswith(".log") or f in outfiles or f in suppfiles)
                    and f not in rp.manifest]
        dirlist = ["SUPP", "OUT"]
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


def cleanup(manifest, rp=None):
    """
    Moves files to SUPP and OUT folders, writes manifest, adds a final
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
            logger.warning("Error sorting files to SUPP/OUT folders: ",
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
    if rp is None:
        elapsedTimeStr = "unknown"
    else:
        elapsedTimeStr = get_elapsed_time_str(timer() - rp.starttime)
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
