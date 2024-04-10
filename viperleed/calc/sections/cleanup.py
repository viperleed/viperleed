"""Cleanup functions, to be used between sections or before/after execution."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2021-06-04'
__license__ = 'GPLv3+'

import logging
import os
from pathlib import Path                                                        # TODO: use everywhere
import re
import shutil
import time
from timeit import default_timer as timer
from zipfile import ZIP_DEFLATED, ZipFile

from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc import LOG_PREFIX
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib.base import copytree_exists_ok,get_elapsed_time_str

_DEFAULT_SUPP = "SUPP"
_DEFAULT_OUT = "OUT"

# files to go in SUPP
_SUPP_FILES = (
    "AUXBEAMS",
    "AUXEXPBEAMS",
    "AUXGEO",
    "AUXLATGEO",
    "AUXNONSTRUCT",
    "BEAMLIST",
    "delta-input",
    "EEASISSS-input.txt",
    "eeasisss-input",
    "EEASISSS-log.txt",
    "muftin.f",
    "Phaseshifts_plots.pdf",
    "POSCAR_bulk_appended",
    "POSCAR_bulk",
    "POSCAR_mincell",
    "POSCAR_oricell",
    "POSCAR_vacuum_corrected",
    "refcalc-FIN",
    "refcalc-PARAM",
    "restrict.f",
    "rfactor-PARAM",
    "rfactor-WEXPEL",
    "search-PARAM",
    "search-rf.info",
    "search.steu",
    "searchpars.info",
    "superpos-CONTRIN",
    "superpos-PARAM",
    )

_SUPP_DIRS = (ORIGINAL_INPUTS_DIR_NAME, "compile_logs")

# files to go in OUT
_OUT_FILES = (
    "Complex_amplitudes_imag.csv",
    "Complex_amplitudes_real.csv"
    "control.chem",
    "Errors_summary.csv",
    "Errors.pdf",
    "Errors.zip",
    "FD_Optimization_beams.pdf",
    "FD_Optimization.csv",
    "FD_Optimization.pdf",
    "FITBEAMS_norm.csv",
    "FITBEAMS.csv",
    "PARAMETERS",
    "PatternInfo.tlm",
    "refcalc-amp.out",
    "Rfactor_analysis_refcalc.pdf",
    "Rfactor_analysis_superpos.pdf",
    "Rfactor_plots_refcalc.pdf",
    "Rfactor_plots_superpos.pdf",
    "SD.TL", "refcalc-fd.out",
    "Search-progress.csv",
    "Search-progress.pdf",
    "Search-report.pdf",
    "superpos-spec.out",
    "THEOBEAMS_norm.csv",
    "THEOBEAMS.csv",
    "THEOBEAMS.pdf",
    )

# Label given to workhistory folders when cleaning up stray remains
# from previous viperleed.calc executions from the work directory
PREVIOUS_LABEL = 'previous'

# output files that can be used as input in future runs - keep during prerun
iofiles = ["control.chem", "refcalc-fd.out", "superpos-spec.out"]

logger = logging.getLogger(__name__)


def prerun_clean(rp, logname=""):
    """Clean up the work directory before viperleed.calc starts.

    Delete workhistory, old executables, and old logfiles.
    Call move_oldruns if required.

    Parameters
    ----------
    rp : Rparams
        The run parameters, needed for move_oldruns.
    logname : str, optional
        Name of the current log file, to be excluded from cleanup.

    Returns
    -------
    None.
    """
    # clean out the workhistory folder, if there is one
    if os.path.isdir(os.path.join(".", DEFAULT_WORK_HISTORY)):
        try:
            shutil.rmtree(os.path.join(".", DEFAULT_WORK_HISTORY))
        except Exception:
            logger.warning(f"Failed to clear {DEFAULT_WORK_HISTORY} folder.")
    # get rid of old POSCAR_OUT, VIBROCC_OUT, PARAMETERS and R_OUT files:
    for d in [".", os.path.join(".", "OUT")]:
        if os.path.isdir(d):
            for s in ["POSCAR_OUT", "VIBROCC_OUT", "PARAMETERS", "R_OUT"]:
                for file in [fn for fn in os.listdir(d) if fn.startswith(s)]:
                    try:
                        os.remove(os.path.join(d, file))
                    except Exception:
                        logger.debug("Failed to delete file {}"
                                     .format(os.path.join(d, file)))
    # clean up old executable files:
    for fn in ["refcalc", "rfactor", "search", "superpos"]:
        p = re.compile(fn+r'-\d{6}-\d{6}')
        for file in [f for f in os.listdir()
                  if len(f) == len(fn) + 14 and p.match(f)]:
            try:
                os.remove(file)
            except Exception:
                logger.debug(f"Failed to delete file {file}")
    # see if there are old logfiles
    old_logs = [f for f in os.listdir() if os.path.isfile(f) and
               f.endswith(".log") and f != logname]
    if len(old_logs) > 0:
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


def organize_workdir(tensor_index, delete_unzipped=False,
                     tensors=True, deltas=True, workdir=Path(),
                     compression_level=2):
    """Reorganize files in workdir into SUPP, OUT, Tensors and Deltas.

    Tensors and Deltas folders are zipped and moved over. All other
    files are copied to appropriate locations in SUPP and OUT.

    Parameters
    ----------
    tensor_index : int
        Which Delta and Tensor files should be considered.
    delete_unzipped : bool, optional
        Whether the original Delta- and Tensor-files should be
        deleted after making the archives. The default is False.
    tensors, deltas : bool, optional
        Whether the Tensor/Delta files contain new information
        and should be saved. The default is True.
    workdir : pathlike, optional
        The path to work folder that contains the files to be
        reorganized. The default is "".
    compression_level : int
        Compression level to be applied for ZIP archives of Tensors and
        Deltas. Default is 2.

    Returns
    -------
    None.
    """
    # outfiles with variable names:
    path = Path(workdir)
    outfiles = set(path / f for f in _OUT_FILES)
    for pattern in ("POSCAR_OUT*", "VIBROCC_OUT*", "R_OUT*"):
        outfiles.update(path.glob(pattern))

    _collect_deltas(tensor_index, path)
    _zip_deltas_and_tensors(delete_unzipped, tensors, deltas, path,
                            compression_level)
    _organize_supp_out(path, outfiles)

def _organize_supp_out(path, outfiles):
    """Helper function for organizing SUPP and OUT directories."""
    # SUPP
    supp_path = path / _DEFAULT_SUPP

    files_to_copy = set(path / f for f in _SUPP_FILES)
    # move directories original_inputs and compile_logs to SUPP
    directories_to_copy = (path / d for d in _SUPP_DIRS)
    # Also add log files into SUPP: skip calc logs (they go to
    # main dir), and compile logs (they go to compile_logs dir)
    logs_to_supp = (f for f in path.glob("*.log")
                    if (not f.name.startswith(LOG_PREFIX)
                        and "compile" not in f.name))
    files_to_copy.update(logs_to_supp)

    _copy_files_and_directories(files_to_copy, directories_to_copy, path,
                                supp_path)

    # OUT
    out_path = path / _DEFAULT_OUT
    _copy_files_and_directories(outfiles, (), path, out_path)


def _copy_files_and_directories(filelist, directory_list, origin, target):
    """Helper function for copying files and directories to SUPP and OUT."""
    folder = target.name  # SUPP or OUT
    # Create the directory
    try:
        target.mkdir(parents=True)
    except FileExistsError:
        pass
    except OSError:
        logger.error(f"Error creating {folder} folder: ", exc_info=True)
        return

    # Copy files and directories
    for file in filelist:
        if not file.is_file():
            continue
        # copies files into SUPP and OUT directories
        try:
            shutil.copy2(file, target / file.name)
        except OSError:
            logger.error(f"Error moving {folder} file {file.name}: ",
                            exc_info=True)

    for _dir in directory_list:
        if not _dir.is_dir():
            continue
        try:
            copytree_exists_ok(_dir, target / _dir.name)
        except OSError:
            logger.error(f"Error moving {folder} directory {_dir.name}: ",
                            exc_info=True)


def _zip_deltas_and_tensors(delete_unzipped, tensors, deltas, path,
                            compression_level):
    # If there are unzipped Tensors or Deltas directories, zip them:
    for folder in ["Tensors", "Deltas"]:
        todo = tensors if folder == "Tensors" else deltas
        origin_base = path / folder
        if not origin_base.is_dir():
            continue
        if not todo and not delete_unzipped:
            continue
        rgx = re.compile(rf"{folder}_[0-9]{{3}}")                               # TODO: maybe we want "three or more" digits, i.e., {{3,}}? Or could we use tensor_index?
        for _dir in origin_base.glob("*"):
            if not _dir.is_dir():
                continue
            match = rgx.match(_dir.name)
            if not match or match.span()[1] != len(folder) + 4:                 # TODO: should this 4 be adjusted to the previous TODO? Unclear what it guards
                continue
            delete = delete_unzipped
            if todo:
                logger.info(f"Packing {_dir.name}.zip...")
                _dir_path = Path(_dir)
                move_to_archive = _dir_path.glob('*')
                arch_name = _dir_path.with_suffix(".zip")
                try:
                    with ZipFile(arch_name, 'a', compression=ZIP_DEFLATED,
                        compresslevel=compression_level) as archive:
                        for fname in move_to_archive:
                            archive.write(fname, fname.relative_to(_dir))
                except OSError:
                    logger.error(f"Error packing {_dir.name}.zip file: ",
                                    exc_info=True)
                    delete = False
            if delete:
                try:
                    shutil.rmtree(_dir)
                except OSError:
                    logger.warning(
                        f"Error deleting unzipped {folder} directory. "
                        "This will increase the size of the work folder, "
                        "but not cause any problems.")


def _collect_deltas(tensor_index, path):
    # Clean up deltas
    deltalist = list(path.glob("DEL_*"))
    if len(deltalist) > 0:
        destination = path / "Deltas" / f"Deltas_{tensor_index:03d}"
        try:
            destination.mkdir(parents=True)
        except FileExistsError:
            pass
        except OSError:
            logger.error(f"Failed to create {destination} folder: ",
                         exc_info=True)
        if destination.exists():
            errors = []
            for delta_file in deltalist:
                try:
                    shutil.move(delta_file, destination / delta_file.name)
                except OSError as err:
                    errors.append(err)
            if errors:
                logger.error(f"Error moving Delta files: {errors}")


def move_oldruns(rp, prerun=False):
    """Copy relevant files to a new 'workhistory' subfolder.

    Files are copied from SUPP, OUT and the list in rp.manifest.
    The main log file is excluded.

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
    None.
    """
    sectionabbrv = {1: "R", 2: "D", 3: "S"}
    try:
        os.makedirs(os.path.join(".", DEFAULT_WORK_HISTORY), exist_ok=True)
    except Exception:
        logger.error(f"Error creating {DEFAULT_WORK_HISTORY} folder: ",
                     exc_info=True)
        raise
    if not prerun:
        rp.manifest.append(DEFAULT_WORK_HISTORY)
    dl = [n for n in os.listdir(DEFAULT_WORK_HISTORY)
          if os.path.isdir(os.path.join(DEFAULT_WORK_HISTORY, n))]
    maxnum = -1
    rgx = re.compile(r't'+'{:03d}'.format(rp.TENSOR_INDEX)+r'.r[0-9]{3}_')             # TODO: would be nicer to use a capture group for the last three digits, used to decide how to number the new folder
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
                              f.endswith(".log") and f.startswith(LOG_PREFIX)
                              and f not in rp.manifest])
        if len(oldlogfiles) > 0:
            oldTimeStamp = oldlogfiles[-1][7:20]
        else:
            oldTimeStamp = "moved-" + rp.timestamp
        dirname = (f"t{rp.TENSOR_INDEX:03d}.r{num:03d}_{PREVIOUS_LABEL}_"
                   + oldTimeStamp)
    else:
        dirname = "t{:03d}.r{:03d}_".format(rp.TENSOR_INDEX, num)
        for ind in rp.runHistory[len(rp.lastOldruns):]:
            if ind in sectionabbrv:
                dirname += sectionabbrv[ind]
        rp.lastOldruns = rp.runHistory[:]
        dirname += "_" + rp.timestamp

    # make workhistory directory
    work_hist_path = Path(".") / DEFAULT_WORK_HISTORY / dirname
    try:
        os.mkdir(work_hist_path)
    except Exception:
        logger.error(f"Error creating {DEFAULT_WORK_HISTORY} subfolder: ",
                     exc_info=True)
        raise
    if not prerun:
        organize_workdir(rp.TENSOR_INDEX, delete_unzipped=False,
                         tensors=False, deltas=False,
                         compression_level=rp.ZIP_COMPRESSION_LEVEL)
        for dp in rp.domainParams:
            organize_workdir(dp.rp.TENSOR_INDEX, delete_unzipped=False,
                             tensors=False, deltas=False,
                             compression_level=rp.ZIP_COMPRESSION_LEVEL)
    if prerun:
        filelist = [f for f in os.listdir() if os.path.isfile(f) and
                    (f.endswith(".log") or f in _OUT_FILES or f in _SUPP_FILES)
                    and f not in rp.manifest and f not in iofiles]
        dirlist = ["SUPP", "OUT"]
    else:
        filelist = [f for f in rp.manifest if os.path.isfile(f) and not
                    (f.startswith(LOG_PREFIX) and f.endswith(".log"))]
        dirlist = [d for d in rp.manifest if os.path.isdir(d) and
                   d not in ["Tensors", "Deltas", DEFAULT_WORK_HISTORY]]
    for f in filelist:
        try:
            if not prerun or f in iofiles:
                shutil.copy2(f, work_hist_path / f)
            else:
                shutil.move(f, work_hist_path / f)
        except Exception:
            logger.warning(f"Error copying {f} to {work_hist_path / f}."
                           " File may get overwritten.")
    for d in dirlist:
        try:
            if not prerun:
                shutil.copytree(d, work_hist_path / d)
            else:
                shutil.move(d, work_hist_path / d)
        except Exception:
            logger.warning(f"Error copying {d} to {work_hist_path / d}."
                           " Files in directory may get overwritten.")
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
        compress_level = 2
    else:
        history = rp.runHistory
        rp.closePdfReportFigs()
        compress_level = rp.ZIP_COMPRESSION_LEVEL
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
            organize_workdir(d["tind"], delete_unzipped=True,
                             tensors=d["newTensors"],
                             deltas=d["newDeltas"], workdir=d["path"],
                             compression_level=compress_level)
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
