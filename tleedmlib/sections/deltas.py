# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Delta Amplitudes
"""

import os
import logging
import shutil
import subprocess
import hashlib
import numpy as np

import viperleed.tleedmlib as tl
import viperleed.tleedmlib.files.iodeltas as io
from viperleed.tleedmlib.files.beams import writeAUXBEAMS
from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS_block
# from viperleed.tleedmlib.files.parameters import updatePARAMETERS
from viperleed.tleedmlib.leedbase import monitoredPool

logger = logging.getLogger("tleedm.deltas")


class DeltaCompileTask():
    """Stores information for a worker to compile a delta file, and keeps
    track of the folder that the compiled file is in afterwards."""

    def __init__(self, param, hash_, index):
        self.param = param
        self.hash = hash_
        self.foldername = "Delta_Compile_{}".format(index)
        self.exename = "delta-{}".format(index)
        self.fortran_comp = ["", ""]
        self.sourcedir = ""  # where the fortran files are
        self.basedir = ""    # where the calculation is based


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
    logger = logging.getLogger("tleedm.deltas")
    home = os.getcwd()
    base = runtask.comptask.basedir
    workname = "calculating_"+runtask.deltaname
    workfolder = os.path.join(base, workname)
    # make folder and go there:
    if os.path.isdir(workfolder):
        logger.warning("Folder "+workname+" already exists. "
                       "Contents may get overwritten.")
    else:
        os.mkdir(workfolder)
    os.chdir(workfolder)
    # get tensor file
    if os.path.isfile(os.path.join(base, "Tensors", runtask.tensorname)):
        try:
            shutil.copy2(os.path.join(base, "Tensors", runtask.tensorname),
                         "AMP")
        except Exception:
            logger.error("Error copying Tensor file: ", exc_info=True)
            return ("Error encountered by DeltaRunTask " + runtask.deltaname
                    + ": Error copying Tensor file.")
    else:
        logger.error("Tensor file not found: " + runtask.tensorname)
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Tensor not found.")
    # get executable
    exename = runtask.comptask.exename
    try:
        shutil.copy2(os.path.join(base, runtask.comptask.foldername, exename),
                     os.path.join(workfolder, exename))
    except Exception:
        logger.error("Error getting delta executable: ", exc_info=True)
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Failed to get delta executable.")
    # run execution
    try:
        with open("delta.log", "w") as log:
            subprocess.run(os.path.join(workfolder, exename),
                           input=runtask.din, encoding="ascii",
                           stdout=log, stderr=log)
    except Exception:
        logger.error("Error while executing delta-amplitudes "
                     "calculation for " + runtask.deltaname + ". Also check "
                     "delta log file.")
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Error during delta execution.")
    # copy delta file out
    try:
        shutil.copy2(os.path.join(workfolder, "DELWV"),
                     os.path.join(base, runtask.deltaname))
    except Exception:
        logger.error("Failed to copy delta output file DELWV"
                     " to main folder as " + runtask.deltaname,
                     exc_info=True)
        return ("Error encountered by DeltaRunTask " + runtask.deltaname
                + ": Failed to copy result file out.")
    # append log
    log = ""
    try:
        with open("delta.log", "r") as rf:
            log = rf.read()
    except Exception:
        logger.warning("Could not read local delta log for "
                       + runtask.deltaname)
    if log != "":
        deltalog = os.path.join(base, runtask.deltalogname)
        try:
            with open(deltalog, "a") as wf:
                wf.write("\n\n### STARTING LOG FOR " + runtask.deltaname
                         + " ###\n\n" + log)
        except Exception:
            logger.warning("Error writing delta log part "
                           + runtask.deltaname + ": ", exc_info=True)
    # clean up
    os.chdir(home)
    try:
        shutil.rmtree(workfolder)
    except Exception:
        logger.warning("Error deleting folder " + workname)
    return ""


def compileDelta(comptask):
    """Function meant to be executed by parallelized workers. Executes a
    DeltaCompileTask."""
    home = os.getcwd()
    workfolder = os.path.join(comptask.basedir, comptask.foldername)
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
    except Exception:
        logger.error("Error writing PARAM file: ", exc_info=True)
        return ("Error encountered by DeltaCompileTask " + comptask.foldername
                + "while trying to write PARAM file.")
    # get Fortran source files
    try:
        libpath = os.path.join(comptask.sourcedir, 'lib')
        libname1 = [f for f in os.listdir(libpath)
                    if f.startswith('lib.tleed')][0]
        shutil.copy2(os.path.join(libpath, libname1), libname1)
        libname2 = [f for f in os.listdir(libpath)
                    if f.startswith('lib.delta')][0]
        shutil.copy2(os.path.join(libpath, libname2), libname2)
        srcpath = os.path.join(comptask.sourcedir, 'src')
        srcname = [f for f in os.listdir(srcpath)
                   if f.startswith('delta')][0]
        shutil.copy2(os.path.join(srcpath, srcname), srcname)
        globalname = "GLOBAL"
        shutil.copy2(os.path.join(srcpath, globalname), globalname)
    except Exception:
        logger.error("Error getting TensErLEED files for "
                     "delta-amplitudes: ", exc_info=True)
        return ("Error encountered by DeltaCompileTask " + comptask.foldername
                + "while trying to fetch fortran source files")
    # compile
    ctasks = [(comptask.fortran_comp[0] + " -o " + oname + " -c",
               fname, comptask.fortran_comp[1]) for (fname, oname)
              in [(srcname, "main.o"), (libname1, "lib.tleed.o"),
                  (libname2, "lib.delta.o")]]
    ctasks.append((comptask.fortran_comp[0] + " -o " + comptask.exename,
                   "main.o lib.tleed.o lib.delta.o",
                   comptask.fortran_comp[1]))
    try:
        tl.leedbase.fortran_compile_batch(ctasks)
    except Exception as e:
        logger.error("Error compiling fortran files: " + str(e))
        return ("Fortran compile error in DeltaCompileTask "
                + comptask.foldername)
    os.chdir(home)
    return ""


def deltas(sl, rp, subdomain=False):
    """Runs the delta-amplitudes calculation."""

    if rp.domainParams:
        deltas_domains(rp)
        return

    # read DISPLACEMENTS block
    if not rp.disp_block_read:
        readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
        rp.disp_block_read = True
    # get Tensors
    if not os.path.isdir(os.path.join(".", "Tensors")):
        logger.error("No Tensors directory found.")
        raise RuntimeError("Tensors not found")
    tl.leedbase.getTensors(rp.TENSOR_INDEX)
    if 1 not in rp.runHistory:
        dn = "Tensors_"+str(rp.TENSOR_INDEX).zfill(3)
        logger.debug(
            "Running without reference calculation, checking "
            "input files in "+dn+" to determine original configuration.")
        tl.leedbase.getTensorOriStates(sl, os.path.join(".", "Tensors", dn))
        sl.restoreOriState(keepDisp=True)
    # if there are old deltas, fetch them
    tl.leedbase.getDeltas(rp.TENSOR_INDEX, required=False)
    dbasic = io.generateDeltaBasic(sl, rp)
    # get AUXBEAMS; if AUXBEAMS is not in work folder, check SUPP folder
    if not os.path.isfile(os.path.join(".", "AUXBEAMS")):
        if os.path.isfile(os.path.join(".", "SUPP", "AUXBEAMS")):
            try:
                shutil.copy2(os.path.join(".", "SUPP", "AUXBEAMS"),
                             "AUXBEAMS")
            except Exception:
                logger.warning("Failed to copy AUXBEAMS from SUPP folder. "
                               "Generating new file...")
                try:
                    writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
                except Exception:
                    logger.error("Exception during writeAUXBEAMS: ")
                    raise
        else:
            try:
                writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist)
            except Exception:
                logger.error("Exception during writeAUXBEAMS: ")
                raise
    try:
        with open("AUXBEAMS", "r") as rf:
            auxbeams = rf.read()
        if auxbeams[-1] != "\n":
            auxbeams += "\n"
    except Exception:
        logger.error("Could not read AUXBEAMS for delta-input")
        raise
    # get PHASESHIFTS
    try:
        with open("PHASESHIFTS", "r") as rf:
            phaseshifts = rf.read()
        if phaseshifts[-1] != "\n":
            phaseshifts += "\n"
    except Exception:
        logger.error("Could not read PHASESHIFTS for delta-input")
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
                if abs(at.disp_vib[el][0]) >= 1e-4:
                    found = True
                    break
        if not found:
            for el in at.disp_geo:
                if np.linalg.norm(at.disp_geo[el][0]) >= 1e-4:
                    found = True
                    break
        if not found:
            occlists = []
            for k in at.disp_occ:
                occlists.append(at.disp_occ[k])
            for i in range(0, len(occlists[0])):
                totalocc = 0.
                for ol in occlists:
                    if len(ol) <= i:
                        break  # error - will pop up again later...
                    else:
                        totalocc += ol[i]
                if totalocc < 1 - 1e-4:
                    found = True
                    break
        if not found:
            attodo.pop(j)
        else:
            j += 1

    vaclist = []    # atoms for which a vacancy delta file is needed
    for at in attodo:
        occlists = []
        for k in at.disp_occ:
            occlists.append(at.disp_occ[k])
        for i in range(0, len(occlists[0])):
            totalocc = 0.
            for ol in occlists:
                if len(ol) <= i:
                    logger.error("Inconsistent occupancy lists for {} "
                                 .format(at))
                    raise ValueError("Inconsistent occupancy lists for {}"
                                     .format(at))
                else:
                    totalocc += ol[i]
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
                if io.checkDelta(df, at, el, rp):
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
        return
    elif countExisting > 0:
        logger.info("{} of {} required Delta-files are already present. "
                    "Generating remaining {} files..."
                    .format(countExisting, len(atElTodo) + countExisting,
                            len(atElTodo)))
    # create log file:
    deltaname = "delta-"+rp.timestamp
    deltalogname = deltaname+".log"
    if not subdomain:
        logger.info(
            "Generating delta files...\n"
            "Delta log will be written to local subfolders, and collected in "
            + deltalogname)
    rp.manifest.append(deltalogname)
    try:
        with open(deltalogname, "w") as wf:
            wf.write("Logs from multiple delta calculations are collected "
                     "here. Their order may not be preserved.\n")
    except Exception:
        logger.warning("Error creating delta log file. This will not "
                       "affect execution, proceeding...")

    # move PARAM file
    if os.path.isfile("PARAM"):
        try:
            os.rename("PARAM", "PARAM-old")
        except Exception:
            try:
                os.remove(os.path.join(".", "PARAM"))
            except Exception:
                logger.warning(
                    "Section Delta-Amplitudes: Cannot rename/remove old PARAM "
                    "file. This might cause the Delta generation to fail!")
    # assemble tasks
    deltaCompTasks = []  # keep track of what versions to compile
    deltaRunTasks = []   # which deltas to run
    tensordir = "Tensors_"+str(rp.TENSOR_INDEX).zfill(3)
    for (at, el) in atElTodo:
        din, din_short, param = io.generateDeltaInput(
            at, el, sl, rp, dbasic, auxbeams, phaseshifts)
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
            except Exception:
                pass
        if nums:
            n = max(nums) + 1
        rt.deltaname = nameBase + "_{}".format(n)
        rt.deltalogname = deltalogname
        at.deltasGenerated.append(rt.deltaname)

    # sort deltasGenerated
    for at in attodo:
        checkEls = list(at.disp_occ.keys())
        if at in vaclist:
            checkEls.append("vac")
        copydel = at.deltasGenerated[:]
        at.deltasGenerated = []
        for el in checkEls:
            at.deltasGenerated.append(
                [df for df in copydel
                 if df.split("_")[-2].lower() == el.lower()][0])
        if len(at.deltasGenerated) != len(copydel):
            logger.error("Failed to sort delta files for {}".format(at))
            raise RuntimeError("Inconsistent delta files")

    # write delta-input file
    dinput = ("""# ABOUT THIS FILE:
# Input for the delta-calculations is collected here. The blocks of data are
# new 'PARAM' files, which are used to recompile the fortran code, and input
# for generation of specific DELTA files. Lines starting with '#' are comments
# on the function of the next block of data.
# In the DELTA file blocks, [AUXBEAMS] and [PHASESHIFTS] denote where the
# entire contents of the AUXBEAMS and PHASESHIFTS files should be inserted.
""")
    for ct in deltaCompTasks:
        dinput += ("\n#### NEW 'PARAM' FILE: ####\n\n" + ct.param + "\n")
        for rt in [t for t in deltaRunTasks if t.comptask == ct]:
            dinput += ("\n#### INPUT for new DELTA file {}: ####\n\n"
                       .format(rt.deltaname) + rt.din_short + "\n")
    try:
        with open("delta-input", "w") as wf:
            wf.write(dinput)
    except Exception:
        logger.warning("Failed to write file 'delta-input'. This will "
                       "not affect TensErLEED execution, proceeding...")

    # if execution is suppressed, stop here
    if rp.SUPPRESS_EXECUTION and not subdomain:
        rp.setHaltingLevel(3)
        return

    # make sure there's a compiler ready:
    if rp.FORTRAN_COMP[0] == "" and not subdomain:
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("No Fortran compiler")
    tlp = tl.leedbase.getTLEEDdir(os.path.abspath(rp.sourcedir),
                                  version=rp.TL_VERSION)
    if not tlp:
        raise RuntimeError("TensErLEED code not found.")
    for ct in deltaCompTasks:
        ct.fortran_comp = rp.FORTRAN_COMP
        ct.sourcedir = tlp
        ct.basedir = os.getcwd()

    if subdomain:   # actual calculations done in deltas_domains
        if len(deltaRunTasks) > 0:
            rp.manifest.append("Deltas")
        return (deltaCompTasks, deltaRunTasks)

    rp.updateCores()
    # compile files
    logger.info("Compiling fortran files...")
    poolsize = min(len(deltaCompTasks), rp.N_CORES)
    monitoredPool(rp, poolsize, compileDelta, deltaCompTasks)
    if rp.STOP:
        return

    # run executions
    logger.info("Running delta calculations...")
    poolsize = min(len(deltaRunTasks), rp.N_CORES)
    monitoredPool(rp, poolsize, runDelta, deltaRunTasks)
    if rp.STOP:
        return
    logger.info("Delta calculations finished.")

    # clean up
    for ct in deltaCompTasks:
        try:
            shutil.rmtree(os.path.join(ct.basedir, ct.foldername))
        except Exception:
            logger.warning("Error deleting delta compile folder "
                           + ct.foldername)
    rp.manifest.append("Deltas")
    return


def deltas_domains(rp):
    """Define and run delta calculations for all domains."""
    home = os.getcwd()
    deltaCompTasks = []
    deltaRunTasks = []
    # get input for all domains
    for dp in rp.domainParams:
        logger.info("Getting input for delta calculations: domain {}"
                    .format(dp.name))
        try:
            os.chdir(dp.workdir)
            r = deltas(dp.sl, dp.rp, subdomain=True)
        except Exception:
            logger.error("Error while creating delta input for domain {}"
                         .format(dp.name))
            raise
        finally:
            os.chdir(home)
        if type(r) == tuple:  # if no deltas need to be calculated returns None
            deltaCompTasks.extend(r[0])
            deltaRunTasks.extend(r[1])
        else:
            raise RuntimeError("Unknown error while creating delta input for "
                               "domain {}".format(dp.name))

    # if execution is suppressed, stop here
    if rp.SUPPRESS_EXECUTION:
        rp.setHaltingLevel(3)
        return

    # make sure there's a compiler ready, and we know the number of cores:
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise RuntimeError("No Fortran compiler")
    for ct in deltaCompTasks:
        ct.fortran_comp = rp.FORTRAN_COMP
    rp.updateCores()  # if number of cores is not defined, try to find it

    # compile files
    if len(deltaCompTasks) > 0:
        logger.info("Compiling fortran files...")
        poolsize = min(len(deltaCompTasks), rp.N_CORES)
        monitoredPool(rp, poolsize, compileDelta, deltaCompTasks)
        if rp.STOP:
            return

    # run executions
    if len(deltaRunTasks) > 0:
        logger.info("Running delta calculations...")
        poolsize = min(len(deltaRunTasks), rp.N_CORES)
        monitoredPool(rp, poolsize, runDelta, deltaRunTasks)
        if rp.STOP:
            return
        logger.info("Delta calculations finished.")

    # clean up
    for ct in deltaCompTasks:
        d = os.path.join(ct.basedir, ct.foldername)
        try:
            shutil.rmtree(d)
        except Exception:
            logger.warning("Error deleting delta compile folder "
                           + os.path.relpath(d))

    return