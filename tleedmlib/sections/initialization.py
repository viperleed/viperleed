# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Initialization
"""

import os
import shutil
import logging
import copy

import tleedmlib as tl
import tleedmlib.beamgen
import tleedmlib.psgen
from tleedmlib.files.poscar import readPOSCAR, writeCONTCAR
from tleedmlib.files.parameters import readPARAMETERS, interpretPARAMETERS
from tleedmlib.files.phaseshifts import readPHASESHIFTS, writePHASESHIFTS
from tleedmlib.files.beams import (readOUTBEAMS, readBEAMLIST, checkEXPBEAMS, 
                                   sortIVBEAMS, writeIVBEAMS)
from tleedmlib.files.patterninfo import writePatternInfo

logger = logging.getLogger("tleedm.initialization")


def initialization(sl, rp, subdomain=False):
    """Runs the initialization. Returns 0 on success."""
    if not subdomain:
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
                    rp.expbeams = readOUTBEAMS(fn, enrange=er)
                    rp.fileLoaded["EXPBEAMS"] = True
                except:
                    logger.error("Error while reading file "+fn, exc_info=True)
    rp.initTheoEnergies()  # may be initialized based on exp. beams
    
    if rp.hasDomains and not subdomain:
        try:
            r = init_domains(rp)
        except:
            raise
        if r != 0:
            return r
        return 0
    
    # check whether _PHASESHIFTS are present & consistent:
    newpsGen, newpsWrite = True, True
                          # new phaseshifts need to be generated/written
    if os.path.isfile(os.path.join(".","_PHASESHIFTS")):
        try:
            (rp.phaseshifts_firstline, rp.phaseshifts,
                 newpsGen, newpsWrite) = readPHASESHIFTS(sl, rp)
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
                        rp.phaseshifts) = tl.psgen.runPhaseshiftGen(sl, rp,
                                               psgensource = rundgrenpath, 
                                               excosource = serneliuspath)
            logger.debug("Finished generating phaseshift data")
        except:
            logger.error("Exception while calling phaseshiftgen: ")
            raise
    if newpsWrite:
        try:
            writePHASESHIFTS(rp.phaseshifts_firstline, rp.phaseshifts)
        except:
            logger.error("Exception during writePHASESHIFTS: ")
            raise
    rp.fileLoaded["PHASESHIFTS"] = True
    rp.updateDerivedParams()
    rp.manifest.append("_PHASESHIFTS")

    # if necessary, run findSymmetry:
    if sl.planegroup == "unknown":
        tl.symmetry.findSymmetry(sl, rp)
        tl.symmetry.enforceSymmetry(sl, rp)
    
    # generate new POSCAR
    tmpslab = copy.deepcopy(sl)
    tmpslab.sortOriginal()
    try:
        writeCONTCAR(tmpslab, filename='POSCAR', comments='all')
    except:
        logger.error("Exception occurred while writing new POSCAR")
        raise
    rp.manifest.append('POSCAR')
    # generate POSCAR_oricell
    tmpslab.revertUnitCell()
    try:
        writeCONTCAR(tmpslab, filename='POSCAR_oricell', comments='nodir')
    except:
        logger.error("Exception occurred while writing POSCAR_oricell, "
                      "execution will continue...")

    # create bulk slab:
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    bsl = sl.bulkslab
    # find minimum in-plane unit cell for bulk:
    logger.info("Checking bulk unit cell...")
    changecell, mincell = bsl.getMinUnitCell(rp)
    if changecell:
        sl.changeBulkCell(rp, mincell)
        bsl = sl.bulkslab
    if not rp.superlattice_defined:
        ws = tl.leedbase.writeWoodsNotation(rp.SUPERLATTICE) 
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
    tl.symmetry.findSymmetry(bsl, rp, bulk=True, output=False)
    bsl.revertUnitCell() # keep origin matched with main slab
    logger.info("Found bulk plane group: "+bsl.foundplanegroup)
    tl.symmetry.findBulkSymmetry(bsl, rp)
    
    # write POSCAR_bulk
    bsl = copy.deepcopy(sl.bulkslab)
    bsl.sortOriginal()
    try:
        writeCONTCAR(bsl, filename='POSCAR_bulk', comments='bulk')
    except:
        logger.error("Exception occurred while writing POSCAR_bulk")
        raise
    
    # generate beamlist
    logger.info("Generating _BEAMLIST...")
    try:
        bgenpath = os.path.join('.', 'source', 'beamgen3.out')
        tl.beamgen.runBeamGen(sl,rp,beamgensource = bgenpath)
        # this does NOT read the resulting file!
    except:
        logger.error("Exception occurred while calling beamgen.")
        raise
    rp.manifest.append("_BEAMLIST")
    try:
        rp.beamlist = readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except:
        logger.error("Error while reading required file _BEAMLIST")
        raise
    
    if not subdomain:
        writePatternInfo(sl, rp)
        
        # if EXPBEAMS was loaded, it hasn't been checked yet - check now
        if rp.fileLoaded["EXPBEAMS"]:
            checkEXPBEAMS(sl, rp)
        # write and sort IVBEAMS
        if not rp.fileLoaded["IVBEAMS"]:
            try:
                rp.ivbeams = writeIVBEAMS(sl, rp)
                rp.ivbeams_sorted = False
                rp.fileLoaded["IVBEAMS"] = True
                rp.manifest.append("IVBEAMS")
            except:
                logger.error("Error while writing IVBEAMS file based on "
                              "EXPBEAMS data.")
                raise
        if not rp.ivbeams_sorted:
            rp.ivbeams = sortIVBEAMS(sl, rp)
            rp.ivbeams_sorted = True
        return 0

def init_domains(rp):
    """Runs an alternative initialization for the domain search. This will 
    include running the 'normal' initialization for each domain."""
    if len(rp.DOMAINS) < 2:
        logger.error("A domain search was defined, but less than two domains "
                     "are defined. Execution will stop.")
        rp.setHaltingLevel(3)
        return 0
    checkFiles = ["POSCAR", "PARAMETERS", "VIBROCC", "_PHASESHIFTS"]
    home = os.getcwd()
    for (name, path) in rp.DOMAINS:
        # determine the target path
        target = os.path.abspath("Domain_"+name)
        dp = tl.DomainParameters(target, home)
        if os.path.isdir(target):
            logger.warning("Folder "+target+" already exists. "
                           "Contents may get overwritten.")
        else:
            os.mkdir(target)
        logger.info("Fetching input files for domain {}".format(name))
        if os.path.isdir(path):
            # check the path for Tensors
            tensorIndex = tl.leedbase.getMaxTensorIndex(path)
            if tensorIndex != 0:
                try:
                    r = tl.leedbase.getTensors(tensorIndex, basedir=path, 
                                               targetdir=target)
                except:
                    tensorIndex = 0
                if r != 0:
                    tensorIndex = 0
            if tensorIndex != 0:
                tensorDir = os.path.join(target, "Tensors", 
                                         "Tensors_"+str(tensorIndex).zfill(3))
                for file in (checkFiles + ["IVBEAMS"]):
                    if os.path.isfile(os.path.join(tensorDir, file)):
                        shutil.copy2(os.path.join(tensorDir, file),
                                     os.path.join(target,file))
                    else:
                        tensorIndex = 0
                        break
            if tensorIndex != 0:
                dp.tensorDir = tensorDir
            else:       # no usable tensors in that dir; get input
                dp.refcalcRequired = True
                for file in checkFiles:
                    if os.path.isfile(os.path.join(path, file)):
                        try:
                            shutil.copy(os.path.join(path,file),
                                        os.path.join(target,file))
                        except:
                            if file != "_PHASESHIFTS":
                                logger.error("Error copying required file {}"
                                        "for domain {} from origin folder {}"
                                        .format(file, name, path))
                                return "Error getting domain input files"
                    elif file != "_PHASESHIFTS":
                        logger.error("Required file {} for domain {} not "
                                     "found in origin folder {}"
                                     .format(file, name, path))
                        return "Error getting domain input files"
        elif os.path.isfile(path):
            if not os.path.isdir(os.path.join(target, "Tensors")):
                os.mkdir(os.path.join(target, "Tensors"))
            tensorDir = os.path.join(target, "Tensors", 
                                     os.path.basename(path)[:-4])
            if not os.path.isdir(os.path.join(target,"Tensors",tensorDir)):
                os.mkdir(os.path.join(target,"Tensors",tensorDir))
            try:
                shutil.unpack_archive(path,tensorDir)
            except:
                logger.error("Failed to unpack Tensors for domain {} from "
                             "file {}".format(name, path))
                return "Error getting domain input files"
            for file in (checkFiles + ["IVBEAMS"]):
                if os.path.isfile(os.path.join(tensorDir, file)):
                    shutil.copy2(os.path.join(tensorDir, file),
                                 os.path.join(target,file))
                else:
                    logger.error("Required file {} for domain {} not found in "
                                 "Tensor directory {}".format(file, name, 
                                                              tensorDir))
                    return "Error getting domain input files"
            dp.tensorDir = tensorDir
        try:
            # initialize for that domain
            os.chdir(target)
            logger.info("Reading input files for domain {}".format(name))
            try:
                dp.sl = readPOSCAR()
                dp.rp = readPARAMETERS()
                interpretPARAMETERS(dp.rp, slab=dp.sl, silent=True)
                dp.sl.fullUpdate(rp)   #gets PARAMETERS data into slab
                dp.rp.fileLoaded["POSCAR"] = True
                if dp.sl.preprocessed:
                    dp.rp.SYMMETRY_FIND_ORI = False
                dp.rp.updateDerivedParams()
            except:
                logger.error("Error loading POSCAR and PARAMETERS for domain "
                             "{}".format(name))
                raise
            logger.info("Running initialization for domain {}".format(name))
            try:
                initialization(dp.sl, dp.rp)
            except:
                logger.error("Error running initialization for domain {}"
                             .format(name))
            rp.domainParams.append(dp)
        except:
            logger.error("Error while initializing domain {}".format(name))
            raise
        finally:
            os.chdir(home)