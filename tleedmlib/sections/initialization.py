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
import numpy as np

import tleedmlib as tl
from tleedmlib.base import angle
from tleedmlib.beamgen import runBeamGen
from tleedmlib.psgen import runPhaseshiftGen
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
    if os.path.isfile("_PHASESHIFTS"):
        try:
            (rp.phaseshifts_firstline, rp.phaseshifts,
                 newpsGen, newpsWrite) = readPHASESHIFTS(sl, rp, 
                                                    ignoreEnRange = subdomain)
        except:
            logger.warning("Found a _PHASESHIFTS file but could not "
                "read it. A new _PHASESHIFTS file will be generated."
                "The exception during read was: ", exc_info=True)
            rp.setHaltingLevel(1)
    if newpsGen:
        try:
            rundgrenpath = os.path.join('source', 'EEASiSSS.x')
            serneliuspath = os.path.join('source', 'seSernelius')
            logger.info("Generating phaseshifts data... ")
            (rp.phaseshifts_firstline, 
                        rp.phaseshifts) = runPhaseshiftGen(sl, rp,
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
    
    # write POSCAR_bulk_appended
    try:
        writeCONTCAR(sl.addBulkLayers(rp)[0], filename='POSCAR_bulk_appended')
    except:
        logger.warning("Exception occurred while writing POSCAR_bulk_appended")
    

    # generate beamlist
    logger.info("Generating _BEAMLIST...")
    try:
        bgenpath = os.path.join('source', 'beamgen3.out')
        runBeamGen(sl,rp,beamgensource = bgenpath)
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
    if rp.fileLoaded["IVBEAMS"] and not rp.ivbeams_sorted:
        rp.ivbeams = sortIVBEAMS(sl, rp)
        rp.ivbeams_sorted = True
    return 0

def init_domains(rp):
    """Runs an alternative initialization for the domain search. This will 
    include running the 'normal' initialization for each domain."""
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
        dp = tl.DomainParameters(target, home, name)
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
                dp.rp.workdir = home
                interpretPARAMETERS(dp.rp, slab=dp.sl, silent=True)
                dp.sl.fullUpdate(dp.rp)   #gets PARAMETERS data into slab
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
                initialization(dp.sl, dp.rp, subdomain=True)
            except:
                logger.error("Error running initialization for domain {}"
                             .format(name))
                raise
            rp.domainParams.append(dp)
        except:
            logger.error("Error while initializing domain {}".format(name))
            raise
        finally:
            os.chdir(home)
    if len(rp.domainParams) < len(rp.DOMAINS):
        return "Failed to read domain parameters"
    # check whether bulk unit cells match
    logger.info("Starting domain consistency check...")
    bulkuc0 = np.transpose(rp.domainParams[0].sl.bulkslab.ucell[:2,:2])
    for dp in rp.domainParams[1:]:
        bulkuc = np.transpose(dp.sl.bulkslab.ucell[:2,:2])
        if np.all(abs(bulkuc-bulkuc0) < 1e-4):
            continue
        # if the unit cells don't match right away, try if rotation matches
        found = False
        for i in [0,1]:
            ang = angle(bulkuc0[i], bulkuc[0])
            rotm = np.array([[np.cos(ang),-np.sin(ang)],
                             [np.sin(ang),np.cos(ang)]])
            rotuc = np.dot(rotm, bulkuc)
            if np.all(abs(bulkuc0-rotuc) < 1e-4):
                logger.info("Bulk unit cells of domain {0} and domain {1} are "
                    "mismatched, but can be matched by rotating domain {1}."
                    .format(rp.domainParams[0].name, dp.name))
                rotm3d = np.array([[0,0,0],[0,0,0],[0,0,1]])
                rotm3d[:2,:2] = rotm
                dp.sl.ucell = np.dot(rotm3d, dp.sl.ucell)
                dp.sl.getCartesianCoordinates()
                dp.sl.bulkslab.ucell = np.dot(rotm3d, dp.sl.bulkslab.ucell)
                dp.sl.bulkslab.getCartesianCoordinates()
                found = True
                break
        if found:
            continue
        logger.error("Bulk unit cells of domain {0} and domain {1} are "
                    "mismatched, and cannot be matched by rotation. Domain "
                    "search cannot proceed. Execution will stop."
                    .format(rp.domainParams[0].name, dp.name))
        rp.setHaltingLevel(3)
        return 0
    logger.debug("Domain bulk unit cells are compatible.")
    uc0 = np.transpose(rp.domainParams[0].sl.ucell[:2,:2])
    largestDomain = rp.domainParams[0]
    allMatched = all([np.all(abs(np.transpose(dp.sl.ucell[:2,:2])-uc0) < 1e-4) 
                      for dp in rp.domainParams[1:]])
    if allMatched:
        logger.debug("Domain surface unit cells are matched.")
    else:
        maxArea = abs(np.linalg.det(uc0))
        for dp in rp.domainParams[1:]:
            uc = np.transpose(dp.sl.ucell[:2,:2])
            if abs(np.linalg.det(uc)) > maxArea:
                maxArea = abs(np.linalg.det(uc))
                largestDomain = dp
        uc0 = np.transpose(largestDomain.sl.ucell[:2,:2])
        for dp in [p for p in rp.domainParams if p != largestDomain]:
            uc = np.transpose(dp.sl.ucell[:2,:2])
            if not np.all(abs(uc-uc0) < 1e-4):
                dp.refcalcRequired = True
                trans = uc0 * np.linalg.inv(uc)
                if np.any(abs(trans - np.round(trans)) > 1e-4):
                    logger.error("Surface unit cell of domain {0} cannot be "
                        "transformed to the largest surface unit cell (domain "
                        "{1}) by an integer transformation. Execution will "
                        "stop. Please supply all domain structures as "
                        "matching supercells.".format(dp.name, 
                                                      largestDomain.name))
                    rp.setHaltingLevel(3)
                    return 0
        logger.info("Domain surface unit cells are mismatched, but can be "
                    "matched by integer transformations.")
    # store some information about the supercell in rp:
    rp.SUPERLATTICE = copy.copy(largestDomain.rp.SUPERLATTICE)
    rp.pseudoSlab = tl.Slab()
    rp.pseudoSlab.ucell = copy.copy(largestDomain.sl.ucell)
    rp.pseudoSlab.bulkslab = tl.Slab()
    rp.pseudoSlab.bulkslab.ucell = copy.copy(largestDomain.sl.bulkslab.ucell)
    # run beamgen for the whole system
    logger.info("Generating _BEAMLIST...")
    try:
        bgenpath = os.path.join('source', 'beamgen3.out')
        runBeamGen(rp.pseudoSlab, rp, beamgensource = bgenpath, domains = True)
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
    # if EXPBEAMS was loaded, it hasn't been checked yet - check now
    if rp.fileLoaded["EXPBEAMS"]:
        checkEXPBEAMS(None, rp, domains=True)
    # write and sort IVBEAMS
    if not rp.fileLoaded["IVBEAMS"]:
        try:
            rp.ivbeams = writeIVBEAMS(None, rp, domains=True)
            rp.ivbeams_sorted = False
            rp.fileLoaded["IVBEAMS"] = True
            rp.manifest.append("IVBEAMS")
        except:
            logger.error("Error while writing IVBEAMS file based on "
                          "EXPBEAMS data.")
            raise
    if not rp.ivbeams_sorted:
        rp.ivbeams = sortIVBEAMS(None, rp, domains=True)
        rp.ivbeams_sorted = True

    rp.updateDerivedParams()
    rp.LMAX = max([dp.rp.LMAX for dp in rp.domainParams])
    for dp in rp.domainParams:
        if dp.refcalcRequired:
            continue
        cmessage = ("Reference calculation required for domain {}: "
                    .format(dp.name))
        # check energies
        if (rp.THEO_ENERGIES[0] < dp.rp.THEO_ENERGIES[0] or
               rp.THEO_ENERGIES[1] > dp.rp.THEO_ENERGIES[1] or
               rp.THEO_ENERGIES[2] != dp.rp.THEO_ENERGIES[2] or
               (rp.THEO_ENERGIES[0] % rp.THEO_ENERGIES[2]
                 != dp.rp.THEO_ENERGIES[0] % dp.rp.THEO_ENERGIES[2])):
            logger.info(cmessage+"Energy range is mismatched.")
            dp.refcalcRequired = True
            dp.rp.THEO_ENERGIES = copy.copy(rp.THEO_ENERGIES)
            continue
        # check LMAX
        if rp.LMAX != dp.rp.LMAX:
            logger.info(cmessage+"LMAX is mismatched.")
            dp.refcalcRequired = True
            dp.rp.LMAX = rp.LMAX
        # check beam incidence
        if rp.THETA != dp.rp.THETA or rp.PHI != dp.rp.PHI:
            logger.info(cmessage+"BEAM_INCIDENCE is mismatched.")
            dp.refcalcRequired = True
            dp.rp.THETA = rp.THETA
            dp.rp.PHI = rp.PHI
        # check IVBEAMS
        if not dp.rp.fileLoaded["IVBEAMS"]:
            logger.info(cmessage+"No IVBEAMS file loaded")
            dp.refcalcRequired = True
            dp.rp.ivbeams = copy.deepcopy(rp.ivbeams)
            continue
        if (len(rp.ivbeams) != len(dp.rp.ivbeams)
                or not all([dp.rp.ivbeams[i].isEqual(rp.ivbeams[i]) 
                            for i in range(0,len(rp.ivbeams))])):
            logger.info(cmessage+"IVBEAMS file mismatched with supercell.")
            dp.refcalcRequired = True
            dp.rp.ivbeams = copy.deepcopy(rp.ivbeams)
            continue
    
    rr = [dp for dp in rp.domainParams if dp.refcalcRequired]
    if rr:
        logger.info("The following domains require new reference "
                    "calculations: "+", ".join([d.name for d in rr]))
        logger.warning("Automatic reference calculations are not supported in "
                       "this version of ViPErLEED. Please manually execute "
                       "appropriate reference calculations, and start a "
                       "domain search based on the Tensors.zip files.")
        rp.setHaltingLevel(3)
