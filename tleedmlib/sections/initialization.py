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

import viperleed.tleedmlib as tl
from viperleed.tleedmlib.base import angle, rotation_matrix
from viperleed.tleedmlib.beamgen import runBeamGen
from viperleed.tleedmlib.psgen import (runPhaseshiftGen, runPhaseshiftGen_old)
from viperleed.tleedmlib.files.poscar import readPOSCAR, writePOSCAR
from viperleed.tleedmlib.files.vibrocc import readVIBROCC
from viperleed.tleedmlib.files.parameters import (
    readPARAMETERS, interpretPARAMETERS, modifyPARAMETERS)
from viperleed.tleedmlib.files.phaseshifts import (
    readPHASESHIFTS, writePHASESHIFTS, plot_phaseshifts)
from viperleed.tleedmlib.files.beams import (
    readOUTBEAMS, readBEAMLIST, checkEXPBEAMS, readIVBEAMS, sortIVBEAMS,
    writeIVBEAMS)
from viperleed.tleedmlib.files.patterninfo import writePatternInfo

logger = logging.getLogger("tleedm.initialization")


def initialization(sl, rp, subdomain=False):
    """Runs the initialization."""
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
                    rp.expbeams = readOUTBEAMS(filename=fn, enrange=er)
                    if len(rp.expbeams) > 0:
                        rp.fileLoaded["EXPBEAMS"] = True
                    else:
                        logger.error("Error reading "+fn+": No data was read.")
                except Exception:
                    logger.error("Error while reading file "+fn, exc_info=True)
    rp.initTheoEnergies()  # may be initialized based on exp. beams

    if (rp.DOMAINS or rp.domainParams) and not subdomain:
        init_domains(rp)
        return

    # check whether PHASESHIFTS are present & consistent:
    newpsGen, newpsWrite = True, True
    # True: new phaseshifts need to be generated/written
    if os.path.isfile("PHASESHIFTS") or os.path.isfile("_PHASESHIFTS"):
        try:
            (rp.phaseshifts_firstline, rp.phaseshifts,
             newpsGen, newpsWrite) = readPHASESHIFTS(sl, rp,
                                                     ignoreEnRange=subdomain)
        except Exception:
            logger.warning(
                "Found a PHASESHIFTS file but could not "
                "read it. A new PHASESHIFTS file will be generated."
                "The exception during read was: ", exc_info=True)
            rp.setHaltingLevel(1)
    if newpsGen:
        try:
            rundgrenpath = os.path.join('tensorleed', 'EEASiSSS.x')
            serneliuspath = os.path.join('tensorleed', 'seSernelius')
            logger.info("Generating phaseshifts data... ")
            if rp.PHASESHIFTS_CALC_OLD:
                (rp.phaseshifts_firstline,
                 rp.phaseshifts) = runPhaseshiftGen_old(sl, rp,
                                                    psgensource=rundgrenpath,
                                                    excosource=serneliuspath)
            else:
                (rp.phaseshifts_firstline,
                 rp.phaseshifts) = runPhaseshiftGen(sl, rp,
                                                        psgensource=rundgrenpath)
            logger.debug("Finished generating phaseshift data")
        except Exception:
            logger.error("Exception while calling phaseshiftgen: ")
            raise
    if newpsWrite:
        try:
            writePHASESHIFTS(rp.phaseshifts_firstline, rp.phaseshifts)
        except Exception:
            logger.error("Exception during writePHASESHIFTS: ")
            raise
    rp.fileLoaded["PHASESHIFTS"] = True
    rp.updateDerivedParams()
    rp.manifest.append("PHASESHIFTS")
    try:
        plot_phaseshifts(sl, rp)
    except Exception:
        logger.warning("Failed to plot phaseshifts", exc_info=rp.LOG_DEBUG)

    # if necessary, run findSymmetry:
    if sl.planegroup == "unknown":
        tl.symmetry.findSymmetry(sl, rp)
        tl.symmetry.enforceSymmetry(sl, rp)

    # check whether the slab unit cell is minimal:
    changecell, mincell = sl.getMinUnitCell(rp)
    transform = np.dot(np.transpose(sl.ucell[:2, :2]),
                       np.linalg.inv(mincell)).round()
    ws = tl.leedbase.writeWoodsNotation(transform)
    if changecell and np.isclose(rp.SYMMETRY_CELL_TRANSFORM,
                                 np.identity(2)).all():
        if ws:
            ws = "= "+ws
        else:
            ws = "M = {} {}, {} {}".format(*[x for y in transform.astype(int)
                                             for x in y])
        ssl = sl.makeSymBaseSlab(rp, transform=transform)
        if subdomain:
            rp.SYMMETRY_CELL_TRANSFORM = transform
            logger.info("Found SYMMETRY_CELL_TRANSFORM "+ws)
            sl.symbaseslab = ssl
            if "M = " not in ws:
                modifyPARAMETERS(rp, "SYMMETRY_CELL_TRANSFORM", new=ws)
            else:
                modifyPARAMETERS(
                    rp, "SYMMETRY_CELL_TRANSFORM",
                    new=("SYMMETRY_CELL_TRANSFORM M = "
                         "{:.0f} {:.0f}, {:.0f} {:.0f}".format(
                                     *[x for y in transform for x in y])),
                    include_left=True)
        else:
            logger.warning(
                "POSCAR unit cell is not minimal (supercell {}). "
                "A minimal POSCAR will be written as POSCAR_mincell for "
                "information. Consider calculating with POSCAR_mincell, "
                "or setting SYMMETRY_CELL_TRANSFORM to conserve "
                "translational symmetry.".format(ws))
            writePOSCAR(ssl, filename="POSCAR_mincell")
            rp.setHaltingLevel(1)
    elif not np.isclose(rp.SYMMETRY_CELL_TRANSFORM, np.identity(2)).all():
        if not np.isclose(rp.SYMMETRY_CELL_TRANSFORM, transform).all():
            logger.warning("SYMMETRY_CELL_TRANSFORM parameter differs from "
                           "automatically detected supercell "+ws)
            rp.setHaltingLevel(1)
        if sl.symbaseslab is None:
            sl.symbaseslab = sl.makeSymBaseSlab(rp)
    if sl.symbaseslab is not None:
        logger.info("A symmetry cell transformation was found. Re-running "
                    "slab symmetry search using base unit cell...")
        tl.symmetry.getSymBaseSymmetry(sl, rp)
        try:
            writePOSCAR(sl.symbaseslab, filename='POSCAR_mincell',
                         comments='all')
        except Exception:
            logger.warning("Exception occurred while writing POSCAR_mincell")

    # generate new POSCAR
    tmpslab = copy.deepcopy(sl)
    tmpslab.sortOriginal()
    try:
        writePOSCAR(tmpslab, filename='POSCAR', comments='all')
    except Exception:
        logger.error("Exception occurred while writing new POSCAR")
        raise
    rp.manifest.append('POSCAR')
    # generate POSCAR_oricell
    tmpslab.revertUnitCell()
    try:
        writePOSCAR(tmpslab, filename='POSCAR_oricell', comments='nodir')
    except Exception:
        logger.error("Exception occurred while writing POSCAR_oricell, "
                     "execution will continue...")

    if rp.BULK_LIKE_BELOW > 0:
        if rp.BULK_REPEAT is not None:
            logger.warning("Both BULK_LIKE_BELOW and BULK_REPEAT are defined."
                           "BULK_LIKE_BELOW will be ignored in favour of the "
                           "explicitly defined bulk repeat vector.")
        else:
            cvec, cuts = sl.detectBulk(rp)
            rp.BULK_REPEAT = cvec
            vec_str = "[{:.5f} {:.5f} {:.5f}]".format(*rp.BULK_REPEAT)
            modifyPARAMETERS(rp, "BULK_REPEAT", vec_str,
                             comment="Automatically detected repeat vector")
            logger.info("Detected bulk repeat vector: " + vec_str)
            layer_cuts = sl.createLayers(rp, bulk_cuts=cuts)
            modifyPARAMETERS(rp, "LAYER_CUTS", " ".join(["{:.4f}".format(f)
                                                         for f in layer_cuts]))
            rp.N_BULK_LAYERS = len(cuts)
            modifyPARAMETERS(rp, "N_BULK_LAYERS", str(len(cuts)))
        modifyPARAMETERS(rp, "BULK_LIKE_BELOW", new="")

    # create bulk slab:
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    bsl = sl.bulkslab

    # try identifying bulk repeat:
    if rp.BULK_REPEAT is None:
        rvec = sl.getBulkRepeat(rp)
        if rvec is not None:
            rp.BULK_REPEAT = rvec
            vec_str = "[{:.5f} {:.5f} {:.5f}]".format(*rp.BULK_REPEAT)
            modifyPARAMETERS(rp, "BULK_REPEAT", vec_str,
                             comment="Automatically detected repeat vector")
            logger.info("Detected bulk repeat vector: " + vec_str)
            # update bulk slab vector
            sl.bulkslab.getCartesianCoordinates()
            sl.bulkslab.ucell[:, 2] = np.copy(rvec)
            sl.bulkslab.collapseCartesianCoordinates()
    if rp.BULK_REPEAT is None:
        # failed to detect repeat vector, use fixed distance instead
        blayers = [lay for lay in sl.layers if lay.isBulk]
        # assume that interlayer vector from bottom non-bulk to top bulk layer
        #   is the same as between bulk units
        # save BULK_REPEAT value for later runs, in case atom above moves
        if rp.N_BULK_LAYERS == 2:
            rp.BULK_REPEAT = (blayers[1].cartbotz
                              - sl.layers[blayers[0].num-1].cartbotz)
        else:
            rp.BULK_REPEAT = (blayers[0].cartbotz
                              - sl.layers[blayers[0].num-1].cartbotz)
        modifyPARAMETERS(rp, "BULK_REPEAT", "{:.5f}".format(rp.BULK_REPEAT),
                         comment=("Automatically detected spacing. "
                                  "Check POSCAR_bulk."))
        logger.warning(
            "The BULK_REPEAT parameter was undefined, which may lead to "
            "unintended changes in the bulk unit cell during optimization if "
            "the lowest non-bulk atom moves.\n"
            "# Automatic detection of a bulk repeat vector failed, possibly "
            "because not enough bulk-like layers were found.\n"
            "# The BULK_REPEAT vector is assumed to be parallel to the POSCAR "
            "c vector. Check POSCAR_bulk and the BULK_REPEAT thickness "
            "written to the PARAMETERS file.")
        rp.checklist.append("Check bulk repeat vector in PARAMETERS")
        rp.setHaltingLevel(2)

    if bsl.planegroup == "unknown":
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
            if ws:
                logger.info("Found SUPERLATTICE = "+ws)
            else:
                logger.info("Found SUPERLATTICE M = {} {}, {} {}".format(
                        *[x for y in rp.SUPERLATTICE.astype(int) for x in y]))

        # bulk plane group detection:
        logger.info("Starting bulk symmetry search...")
        tl.symmetry.findSymmetry(bsl, rp, bulk=True, output=False)
        bsl.revertUnitCell()  # keep origin matched with main slab
        logger.info("Found bulk plane group: "+bsl.foundplanegroup)
        tl.symmetry.findBulkSymmetry(bsl, rp)

        # write POSCAR_bulk
        bsl = copy.deepcopy(sl.bulkslab)
        bsl.sortOriginal()
        try:
            writePOSCAR(bsl, filename='POSCAR_bulk', comments='bulk')
        except Exception:
            logger.error("Exception occurred while writing POSCAR_bulk")
            raise

    # write POSCAR_bulk_appended
    n = 1
    if len(bsl.sublayers) <= len(bsl.elements)*2:
        n += 1
        if len(bsl.sublayers) <= len(bsl.elements):
            n += 1
    try:
        writePOSCAR(sl.addBulkLayers(rp, n=n)[0],
                     filename='POSCAR_bulk_appended')
    except Exception:
        logger.warning("Exception occurred while writing POSCAR_bulk_appended")

    # generate beamlist
    logger.info("Generating BEAMLIST...")
    try:
        bgenpath = os.path.join('tensorleed', 'beamgen3.out')
        runBeamGen(sl, rp, beamgensource=bgenpath)
        # this does NOT read the resulting file!
    except Exception:
        logger.error("Exception occurred while calling beamgen.")
        raise
    try:
        rp.beamlist = readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except Exception:
        logger.error("Error while reading required file BEAMLIST")
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
            except Exception:
                logger.error("Error while writing IVBEAMS file based on "
                             "EXPBEAMS data.")
                raise
    if rp.fileLoaded["IVBEAMS"] and not rp.ivbeams_sorted:
        rp.ivbeams = sortIVBEAMS(sl, rp)
        rp.ivbeams_sorted = True

    # Create directory compile_logs in which logs from compilation will be saved
    make_compile_logs_dir(rp, logger)

    # At the end of initialization preserve the original input files
    preserve_original_input(rp, logger)
    return


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
                rp.expbeams = readOUTBEAMS(filename=fn, enrange=er)
                if len(rp.expbeams) > 0:
                    rp.fileLoaded["EXPBEAMS"] = True
                else:
                    logger.error("Error reading "+fn+": No data was read.")
            except Exception:
                logger.error("Error while reading file "+fn, exc_info=True)
    rp.initTheoEnergies()  # may be initialized based on exp. beams
    if len(rp.DOMAINS) < 2:
        logger.error("A domain search was defined, but less than two domains "
                     "are defined. Execution will stop.")
        rp.setHaltingLevel(3)
        return
    checkFiles = ["POSCAR", "PARAMETERS", "VIBROCC", "PHASESHIFTS"]
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
                    tl.leedbase.getTensors(tensorIndex, basedir=path,
                                           targetdir=target)
                except Exception as e:
                    tensorIndex = 0
                    logger.warning("Error fetching Tensors: " + str(e))
            if tensorIndex != 0:
                tensorDir = os.path.join(target, "Tensors",
                                         "Tensors_"+str(tensorIndex).zfill(3))
                for file in (checkFiles + ["IVBEAMS"]):
                    if os.path.isfile(os.path.join(tensorDir, file)):
                        shutil.copy2(os.path.join(tensorDir, file),
                                     os.path.join(target, file))
                    else:
                        logger.warning("Input file {} is missing in Tensors "
                                       "directory. A new reference "
                                       "calculation is required.".format(file))
                        tensorIndex = 0
                        break
            if tensorIndex != 0:
                dp.tensorDir = tensorDir
            else:       # no usable tensors in that dir; get input
                dp.refcalcRequired = True
                logger.info("No previous Tensors found, reference calculation "
                            "is required.")
                for file in checkFiles:
                    if os.path.isfile(os.path.join(path, file)):
                        try:
                            shutil.copy(os.path.join(path, file),
                                        os.path.join(target, file))
                        except Exception:
                            if file != "PHASESHIFTS":
                                logger.error(
                                    "Error copying required file {} for "
                                    "domain {} from origin folder {}"
                                    .format(file, name, path))
                                raise RuntimeError("Error getting domain "
                                                   "input files")
                    elif file != "PHASESHIFTS":
                        logger.error("Required file {} for domain {} not "
                                     "found in origin folder {}"
                                     .format(file, name, path))
                        raise RuntimeError("Error getting domain input files")
        elif os.path.isfile(path):
            try:
                tensorIndex = tl.leedbase.getMaxTensorIndex(target)
            except Exception:
                tensorIndex = 0
            tensorDir = os.path.join(target, "Tensors",
                                     "Tensors_"+str(tensorIndex+1).zfill(3))
            try:
                os.makedirs(os.path.join(target, "Tensors", tensorDir),
                            exist_ok=True)
            except Exception:
                raise
            try:
                shutil.unpack_archive(path, tensorDir)
            except Exception:
                logger.error("Failed to unpack Tensors for domain {} from "
                             "file {}".format(name, path))
                raise RuntimeError("Error getting domain input files")
            for file in (checkFiles + ["IVBEAMS"]):
                if os.path.isfile(os.path.join(tensorDir, file)):
                    shutil.copy2(os.path.join(tensorDir, file),
                                 os.path.join(target, file))
                else:
                    logger.error("Required file {} for domain {} not found in "
                                 "Tensor directory {}".format(file, name,
                                                              tensorDir))
                    raise RuntimeError("Error getting domain input files")
            dp.tensorDir = tensorDir
        try:
            # initialize for that domain
            os.chdir(target)
            logger.info("Reading input files for domain {}".format(name))
            try:
                dp.sl = readPOSCAR()
                dp.rp = readPARAMETERS()
                dp.rp.workdir = home
                dp.rp.sourcedir = rp.sourcedir
                dp.rp.timestamp = rp.timestamp
                interpretPARAMETERS(dp.rp, slab=dp.sl, silent=True)
                dp.sl.fullUpdate(dp.rp)   # gets PARAMETERS data into slab
                dp.rp.fileLoaded["POSCAR"] = True
                if dp.sl.preprocessed:
                    dp.rp.SYMMETRY_FIND_ORI = False
                dp.rp.updateDerivedParams()
                try:
                    readVIBROCC(dp.rp, dp.sl)
                    dp.rp.fileLoaded["VIBROCC"] = True
                except Exception:
                    logger.error("Error while reading required file VIBROCC")
                    raise
                dp.sl.fullUpdate(dp.rp)
                try:
                    dp.rp.ivbeams = readIVBEAMS()
                    dp.rp.ivbeams_sorted = False
                    dp.rp.fileLoaded["IVBEAMS"] = True
                except FileNotFoundError:
                    pass
                except Exception:
                    logger.error("Error while reading IVBEAMS for domain {}"
                                 .format(name))
            except Exception:
                logger.error("Error loading POSCAR and PARAMETERS for domain "
                             "{}".format(name))
                raise
            logger.info("Running initialization for domain {}".format(name))
            try:
                initialization(dp.sl, dp.rp, subdomain=True)
            except Exception:
                logger.error("Error running initialization for domain {}"
                             .format(name))
                raise
            rp.domainParams.append(dp)
        except Exception:
            logger.error("Error while initializing domain {}".format(name))
            raise
        finally:
            os.chdir(home)
    if len(rp.domainParams) < len(rp.DOMAINS):
        raise RuntimeError("Failed to read domain parameters")
    # check whether bulk unit cells match
    logger.info("Starting domain consistency check...")
    bulkuc0 = rp.domainParams[0].sl.bulkslab.ucell[:2, :2].T
    eps = 1e-4
    for dp in rp.domainParams[1:]:
        bulkuc = dp.sl.bulkslab.ucell[:2, :2].T
        if np.all(abs(bulkuc-bulkuc0) < eps):
            continue
        # if the unit cells don't match right away, try if rotation matches
        if (all([abs(np.linalg.norm(bulkuc0[i]) - np.linalg.norm(bulkuc[i]))
                 < eps for i in range(0, 2)])
                and abs(angle(bulkuc[0], bulkuc[1])
                        - angle(bulkuc0[0], bulkuc0[1])) < eps):
            logger.info("Bulk unit cells of domain {0} and domain {1} are "
                        "mismatched, but can be matched by rotating domain "
                        "{1}.".format(rp.domainParams[0].name, dp.name))
            ang = angle(bulkuc0[0], bulkuc[0])
            rotm = rotation_matrix(ang, dim=3)
            rotm_t = np.transpose(rotm)
            # dp.sl.ucell = np.transpose(np.dot(np.transpose(dp.sl.ucell),
            #                                   rotuc))
            dp.sl.ucell = np.dot(rotm_t, dp.sl.ucell)
            dp.sl.getCartesianCoordinates()
            dp.sl.bulkslab.ucell = np.dot(rotm_t, dp.sl.bulkslab.ucell)
            dp.sl.bulkslab.getCartesianCoordinates()
        else:
            logger.error("Bulk unit cells of domain {0} and domain {1} are "
                         "mismatched, and cannot be matched by rotation. "
                         "Domain search cannot proceed. Execution will stop."
                         .format(rp.domainParams[0].name, dp.name))
            rp.setHaltingLevel(3)
            return
    logger.debug("Domain bulk unit cells are compatible.")
    uc0 = rp.domainParams[0].sl.ucell[:2, :2].T
    largestDomain = rp.domainParams[0]
    allMatched = all([np.all(abs(dp.sl.ucell[:2, :2].T - uc0) < 1e-4)
                      for dp in rp.domainParams[1:]])
    supercellRequired = []
    if allMatched:
        logger.debug("Domain surface unit cells are matched.")
    else:
        maxArea = abs(np.linalg.det(uc0))
        for dp in rp.domainParams[1:]:
            uc = dp.sl.ucell[:2, :2].T
            if abs(np.linalg.det(uc)) > maxArea:
                maxArea = abs(np.linalg.det(uc))
                largestDomain = dp
        uc0 = largestDomain.sl.ucell[:2, :2].T
        for dp in [p for p in rp.domainParams if p != largestDomain]:
            uc = dp.sl.ucell[:2, :2].T
            if not np.all(abs(uc-uc0) < 1e-4):
                dp.refcalcRequired = True
                trans = np.dot(uc0, np.linalg.inv(uc))
                if np.any(abs(trans - np.round(trans)) > 1e-4):
                    logger.error(
                        "Surface unit cell of domain {0} cannot be "
                        "transformed to the largest surface unit cell (domain "
                        "{1}) by an integer transformation. Execution will "
                        "stop. Please supply all domain structures as "
                        "matching supercells.".format(dp.name,
                                                      largestDomain.name))
                    rp.setHaltingLevel(3)
                    return
                else:
                    supercellRequired.append(dp)
                    oldslab = dp.sl
                    dp.sl = dp.sl.makeSupercell(np.round(trans))
                    dp.rp.SUPERLATTICE = copy.copy(largestDomain
                                                   .rp.SUPERLATTICE)
                    dp.sl.symbaseslab = oldslab
                    dp.rp.SYMMETRY_CELL_TRANSFORM = trans
                    ws = tl.leedbase.writeWoodsNotation(trans)
                    if ws:
                        modifyPARAMETERS(dp.rp, "SYMMETRY_CELL_TRANSFORM",
                                         new=ws, path=dp.workdir)
                    else:
                        modifyPARAMETERS(
                            dp.rp, "SYMMETRY_CELL_TRANSFORM",
                            new=("SYMMETRY_CELL_TRANSFORM M = "
                                 "{:.0f} {:.0f}, {:.0f} {:.0f}".format(
                                       *[x for y in trans for x in y])),
                            path=dp.workdir, include_left=True)
        logger.info("Domain surface unit cells are mismatched, but can be "
                    "matched by integer transformations.")
    # store some information about the supercell in rp:
    rp.SUPERLATTICE = copy.copy(largestDomain.rp.SUPERLATTICE)
    rp.pseudoSlab = tl.Slab()
    rp.pseudoSlab.ucell = copy.copy(largestDomain.sl.ucell)
    rp.pseudoSlab.bulkslab = tl.Slab()
    rp.pseudoSlab.bulkslab.ucell = copy.copy(largestDomain.sl.bulkslab.ucell)
    # run beamgen for the whole system
    logger.info("Generating BEAMLIST...")
    try:
        bgenpath = os.path.join('tensorleed', 'beamgen3.out')
        runBeamGen(rp.pseudoSlab, rp, beamgensource=bgenpath, domains=True)
        # this does NOT read the resulting file!
    except Exception:
        logger.error("Exception occurred while calling beamgen.")
        raise
    try:
        rp.beamlist = readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except Exception:
        logger.error("Error while reading required file BEAMLIST")
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
        except Exception:
            logger.error("Error while writing IVBEAMS file based on "
                         "EXPBEAMS data.")
            raise
    if not rp.ivbeams_sorted:
        rp.ivbeams = sortIVBEAMS(None, rp)
        rp.ivbeams_sorted = True

    rp.updateDerivedParams()
    if rp.LMAX[1] == 0:
        rp.LMAX[1] = max([dp.rp.LMAX[1] for dp in rp.domainParams])
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
            continue
        # check LMAX - should be obsolete since TensErLEED version 1.6
        if rp.LMAX[1] != dp.rp.LMAX[1] and rp.TL_VERSION <= 1.6:
            logger.info(cmessage+"LMAX is mismatched.")
            dp.refcalcRequired = True
        # check beam incidence
        if rp.THETA != dp.rp.THETA or rp.PHI != dp.rp.PHI:
            logger.info(cmessage+"BEAM_INCIDENCE is mismatched.")
            dp.refcalcRequired = True
        # check IVBEAMS
        if not dp.rp.fileLoaded["IVBEAMS"]:
            logger.info(cmessage+"No IVBEAMS file loaded")
            dp.refcalcRequired = True
            continue
        if (len(rp.ivbeams) != len(dp.rp.ivbeams)
                or not all([dp.rp.ivbeams[i].isEqual(rp.ivbeams[i])
                            for i in range(0, len(rp.ivbeams))])):
            logger.info(cmessage+"IVBEAMS file mismatched with supercell.")
            dp.refcalcRequired = True
            continue

    rr = [dp for dp in rp.domainParams if dp.refcalcRequired]
    if rr:
        logger.info("The following domains require new reference "
                    "calculations: "+", ".join([d.name for d in rr]))
        for dp in rp.domainParams:
            for var in ["THEO_ENERGIES", "THETA", "PHI", "N_CORES", "ivbeams",
                        "sourcedir"]:
                setattr(dp.rp, var, copy.deepcopy(getattr(rp, var)))
            if rp.TL_VERSION <= 1.6:  # not required since TensErLEED v1.61
                dp.rp.LMAX[1] = rp.LMAX[1]

    # repeat initialization for all slabs that require a supercell
    for dp in supercellRequired:
        logger.info("Re-running intialization with supercell slab for domain "
                    "{}".format(dp.name))
        try:
            os.chdir(dp.workdir)
            dp.sl.resetSymmetry()
            dp.rp.SYMMETRY_FIND_ORI = True
            initialization(dp.sl, dp.rp, subdomain=True)
        except Exception:
            logger.error("Error while re-initializing domain {}".format(name))
            raise
        finally:
            os.chdir(home)

    if 4 not in rp.RUN and 1 not in rp.RUN and rr:
        logger.error(
            "Some domains require new reference calculations before "
            "a domain search can be executed. Please either manually "
            "execute appropriate reference calculations, or set RUN = 4")
        rp.setHaltingLevel(3)
        return

    while 4 in rp.RUN:
        if rr:
            rp.RUN.insert(rp.RUN.index(4), 1)
        rp.RUN.insert(rp.RUN.index(4), 2)
        rp.RUN.insert(rp.RUN.index(4), 3)
        rp.RUN.remove(4)
    return


def preserve_original_input(rp, init_logger, path=""):
    """
    Creates directory original_inputs and copies input files there for preservation.
    """
    # Folder name "original_inputs" hardcoded. If changed, change also in cleanup.py !
    folder_name = "original_inputs"
    if not path:
        path = "."

    # makes orig_inputs directory
    try:
        orig_inputs_path = os.path.join(path, folder_name)
        os.makedirs(orig_inputs_path, exist_ok=True)
    except Exception:
        logger.warning("Could not create directory {}".format(folder_name))
        rp.setHaltingLevel(1)

    # copy all files to orig_inputs that were used as original input
    # !!! TODO: rp.fileLoaded contains only files that were needed so far,
    #  and EXPBEAMS is missing the .csv extension -> probably better to list
    #  all input files explicitly here
    for file in rp.fileLoaded:
        if rp.fileLoaded[file]:
            # copy to original input
            try:
                if os.path.isfile(file):
                    shutil.copy2(os.path.join(path, file), orig_inputs_path)
                else:
                    raise FileNotFoundError
            except Exception:
                init_logger.warning("Could not copy file {} to ".format(file)
                                    + folder_name)
                rp.setHaltingLevel(1)
    return

def make_compile_logs_dir(rp, init_logger, path=""):
    """
    Creates directory compile_logs in which logs from compilation will be saved.
    """
    folder_name = "compile_logs"
    if not path:
        path = "."

    # makes compile_logs directory
    try:
        os.makedirs(os.path.join(path, folder_name), exist_ok=True)
    except Exception:
        logger.warning("Could not create directory {}".format(folder_name))
        rp.setHaltingLevel(1)
    return