# -*- coding: utf-8 -*-
"""
Created on Aug 11 2020

@author: Florian Kraushofer
@author: Alexander M. Imre
@author: Michele Riva

Tensor LEED Manager section Initialization
"""

import copy
import logging
import os
from pathlib import Path
import shutil
from zipfile import ZipFile

import numpy as np

from viperleed.tleedmlib import leedbase
from viperleed.tleedmlib import symmetry as tl_symmetry
from viperleed.tleedmlib.base import angle, rotation_matrix
from viperleed.tleedmlib.beamgen import calc_and_write_beamlist
from viperleed.tleedmlib.classes.slab import Slab
from viperleed.tleedmlib.classes.rparams import DomainParameters
from viperleed.tleedmlib.files import beams as tl_beams, parameters
from viperleed.tleedmlib.files import patterninfo, phaseshifts, poscar, vibrocc
from viperleed.tleedmlib.files.woods_notation import writeWoodsNotation
from viperleed.tleedmlib.psgen import runPhaseshiftGen, runPhaseshiftGen_old
from viperleed.tleedmlib.sections._sections import (ALL_INPUT_FILES,
                                                    EXPBEAMS_NAMES)


logger = logging.getLogger("tleedm.initialization")

ORIGINAL_INPUTS_DIR_NAME = 'original_inputs'


def initialization(sl, rp, subdomain=False):
    """Runs the initialization."""
    if not subdomain:
        rp.try_loading_expbeams_file()
    rp.initTheoEnergies()  # may be initialized based on exp. beams

    if (rp.DOMAINS or rp.domainParams) and not subdomain:
        init_domains(rp)
        return

    # if necessary, run findSymmetry:
    if sl.planegroup == "unknown":
        tl_symmetry.findSymmetry(sl, rp)
        tl_symmetry.enforceSymmetry(sl, rp)

    # check whether the slab unit cell is minimal:
    changecell, mincell = sl.getMinUnitCell(rp)
    transform = np.dot(np.transpose(sl.ucell[:2, :2]),
                       np.linalg.inv(mincell)).round()
    ws = writeWoodsNotation(transform)
    if ws:
        ws = f"= {ws}"
    else:
        ws = "M = {} {}, {} {}".format(*transform.astype(int).ravel())
    if changecell and np.isclose(rp.SYMMETRY_CELL_TRANSFORM,
                                 np.identity(2)).all():
        ssl = sl.makeSymBaseSlab(rp, transform=transform)
        if subdomain:
            rp.SYMMETRY_CELL_TRANSFORM = transform
            logger.info(f"Found SYMMETRY_CELL_TRANSFORM {ws}")
            sl.symbaseslab = ssl
            parameters.modify(rp, "SYMMETRY_CELL_TRANSFORM")                    #TODO: there should probably a comment in the PARAMETERS file?
        else:
            logger.warning(
                f"POSCAR unit cell is not minimal (supercell {ws}). "
                "A minimal POSCAR will be written as POSCAR_mincell for "
                "information. Consider calculating with POSCAR_mincell, "
                "or setting SYMMETRY_CELL_TRANSFORM to conserve "
                "translational symmetry."
                )
            poscar.write(ssl, filename='POSCAR_mincell')
            rp.setHaltingLevel(1)
    elif not np.isclose(rp.SYMMETRY_CELL_TRANSFORM, np.identity(2)).all():
        if not np.isclose(rp.SYMMETRY_CELL_TRANSFORM, transform).all():
            logger.warning("SYMMETRY_CELL_TRANSFORM parameter differs from "
                           f"automatically detected supercell {ws}")
            rp.setHaltingLevel(1)
        if sl.symbaseslab is None:
            sl.symbaseslab = sl.makeSymBaseSlab(rp)
    if sl.symbaseslab is not None:
        logger.info("A symmetry cell transformation was found. Re-running "
                    "slab symmetry search using base unit cell...")
        tl_symmetry.getSymBaseSymmetry(sl, rp)
        try:
            poscar.write(sl.symbaseslab, filename='POSCAR_mincell',
                         comments='all')
        except Exception:
            logger.warning("Exception occurred while writing POSCAR_mincell")

    # generate new POSCAR
    tmpslab = copy.deepcopy(sl)
    tmpslab.sortOriginal()
    try:
        poscar.write(tmpslab, filename='POSCAR', comments='all')
    except Exception:
        logger.error("Exception occurred while writing new POSCAR")
        raise
    rp.manifest.append('POSCAR')
    # generate POSCAR_oricell
    tmpslab.revertUnitCell()
    try:
        poscar.write(tmpslab, filename='POSCAR_oricell', comments='nodir')
    except Exception:
        logger.error("Exception occurred while writing POSCAR_oricell, "
                     "execution will continue...")

    if rp.BULK_LIKE_BELOW > 0:
        if rp.BULK_REPEAT is not None:
            logger.warning("Both BULK_LIKE_BELOW and BULK_REPEAT are defined."
                           "BULK_LIKE_BELOW will be ignored in favour of the "
                           "explicitly defined bulk repeat vector.")
        else:
            # The modifications to the PARAMETERS file below are currently not in the docs. Should we add that?
            cvec, cuts = sl.detectBulk(rp)
            rp.BULK_REPEAT = cvec
            vec_str = parameters.modify(
                rp, "BULK_REPEAT",
                comment="Automatically detected repeat vector"
                )
            logger.info(f"Detected bulk repeat vector: {vec_str}")
            new_cuts = sl.createLayers(rp, bulk_cuts=cuts)
            rp.LAYER_CUTS.update_from_sequence(new_cuts)
            parameters.modify(rp, "LAYER_CUTS")
            rp.N_BULK_LAYERS = len(cuts)
            parameters.modify(rp, "N_BULK_LAYERS")
        parameters.comment_out(rp, "BULK_LIKE_BELOW")

    # create bulk slab:
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    bsl = sl.bulkslab

    # try identifying bulk repeat:
    if rp.BULK_REPEAT is None:
        rvec = sl.getBulkRepeat(rp)
        if rvec is not None:
            rp.BULK_REPEAT = rvec
            vec_str = parameters.modify(
                rp, "BULK_REPEAT",
                comment="Automatically detected repeat vector"
                )
            logger.info(f"Detected bulk repeat vector: {vec_str}")
            # update bulk slab vector
            sl.bulkslab.getCartesianCoordinates()
            sl.bulkslab.ucell[:, 2] = np.copy(rvec)
            sl.bulkslab.collapseCartesianCoordinates()
    if rp.BULK_REPEAT is None:
        # failed to detect repeat vector, use fixed distance instead
        blayers = [lay for lay in sl.layers if lay.is_bulk]
        # assume that interlayer vector from bottom non-bulk to
        # top bulk layer is the same as between bulk units
        # save BULK_REPEAT value for later runs, in case atom above moves
        if rp.N_BULK_LAYERS == 2:
            rp.BULK_REPEAT = (blayers[1].cartbotz
                              - sl.layers[blayers[0].num-1].cartbotz)
        else:
            rp.BULK_REPEAT = (blayers[0].cartbotz
                              - sl.layers[blayers[0].num-1].cartbotz)
        parameters.modify(
            rp, "BULK_REPEAT",
            comment="Automatically detected spacing. Check POSCAR_bulk."
            )
        logger.warning(
            "The BULK_REPEAT parameter was undefined, which may lead to "
            "unintended changes in the bulk unit cell during optimization if "
            "the lowest non-bulk atom moves.\n"
            "# Automatic detection of a bulk repeat vector failed, possibly "
            "because not enough bulk-like layers were found.\n"
            "# The BULK_REPEAT vector is assumed to be parallel to the POSCAR "
            "c vector. Check POSCAR_bulk and the BULK_REPEAT thickness "
            "written to the PARAMETERS file."
            )
        rp.checklist.append("Check bulk repeat vector in PARAMETERS")
        rp.setHaltingLevel(2)

    if bsl.planegroup == "unknown":
        # find minimum in-plane unit cell for bulk:
        logger.info("Checking bulk unit cell...")
        changecell, mincell = bsl.getMinUnitCell(rp, warn_convention=True)
        if changecell:
            sl.changeBulkCell(rp, mincell)
            bsl = sl.bulkslab
        if not rp.superlattice_defined:
            ws = writeWoodsNotation(rp.SUPERLATTICE)                   # TODO: replace writeWoodsNotation with guilib functions
            superlattice = rp.SUPERLATTICE.astype(int)
            if ws:
                info = f"= {ws}"
            else:
                info = "M = {} {}, {} {}".format(*superlattice.ravel())
            logger.info(f"Found SUPERLATTICE {info}")

        # bulk plane group detection:
        logger.info("Starting bulk symmetry search...")
        tl_symmetry.findSymmetry(bsl, rp, bulk=True, output=False)
        bsl.revertUnitCell()  # keep origin matched with main slab
        logger.info(f"Found bulk plane group: {bsl.foundplanegroup}")
        tl_symmetry.findBulkSymmetry(bsl, rp)

        # write POSCAR_bulk
        bsl = copy.deepcopy(sl.bulkslab)
        bsl.sortOriginal()
        try:
            poscar.write(bsl, filename='POSCAR_bulk', comments='bulk')
        except Exception:
            logger.error("Exception occurred while writing POSCAR_bulk")
            raise

    # write POSCAR_bulk_appended
    n = 1
    if len(bsl.sublayers) <= len(bsl.elements)*2:
        # For bulk slabs with very few layers, add a few more
        # repetitions to help the user see better the bulk part
        n += 1
        if len(bsl.sublayers) <= len(bsl.elements):
            n += 1
    try:
        poscar.write(sl.addBulkLayers(rp, n=n)[0],
                     filename='POSCAR_bulk_appended')
    except Exception:
        logger.warning("Exception occurred while writing POSCAR_bulk_appended")

    # Check for an ambiguous angle phi
    _check_and_warn_ambiguous_phi(sl, rp, angle_eps=0.1)

    # Check that layer cuts are not too close together
    _check_and_warn_layer_cuts(rp, sl)

    # check whether PHASESHIFTS are present & consistent:
    newpsGen, newpsWrite = True, True
    # True: new phaseshifts need to be generated/written
    if os.path.isfile("PHASESHIFTS") or os.path.isfile("_PHASESHIFTS"):
        try:
            (rp.phaseshifts_firstline,
             rp.phaseshifts,
             newpsGen,
             newpsWrite) = phaseshifts.readPHASESHIFTS(sl, rp,
                                                       ignoreEnRange=subdomain)
        except Exception:
            logger.warning(
                "Found a PHASESHIFTS file but could not "
                "read it. A new PHASESHIFTS file will be generated."
                "The exception during read was: ", exc_info=True)
            rp.setHaltingLevel(1)
    if rp.TL_VERSION >= 1.71 and rp.V0_REAL != 'default':
        # check V0_REAL - may have to replace firstline
        if type(rp.V0_REAL) != list:
            logger.warning("Parameter V0_REAL currently does not support "
                           "values other than Rundgren-type potentials."
                           "Input for V0_REAL will be ignored.")
            rp.setHaltingLevel(2)
            newpsGen, newpsWrite = True, True
        else:
            llist = rp.phaseshifts_firstline.split()
            try:
                flist = [float(llist[i+1]) for i in range(4)]
            except (ValueError, IndexError):
                firstline_write = True
            else:
                firstline_write = any(abs(flist[i] - rp.V0_REAL[i]) >= 1e-2
                                      for i in range(4))
            if firstline_write:
                newpsWrite = True
                rp.phaseshifts_firstline = (
                    rp.phaseshifts_firstline[:4]
                    + "{:8.2f}{:8.2f}{:8.2f}{:8.2f}".format(*rp.V0_REAL)
                    + rp.phaseshifts_firstline[36:]
                    )
    if newpsGen:
        rundgrenpath = 'EEASiSSS.x'
        serneliuspath = 'seSernelius'
        logger.info("Generating phaseshifts data... ")
        ps_gen, kwargs = runPhaseshiftGen, {}
        if rp.PHASESHIFTS_CALC_OLD:
            ps_gen = runPhaseshiftGen_old
            kwargs = {"excosource": serneliuspath}
        try:
            (rp.phaseshifts_firstline,
             rp.phaseshifts) = ps_gen(sl, rp, psgensource=rundgrenpath,
                                      **kwargs)
        except Exception:
            logger.error("Exception while calling phaseshiftgen: ")
            raise
        else:
            logger.debug("Finished generating phaseshift data")
    if newpsWrite:
        try:
            phaseshifts.writePHASESHIFTS(rp.phaseshifts_firstline,
                                         rp.phaseshifts)
        except Exception:
            logger.error("Exception during writePHASESHIFTS: ")
            raise
    rp.fileLoaded["PHASESHIFTS"] = True
    rp.updateDerivedParams()
    rp.manifest.append("PHASESHIFTS")
    try:
        phaseshifts.plot_phaseshifts(sl, rp)
    except Exception:
        logger.warning("Failed to plot phaseshifts", exc_info=rp.is_debug_mode)

    # generate beamlist
    logger.info("Generating BEAMLIST...")
    calc_and_write_beamlist(sl, rp, beamlist_name="BEAMLIST")

    try:
        rp.beamlist = tl_beams.readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except Exception:
        logger.error("Error while reading required file BEAMLIST")
        raise

    if not subdomain:
        patterninfo.writePatternInfo(sl, rp)

        # if EXPBEAMS was loaded, it hasn't been checked yet - check now
        if rp.fileLoaded["EXPBEAMS"]:
            tl_beams.checkEXPBEAMS(sl, rp)
        # write and sort IVBEAMS
        if not rp.fileLoaded["IVBEAMS"]:
            try:
                rp.ivbeams = tl_beams.writeIVBEAMS(sl, rp)
                rp.ivbeams_sorted = False
                rp.fileLoaded["IVBEAMS"] = True
                rp.manifest.append("IVBEAMS")
            except Exception:
                logger.error("Error while writing IVBEAMS file based on "
                             "EXPBEAMS data.")
                raise
    if rp.fileLoaded["IVBEAMS"] and not rp.ivbeams_sorted:
        rp.ivbeams = tl_beams.sortIVBEAMS(sl, rp)
        rp.ivbeams_sorted = True

    # Create directory compile_logs in which logs from compilation will be saved
    make_compile_logs_dir(rp)

    # At the end of initialization preserve the original input files
    preserve_original_input(rp, logger)
    return


def init_domains(rp):
    """Runs an alternative initialization for the domain search. This will
    include running the 'normal' initialization for each domain."""
    rp.try_loading_expbeams_file()
    rp.initTheoEnergies()  # may be initialized based on exp. beams
    if len(rp.DOMAINS) < 2:
        logger.error("A domain search was defined, but less than two domains "
                     "are defined. Execution will stop.")
        rp.setHaltingLevel(3)
        return
    checkFiles = ["POSCAR", "PARAMETERS", "VIBROCC", "PHASESHIFTS"]
    home = Path.cwd()
    for name, path in rp.DOMAINS.items():
        # determine the target path
        target = Path(f"Domain_{name}").resolve()
        dp = DomainParameters(target, home, name)
        if target.is_dir():
            logger.warning(f"Folder {target} already exists. "
                           "Contents may get overwritten.")
        else:
            target.mkdir()
        logger.info(f"Fetching input files for domain {name}")
        if os.path.isdir(path):
            # check the path for Tensors
            tensorIndex = leedbase.getMaxTensorIndex(path)
            if tensorIndex != 0:
                try:
                    leedbase.getTensors(tensorIndex, basedir=path,
                                        targetdir=target)
                except Exception as exc:
                    tensorIndex = 0
                    logger.warning(f"Error fetching Tensors: {exc}")
            if tensorIndex != 0:
                tensorDir = target / "Tensors" / f"Tensors_{tensorIndex:03d}"
                for file in (checkFiles + ["IVBEAMS"]):
                    if os.path.isfile(tensorDir / file):
                        shutil.copy2(tensorDir / file, target / file)
                    else:
                        logger.warning(f"Input file {file} is missing in "
                                       "Tensors directory. A new reference "
                                       "calculation is required.")
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
                                    f"Error copying required file {file} for "
                                    f"domain {name} from origin folder {path}"
                                    )
                                raise RuntimeError("Error getting domain "
                                                   "input files")
                    elif file != "PHASESHIFTS":
                        logger.error(f"Required file {file} for domain {name} "
                                     f"not found in origin folder {path}")
                        raise RuntimeError("Error getting domain input files")
        elif os.path.isfile(path):
            try:
                tensorIndex = leedbase.getMaxTensorIndex(target)
            except Exception:
                tensorIndex = 0
            tensorDir = target / "Tensors" / f"Tensors_{tensorIndex + 1:03d}"
            try:
                tensorDir.mkdir(parents=True, exist_ok=True)
            except Exception:
                raise
            try:
                with ZipFile(path, 'r') as archive:
                    archive.extractall(tensorDir)                               # TODO: maybe it would be nicer to read directly from the zip file
            except Exception:
                logger.error("Failed to unpack Tensors for domain "
                             f"{name} from file {path}")
                raise RuntimeError("Error getting domain input files")
            for file in (checkFiles + ["IVBEAMS"]):
                if (tensorDir / file).is_file():
                    shutil.copy2(tensorDir / file, target / file)
                else:
                    logger.error(f"Required file {file} for domain {name} not "
                                 f"found in Tensor directory {tensorDir}")
                    raise RuntimeError("Error getting domain input files")
            dp.tensorDir = tensorDir
        try:
            # initialize for that domain
            os.chdir(target)
            logger.info(f"Reading input files for domain {name}")
            try:
                dp.sl = poscar.read()
                dp.rp = parameters.read()                                       # NB: if we are running from stored Tensors, then these parameters will be stored versions, not current PARAMETERS from Domain directory
                dp.rp.workdir = home
                dp.rp.source_dir = rp.source_dir
                dp.rp.timestamp = rp.timestamp
                interpret_domain_params_silent = rp.LOG_LEVEL > logging.DEBUG
                parameters.interpret(dp.rp, slab=dp.sl,
                                     silent=interpret_domain_params_silent)
                dp.sl.fullUpdate(dp.rp)   # gets PARAMETERS data into slab
                dp.rp.fileLoaded["POSCAR"] = True
                dp.rp.updateDerivedParams()
                try:
                    vibrocc.readVIBROCC(dp.rp, dp.sl)
                    dp.rp.fileLoaded["VIBROCC"] = True
                except Exception:
                    logger.error("Error while reading required file VIBROCC")
                    raise
                dp.sl.fullUpdate(dp.rp)
                try:
                    dp.rp.ivbeams = tl_beams.readIVBEAMS()
                except FileNotFoundError:
                    pass
                except Exception:
                    logger.error(
                        f"Error while reading IVBEAMS for domain {name}"
                        )
                else:
                    dp.rp.ivbeams_sorted = False
                    dp.rp.fileLoaded["IVBEAMS"] = True
            except Exception:
                logger.error("Error loading POSCAR and "
                             f"PARAMETERS for domain {name}")
                raise
            logger.info(f"Running initialization for domain {name}")
            try:
                initialization(dp.sl, dp.rp, subdomain=True)
            except Exception:
                logger.error(f"Error running initialization for domain {name}")
                raise
            rp.domainParams.append(dp)
        except Exception:
            logger.error(f"Error while initializing domain {name}")
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
        if (all(abs(np.linalg.norm(bulkuc0[i]) - np.linalg.norm(bulkuc[i]))
                < eps for i in range(0, 2))
                and abs(angle(bulkuc[0], bulkuc[1])
                        - angle(bulkuc0[0], bulkuc0[1])) < eps):
            logger.info(f"Bulk unit cells of domain {rp.domainParams[0].name} "
                        f"and domain {dp.name} are mismatched, but can be "
                        f"matched by rotating domain {dp.name}.")
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
            logger.error(f"Bulk unit cells of domain {rp.domainParams[0].name}"
                         f" and domain {dp.name} are mismatched, and cannot be"
                         "matched by rotation. Domain search cannot proceed. "
                         "Execution will stop.")
            rp.setHaltingLevel(3)
            return
    logger.debug("Domain bulk unit cells are compatible.")
    uc0 = rp.domainParams[0].sl.ucell[:2, :2].T
    largestDomain = rp.domainParams[0]
    allMatched = all(np.all(abs(dp.sl.ucell[:2, :2].T - uc0) < 1e-4)
                     for dp in rp.domainParams[1:])
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
                        f"Surface unit cell of domain {dp.name} cannot be "
                        "transformed to the largest surface unit cell (domain "
                        f"{largestDomain.name}) by an integer transformation. "
                        "Execution will stop. Please supply all domain "
                        "structures as matching supercells."
                        )
                    rp.setHaltingLevel(3)
                    return
                else:
                    supercellRequired.append(dp)
                    oldslab = dp.sl
                    dp.sl = dp.sl.makeSupercell(np.round(trans))
                    dp.rp.SUPERLATTICE = largestDomain.rp.SUPERLATTICE.copy()
                    dp.sl.symbaseslab = oldslab
                    dp.rp.SYMMETRY_CELL_TRANSFORM = trans
                    parameters.modify(dp.rp, "SYMMETRY_CELL_TRANSFORM",
                                      path=dp.workdir)
        logger.info("Domain surface unit cells are mismatched, but can be "
                    "matched by integer transformations.")
    # store some information about the supercell in rp:
    rp.pseudoSlab = Slab()
    rp.pseudoSlab.ucell = largestDomain.sl.ucell.copy()
    rp.pseudoSlab.bulkslab = Slab()
    rp.pseudoSlab.bulkslab.ucell = largestDomain.sl.bulkslab.ucell.copy()
    # run beamgen for the whole system
    logger.info("Generating BEAMLIST...")
    calc_and_write_beamlist(copy.deepcopy(largestDomain.sl),
                      rp,
                      domains=True,
                      beamlist_name='BEAMLIST')
    try:
        rp.beamlist = tl_beams.readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except Exception:
        logger.error("Error while reading required file BEAMLIST")
        raise
    # if EXPBEAMS was loaded, it hasn't been checked yet - check now
    if rp.fileLoaded["EXPBEAMS"]:
        tl_beams.checkEXPBEAMS(None, rp, domains=True)
    # write and sort IVBEAMS
    if not rp.fileLoaded["IVBEAMS"]:
        try:
            rp.ivbeams = tl_beams.writeIVBEAMS(None, rp, domains=True)
            rp.ivbeams_sorted = False
            rp.fileLoaded["IVBEAMS"] = True
            rp.manifest.append("IVBEAMS")
        except Exception:
            logger.error("Error while writing IVBEAMS file based on "
                         "EXPBEAMS data.")
            raise
    if not rp.ivbeams_sorted:
        rp.ivbeams = tl_beams.sortIVBEAMS(None, rp)
        rp.ivbeams_sorted = True

    rp.updateDerivedParams()  # Also sets LMAX
    if not rp.LMAX.has_max:
        rp.LMAX.max = max(dp.rp.LMAX.max for dp in rp.domainParams)
    for dp in rp.domainParams:
        if dp.refcalcRequired:
            continue
        cmessage = f"Reference calculation required for domain {dp.name}: "
        # check energies
        if not dp.rp.THEO_ENERGIES.contains(rp.THEO_ENERGIES):
            logger.info("%sEnergy range is mismatched.", cmessage)
            dp.refcalcRequired = True
            continue
        # check LMAX - should be obsolete since TensErLEED version 1.6
        if rp.TL_VERSION <= 1.6 and rp.LMAX.max != dp.rp.LMAX.max:
            logger.info("%sLMAX is mismatched.", cmessage)
            dp.refcalcRequired = True
        # check beam incidence
        if rp.THETA != dp.rp.THETA or rp.PHI != dp.rp.PHI:
            logger.info("%sBEAM_INCIDENCE is mismatched.", cmessage)
            dp.refcalcRequired = True
        # check IVBEAMS
        if not dp.rp.fileLoaded["IVBEAMS"]:
            logger.info("%sNo IVBEAMS file loaded", cmessage)
            dp.refcalcRequired = True
            continue
        if (len(rp.ivbeams) != len(dp.rp.ivbeams)
                or not all(dp.rp.ivbeams[i].isEqual(rp.ivbeams[i])
                           for i in range(len(rp.ivbeams)))):
            logger.info("%sIVBEAMS file mismatched with supercell.", cmessage)
            dp.refcalcRequired = True
            continue

    rr = [dp for dp in rp.domainParams if dp.refcalcRequired]
    if rr:
        logger.info("The following domains require new reference "
                    f"calculations: {', '.join(d.name for d in rr)}")
        for dp in rp.domainParams:
            for var in ["THEO_ENERGIES", "THETA", "PHI", "N_CORES", "ivbeams",
                        "sourcedir"]:
                setattr(dp.rp, var, copy.deepcopy(getattr(rp, var)))
            if rp.TL_VERSION <= 1.6:  # not required since TensErLEED v1.61
                dp.rp.LMAX.max = rp.LMAX.max

    # repeat initialization for all slabs that require a supercell
    for dp in supercellRequired:
        logger.info("Re-running initialization with "
                    f"supercell slab for domain {dp.name}")
        try:
            os.chdir(dp.workdir)
            dp.sl.resetSymmetry()
            dp.rp.SYMMETRY_FIND_ORI = True
            initialization(dp.sl, dp.rp, subdomain=True)
        except Exception:
            logger.error(f"Error while re-initializing domain {dp.name}")
            raise
        finally:
            os.chdir(home)

    if 4 not in rp.RUN and 1 not in rp.RUN and rr:
        logger.error(
            "Some domains require new reference calculations before "
            "a domain search can be executed. Please either manually "
            "execute appropriate reference calculations, or set RUN = 4"
            )
        rp.setHaltingLevel(3)
        return

    while 4 in rp.RUN:
        if rr:
            rp.RUN.insert(rp.RUN.index(4), 1)
        rp.RUN.insert(rp.RUN.index(4), 2)
        rp.RUN.insert(rp.RUN.index(4), 3)
        rp.RUN.remove(4)


def preserve_original_input(rp, init_logger, path=""):
    """Create original_inputs directory and copy there the input files."""
    path = Path(path)
    orig_inputs_path = path / ORIGINAL_INPUTS_DIR_NAME
    try:
        orig_inputs_path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise RuntimeError(f"Could not create directory "
                           f"{ORIGINAL_INPUTS_DIR_NAME}. "
                           "Check disk permissions.") from exc

    # We will copy all files that have potentially been used as
    # inputs. Make sure the correct version of EXPBEAMS is stored
    files_to_preserve = ALL_INPUT_FILES.copy()
    files_to_preserve.remove('EXPBEAMS')
    if rp.expbeams_file_name:
        files_to_preserve.add(rp.expbeams_file_name)

    # copy all files to orig_inputs that were used as original input
    for file in files_to_preserve:
        # save under name EXPBEAMS.csv
        file_path = path / file
        if not file_path.is_file():
            init_logger.warning(f"Could not find file {file}. "
                                    "It will not be stored in "
                                    f"{ORIGINAL_INPUTS_DIR_NAME}.")
            rp.setHaltingLevel(1)
            return
        try:
            shutil.copy2(file_path, orig_inputs_path)
        except OSError:
            init_logger.warning(f"Could not copy file {file} to "
                                f"{ORIGINAL_INPUTS_DIR_NAME}.")
            rp.setHaltingLevel(1)


def make_compile_logs_dir(rp):
    """Create compile_logs directory where compilation logs are saved."""
    # put into rp
    rp.compile_logs_dir = Path(rp.workdir) / "compile_logs"

    # makes compile_logs directory
    try:
        rp.compile_logs_dir.mkdir(exist_ok=True)
    except OSError:
        logger.warning(f"Could not create directory {rp.compile_logs_dir}")
        rp.setHaltingLevel(1)


def _check_and_warn_ambiguous_phi(sl, rp, angle_eps=0.1):
    """Check if phi is ambiguous and warn if so."""
    angle_between_first_uc_vec_and_x = sl.angle_between_ucell_and_coord_sys
    if angle_between_first_uc_vec_and_x > angle_eps and rp.THETA > angle_eps:
        logger.info(
            f"Detected non-zero angle theta ({rp.THETA:.2f})° and"
            f"an angle of {angle_between_first_uc_vec_and_x:.2f}° "
            "between the first unit cell vector and the x direction of "
            "the coordinate system in the POSCAR file.\n"
            "Make sure the angle phi is interpreted correctly: "
            f"Phi is {rp.PHI:.2f}° from x, which is "
            f"{(rp.PHI+ angle_between_first_uc_vec_and_x):.2f}° from a.\n"
            "See the ViPErLEED documentation for the parameter BEAM_INCDIDENCE "
            "for details."
            )

def _check_and_warn_layer_cuts(rpars, slab):
    """Check if layer cuts are too close together and warn if so."""
    layer_cuts = rpars.LAYER_CUTS
    min_spacing = slab.getMinLayerSpacing()
    if min_spacing < 1.0:
        logger.warning(
            f"Layer cuts are very close together. The minimum spacing "
            f"between layers is {min_spacing:.2f} Å. This may lead to "
            "covergence issues in the reference calculation. Check the "
            "LAYERS_CUTS parameter in the PARAMETERS file."
            )
