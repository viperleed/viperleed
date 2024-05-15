"""Module initialization of viperleed.calc.sections.

ViPErLEED calculation section INITIALIZATION.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path
import shutil
from zipfile import ZipFile

import numpy as np

from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc import symmetry
from viperleed.calc.classes.rparams import DomainParameters
from viperleed.calc.classes.slab import AlreadyMinimalError
from viperleed.calc.classes.slab import BulkSlab
from viperleed.calc.classes.slab import NoBulkRepeatError
from viperleed.calc.classes.slab import NoVacuumError
from viperleed.calc.classes.slab import Slab
from viperleed.calc.classes.slab import VacuumError
from viperleed.calc.classes.slab import WrongVacuumPositionError
from viperleed.calc.files import beams as iobeams
from viperleed.calc.files import iotensors
from viperleed.calc.files import parameters
from viperleed.calc.files import patterninfo
from viperleed.calc.files import phaseshifts
from viperleed.calc.files import poscar
from viperleed.calc.files import vibrocc
from viperleed.calc.files.beamgen import calc_and_write_beamlist
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.base import NonIntegerMatrixError
from viperleed.calc.lib.base import angle, rotation_matrix
from viperleed.calc.lib.woods_notation import writeWoodsNotation
from viperleed.calc.psgen import runPhaseshiftGen, runPhaseshiftGen_old
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.calc_section import EXPBEAMS_NAMES

logger = logging.getLogger(__name__)


OPTIONAL_INPUT_FILES = ('BEAMLIST')

def initialization(sl, rp, subdomain=False):
    """Runs the initialization."""
    if not subdomain:
        rp.try_loading_expbeams_file()
    rp.initTheoEnergies()  # may be initialized based on exp. beams

    # preserve unmodified input files to original_inputs directory
    if not subdomain:
        rp.inputs_dir = rp.workdir / ORIGINAL_INPUTS_DIR_NAME
        _preserve_original_input(rp, origin_dir=Path.cwd().resolve())

    if (rp.DOMAINS or rp.domainParams) and not subdomain:
        init_domains(rp)
        return

    _check_slab_duplicates_and_vacuum(sl, rp)

    # if necessary, run findSymmetry:
    if sl.planegroup == "unknown":
        symmetry.findSymmetry(sl, rp)
        symmetry.enforceSymmetry(sl, rp)

    # check whether the slab unit cell is minimal:
    try:
        mincell = sl.get_minimal_ab_cell(rp.SYMMETRY_EPS, rp.SYMMETRY_EPS.z)
    except AlreadyMinimalError:
        transform = np.identity(2)
        reducible = False
    else:
        transform = np.dot(sl.ab_cell.T, np.linalg.inv(mincell)).round()
        reducible = True
    ws = writeWoodsNotation(transform)
    if ws:
        ws = f"= {ws}"
    else:
        ws = "M = {} {}, {} {}".format(*transform.astype(int).ravel())
    if reducible and np.isclose(rp.SYMMETRY_CELL_TRANSFORM,
                                np.identity(2)).all():
        ssl = sl.make_subcell(rp, transform=transform)
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
            try:
                sl.symbaseslab = sl.make_subcell(rp,
                                                 rp.SYMMETRY_CELL_TRANSFORM)
            except ValueError as exc:
                logger.error(f'Invalid SYMMETRY_CELL_TRANSFORM. {exc}')
                raise exc
    if sl.symbaseslab is not None:
        logger.info("A symmetry cell transformation was found. Re-running "
                    "slab symmetry search using base unit cell...")
        symmetry.getSymBaseSymmetry(sl, rp)
        try:
            poscar.write(sl.symbaseslab, filename='POSCAR_mincell',
                         comments='all')
        except OSError:
            logger.warning("Exception occurred while writing POSCAR_mincell")

    # generate new POSCAR
    tmpslab = copy.deepcopy(sl)
    tmpslab.sort_original()
    try:
        poscar.write(tmpslab,
                     filename=rp.workdir/'POSCAR_OUT', comments='all')
    except Exception:
        logger.error("Exception occurred while writing new POSCAR")
        raise
    # generate POSCAR_oricell
    tmpslab.revert_unit_cell()
    try:
        poscar.write(tmpslab, filename='POSCAR_oricell', comments='nodir')
    except OSError:
        logger.error("Exception occurred while writing POSCAR_oricell, "
                     "execution will continue...")

    if rp.BULK_LIKE_BELOW > 0:
        if rp.BULK_REPEAT is not None:
            logger.warning('Both BULK_LIKE_BELOW and BULK_REPEAT are defined.'
                           'BULK_LIKE_BELOW will be ignored in favour of the '
                           'explicitly defined bulk repeat vector.')
        else:
            # If successful, the next call updates the relevant rpars
            # attributes: SUPERLATTICE (not written to file),
            # BULK_REPEAT, LAYER_CUTS, and N_BULK_LAYERS (written
            # to file below). The slab will have .layers and the
            # freshly detected .bulkslab.
            sl.detect_bulk(rp)
            vec_str = parameters.modify(
                rp, 'BULK_REPEAT',
                comment='Automatically detected repeat vector'
                )
            parameters.modify(rp, 'LAYER_CUTS')
            parameters.modify(rp, 'N_BULK_LAYERS')
            logger.info(f'Detected bulk repeat vector: {vec_str}')
        parameters.comment_out(rp, 'BULK_LIKE_BELOW')

    # create bulk slab:
    if sl.bulkslab is None:
        sl.make_bulk_slab(rp)
    bsl = sl.bulkslab

    # try identifying bulk repeat:
    if rp.BULK_REPEAT is None:
        try:
            rp.BULK_REPEAT = sl.identify_bulk_repeat(rp.SYMMETRY_EPS,
                                                     rp.SYMMETRY_EPS.z)
        except NoBulkRepeatError:
            pass
        else:
            vec_str = parameters.modify(
                rp, 'BULK_REPEAT',
                comment='Automatically detected repeat vector'
                )
            logger.info(f'Detected bulk repeat vector: {vec_str}')
            # update bulk slab vector
            sl.bulkslab.update_cartesian_from_fractional()
            sl.bulkslab.c_vector[:] = rp.BULK_REPEAT
            sl.bulkslab.collapse_cartesian_coordinates()
    if rp.BULK_REPEAT is None:
        # Failed to detect repeat vector, use fixed distance instead.
        # Assume that interlayer vector from bottom non-bulk to top
        # bulk layer is the same as between bulk units. Save
        # BULK_REPEAT value for later runs, in case atom above moves
        rp.BULK_REPEAT = sl.get_bulk_repeat(rp, only_z_distance=True)
        parameters.modify(
            rp, 'BULK_REPEAT',
            comment='Automatically detected spacing. Check POSCAR_bulk.'
            )
        logger.warning(
            'The BULK_REPEAT parameter was undefined, which may lead to '
            'unintended changes in the bulk unit cell during optimization if '
            'the lowest non-bulk atom moves.\n'
            '# Automatic detection of a bulk repeat vector failed, possibly '
            'because not enough bulk-like layers were found.\n'
            '# The BULK_REPEAT vector is assumed to be parallel to the POSCAR '
            'c vector. Check POSCAR_bulk and the BULK_REPEAT thickness '
            'written to the PARAMETERS file.'
            )
        rp.checklist.append('Check bulk repeat vector in PARAMETERS')
        rp.setHaltingLevel(2)

    if bsl.planegroup == "unknown":
        # find minimum in-plane unit cell for bulk:
        logger.info('Checking bulk unit cell...')
        try:
            sl.ensure_minimal_bulk_ab_cell(rp, warn_convention=True)
        except NonIntegerMatrixError as exc:
            logger.error(
                f'{exc}\n# User-defined SUPERLATTICE will be used instead.'
                )
            rp.setHaltingLevel(2)
        except parameters.errors.InconsistentParameterError as exc:
            # SUPERLATTICE does not match auto-detected one
            logger.warning(exc)
            rp.setHaltingLevel(2)

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
        symmetry.findSymmetry(bsl, rp, bulk=True, output=False)
        bsl.revert_unit_cell()  # keep origin matched with main slab
        logger.info(f"Found bulk plane group: {bsl.foundplanegroup}")
        symmetry.findBulkSymmetry(bsl, rp)

        # write POSCAR_bulk
        bsl = copy.deepcopy(sl.bulkslab)
        bsl.sort_original()
        try:
            poscar.write(bsl, filename='POSCAR_bulk', comments='bulk')
        except Exception:
            logger.error("Exception occurred while writing POSCAR_bulk")
            raise

    # write POSCAR_bulk_appended
    n_cells = 1
    if bsl.n_sublayers <= len(bsl.elements)*2:
        # For bulk slabs with very few layers, add a few more
        # repetitions to help the user see better the bulk part
        n_cells += 1
        if bsl.n_sublayers <= len(bsl.elements):
            n_cells += 1
    try:
        poscar.write(sl.with_extra_bulk_units(rp, n_cells)[0],
                     filename='POSCAR_bulk_appended')
    except OSError:
        logger.warning('Exception occurred while writing POSCAR_bulk_appended')

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
        # Check for old executable. Used to be called EEASiSSS.x
        if (not Path(rp.source_dir / "eeasisss").is_file()
                and Path(rp.source_dir / "EEASiSSS.x").is_file()):
            rundgrenpath = 'EEASiSSS.x'
        else:
            # let psgen catch the error if neither executable is found
            rundgrenpath = 'eeasisss'
        serneliuspath = 'seSernelius'
        atom_density_path = 'atom_density_files'
        logger.info("Generating phaseshifts data... ")
        ps_gen, kwargs = runPhaseshiftGen, {}
        if rp.PHASESHIFTS_CALC_OLD:
            ps_gen = runPhaseshiftGen_old
            kwargs = {"excosource": serneliuspath,
                      "atdenssource": atom_density_path}
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
        rp.beamlist = iobeams.readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except Exception:
        logger.error("Error while reading required file BEAMLIST")
        raise

    if not subdomain:
        patterninfo.writePatternInfo(sl, rp)

        # if EXPBEAMS was loaded, it hasn't been checked yet - check now
        if rp.fileLoaded["EXPBEAMS"]:
            iobeams.checkEXPBEAMS(sl, rp)
        # write and sort IVBEAMS
        if not rp.fileLoaded["IVBEAMS"]:
            try:
                rp.ivbeams = iobeams.writeIVBEAMS(sl, rp)
                rp.ivbeams_sorted = False
                rp.fileLoaded["IVBEAMS"] = True
                rp.manifest.append("IVBEAMS")
            except Exception:
                logger.error("Error while writing IVBEAMS file based on "
                             "EXPBEAMS data.")
                raise
    if rp.fileLoaded["IVBEAMS"] and not rp.ivbeams_sorted:
        rp.ivbeams = iobeams.sortIVBEAMS(sl, rp)
        rp.ivbeams_sorted = True

    # Create directory compile_logs in which logs from compilation will be saved
    make_compile_logs_dir(rp)

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
                    iotensors.getTensors(tensorIndex, base_dir=path,
                                         target_dir=target)
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
                dp.rp.workdir = home
                dp.rp.source_dir = rp.source_dir
                dp.rp.timestamp = rp.timestamp

                # run _preserve_original_input separately for each domain
                dp.rp.inputs_dir = dp.rp.workdir / ORIGINAL_INPUTS_DIR_NAME / name
                _preserve_original_input(dp.rp, origin_dir=Path.cwd().resolve())

                dp.sl = poscar.read(dp.rp.inputs_dir / "POSCAR")
                dp.rp = parameters.read()                                       # NB: if we are running from stored Tensors, then these parameters will be stored versions, not current PARAMETERS from Domain directory
                warn_if_slab_has_atoms_in_multiple_c_cells(dp.sl, dp.rp, name)
                dp.rp.inputs_dir = Path.cwd().resolve()

                parameters.interpret(dp.rp, slab=dp.sl,
                                     silent=rp.LOG_LEVEL > logging.DEBUG)
                dp.sl.full_update(dp.rp)
                dp.rp.fileLoaded["POSCAR"] = True
                dp.rp.updateDerivedParams()
                try:
                    vibrocc.readVIBROCC(dp.rp, dp.sl)
                    dp.rp.fileLoaded["VIBROCC"] = True
                except Exception:
                    logger.error("Error while reading required file VIBROCC")
                    raise
                dp.sl.full_update(dp.rp)
                try:
                    dp.rp.ivbeams = iobeams.readIVBEAMS()
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
    bulkuc0 = rp.domainParams[0].sl.bulkslab.ab_cell.T
    eps = 1e-4
    for dp in rp.domainParams[1:]:
        bulkuc = dp.sl.bulkslab.ab_cell.T
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
            ang = angle(bulkuc[0], bulkuc0[0])
            dp.sl.apply_matrix_transformation(rotation_matrix(ang, dim=3))      # TODO: this changes the coordinate frame. We need to modify BEAM_INCIDENCE! Issue #69, PR #73
        else:
            logger.error(f"Bulk unit cells of domain {rp.domainParams[0].name}"
                         f" and domain {dp.name} are mismatched, and cannot be"
                         "matched by rotation. Domain search cannot proceed. "
                         "Execution will stop.")
            rp.setHaltingLevel(3)
            return
    logger.debug("Domain bulk unit cells are compatible.")
    uc0 = rp.domainParams[0].sl.ab_cell.T
    largestDomain = rp.domainParams[0]
    allMatched = all(np.all(abs(dp.sl.ab_cell.T - uc0) < 1e-4)
                     for dp in rp.domainParams[1:])
    supercellRequired = []
    if allMatched:
        logger.debug("Domain surface unit cells are matched.")
    else:
        maxArea = abs(np.linalg.det(uc0))
        for dp in rp.domainParams[1:]:
            uc = dp.sl.ab_cell.T
            if abs(np.linalg.det(uc)) > maxArea:
                maxArea = abs(np.linalg.det(uc))
                largestDomain = dp
        uc0 = largestDomain.sl.ab_cell.T
        for dp in [p for p in rp.domainParams if p != largestDomain]:
            uc = dp.sl.ab_cell.T
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
                    dp.sl = dp.sl.make_supercell(np.round(trans))
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
    rp.pseudoSlab.bulkslab = BulkSlab()
    rp.pseudoSlab.bulkslab.ucell = largestDomain.sl.bulkslab.ucell.copy()
    # run beamgen for the whole system
    logger.info("Generating BEAMLIST...")
    calc_and_write_beamlist(copy.deepcopy(largestDomain.sl),
                      rp,
                      domains=True,
                      beamlist_name='BEAMLIST')
    try:
        rp.beamlist = iobeams.readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except Exception:
        logger.error("Error while reading required file BEAMLIST")
        raise
    # if EXPBEAMS was loaded, it hasn't been checked yet - check now
    if rp.fileLoaded["EXPBEAMS"]:
        iobeams.checkEXPBEAMS(None, rp, domains=True)
    # write and sort IVBEAMS
    if not rp.fileLoaded["IVBEAMS"]:
        try:
            rp.ivbeams = iobeams.writeIVBEAMS(None, rp, domains=True)
            rp.ivbeams_sorted = False
            rp.fileLoaded["IVBEAMS"] = True
            rp.manifest.append("IVBEAMS")
        except Exception:
            logger.error("Error while writing IVBEAMS file based on "
                         "EXPBEAMS data.")
            raise
    if not rp.ivbeams_sorted:
        rp.ivbeams = iobeams.sortIVBEAMS(None, rp)
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
                        "source_dir"]:
                setattr(dp.rp, var, copy.deepcopy(getattr(rp, var)))
            if rp.TL_VERSION <= 1.6:  # not required since TensErLEED v1.61
                dp.rp.LMAX.max = rp.LMAX.max

    # repeat initialization for all slabs that require a supercell
    for dp in supercellRequired:
        logger.info("Re-running initialization with "
                    f"supercell slab for domain {dp.name}")
        try:
            os.chdir(dp.workdir)
            dp.sl.clear_symmetry_and_ucell_history()
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


def _preserve_original_input(rp, origin_dir):
    """Create original_inputs directory and copies the input files there."""

    try:
        rp.inputs_dir.mkdir(parents=True, exist_ok=True)
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
        file_path = origin_dir / file
        if not file_path.is_file():
            if file in OPTIONAL_INPUT_FILES:
                continue
            logger.warning(f"Could not find file {file}. "
                            "It will not be stored in "
                            f"{ORIGINAL_INPUTS_DIR_NAME}.")
            rp.setHaltingLevel(1)
            return
        try:
            shutil.copy2(file_path, rp.inputs_dir)
        except OSError:
            logger.warning(f"Could not copy file {file} to "
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


def warn_if_slab_has_atoms_in_multiple_c_cells(slab, rpars, domain_name=''):
    """Log a WARNING if slab's atoms do not all belong to the same cell.

    It only makes sense to use this function right after `slab` has
    been loaded (e.g., via poscar.read or .from_ase) and before any
    of the c-collapsing functions are called. These include:
        .collapse_cartesian_coordinates
        .collapse_fractional_coordinates
        .full_update

    Parameters
    ----------
    slab : SurfaceSlab
        The slab to be checked.
    rpars : Rparams
        The current PARAMETERS object. Used for logging purposes only.
    domain_name : str, optional
        The name of the structural domain to which slab belongs.
        Used only for logging purposes. Default is an empty string.

    Returns
    -------
    None.
    """
    _msg = 'POSCAR file ' + f'of domain {domain_name} ' if domain_name else ''
    _msg += ('has some atoms outside the base unit cell along the third unit '
             'vector (i.e., they have fractional c coordinates smaller/larger '
             'than 0/1). These atoms will be back-folded into the base unit '
             'cell assuming periodicity along c. This will impact layer '
             'creation and may modify the thickness of the vacuum gap '
             'detected. Check POSCAR to ensure correct layer assignment.')
    if slab.has_atoms_in_multiple_c_cells():
        logger.warning(_msg)
        rpars.setHaltingLevel(1)


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
    min_spacing = slab.smallest_interlayer_gap
    if min_spacing < 1.0:
        logger.warning(
            f"Layer cuts are very close together. The minimum spacing "
            f"between layers is {min_spacing:.2f} Å. This may lead to "
            "covergence issues in the reference calculation. Check the "
            "LAYERS_CUTS parameter in the PARAMETERS file."
            )


def _check_slab_duplicates_and_vacuum(slab, rpars):
    """Complain if slab has duplicate atoms or an inadequate vacuum gap."""
    # Make sure that there are no duplicate atoms
    slab.check_atom_collisions(rpars.SYMMETRY_EPS)

    # Make sure that there's enough vacuum. Stop in complex situations
    try:
        slab.check_vacuum_gap()
    except VacuumError as exc:
        exc.fixed_slab.sort_original()  # Rather than by z
        poscar.write(exc.fixed_slab, 'POSCAR_vacuum_corrected')
        exc_type = type(exc)
        _msg = exc.message
        _msg += '. This may cause problems with layer assignment! '
        _msg += 'You can find a POSCAR_vacuum_corrected file in SUPP with '
        _msg += ('the correct position of vacuum. '
                 if exc_type is WrongVacuumPositionError
                 else 'a large enough vacuum gap. ')
        _msg += (
            'If you intend to use it for a new run, be careful to edit any '
            'PARAMETERS expressed as fractions of the c unit vector (e.g., '
            'LAYERS_CUTS, BULK_LIKE_BELOW, BULK_REPEAT, ...)'
            )
        # Notice the straight type check rather than isinstance, as
        # NoVacuumError is a subclass of NotEnoughVacuumError, which
        # is, instead, the acceptable case in which we warn below
        if exc_type in (NoVacuumError, WrongVacuumPositionError):
            raise exc_type(_msg, exc.fixed_slab) from None
        rpars.setHaltingLevel(1)
        logger.warning(_msg)
    if not slab.layers:
        # May have been cleared by shifting slab away from c==0
        slab.create_layers(rpars)
