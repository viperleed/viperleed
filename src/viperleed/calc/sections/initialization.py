"""Module initialization of viperleed.calc.sections.

ViPErLEED calculation section INITIALIZATION.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path

import numpy as np

from viperleed.calc import symmetry
from viperleed.calc.classes.rparams.domain_params import DomainParameters
from viperleed.calc.classes.slab import AlreadyMinimalError
from viperleed.calc.classes.slab import BulkSlab
from viperleed.calc.classes.slab import NoBulkRepeatError
from viperleed.calc.classes.slab import NoVacuumError
from viperleed.calc.classes.slab import Slab
from viperleed.calc.classes.slab import VacuumError
from viperleed.calc.classes.slab import WrongVacuumPositionError
from viperleed.calc.constants import DEFAULT_DOMAIN_FOLDER_PREFIX
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.files import beams as iobeams
from viperleed.calc.files import parameters
from viperleed.calc.files import experiment_symmetry
from viperleed.calc.files import phaseshifts
from viperleed.calc.files import poscar
from viperleed.calc.files import vibrocc
from viperleed.calc.files.beamgen import calc_and_write_beamlist
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.math_utils import angle
from viperleed.calc.lib.matrix import NonIntegerMatrixError
from viperleed.calc.lib.matrix import rotation_matrix
from viperleed.calc.lib.version import Version
from viperleed.calc.lib.woods_notation import writeWoodsNotation
from viperleed.calc.psgen import runPhaseshiftGen, runPhaseshiftGen_old
from viperleed.calc.sections.cleanup import preserve_original_inputs

logger = logging.getLogger(__name__)


def initialization(sl, rp, subdomain=False):
    """Runs the initialization."""
    if not subdomain:
        rp.try_loading_expbeams_file()
    rp.initTheoEnergies()  # may be initialized based on exp. beams

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
            parameters.modify(rp, "SYMMETRY_CELL_TRANSFORM")                    # TODO: there should probably a comment in the PARAMETERS file?
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
        poscar.write(tmpslab, filename='POSCAR', comments='all')
    except Exception:
        logger.error("Exception occurred while writing new POSCAR")
        raise
    rp.files_to_out.add('POSCAR')
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
            bulk_repeat = parameters.modify(
                rp, 'BULK_REPEAT',
                comment='Automatically detected repeat vector',
                )
            parameters.modify(rp, 'LAYER_CUTS')
            parameters.modify(rp, 'N_BULK_LAYERS')
            logger.info(
                f'Detected bulk repeat vector: {bulk_repeat.fmt_value}'
                )
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
            bulk_repeat = parameters.modify(
                rp, 'BULK_REPEAT',
                comment='Automatically detected repeat vector'
                )
            logger.info('Detected bulk repeat vector: %s',
                        bulk_repeat.fmt_value)
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
            ws = writeWoodsNotation(rp.SUPERLATTICE)                   # TODO: replace writeWoodsNotation with gui functions
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
    if rp.TL_VERSION >= Version('1.7.1') and rp.V0_REAL != 'default':
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
        tensorleed = rp.paths.tensorleed
        # Check for old executable. Used to be called EEASiSSS.x
        if (not (tensorleed / 'eeasisss').is_file()
                and (tensorleed / 'EEASiSSS.x').is_file()):
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
    rp.manifest.add("PHASESHIFTS")
    try:
        phaseshifts.plot_phaseshifts(sl, rp)
    except Exception:
        logger.warning("Failed to plot phaseshifts", exc_info=rp.is_debug_mode)

    # generate beamlist
    logger.info("Generating BEAMLIST...")                                       # TODO: this bit is largely repeated in init_domains
    calc_and_write_beamlist(sl, rp, beamlist_name="BEAMLIST")

    try:
        rp.beamlist = iobeams.readBEAMLIST()
        rp.fileLoaded["BEAMLIST"] = True
    except Exception:
        logger.error("Error while reading required file BEAMLIST")
        raise

    if not subdomain:
        try:
            experiment_symmetry.write(sl, rp)
        except (OSError, ValueError):
            # OSError: failed to write file. It's not that critical,
            # so we can probably go ahead. We logged the problem
            # already. If it is a more fundamental issue it will
            # pop up when we try to do more file-system operations.
            # ValueError: SUPERLATTICE not integer. Probably
            # we complain already somewhere else. Surely in
            # iobeams.writeIVBEAMS, likely already earlier
            # when we work on the slab.
            pass

        # if EXPBEAMS was loaded, it hasn't been checked yet - check now
        if rp.fileLoaded["EXPBEAMS"]:
            iobeams.checkEXPBEAMS(sl, rp)
        # write and sort IVBEAMS
        if not rp.fileLoaded["IVBEAMS"]:
            try:
                rp.ivbeams = iobeams.writeIVBEAMS(sl, rp)
                rp.ivbeams_sorted = False
                rp.fileLoaded["IVBEAMS"] = True
                rp.manifest.add("IVBEAMS")
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
    # Make sure we always have unique work folders for all domains
    nr_unique_paths = len(set(p for p, _ in rp.DOMAINS.values()))
    must_use_auto_name = nr_unique_paths < len(rp.DOMAINS)
    for name, (path, user_given) in rp.DOMAINS.items():
        workdir = _make_domain_workdir(name,
                                       path,
                                       rp.paths.home,
                                       must_use_auto_name)
        domain = DomainParameters(workdir, name)
        domain.collect_input_files(path)
        mod_path = Path(domain.workdir.name)
        mod_value = parameters.modify(rp, 'DOMAIN', f'./{mod_path}',
                                      original=user_given)
        rp.DOMAINS[name] = (mod_path, mod_value.to_assignment())
        with execute_in_dir(domain.workdir):
            try:  # Initialize for that domain
                _run_initialization_for_domain(domain, rp)
            except Exception:
                logger.error(f'Error while initializing {domain}')
                raise
    if len(rp.domainParams) < len(rp.DOMAINS):
        raise RuntimeError("Failed to read domain parameters")
    # check whether bulk unit cells match
    logger.info("Starting domain consistency check...")
    bulkuc0 = rp.domainParams[0].slab.bulkslab.ab_cell.T
    eps = 1e-4
    for dp in rp.domainParams[1:]:
        bulkuc = dp.slab.bulkslab.ab_cell.T
        if np.all(abs(bulkuc-bulkuc0) < eps):
            continue
        # if the unit cells don't match right away, try if rotation matches
        if (all(abs(np.linalg.norm(bulkuc0[i]) - np.linalg.norm(bulkuc[i]))
                < eps for i in range(0, 2))
                and abs(angle(*bulkuc) - angle(*bulkuc0)) < eps):
            logger.info(f"Bulk unit cells of {rp.domainParams[0]} "
                        f"and {dp} are mismatched, but can be "
                        f"matched by rotating {dp}.")
            ang = angle(bulkuc[0], bulkuc0[0])
            dp.slab.apply_matrix_transformation(rotation_matrix(ang, dim=3))    # TODO: this changes the coordinate frame. We need to modify BEAM_INCIDENCE! Issue #69, PR #73
        else:
            logger.error(f'Bulk unit cells of {rp.domainParams[0]} '
                         f'and {dp} are mismatched, and cannot be '
                         'matched by rotation. Domain search cannot '
                         'proceed. Execution will stop.')
            rp.setHaltingLevel(3)
            return
    logger.debug("Domain bulk unit cells are compatible.")
    uc0 = rp.domainParams[0].slab.ab_cell.T
    largestDomain = rp.domainParams[0]
    allMatched = all(np.all(abs(dp.slab.ab_cell.T - uc0) < 1e-4)
                     for dp in rp.domainParams[1:])
    supercellRequired = []
    if allMatched:
        logger.debug("Domain surface unit cells are matched.")
    else:
        maxArea = abs(np.linalg.det(uc0))
        for dp in rp.domainParams[1:]:
            uc = dp.slab.ab_cell.T
            if abs(np.linalg.det(uc)) > maxArea:
                maxArea = abs(np.linalg.det(uc))
                largestDomain = dp
        uc0 = largestDomain.slab.ab_cell.T
        for dp in [p for p in rp.domainParams if p != largestDomain]:
            uc = dp.slab.ab_cell.T
            if not np.all(abs(uc-uc0) < 1e-4):
                dp.refcalc_required = True
                trans = np.dot(uc0, np.linalg.inv(uc))
                if np.any(abs(trans - np.round(trans)) > 1e-4):
                    logger.error(
                        f'Surface unit cell of {dp} cannot be transformed '
                        f'to the largest surface unit cell ({largestDomain}) '
                        'by an integer transformation. Execution will stop. '
                        'Please supply all domain structures as matching '
                        'supercells.'
                        )
                    rp.setHaltingLevel(3)
                    return
                else:
                    supercellRequired.append(dp)
                    oldslab = dp.slab
                    dp.slab = dp.slab.make_supercell(np.round(trans))
                    dp.rpars.SUPERLATTICE = (
                        largestDomain.rpars.SUPERLATTICE.copy()
                        )
                    dp.slab.symbaseslab = oldslab
                    dp.rpars.SYMMETRY_CELL_TRANSFORM = trans
                    with execute_in_dir(dp.workdir):
                        parameters.modify(dp.rpars, 'SYMMETRY_CELL_TRANSFORM')
        logger.info("Domain surface unit cells are mismatched, but can be "
                    "matched by integer transformations.")
    # store some information about the supercell in rp:
    rp.pseudoSlab = Slab()
    rp.pseudoSlab.ucell = largestDomain.slab.ucell.copy()
    rp.pseudoSlab.bulkslab = BulkSlab()
    rp.pseudoSlab.bulkslab.ucell = largestDomain.slab.bulkslab.ucell.copy()
    # run beamgen for the whole system
    logger.info("Generating BEAMLIST...")                                       # TODO: this bit is largely repeated in the end of initialization
    calc_and_write_beamlist(copy.deepcopy(largestDomain.slab),
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
            rp.manifest.add("IVBEAMS")
        except Exception:
            logger.error("Error while writing IVBEAMS file based on "
                         "EXPBEAMS data.")
            raise
    if not rp.ivbeams_sorted:
        rp.ivbeams = iobeams.sortIVBEAMS(None, rp)
        rp.ivbeams_sorted = True

    rp.updateDerivedParams()  # Also sets LMAX
    if not rp.LMAX.has_max:
        rp.LMAX.max = max(dp.rpars.LMAX.max for dp in rp.domainParams)
    for dp in rp.domainParams:
        if dp.refcalc_required:
            continue
        cmessage = f'Reference calculation required for {dp}: '
        # check energies
        if not dp.rpars.THEO_ENERGIES.contains(rp.THEO_ENERGIES):
            logger.info("%sEnergy range is mismatched.", cmessage)
            dp.refcalc_required = True
            continue
        # check LMAX - should be obsolete since TensErLEED version 1.6
        if (rp.TL_VERSION <= Version('1.6.0')
                and rp.LMAX.max != dp.rpars.LMAX.max):
            logger.info("%sLMAX is mismatched.", cmessage)
            dp.refcalc_required = True
        # check beam incidence
        if rp.THETA != dp.rpars.THETA or rp.PHI != dp.rpars.PHI:
            logger.info("%sBEAM_INCIDENCE is mismatched.", cmessage)
            dp.refcalc_required = True
        # check IVBEAMS
        if not dp.rpars.fileLoaded["IVBEAMS"]:
            logger.info("%sNo IVBEAMS file loaded", cmessage)
            dp.refcalc_required = True
            continue
        if (len(rp.ivbeams) != len(dp.rpars.ivbeams)
                or not all(dp.rpars.ivbeams[i].isEqual(rp.ivbeams[i])
                           for i in range(len(rp.ivbeams)))):
            logger.info("%sIVBEAMS file mismatched with supercell.", cmessage)
            dp.refcalc_required = True
            continue

    rr = [dp for dp in rp.domainParams if dp.refcalc_required]
    if rr:
        logger.info("The following domains require new reference "
                    f"calculations: {', '.join(d.name for d in rr)}")
        inherited = (
            'THEO_ENERGIES',
            'THETA',
            'PHI',
            'N_CORES',
            'ivbeams',
            )
        for dp in rp.domainParams:
            dp.rpars.inherit_from(rp, *inherited, override=True)
            if rp.TL_VERSION <= Version('1.6.0'):  # not required since TensErLEED v1.61
                dp.rpars.LMAX.max = rp.LMAX.max

    # repeat initialization for all slabs that require a supercell
    for dp in supercellRequired:
        logger.info(f'Re-running initialization with supercell slab for {dp}')
        dp.slab.clear_symmetry_and_ucell_history()
        dp.rpars.SYMMETRY_FIND_ORI = True
        with execute_in_dir(dp.workdir):
            try:
                initialization(dp.slab, dp.rpars, subdomain=True)
            except Exception:
                logger.error(f'Error while re-initializing {dp}')
                raise

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


def make_compile_logs_dir(rpars):
    """Create compile_logs directory where compilation logs are saved."""
    directory = rpars.paths.compile_logs
    try:
        directory.mkdir(exist_ok=True)
    except OSError:
        logger.warning(f'Could not create directory {directory}')
        rpars.setHaltingLevel(1)


def warn_if_slab_has_atoms_in_multiple_c_cells(slab, rpars, domain=None):
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
    domain : DomainParameters, optional
        The structural domain to which `slab` belongs. Used only for
        logging purposes. Default is None.

    Returns
    -------
    None.
    """
    _msg = 'POSCAR file ' + f'of {domain} ' if domain else ''
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
            f"Detected non-zero angle theta ({rp.THETA:.2f})° and "
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
    min_spacing = slab.smallest_interlayer_gap
    if min_spacing < 1.0:
        logger.warning(
            f"Layer cuts are very close together. The minimum spacing "
            f"between layers is {min_spacing:.2f} Å. This may lead to "
            "convergence issues in the reference calculation. Check the "
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
        _msg += (
            '. This may cause problems with layer assignment! You can '
            f'find a POSCAR_vacuum_corrected file in {DEFAULT_SUPP} with '
            )
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


def _make_domain_workdir(name, src, calc_started_at, must_use_auto_name):
    """Create a work directory (as a subfolder of CWD) for a domain.

    The work directory is named 'Domain_<name>', unless
    `must_use_auto_name` is False and the inputs for this
    domain come from a direct subfolder of the path at
    which viperleed.calc was originally started.

    Parameters
    ----------
    name : str
        The user-given (or auto-generated) name of the domain for which
        the work folder should be created. Typically the left-side flag
        in a DOMAIN assignment.
    src : Path
        The path from where the domain input files should be collected.
    calc_started_at : Path or None
        The path in which viperleed.calc was started.
    must_use_auto_name : bool
        Whether 'Domain_<name>' should be used irrespective of the
        value of `src`. This is used to ensure that there always
        is one unique path for each domain.

    Returns
    -------
    workdir : Path
        Absolute path to the work directory that was created.
    """
    should_use_src = (
       not must_use_auto_name
       and calc_started_at is not None
       and src.is_dir()
       and src.parent == calc_started_at
       )
    workdir = Path(src.name if should_use_src
                   else f'{DEFAULT_DOMAIN_FOLDER_PREFIX}{name}')
    try:
        workdir.mkdir()
    except FileExistsError:
        logger.warning(f'Folder {workdir} already exists. '
                       'Contents may get overwritten.')
    return workdir.resolve()


def _read_inputs_for_domain(domain, main_rpars):
    """Read input files for a single domain.

    Parameters
    ----------
    domain : DomainParameters
        The parameters for this domain. At the end of this call,
        domain.rpars and domain.slab are updated to contain the
        Rparams and Slab objects read (and updated) from file.
        domain.rpars is also up to date with respect to other
        input files loaded.
    main_rpars : Rparams
        The run parameters of the **main** viperleed.calc run.

    Raises
    ------
    Exception
        If reading VIBROCC or an existing IVBEAMS fails.
    """
    # NB: if we are running from stored Tensors, then
    # this PARAMETERS will be a copy of the one stored
    # in the Tensors, not the one the user may have given
    # in the current Domain directory (it is overwritten
    # when fetching files in init_domains).
    domain.rpars = rpars = parameters.read()

    # Inherit some values from the main PARAMETERS
    inherited = (
        'paths',
        'timestamp',
        'ZIP_COMPRESSION_LEVEL',
        )
    rpars.inherit_from(main_rpars, *inherited)

    # Store input files for each domain, BEFORE any edit
    preserve_original_inputs(rpars)

    domain.slab = slab = poscar.read()
    warn_if_slab_has_atoms_in_multiple_c_cells(slab, rpars, domain)

    silent = rpars.LOG_LEVEL > logging.DEBUG
    parameters.interpret(rpars, slab=slab, silent=silent)
    slab.full_update(rpars)
    rpars.fileLoaded['POSCAR'] = True
    rpars.updateDerivedParams()
    try:
        vibrocc.readVIBROCC(rpars, slab)
    except Exception:
        logger.error('Error while reading required file VIBROCC')
        raise
    rpars.fileLoaded['VIBROCC'] = True
    slab.full_update(rpars)
    try:
        rpars.ivbeams = iobeams.readIVBEAMS()
    except FileNotFoundError:
        pass
    except Exception:
        logger.error(f'Error while reading IVBEAMS for {domain}')
    else:
        rpars.ivbeams_sorted = False
        rpars.fileLoaded['IVBEAMS'] = True


def _run_initialization_for_domain(domain, main_rpars):
    """Execute initialization for a specific domain.

    Parameters
    ----------
    domain : DomainParameters
        The parameters for the domain to be initialized.
    main_rpars : Rparams
        The run parameters of the **main** viperleed.calc run.

    Raises
    ------
    Exception
        If reading input files or executing the initialization fail.
    """
    logger.info(f'Reading input files for {domain}')
    try:
        _read_inputs_for_domain(domain, main_rpars)
    except Exception:
        logger.error(f'Error loading input files for {domain}')
        raise

    logger.info(f'Running initialization for {domain}')
    try:
        initialization(domain.slab, domain.rpars, subdomain=True)
    except Exception:
        logger.error(f'Error running initialization for {domain}')
        raise

    main_rpars.domainParams.append(domain)
