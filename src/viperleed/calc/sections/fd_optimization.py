"""Section Full-dynamic Optimization."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-10-22'
__license__ = 'GPLv3+'

import copy
import logging
from pathlib import Path
import shutil

import numpy as np
from numpy.polynomial import Polynomial

from viperleed.calc.files import iofdopt
from viperleed.calc.files import parameters
from viperleed.calc.files import poscar
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.initialization import make_compile_logs_dir
from viperleed.calc.sections.refcalc import refcalc as section_refcalc
from viperleed.calc.sections.rfactor import rfactor as section_rfactor


_FD_COMMON_INPUT_FILES = (
    # Input files common to all the full-dynamic calculation runs.
    # Collected from the main work directory into the subfolders
    # for each parameter variation.
    'BEAMLIST',
    'PHASESHIFTS',
    )
logger = logging.getLogger(__name__)


class FullDynamicCalculationError(Exception):
    """Base class for exceptions raised by the full dynamic calculation."""


class FullDynamicOptimizationOutOfBoundsError(FullDynamicCalculationError):
    """Raised when the optimization is out of bounds."""


def get_fd_r(slab, rpars, workdir):
    """Run reference and R-factor calculations.

    Parameters
    ----------
    slab : Slab
        Object containing atomic configuration
    rpars : Rparams
        Object containing run parameters
    workdir : Path
        Path in which the calculations should be executed
        before returning to the current directory.

    Returns
    -------
    rfactor : float
        The overall R-factor obtained for the (`slab`, `rpars`)
        combination.
    rfactor_per_beam : list
        R-factor for each individual beam.
    """
    rpars.TENSOR_OUTPUT = [0]
    # internally transform theta, phi to within range
    if rpars.THETA < 0:
        rpars.THETA = abs(rpars.THETA)
        rpars.PHI += 180
    rpars.THETA = min(rpars.THETA, 90.)
    rpars.PHI = rpars.PHI % 360
    main_work = Path.cwd()
    with execute_in_dir(workdir, mkdir=True):
        make_compile_logs_dir(rpars)
        for input_file in _FD_COMMON_INPUT_FILES:
            try:
                shutil.copy2(main_work / input_file, workdir)
            except OSError:
                logger.error(f'Error copying {input_file} to subfolder.')
                raise
        logger.info('Starting full-dynamic calculation')
        try:
            section_refcalc(slab, rpars, parent_dir=main_work)
        except Exception:
            logger.error('Error running reference calculation')
            raise
        logger.info('Starting R-factor calculation...')
        try:
            rfaclist = section_rfactor(slab, rpars, 11)
        except Exception:
            logger.error('Error running R-factor calculation')
            raise
        return rpars.last_R, rfaclist


def apply_scaling(sl, rp, which, scale):
    scaling = (scale if char in which else 1 for char in 'abc')
    sl.update_fractional_from_cartesian()                                       # TODO: needed?
    sl.bulkslab.update_fractional_from_cartesian()                              # TODO: needed?
    sl.apply_scaling(*scaling)
    rp.BULK_REPEAT = sl.bulkslab.get_bulk_repeat(rp)


def fd_optimization(sl, rp):
    """
    Runs multiple consecutive reference calculations and r-factor calculations
    optimizing a parameter, as defined by OPTIMIZE in PARAMETERS.

    Parameters
    ----------
    sl : Slab
        Object containing atomic configuration
    rp : Rparams
        Object containing run parameters

    Returns
    -------
    None.
    """

    which = rp.OPTIMIZE["which"]
    if which == "none":
        logger.error("No parameter is defined for optimization. Set OPTIMIZE "
                     "in the PARAMETERS file to define behaviour.")
        rp.setHaltingLevel(2)
        return
    logger.info(f"Starting optimization of {rp.OPTIMIZE['which']}")
    # make sure there's a compiler ready, and we know the number of cores:
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise
    rp.updateCores()  # if number of cores is not defined, try to find it
    # check whether step is set; if not, choose default value
    if rp.OPTIMIZE["step"] == 0.:
        if which in ["v0i", "theta"]:
            rp.OPTIMIZE["step"] = 0.5
        elif which == "phi":
            rp.OPTIMIZE["step"] = 2.
        else:   # unit cell size
            rp.OPTIMIZE["step"] = 0.01
        logger.debug("Initial step size undefined, defaulting to {}"
                     .format(rp.OPTIMIZE["step"]))
    if rp.OPTIMIZE["maxstep"] == 0.:
        rp.OPTIMIZE["maxstep"] = 3 * rp.OPTIMIZE["step"]
    if rp.OPTIMIZE["convergence"] == 0.:
        rp.OPTIMIZE["convergence"] = 0.1 * rp.OPTIMIZE["step"]
    rp.OPTIMIZE["maxstep"] = abs(rp.OPTIMIZE["maxstep"])  # always positive
    rp.OPTIMIZE["convergence"] = abs(rp.OPTIMIZE["convergence"])
    # loop over required configurations, copying sl and rp every time
    known_points = np.empty((0, 2))
    # known_points will contain columns [x, r], with x the parameter value
    rfactor_lists = []    # for each point, store r-factors per beam
    tmpdirs = []
    if which == "v0i":
        x0 = rp.V0_IMAG
    elif which == "theta":
        x0 = rp.THETA
    elif which == "phi":
        x0 = rp.PHI
    else:
        x0 = 1.   # geometry: x is a scaling factor



    # optimization loop
    curvature_fail = 0
    parabola = None
    best_explicit_r = None
    while True:
        if rp.STOP:
            break
        # Make sure that there are no duplicate atoms. This is more
        # a safeguard for the future, if we allow also full-dynamic
        # optimization of atoms.
        sl.check_atom_collisions(rp.SYMMETRY_EPS)
        if len(known_points) == 0:
            x = x0   # x will be the value of the parameter under variation
        elif len(known_points) == 1:
            x = x0 + rp.OPTIMIZE["step"]   # first step is given
        elif len(known_points) == 2:  # get third point
            if known_points[1, 1] < known_points[0, 1]:
                x = x0 + 2*rp.OPTIMIZE["step"]
            else:
                x = x0 - rp.OPTIMIZE["step"]
        else:
            parabola = Polynomial.fit(known_points[:, 0], known_points[:, 1],
                                      2)
            coefs = parabola.convert(domain=[-1, 1]).coef
            current_scope = (min(known_points[:, 0]), max(known_points[:, 0]))
            # first check curvature
            if coefs[2] <= 0:  # TODO: some epsilon would be better
                # no good slope - keep adding points on either side
                curvature_fail += 1
                if (len([v for v in known_points[:, 0] if v > x0])
                        > len([v for v in known_points[:, 0] if v < x0])):
                    x = current_scope[0] - (abs(rp.OPTIMIZE["step"])
                                            * curvature_fail)
                else:
                    x = current_scope[1] + (abs(rp.OPTIMIZE["step"])
                                            * curvature_fail)
                logger.info("Current fit parabola is not concave, no minimum "
                            "predicted. Adding data point at {:.4f}".format(x))
            else:
                # try to use minimum
                new_min = -0.5*coefs[1] / coefs[2]
                x = new_min
                # if at the edge, better to step a bit farther again
                if x < current_scope[0]:
                    x = min(x, current_scope[0] - abs(rp.OPTIMIZE["step"]))
                elif x > current_scope[1]:
                    x = max(x, current_scope[1] + abs(rp.OPTIMIZE["step"]))
                # limit to maxstep
                x = min(x, current_scope[1] + rp.OPTIMIZE["maxstep"])
                x = max(x, current_scope[0] - rp.OPTIMIZE["maxstep"])
                if which == "v0i":
                    x = max(0, x)
                if which in "abc":
                    x = max(0.1, x)   # shouldn't happen, just in case
                # check whether we're close to point that is already known
                if any(abs(v - x) < rp.OPTIMIZE["convergence"]
                       for v in known_points[:, 0]):
                    left = [v for v in known_points[:, 0] if v < x]
                    right = [v for v in known_points[:, 0] if v > x]
                    if (not left or
                        all(abs(v - x) < rp.OPTIMIZE["convergence"]
                            for v in left)):     # take step to left
                        x = current_scope[0] - abs(rp.OPTIMIZE["step"])
                    elif (not right or
                          all(abs(v - x) < rp.OPTIMIZE["convergence"]
                              for v in right)):  # take step to right
                        x = current_scope[1] + abs(rp.OPTIMIZE["step"])
                    elif (all(abs(v - x) < rp.OPTIMIZE["convergence"]
                              for v in (max(left), min(right)))
                          and len(known_points) < rp.OPTIMIZE["minpoints"]):
                        # points are very close already, default to adding
                        #   points outside of scope
                        if (len(right) > len(left)):
                            x = current_scope[0] - abs(rp.OPTIMIZE["step"])
                        else:
                            x = current_scope[1] + abs(rp.OPTIMIZE["step"])
                    else:        # take midpoint
                        x = (max(left) + min(right)) / 2
                    # stop if minpoints and convergence
                    if len(known_points) >= rp.OPTIMIZE["minpoints"]:
                        if not (current_scope[0] < x < current_scope[1]):
                            logger.info(
                                "Minimum has converged, but is at the edge of "
                                "the tested range. Proceeding with "
                                "optimization.")
                        else:
                            lowest_dist = min([abs(v - new_min) for v
                                               in known_points[:, 0]])
                            closest_point = [
                                v for v in known_points[:, 0]
                                if abs(v - new_min) == lowest_dist][0]
                            logger.info(
                                "Stopping optimization: Predicted new minimum "
                                f"{new_min:.4f} is within convergence limits to "
                                f"known point {closest_point:.4f}.")
                            break  # within convergence limit
                logger.info("Currently predicting minimum at {:.4f} with R = "
                            "{:.4f}, adding data point at {:.4f}"
                            .format(new_min, parabola(new_min), x))

        # write out results
        if len(known_points) != 0:
            iofdopt.write_fd_opt_csv(known_points, which)
        if len(known_points) > 2:
            iofdopt.write_fd_opt_pdf(known_points, which, parabola=parabola)

        # create test objects tsl, trp and set parameters
        tsl = copy.deepcopy(sl)
        trp = copy.deepcopy(rp)
        if which == "v0i":
            x = max(0, x)
            trp.V0_IMAG = x
        elif which == "theta":
            trp.THETA = x
        elif which == "phi":
            trp.PHI = x
        else:       # geometry: x is a scaling factor for the unit cell
            x = max(0.1, x)
            apply_scaling(tsl, trp, which, x)

        # create subfolder and calculate there
        if isinstance(x, int):
            dname = f"{which}_{x}"
        else:
            dname = f"{which}_{x:.4f}"
        workdir = Path(dname).resolve()
        tmpdirs.append(workdir)
        logger.info(f"STARTING CALCULATION AT {which} = {x:.4f}")
        r, rfaclist = get_fd_r(tsl, trp, workdir)
        known_points = np.append(known_points, np.array([[x, r]]), 0)
        rfactor_lists.append(rfaclist)
        if not best_explicit_r or r < best_explicit_r[0]:
            best_explicit_r = trp.stored_R["refcalc"]    # has int & frac components
        # decide how to proceed
        if len(known_points) >= rp.OPTIMIZE["maxpoints"]:
            logger.warning(
                "Stopping optimization because the maximum number of data "
                "points has been reached. This may indicate that optimization "
                "did not fully converge.")
            break

    # optimization loop finished; re-fit points, analyze:
    parabola = Polynomial.fit(known_points[:, 0], known_points[:, 1], 2)
    coefs = parabola.convert(domain=[-1, 1]).coef
    new_min = -0.5*coefs[1] / coefs[2]
    logger.info(f"Optimization of {which}: Predicted minimum at "
                f"{new_min:.4f}, R = {parabola(new_min):.4f}")
    current_best = known_points[np.argmin(known_points, 0)[1]]
    rp.stored_R["refcalc"] = best_explicit_r
    if (round(parabola(new_min), 4) > round(current_best[1], 4)
            and current_best[0] != known_points[-1, 0]):
        logger.warning(
            "Fit prediction at {:.4f} is worse than explicit calculation "
            "result at {:.4f} (R = {:.4f}). Outputting best known result "
            "instead of parabolic prediction.".format(new_min, *current_best))
        new_min = current_best[0]

    # output analysis
    iofdopt.write_fd_opt_csv(known_points, which)
    if len(known_points) > 2:
        iofdopt.write_fd_opt_pdf(known_points, which, parabola=parabola)

    # output modified files
    comment = "Found by full-dynamic optimization"
    if which == "v0i":
        rp.V0_IMAG = new_min
        parameters.modify(rp, "V0_IMAG", comment=comment)
    elif which in ("theta", "phi"):
        setattr(rp, which.upper(), new_min)
        if rp.THETA < 0:
            rp.THETA = abs(rp.THETA)
            rp.PHI += 180
        rp.PHI = rp.PHI % 360
        parameters.modify(rp, "BEAM_INCIDENCE", comment=comment)
    else:       # geometry: x is a scaling factor for the unit cell
        apply_scaling(sl, rp, which, new_min)
        if not isinstance(rp.BULK_REPEAT, float) or "c" in which:
            parameters.modify(rp, "BULK_REPEAT", comment=comment)
        poscar.write(sl, filename='POSCAR', comments='all')
        rp.files_to_out.add('POSCAR')

    # fetch I(V) data from all, plot together
    best_rfactors = rfactor_lists[np.argmin(known_points[:, 1])]
    try:
        iofdopt.write_fd_opt_beams_pdf(rp, known_points, which, tmpdirs,
                                       best_rfactors)
    except Exception as exc:
        logger.warning(f"Failed to plot I(V) curves from optimization: {exc}")

    # clean up tmpdirs
    # for path in tmpdirs:
    #     try:
    #         shutil.rmtree(path)
    #     except Exception as e:
    #         logger.warning("Failed to delete temporary directory {}: "
    #                        .format(path) + str(e))
