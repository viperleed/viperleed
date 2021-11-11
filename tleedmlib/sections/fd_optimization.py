# -*- coding: utf-8 -*-

"""
Created on Oct 22 2021

@author: Florian Kraushofer

Tensor LEED Manager section Full-dynamic Optimization
"""

import os
import logging
import copy
import shutil
import numpy as np
from numpy.polynomial import Polynomial

import viperleed.tleedmlib as tl
import viperleed.tleedmlib.files.iofdopt as io
import viperleed.tleedmlib.psgen as psgen
from viperleed.tleedmlib.files.parameters import modifyPARAMETERS


logger = logging.getLogger("tleedm.fdopt")


def get_fd_r(sl, rp, work_dir=".", home_dir=""):
    """
    Runs reference calculation and r-factor calculation, returns R.

    Parameters
    ----------
    sl : Slab
        Object containing atomic configuration
    rp : Rparams
        Object containing run parameters
    work_dir : str, optional
        The directory to execute the calculations in. The default is ".".
    home_dir : str, optional
        The directory to return to after finishing. The default is the current
        working directory.

    Returns
    -------
    rfactor : float
        The r-factor obtained for the sl, rp combination

    """
    if not home_dir:
        home_dir = os.getcwd()
    rp.TENSOR_OUTPUT = [0]
    rp.workdir = work_dir
    # internally transform theta, phi to within range
    if rp.THETA < 0:
        rp.THETA = abs(rp.THETA)
        rp.PHI += 180
    rp.THETA = min(rp.THETA, 90.)
    rp.PHI = rp.PHI % 360
    # create work directoy, go there, execute
    try:
        if work_dir != ".":
            os.makedirs(work_dir, exist_ok=True)
            os.chdir(work_dir)
        for fn in ["BEAMLIST", "PHASESHIFTS"]:
            try:
                shutil.copy2(os.path.join(home_dir, fn), fn)
            except Exception:
                logger.error("Error copying " + fn + " to subfolder.")
                raise
        try:
            logger.info("Starting full-dynamic calculation")
            tl.sections.refcalc(sl, rp, parent_dir=home_dir)
        except Exception:
            logger.error("Error running reference calculation")
            raise
        try:
            logger.info("Starting r-factor calculation")
            rfaclist = tl.sections.rfactor(sl, rp, 11)
        except Exception:
            logger.error("Error running rfactor calculation")
            raise
    finally:
        os.chdir(home_dir)
    return rp.last_R, rfaclist


def apply_scaling(sl, rp, which, scale):
    m = np.eye(3)
    if "a" in which:
        m[0, 0] *= scale
    if "b" in which:
        m[1, 1] *= scale
    if "c" in which:
        m[2, 2] *= scale
    sl.getFractionalCoordinates()
    sl.ucell = np.dot(sl.ucell, m)
    sl.getCartesianCoordinates(updateOrigin=True)
    sl.bulkslab.getFractionalCoordinates()
    sl.bulkslab.ucell = np.dot(sl.bulkslab.ucell, m)
    sl.bulkslab.getCartesianCoordinates()
    if type(rp.BULK_REPEAT) == float:
        rp.BULK_REPEAT *= scale
    elif rp.BULK_REPEAT is not None:
        rp.BULK_REPEAT = np.dot(rp.BULK_REPEAT, m)


def get_fd_phaseshifts(sl, rp, S_ovl):
    rundgrenpath = os.path.join('tensorleed', 'EEASiSSS.x')
    (first_line, phaseshifts) = psgen.runPhaseshiftGen(tsl, trp, psgensource=rundgrenpath)

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
    logger.info("Starting optimization of {}".format(rp.OPTIMIZE["which"]))
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
        elif which == "S_ovl":
            rp.OPTIMIZE["step"] = 0.05
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
    elif which == "s_ovl":
        x0 = rp.S_OVL
    else:
        x0 = 1.   # geometry: x is a scaling factor

    # optimization loop
    curvature_fail = 0
    parabola = None
    while True:
        if rp.STOP:
            break
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
                if which in ["a", "b", "ab", "c", "abc"]:
                    x = max(0.1, x)   # shouldn't happen, just in case
                if which == "S_ovl":
                    x = min(1, x)
                    x = max(0, x)
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
                                "{:.4f} is within convergence limits to known "
                                "point {:.4f}".format(new_min, closest_point))
                            break  # within convergence limit
                logger.info("Currently predicting minimum at {:.4f} with R = "
                            "{:.4f}, adding data point at {:.4f}"
                            .format(new_min, parabola(new_min), x))

        # write out results
        if len(known_points) != 0:
            io.write_fd_opt_csv(known_points, which)
        if len(known_points) > 2:
            io.write_fd_opt_pdf(known_points, which, parabola=parabola)

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
        elif which == "S_ovl":
            trp.S_OVL = x
            get_fd_phaseshifts(sl, rp, x)
        else:       # geometry: x is a scaling factor for the unit cell
            x = max(0.1, x)
            apply_scaling(tsl, trp, which, x)

        # create subfolder and calculate there
        if type(x) == int:
            dname = which + "_{}".format(x)
        else:
            dname = which + "_{:.4f}".format(x)
        workdir = os.path.join(rp.workdir, dname)
        tmpdirs.append(workdir)
        logger.info("STARTING CALCULATION AT {} = {:.4f}".format(which, x))
        r, rfaclist = get_fd_r(tsl, trp, work_dir=workdir, home_dir=rp.workdir)
        known_points = np.append(known_points, np.array([[x, r]]), 0)
        rfactor_lists.append(rfaclist)
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
    logger.info("Optimization of {}: Predicted minimum at {:.4f}, R = {:.4f}"
                .format(which, new_min, parabola(new_min)))
    current_best = known_points[np.argmin(known_points, 0)[1]]
    if (round(parabola(new_min), 4) > round(current_best[1], 4)
            and current_best[0] != known_points[-1, 0]):
        logger.warning(
            "Fit prediction at {:.4f} is worse than explicit calculation "
            "result at {:.4f} (R = {:.4f}). Outputting best known result "
            "instead of parabolic prediction.".format(new_min, *current_best))
        new_min = current_best[0]

    # output analysis
    io.write_fd_opt_csv(known_points, which)
    if len(known_points) > 2:
        io.write_fd_opt_pdf(known_points, which, parabola=parabola)

    # output modified files
    comment = "Found by full-dynamic optimization"
    if which == "v0i":
        rp.V0_IMAG = new_min
        modifyPARAMETERS(rp, "V0_IMAG", new="{:.4f}".format(new_min),
                         comment=comment)
    elif which in ("theta", "phi"):
        if which == "theta":
            rp.THETA = new_min
        else:
            rp.PHI = new_min
        if rp.THETA < 0:
            rp.THETA = abs(rp.THETA)
            rp.PHI += 180
        rp.PHI = rp.PHI % 360
        modifyPARAMETERS(rp, "BEAM_INCIDENCE", new=("THETA {:.4f}, PHI {:.4f}"
                                                    .format(rp.THETA, rp.PHI)),
                         comment=comment)
    elif which == "S_ovl":
        rp.S_ovl = new_min
        modifyPARAMETERS(rp, "S_OVL", "{:.4f}".format(new_min), comment=comment)
    else:       # geometry: x is a scaling factor for the unit cell
        apply_scaling(sl, rp, which, new_min)
        if type(rp.BULK_REPEAT) != float or "c" in which:
            vec_str = "[{:.5f} {:.5f} {:.5f}]".format(*rp.BULK_REPEAT)
            modifyPARAMETERS(rp, "BULK_REPEAT", new=vec_str, comment=comment)
        tl.files.poscar.writePOSCAR(
            sl, filename="POSCAR_OUT_" + rp.timestamp, comments="all")

    # fetch I(V) data from all, plot together
    best_rfactors = rfactor_lists[np.argmin(known_points[:, 1])]
    try:
        io.write_fd_opt_beams_pdf(rp, known_points, which, tmpdirs,
                                  best_rfactors)
    except Exception as e:
        logger.warning("Failed to plot I(V) curves from optimization: "
                       + str(e))

    # clean up tmpdirs
    # for path in tmpdirs:
    #     try:
    #         shutil.rmtree(path)
    #     except Exception as e:
    #         logger.warning("Failed to delete temporary directory {}: "
    #                        .format(path) + str(e))
    return
