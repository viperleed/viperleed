# -*- coding: utf-8 -*-

"""
Created on Oct 22 2021

@author: Florian Kraushofer

Tensor LEED Manager section Full-dynamic Optimization
"""
from abc import ABC, abstractmethod
import copy
import logging
import os
from pathlib import Path
import shutil

import numpy as np
from numpy.polynomial import Polynomial

from viperleed.tleedmlib import psgen
from viperleed.tleedmlib.files import iofdopt as tl_io
from viperleed.tleedmlib.files.parameters import modifyPARAMETERS
from viperleed.tleedmlib.files.poscar import writePOSCAR
from viperleed.tleedmlib.sections.refcalc import refcalc as section_refcalc
from viperleed.tleedmlib.sections.rfactor import rfactor as section_rfactor


logger = logging.getLogger("tleedm.fdopt")


class FullDynamicCalculationError(Exception):
    """Base class for exceptions raised by the full dynamic calculation."""

    def __init__(self, message):
        super().__init__(message)

class FullDynamicOptimizationOutOfBoundsError(FullDynamicCalculationError):
    """Raised when the optimization is out of bounds."""
    def __init__(self, message):
        super().__init__(message)

def get_fd_r(sl, rp, work_dir=Path(), home_dir=Path()):
    """
    Runs reference calculation and r-factor calculation, returns R.

    Parameters
    ----------
    sl : Slab
        Object containing atomic configuration
    rp : Rparams
        Object containing run parameters
    work_dir : pathlike, optional
        The directory to execute the calculations in. The default is ".".
    home_dir : pathlike, optional
        The directory to return to after finishing. The default is the current
        working directory.

    Returns
    -------
    rfactor : float
        The r-factor obtained for the sl, rp combination
    """
    rp.TENSOR_OUTPUT = [0]
    rp.workdir = Path(work_dir)
    # internally transform theta, phi to within range
    if rp.THETA < 0:
        rp.THETA = abs(rp.THETA)
        rp.PHI += 180
    rp.THETA = min(rp.THETA, 90.)
    rp.PHI = rp.PHI % 360
    # create work directory, go there, execute
    try:
        if work_dir != ".":
            os.makedirs(work_dir, exist_ok=True)
            os.chdir(work_dir)
        for fn in ["BEAMLIST", "PHASESHIFTS"]:
            try:
                shutil.copy2(home_dir / fn, fn)
            except Exception:
                logger.error(f"Error copying {fn} to subfolder.")
                raise
        logger.info("Starting full-dynamic calculation")
        try:
            section_refcalc(sl, rp, parent_dir=home_dir)
        except Exception:
            logger.error("Error running reference calculation")
            raise
        logger.info("Starting R-factor calculation...")
        try:
            rfaclist = section_rfactor(sl, rp, 11)
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
            tl_io.write_fd_opt_csv(known_points, which)
        if len(known_points) > 2:
            tl_io.write_fd_opt_pdf(known_points, which, parabola=parabola)

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
        workdir = rp.workdir / dname
        tmpdirs.append(workdir)
        logger.info(f"STARTING CALCULATION AT {which} = {x:.4f}")
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
    logger.info(f"Optimization of {which}: Predicted minimum at "
                f"{new_min:.4f}, R = {parabola(new_min):.4f}")
    current_best = known_points[np.argmin(known_points, 0)[1]]
    if (round(parabola(new_min), 4) > round(current_best[1], 4)
            and current_best[0] != known_points[-1, 0]):
        logger.warning(
            "Fit prediction at {:.4f} is worse than explicit calculation "
            "result at {:.4f} (R = {:.4f}). Outputting best known result "
            "instead of parabolic prediction.".format(new_min, *current_best))
        new_min = current_best[0]

    # output analysis
    tl_io.write_fd_opt_csv(known_points, which)
    if len(known_points) > 2:
        tl_io.write_fd_opt_pdf(known_points, which, parabola=parabola)

    # output modified files
    comment = "Found by full-dynamic optimization"
    if which == "v0i":
        rp.V0_IMAG = new_min
        modifyPARAMETERS(rp, "V0_IMAG", new=f"{new_min:.4f}", comment=comment)
    elif which in ("theta", "phi"):
        setattr(rp, which.upper(), new_min)
        if rp.THETA < 0:
            rp.THETA = abs(rp.THETA)
            rp.PHI += 180
        rp.PHI = rp.PHI % 360
        modifyPARAMETERS(rp, "BEAM_INCIDENCE",
                         new=f"THETA {rp.THETA:.4f}, PHI {rp.PHI:.4f}",
                         comment=comment)
    else:       # geometry: x is a scaling factor for the unit cell
        apply_scaling(sl, rp, which, new_min)
        if not isinstance(rp.BULK_REPEAT, float) or "c" in which:
            vec_str = "[{:.5f} {:.5f} {:.5f}]".format(*rp.BULK_REPEAT)
            modifyPARAMETERS(rp, "BULK_REPEAT", new=vec_str, comment=comment)
        writePOSCAR(sl, filename=f"POSCAR_OUT_{rp.timestamp}", comments="all")

    # fetch I(V) data from all, plot together
    best_rfactors = rfactor_lists[np.argmin(known_points[:, 1])]
    try:
        tl_io.write_fd_opt_beams_pdf(rp, known_points, which, tmpdirs,
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


class FDParameter():
    """Base class for parameters accessible via full dynamic calculation."""

    def __init__(self, name, bounds, transform_func):
        self.name = name
        self.bounds = bounds
        self.transform_func = transform_func

    def apply(self, rparams, slab, val):
        """Apply the parameter to the given rparams and slab."""
        self.transform_func(rparams, slab, val)

fd_param_v0i = FDParameter("v0i", (0, np.inf), lambda r, s, v: setattr(r, "V0_IMAG", v))
fd_param_theta = FDParameter("theta", (-90, 90), lambda r, s, v: setattr(r, "THETA", v))



class OneDimensionalFDOptimizer(ABC):                                           # TODO: extend to n dimensions
    """Base class for full-dynamic optimization."""

    def __init__(self, eval_func, x, fd_parameter, convergence, R=None):
        self.eval_func = eval_func # callback function to evaluate
        self.x = x # points to evaluate initially
        self.convergence = convergence # convergence criterion
        self.param_name = fd_parameter.name # name of parameter to optimize

        # bounds for x
        self.bounds = (fd_parameter.bounds if fd_parameter.bounds is not None
                       else (-np.inf, np.inf))

        if not all(self._in_bounds(x) for x in self.x):
            raise FullDynamicOptimizationOutOfBoundsError(
                "Initial points must be within bounds"
            )

        if not (isinstance(R, list) or R is None):
            raise TypeError("R must be a list or None")

        # results of evaluation (R-factors)
        self.R = R if R is not None else [self.eval_func(x) for x in self.x]

        self.evaluate_missing_points()

    @property
    def x_R_min(self):
        if any(R is None for R in self.R):
            raise ValueError("Cannot find minimum if not requested points "
                             "have been evaluated")
        return min(zip(self.x, self.R), key=lambda x: x[1])

    @property
    def x_R_sorted(self):
        _x_R_sorted = sorted(zip(self.x, self.R), key=lambda x: x[0])
        _x_sorted = [x for x, _ in _x_R_sorted]
        _R_sorted = [R for _, R in _x_R_sorted]
        return _x_sorted, _R_sorted


    def closest_sample_point(self, compare_to):
        return min(self.x, key=lambda x: abs(x - compare_to))

    def _in_bounds(self, x_val):
        return self.bounds[0] <= x_val <= self.bounds[1]

    def evaluate_missing_points(self):
        for index, (x, R) in enumerate(zip(self.x, self.R)):
            if R is None:
                self.R[index] = self.eval_func(x)

    @abstractmethod
    def optimize(self):
        pass

    @abstractmethod
    def finalize(self):
        """Perform final steps after optimization has converged."""
        pass

    @abstractmethod
    def updated_intermediate_output(self):
        pass



class SingleParameterParabolaFit(OneDimensionalFDOptimizer):
    def __init__(self, x, eval_func, fd_parameter, convergence, initial_step, min_points, max_points, R=None, bounds= None, max_step=None):
        super().__init__(eval_func, x, fd_parameter, convergence, R, bounds)
        self.min_points = min_points
        self.max_points = max_points
        self.initial_step = initial_step
        self.max_step = max_step if max_step is not None else np.inf
        self.curvature_fails = 0

        if len(self.x) == 0:
            raise RuntimeError("Must provide at least one evaluation point to "
                               "start.")

    def optimize(self):
        while len(self.x) < self.max_points:
            logger.log(5, "Fitting parabola to points:")
            logger.log(5, f"x: [{' ,'.join(f'{x:.4f}' for x in self.x)}]")
            logger.log(5, f"R: [{' ,'.join(f'{R:.4f}' for R in self.R)}]")

            next_x = self.next_point()

            if next_x is None:
                logger.info(
                    "Stopping optimization: Predicted new minimum "
                    f"{self._predicted_min_x:.4f} is within convergence limits "
                    "to known point "
                    f"{self.closest_sample_point(self._predicted_min_x):.4f}."
                    )
                return

            if len(self.x) >= self.min_points:
                logger.debug(
                    "Currently predicting minimum at "
                    f"{self._predicted_min_x:.4f} "
                    f"with R = {self._predicted_R:.4f}, adding data point "
                    f"at {next_x:.4f}."
                )

            next_R = self.eval_func(next_x)
            self.x.append(next_x)
            self.R.append(next_R)

            # update plot and csv
            self.updated_intermediate_output()

        # max points reached
        logger.warning(
            "Stopping optimization because the maximum number of data "
            "points has been reached. This may indicate that optimization "
            "did not fully converge.")
        return

    @property
    def _predicted_min_x(self):
        self._fit_parabola()
        return -0.5*self._parabola_linear / self._parabola_curvature

    @property
    def _predicted_R(self):
        return self.parabola(self._predicted_min_x)  # implicit _fit_parabola

    @property
    def _parabola_curvature(self):
        self._fit_parabola()
        return self.coefficients[2]

    @property
    def _parabola_linear(self):
        self._fit_parabola()
        return self.coefficients[1]

    def _fit_parabola(self):
        if len(self.x) < 3:
            raise RuntimeError("Cannot fit parabola until at least 3 points "
                               "have been evaluated.")
        self.parabola = Polynomial.fit(self.x, self.R, 2)
        self.coefficients = self.parabola.convert(domain=[-1, 1]).coef

    def next_point(self):
        if any([R is None for R in self.R]):
            raise RuntimeError("Cannot predict next point until all points "
                               "have been evaluated.")
        if len(self.x) == 1:
            return self.x[0] + self.initial_step
        elif len(self.x) == 2:
            if (self.R[0] < self.R[0]
                and (self.x[0] + 2*self.initial_step) > self.convergence):
                return self.x[0] + 2*self.initial_step
            else:
                return self.x[0] - self.initial_step

        # Otherwise we have enough points to fit a parabola

        # if curvature is convex, add point outside current scope
        if self._parabola_curvature < 0:
            self.curvature_fails += 1
            next_point_convex = self._point_outside()
            logger.info("Current fit parabola is not concave, no minimum "
                    f"predicted. Adding data point at {next_point_convex:.4f}")
            return next_point_convex

        # if predicted minimum is outside scope, go a bit further
        if self._predicted_min_x <= min(self.x):
            return min(self.x) - (abs(self.initial_step))
        elif self._predicted_min_x >= max(self.x):
            return max(self.x) + (abs(self.initial_step))
        else:
            min_x = self._predicted_min_x

        # have we tried this point before?
        if not any([abs(x - min_x) < self.convergence for x in self.x]):
            return min_x

        # next, ensure we have at least one point either side of the minimum
        left_of_min = tuple(x < min_x for x in self.x)
        all_left = all(left_of_min)
        if all_left or all(abs(x-min_x) < self.convergence for x in self.x):
            return min_x + abs(self.initial_step)
        right_of_min = tuple(x > min_x for x in self.x)
        all_right = all(right_of_min)
        if all_right or all(abs(x-min_x) < self.convergence for x in self.x):
            return min_x - abs(self.initial_step)

        # check if we have converged
        if self._has_converged():
            return None # we're done

        # check midpoint of current points left & right of predicted minimum
        midpoint = (max(left_of_min) + min(right_of_min)) / 2
        if abs(midpoint - min_x) > self.convergence:
            return midpoint

        # otherwise, as a last resort, add a point outside the current scope
        return self._point_outside()

    def _point_outside(self):
        if (len([v for v in self.x if v > np.median(self.x)])
                > len([v for v in self.x if v < np.median(self.x)])):
            next_point = min(self.x) - (abs(self.initial_step)
                                             * self.curvature_fails)
        else:
            next_point = max(self.x) + (abs(self.initial_step)
                                             * self.curvature_fails)
        # limit to max_step
        next_point = min(next_point, max(self.x) + self.max_step)
        next_point = max(next_point, min(self.x) - self.max_step)

        # check if point is in bounds, otherwise cap it
        next_point = min(next_point, self.bounds[1])
        next_point = max(next_point, self.bounds[0])
        return next_point

    def _has_converged(self):
        if any([R is None for R in self.R]):
            raise RuntimeError("Cannot predict check for convergence until "
                               "all points are evaluated.")

        if len(self.x) < self.min_points:
            return False

        sorted_x, sorted_R = self.x_R_sorted
        min_id = np.argmin(sorted_R)
        min_x = sorted_x[min_id]

        # check if there is a point within convergence range of the minimum
        if any([abs(x - min_x) < self.convergence for x in self.x]):
            return True

        return False

    def finalize(self):
        # check if predicted minimum is worse than any of the points
        if self._predicted_R > min(self.R):
            logger.info(
                f"Fit prediction at {self._predicted_min_x:.4f} is worse than "
                f"explicit calculation at {min(self.x):.4f}. Outputting "
                f"best known result instead of parabolic prediction."
            )
            return self._predicted_min_x
        # else return the predicted minimum
        return self._predicted_min_x, self._predicted_R

    def updated_intermediate_output(self):
        if len(self.x) > 0:
            tl_io.write_fd_opt_csv(np.array([self.x, self.R]).T, which=self.param_name)
        if len(self.x) >= 2:
            tl_io.write_fd_opt_plot(np.array([self.x, self.R]).T, which=self.param_name)


def evaluate_fd_calculation(rp, sl, fd_parameters, parameter_vals):
    # make temporary Rparams and slab
    test_sl = copy.deepcopy(sl)
    test_rp = copy.deepcopy(rp)

    # make temporary workdir
    name = _temp_dir_name(fd_parameters, parameter_vals)
    temp_dir = _fd_dir(name, rp)

    # apply parameters
    for param, val in zip(fd_parameters, parameter_vals):
        param.apply(test_rp, test_sl, val)

    # run calculation
    r, _ = get_fd_r(test_sl, test_rp,
                           work_dir=temp_dir, home_dir=rp.workdir)

    # TODO: update csv and plot
    return r



def _fd_dir(name, rp):
    name  = "fd_test"
    fd_work_dir = rp.workdir / "fd_calc" / name
    return fd_work_dir


def _temp_dir_name(parameters, parameter_vals):
    return "_".join([f"{p.name}={v}" for p, v in zip(parameters, parameter_vals)])

