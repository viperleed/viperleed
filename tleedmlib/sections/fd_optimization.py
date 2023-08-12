# -*- coding: utf-8 -*-

"""
Created on Oct 22 2021

@author: Alexander M. Imre
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
import scipy.optimize

from viperleed.tleedmlib import psgen
from viperleed.tleedmlib.files import iofdopt as tl_io
from viperleed.tleedmlib.files.parameters import modifyPARAMETERS
from viperleed.tleedmlib.files.poscar import writePOSCAR
from viperleed.tleedmlib.sections.refcalc import refcalc as section_refcalc
from viperleed.tleedmlib.sections.rfactor import rfactor as section_rfactor


logger = logging.getLogger("tleedm.fdopt")

FD_PARAMETERS = {
    'v0i': {
        'bounds': (0,np.inf),
        'eval': lambda r, s, v: setattr(r, "V0_IMAG", v),
        'x0': lambda rp: rp.V0_IMAG,
    },
    'theta': {
        'bounds': (0,90),
        'eval': lambda r, s, v: setattr(r, "THETA", v),
        'x0': lambda rp: rp.THETA,
    },
    'phi': {
        'bounds': (0,360),
        'eval': lambda r, s, v: setattr(r, "PHI", v),
        'x0': lambda rp: rp.PHI,
    },
}

AVAILABLE_MINIMIZERS = (
    'Nelder-Mead',
    'Powell',
    'Newton-CG',
)

for scaling in ('a', 'b', 'c', 'ab', 'bc', 'abc'): # scaling of lattice vectors
    FD_PARAMETERS[scaling] = {
        'bounds': (0.1, 10),
        'eval': lambda r, s, v: apply_scaling(s, r, scaling, v),
        'x0': lambda rp: 1,
    }


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

    # TODO: update theses to use the new parameters

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

    # get FDParameter object
    fd_param = FDParameter(
        which,
        FD_PARAMETERS[which]["x0"](rp),
        FD_PARAMETERS[which]["bounds"],
        FD_PARAMETERS[which]["eval"],
    )

    fd_evaluate = lambda x: evaluate_fd_calculation(rp, sl, fd_param, x)
    optimizer = SingleParameterParabolaFit(
        x0=[fd_param.start_value,],
        eval_func=fd_evaluate,
        fd_parameter=fd_param,
        convergence=rp.OPTIMIZE["convergence"],
        step=rp.OPTIMIZE["step"],
        max_points=rp.OPTIMIZE["maxpoints"],
        min_points=rp.OPTIMIZE["minpoints"],
    )

    # perform optimization
    optimizer.optimize()
    opt_x, opt_R = optimizer.finalize()


    logger.info(
        f"Optimization finished. Best {which} = {opt_x:.4f} with "
        f"R = {opt_R:.4f}")

    # update PARAMETERS and POSCAR_OUT
    store_fd_param_to_file(sl, rp, fd_param.name, opt_x)

    # fetch I(V) from all and plot together
    known_points = np.array([optimizer.x, optimizer.R]).T
    try:
        tl_io.write_fd_opt_beams_pdf(rp, known_points,
                                     fd_param.name,
                                     optimizer.eval_dirs,
                                     optimizer.x_R_min[2])  # best R per beam
    except Exception as exc:
        logger.warning(f"Failed to plot I(V) curves from optimization: {exc}")



def store_fd_param_to_file(sl, rp, which, new_min):
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

    # clean up tmpdirs
    # for path in tmpdirs:
    #     try:
    #         shutil.rmtree(path)
    #     except Exception as e:
    #         logger.warning("Failed to delete temporary directory {}: "
    #                        .format(path) + str(e))


class FDParameter():
    """Base class for parameters accessible via full dynamic calculation."""

    def __init__(self, name, start_value, bounds, transform_func):
        self.name = name
        self.start_value = start_value
        # make sure start value is within bounds
        if start_value < bounds[0] or start_value > bounds[1]:
            raise ValueError(f"Start value for {start_value} is not within "
                             "bounds {bounds}")
        self.bounds = bounds

        self.transform_func = transform_func

    def apply(self, rparams, slab, val):
        """Apply the parameter to the given rparams and slab."""
        self.transform_func(rparams, slab, val)


class OneDimensionalFDOptimizer(ABC):                                           # TODO: extend to n dimensions
    """Base class for full-dynamic optimization."""

    def __init__(self, eval_func, x0, fd_parameter):
        self.eval_func = eval_func # callback function to evaluate
        initial_x = x0 if isinstance(x0, list) else [x0]
        self.x = [] # points to evaluate initially
        self.R, self.R_per_beam = [], [] # results of evaluations
        self.eval_dirs = [] # directories where evaluations are preformed
        self.param_name = fd_parameter.name # name of parameter to optimize

        # bounds for x
        self.bounds = (fd_parameter.bounds if fd_parameter.bounds is not None
                       else (-np.inf, np.inf))

        if not all(self._in_bounds(x) for x in self.x):
            raise FullDynamicOptimizationOutOfBoundsError(
                "Initial points must be within bounds"
            )

        # evaluate all requested initial points
        self._evaluate_initials(x0)

    @property
    def x_R_min(self):
        if any(R is None for R in self.R):
            raise ValueError("Cannot find minimum if not requested points "
                             "have been evaluated")
        return min(zip(self.x, self.R, self.R_per_beam), key=lambda x: x[1])

    @property
    def x_R_sorted(self):
        _x_R_sorted = sorted(zip(self.x, self.R), key=lambda x: x[0])
        _x_sorted = [x for x, _ in _x_R_sorted]
        _R_sorted = [R for _, R in _x_R_sorted]
        return _x_sorted, _R_sorted

    def closest_sample_point(self, compare_to):
        return min(self.x, key=lambda x: abs(x - compare_to))

    def evaluate(self, x_val):
        new_R, R_per_beam, eval_dir = self.eval_func(x_val)
        self.x.append(float(x_val))
        self.R.append(new_R)
        self.R_per_beam.append(R_per_beam)
        self.eval_dirs.append(eval_dir)
        return new_R

    def _evaluate_initials(self, initial_x):
        for x_val in initial_x:
            self.evaluate(x_val)

    def _in_bounds(self, x_val):
        return self.bounds[0] <= x_val <= self.bounds[1]

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
    def __init__(self, x0, eval_func, fd_parameter, convergence, step, min_points, max_points, max_step=None):
        super().__init__(eval_func, x0, fd_parameter)
        self.convergence = convergence
        self.min_points = min_points
        self.max_points = max_points
        self.step = step
        self.max_step = max_step if max_step is not None else np.inf
        self.curvature_fails = 0

        if len(self.x) == 0:
            raise RuntimeError("Must provide at least one evaluation point to "
                               "start.")

    def optimize(self):
        while len(self.x) < self.max_points:
            logger.log(5, "Fitting parabola to points:")
            logger.log(5, f"x: [{', '.join(f'{x:.4f}' for x in self.x)}]")
            logger.log(5, f"R: [{', '.join(f'{R:.4f}' for R in self.R)}]")

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

            self.evaluate(next_x)  # updates self.x, self.R

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
            return self.x[0] + self.step
        elif len(self.x) == 2:
            if (self.R[0] < self.R[0]
                and (self.x[0] + 2*self.step) > self.convergence):
                return self.x[0] + 2*self.step
            else:
                return self.x[0] - self.step

        # Otherwise we have enough points to fit a parabola

        # if curvature is convex, add point outside current scope
        if self._parabola_curvature < 0:
            self.curvature_fails += 1
            next_point_convex = self._point_outside()
            logger.info("Current fit parabola is not concave, no minimum "
                    f"predicted. Adding data point at {next_point_convex:.4f}")
            return next_point_convex

        # if predicted minimum is outside scope, potentially go a bit further
        if self._predicted_min_x <= min(self.x):
            return min(self._predicted_min_x,
                       min(self.x) - (abs(self.step)))
        elif self._predicted_min_x >= max(self.x):
            return max(self._predicted_min_x,
                       max(self.x) + (abs(self.step)))
        else:
            min_x = self._predicted_min_x

        # have we tried this point before?
        if not any([abs(x - min_x) < self.convergence for x in self.x]):
            return min_x

        # next, ensure we have at least one point either side of the minimum
        left_of_min = tuple(x < min_x for x in self.x)
        all_left = all(left_of_min)
        if all_left or all(abs(x-min_x) < self.convergence for x in self.x):
            return min_x + abs(self.step)
        right_of_min = tuple(x > min_x for x in self.x)
        all_right = all(right_of_min)
        if all_right or all(abs(x-min_x) < self.convergence for x in self.x):
            return min_x - abs(self.step)

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
            next_point = min(self.x) - (abs(self.step)
                                             * self.curvature_fails)
        else:
            next_point = max(self.x) + (abs(self.step)
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
        # plot for a final time
        known_points = np.array([self.x, self.R]).T
        tl_io.write_fd_opt_csv(known_points, "v0i")
        tl_io.write_fd_opt_pdf(known_points, "v0i", parabola=self.parabola)

        # check if predicted minimum is worse than any of the points
        if self._predicted_R > min(self.R):
            logger.info(
                f"Fit prediction at {self._predicted_min_x:.4f} is worse than "
                f"explicit calculation at {min(self.x):.4f}. Outputting "
                f"best known result instead of parabolic prediction."
            )
            return min(self.x), min(self.R)

        # else return the predicted minimum
        return self._predicted_min_x, self._predicted_R

    def updated_intermediate_output(self):
        if len(self.x) > 0:
            tl_io.write_fd_opt_csv(np.array([self.x, self.R]).T, which=self.param_name)
        if len(self.x) >= 2:
            tl_io.write_fd_opt_pdf(np.array([self.x, self.R]).T, which=self.param_name)


class SingleParameterBruteForceOptimiser(OneDimensionalFDOptimizer):

    def __init__(self, eval_func, fd_parameter, min_val, max_val, steps):
        super().__init__(eval_func, [], fd_parameter)
        self.min_val = min_val
        self.max_val = max_val
        self.steps = steps

        # make sure requested min/max are within bounds
        if not (self._in_bounds(min_val) and self._in_bounds(max_val)):
            raise ValueError("Requested min/max values are outside of "
                             f"parameter bounds for {fd_parameter.name}.")

    def optimize(self):
        self.res = scipy.optimize.brute(
            func=self.evaluate,
            ranges=((self.min_val, self.max_val),),
            Ns=self.steps,
            full_output=True,
            finish=None,
            disp=logger.level <= logging.DEBUG,
            workers=1
        )


    def finalize(self):
        # write csv output
        known_points = np.array([self.x, self.R]).T
        tl_io.write_fd_opt_csv(known_points, "v0i")
        x_opt, R_opt, _, _ = self.res
        return x_opt, R_opt

    def updated_intermediate_output(self):
        pass


class SingleParameterMinimizer(OneDimensionalFDOptimizer):                      # TODO: can implement a callback method in the minimizer
    def __init__(self, eval_func, x0, fd_parameter, minimizer_method, tol=None):
        super().__init__(eval_func, [], fd_parameter)

        if minimizer_method not in AVAILABLE_MINIMIZERS:
            raise ValueError(f"Minimizer method {minimizer_method} not "
                             f"available. Available methods are: "
                             f"{AVAILABLE_MINIMIZERS}")
        self.minimizer_method = minimizer_method
        self.x0 = float(x0)
        self.tol = tol

    def optimize(self):
        self.res = scipy.optimize.minimize(
            fun=self.evaluate,
            x0=self.x0,
            method=self.minimizer_method,
            bounds=(self.bounds,),
            options={"disp": logger.level <= logging.DEBUG},
            tol=self.tol
        )

    def finalize(self):
        if not self.res.success:
            logger.warning("Minimization failed with message: "
                           f"{self.res.message}")
        known_points = np.array([self.x, self.R]).T
        logger.info(f"Minimizer {self.minimizer_method} finished successfully "
                    f"in {self.res.nit} iterations.")
        tl_io.write_fd_opt_csv(known_points, self.param_name)                   # TODO: plot

        x_opt, R_opt, _ = self.x_R_min
        return x_opt, R_opt


def evaluate_fd_calculation(rp, sl, fd_parameters, parameter_vals):
    # make temporary Rparams and slab
    test_sl = copy.deepcopy(sl)
    test_rp = copy.deepcopy(rp)

    # make list if only single parameter
    if not isinstance(fd_parameters, list):
        _fd_parameters = [fd_parameters,]
        _parameter_vals = [parameter_vals,]
    else:
        _fd_parameters = fd_parameters
        _parameter_vals = parameter_vals

    # make temporary workdir
    name = _temp_dir_name(_fd_parameters, _parameter_vals)
    temp_dir = _fd_dir(name, rp)

    # apply parameters
    for param, val in zip(_fd_parameters, _parameter_vals):
        param.apply(test_rp, test_sl, val)

    # run calculation
    r, r_per_beam = get_fd_r(test_sl, test_rp,
                           work_dir=temp_dir, home_dir=rp.workdir)

    # TODO: update csv and plot
    return r, r_per_beam, temp_dir



def _fd_dir(name, rp):
    fd_work_dir = rp.workdir / "fd_calc" / name
    return fd_work_dir


def _temp_dir_name(parameters, parameter_vals):
    return "_".join([f"{p.name}={v}" for p, v in zip(parameters, parameter_vals)])

