# -*- coding: utf-8 -*-
"""
Module fd_optimizers of viperleed.tleedmlib.classes.fd_optimizer.

Created on 2023-12-19

@author: Alexander M. Imre (@amimre)

This module contains the classes that perform the actual optimization of
the full-dynamic calculation. The base class is FDOptimizer, which is
abstract and cannot be instantiated. The subclasses implement the
different optimization strategies for one- and multi-dimensional
optimization.
"""
from abc import ABC, abstractmethod
import copy
import csv
from collections.abc import Iterable
import shutil
import os
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
import scipy.optimize
import logging

from .fd_parameter import FDParameter

from viperleed.tleedmlib.classes.fd_optimizer.fd_parameter import apply_scaling
from viperleed.tleedmlib.classes.r_error import get_zero_crossing, get_n_zero_crossings
from viperleed.tleedmlib.files.iorfactor import read_rfactor_columns
from viperleed.tleedmlib.files.ivplot import plot_iv
from viperleed.tleedmlib.files.ioerrorcalc import plot_r_plus_var_r, draw_error
from viperleed.tleedmlib.sections.refcalc import refcalc as section_refcalc
from viperleed.tleedmlib.sections.rfactor import rfactor as section_rfactor
from viperleed.tleedmlib.files import parameters, poscar


# TODO: move to io_fd_optimization ?
try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')                                                       # TODO: check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    plt.style.use('viperleed.tleedm')
    from matplotlib import cm
except Exception:
    _CAN_PLOT = False
else:
    _CAN_PLOT = True

logger = logging.getLogger("tleedm.fdopt")

_DEFAULT_MINIMIZER_TOL = 1e-4

AVAILABLE_MINIMIZERS = (
    'nelder-mead',
    'powell',
    'cg',
    'bfgs',
    'l-bfgs-b',
    'cobyla',
    'TNC'
)

PARABOLA_SYNONYMS = (
    'parabola',
    'parabolic',
)

ERROR_SYNONYMS = (
    'error',
    'errors',
    'brute-force',
)

SCIPY_SYNONYMS = (
    'scipy',
    'scipy.optimize',
    'scipy.optimize.minimize',
)

AVAILABLE_METHODS = (
    *PARABOLA_SYNONYMS,
    *ERROR_SYNONYMS,
    *SCIPY_SYNONYMS
)


#TODO: should probably be moved to a more appropriate location
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
            rfaclist = section_rfactor(sl, rp, 11)                              # TODO: disable plotting
        except Exception:
            logger.error("Error running rfactor calculation")
            raise
    finally:
        os.chdir(home_dir)
    return rp.last_R, rfaclist


class FDOptimizer(ABC):
    def __init__(self, fd_parameters, sl, rp) -> None:
        self.fd_parameters = (
            [fd_parameters, ] if isinstance(fd_parameters, FDParameter)
            else fd_parameters
        )
        self.slab = copy.deepcopy(sl)
        self.rparams = copy.deepcopy(rp)
        self.points = {}

    @property
    def x(self):
        return tuple(self.points.keys())

    @property
    def R(self):
        return tuple(res[0] for res in self.points.values())

    @property
    def R_per_beam(self):
        return tuple(res[1] for res in self.points.values())

    @property
    def theo_spectra(self):
        return tuple(res[2] for res in self.points.values())

    @property
    def exp_spectra(self):
        return tuple(res[3] for res in self.points.values())


    def evaluate(self, x, keep_dirs=False):
        # must be tuple for hashing
        return self._evaluate(tuple(x), keep_dirs=keep_dirs)

    def _evaluate(self, x, keep_dirs=False):
        # check if we have already evaluated the point
        if x in self.points.keys():
            logger.debug("Point already evaluated. Returning cached value.")
            return self.points[x][0]  # return R only

        tmp_slab, tmp_rparams = self._apply_params(x)

        # silence logger up to level error
        logging.disable(max(logger.getEffectiveLevel(), logging.ERROR))
        # run reference calculation
        r, r_per_beam = get_fd_r(tmp_slab,
                                 tmp_rparams,
                                 work_dir=self._eval_dir(x),
                                 home_dir=self.rparams.workdir)
        # reset logger level
        logging.disable(logging.NOTSET)

        # read spectra
        try:
            theo_spec, exp_spec = read_rfactor_columns(self._eval_dir(x))
        except:
            logger.warning("Failed to read spectrum data.")                     # TODO: longer error message
            theo_spec, exp_spec = None, None

        self.points[x] = r, r_per_beam, theo_spec, exp_spec                     # TODO: Don't store exp_spec for all params

        if not keep_dirs:
            shutil.rmtree(self._eval_dir(x))

        return r


    def _apply_params(self, x):
        """Applies a FDParamter to Rparams and Slab

        Parameters
        ----------
        x : float
            Value of the parameter

        Returns
        -------
        Slab
            Slab object with FD parameter value applied
        Rpars
            Rparams object with FD parameter value applied
        """
        tmp_slab = copy.deepcopy(self.slab)
        tmp_rparams = copy.deepcopy(self.rparams)
        _x = np.array(x)
        for param, val in zip(self.fd_parameters, x):
            param.apply(tmp_rparams, tmp_slab, val)
        return tmp_slab, tmp_rparams


    def _eval_dir(self, x,):
        name = "_".join([f"{p.name}={v}" for p, v in
                         zip(self.fd_parameters, x)])
        return self.rparams.workdir / "fd_calc" / name

    @abstractmethod
    def optimize(self, x0):
        """Perform optimization."""
        _x0 = np.array(x0)
        if _x0.shape[0] != len(self.fd_parameters):
            raise ValueError("x0 must have the same length as fd_parameters")
        raise NotImplementedError

    def write_csv(self, file_path=None, delimiter=","):
        _file_path = (file_path if file_path is not None
                      else self.rparams.workdir / "FD_Optimization.csv")

        if len(self.points) == 0:
            logger.warning(f"Writing {_file_path.name} without "
                           "evaluated points.")
        titles = [param.name for param in self.fd_parameters]
        # add "scaling" to name if it's unit cell scaling
        titles = [name if name not in ("a", "b", "c", "ab", "abc")
                else f"{name} scaling"
                for name in titles]
        titles.append("R")

        try:
            # write to file with csv module
            with open(_file_path, 'w') as csv_file:
                writer = csv.writer(csv_file, delimiter=delimiter)
                writer.writerow(titles)
                for x, (r, *_) in self.points.items():
                    param_values = list(x)
                    param_values.append(r)
                    writer.writerow(param_values)
        except Exception as err:
            logger.warning(f"Failed to write {_file_path.name}: {err}")
            raise err# TODO do not raise here

    def write_beams_pdf(self, file_path=None):
        _file_path = (file_path if file_path is not None
                      else self.rparams.workdir / "FD_Optimization_beams.pdf")

        global _CAN_PLOT
        if not _CAN_PLOT:
            logger.debug("Necessary modules for plotting not found. Skipping "
                        "error plotting.")
            return

        # sort evaluated points by R factor
        by_best_r = np.argsort(self.R)[::-1]  # ascending order
        sorted_params = np.array(self.x)[by_best_r]
        exp_spectrum = np.array(self.exp_spectra, dtype=object)[by_best_r][0]

        if exp_spectrum is None:
            logger.warning("Failed to read in experimental spectra from best "
                        "full-dynamic run. Output of collected I(V) spectra "
                        "will be skipped.")
            return

        # generate list of legends for every spectrum
        legends = [
            ", ".join(
                f"{param.name} = {val}" for param, val
                in zip(self.fd_parameters, 
                       values if isinstance(values, Iterable) else (values,))
            )
            for values in sorted_params
        ]
        legends += ["Experiment"]
        # TODO: below could use some refactoring
        best_r_factors = np.array(self.R_per_beam)[by_best_r][0]
        annotations = [f"R = {r:.4f}" for r in best_r_factors]
        label_style = "overbar" if self.rparams.PLOT_IV["overbar"] else "minus"
        label_width = max([beam.getLabel(style=label_style)[1]
                           for beam in self.rparams.expbeams])
        labels = [beam.getLabel(lwidth=label_width, style=label_style)[0]
                for beam in self.rparams.expbeams]
        formatting = copy.deepcopy(self.rparams.PLOT_IV)  # use to set colors
        formatting['colors'] = (
            list(cm.get_cmap('viridis', len(self.points)).colors)
            + [np.array([0, 0, 0, 1])])
        formatting['linewidths'] = [0.5,] * len(self.points) + [1.,]
        formatting['linewidths'][np.argmin(self.R)] = 1.

        try:
            plot_iv(self.theo_spectra + (exp_spectrum,), _file_path,
                    labels=labels, annotations=annotations, legends=legends,
                    formatting=formatting)
        except Exception:
            logger.warning("Error plotting collected I(V) curves.",
                        exc_info=True)


class OneDimensionalFDOptimizer(FDOptimizer):
    """Base class for full-dynamic optimization."""

    def __init__(self, fd_parameters, sl, rp):
        if len(fd_parameters) != 1:
            raise ValueError("OneDimensionalFDOptimizer can only optimize "
                             "one parameter. Multiple parameters given.")
        fd_parameter = fd_parameters[0]
        super().__init__(fd_parameter, sl, rp)
        self.eval_dirs = [] # directories where evaluations are preformed
        self.param_name = fd_parameter.name # name of parameter to optimize

        # bounds for x
        self.bounds = (fd_parameter.bounds if fd_parameter.bounds is not None
                       else (-np.inf, np.inf))

        if not all(self._in_bounds(x) for x in self.x):
            raise FullDynamicOptimizationOutOfBoundsError(
                "Initial points must be within bounds"
            )

    # overwrite x from parent
    @property
    def x(self):
        return tuple((param[0] for param in self.points.keys()))

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
        return super().evaluate(x_val, keep_dirs=True)


    def _in_bounds(self, x_val):
        return self.bounds[0] <= x_val <= self.bounds[1]

    @abstractmethod
    def optimize(self, x0):
        """Perform optimization."""
        raise NotImplementedError

    @abstractmethod
    def finalize(self):
        """Perform final steps after optimization has converged."""
        pass

    @abstractmethod
    def updated_intermediate_output(self):
        pass

    def write_opt_pdf(self, file_path=None):
        global _CAN_PLOT
        if not _CAN_PLOT:
            logger.debug("Necessary modules for plotting not found. Skipping "
                        "error plotting.")
            return

        _file_path = (file_path if file_path is not None
                      else self.rparams.workdir / "FD_Optimization.pdf")

        title = self.fd_parameters[0].name
        if title in ["a", "b", "c", "ab", "abc"]:
            title += " scaling"

        fig = plt.figure(figsize=(5.8, 4.1))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel(title)
        ax.set_ylabel('Pendry R-factor')
        fig.tight_layout()
        self._annotate_opt_pdf(fig, ax)

        # write
        try:
            pdf = PdfPages(_file_path)
            pdf.savefig(fig)
        except PermissionError:
            logger.warning("Failed to write to " + _file_path.name
                        + ": Permission denied.")
        except KeyboardInterrupt:
            raise
        except Exception:
            logger.warning("Failed to write to "+file_path.name, exc_info=True)
        finally:
            try:
                pdf.close()
            except Exception:
                pass
        try:
            plt.close(fig)
        except Exception:
            pass

    @abstractmethod
    def _annotate_opt_pdf(self, fig, ax):
        ax.plot(self.x, self.R, 'o', c='darkslategray')

    def write_fd_params_to_file(self):
        store_fd_param_to_file(
            self.slab, self.rparams, self.param_name, self.x_R_min[0]
        )


class SingleParameterParabolaFit(OneDimensionalFDOptimizer):

    required_settings = {
        'convergence': float,
        'step': float,
    }
    optional_settings = {
        'min_points': int,
        'max_points': int,
        'max_step': float,
    }

    def __init__(self, fd_parameter, sl, rp, convergence, step, min_points, max_points, max_step=None):
        super().__init__(fd_parameter, sl, rp)
        self.convergence = convergence
        self.min_points = min_points
        self.max_points = max_points
        self.step = step
        self.max_step = max_step if max_step is not None else np.inf
        self.curvature_fails = 0


    def optimize(self, x0):
        if x0 is None:
            raise RuntimeError("Must provide at least one evaluation point to "
                               "start.")
        if not self._in_bounds(x0):
            raise ValueError("Initial point must be within bounds")
        self.evaluate((x0,))

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

            logger.info(f"Adding data point at {self.param_name} = {next_x:.4f}")
            self.evaluate((next_x,))  # updates self.x, self.R

            # update plot and csv
            self.updated_intermediate_output()

        # max points reached
        logger.warning(
            "Stopping optimization because the maximum number of data "
            "points has been reached. This may indicate that the optimization "
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
        x_without_min = [x for x in self.x if x != min_x]
        if any([abs(x - min_x) < self.convergence for x in x_without_min]):
            return True

        return False

    def finalize(self):
        # plot for a final time
        known_points = np.array([self.x, self.R]).T
        self.write_csv()
        self.write_opt_pdf()

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
            self.write_csv()
        if len(self.x) >= 3:
            self.write_opt_pdf()

    def _annotate_opt_pdf(self, fig, ax):
        if len(self.x) >= 2:
            self._fit_parabola()
            ax.plot(self.parabola.linspace(), lw=2)
            annotation_pos = ((min(self.x) + max(self.x))/2,
                              max(self.R) - 0.1*(max(self.R) - min(self.R)))
            if self._parabola_curvature > 0:
                ax.annotate(
                    f"Minimum at {self._predicted_min_x:.4f}\n"
                    f"R = {self._predicted_R:.4f}",
                    annotation_pos,
                    fontsize=10,
                    ha="center")
        super()._annotate_opt_pdf(fig, ax)


class SingleParameterBruteForceOptimizer(OneDimensionalFDOptimizer):

    required_settings = {
        'min_val': float,
        'max_val': float,
        'steps': int,
    }
    optional_settings = {
        'error': bool,
    }

    def __init__(self, fd_parameter, sl, rp, min_val, max_val, steps, error=True):
        super().__init__(fd_parameter, sl, rp)
        if None in (min_val, max_val, steps):
            raise ValueError("min_val, max_val, and steps must all be "
                             "provided for brute force optimizer.")
        self.min_val = min_val
        self.max_val = max_val
        self.steps = steps
        self.error = error

        # make sure requested min/max are within bounds
        if not (self._in_bounds(min_val) and self._in_bounds(max_val)):
            raise ValueError("Requested min/max values are outside of "
                             f"parameter bounds for {fd_parameter.name}.")

    def optimize(self, x0):
        if x0 is not None:
            logger.debug("Brute force optimizer does not use an initial "
                         "guess. Ignoring provided x0.")
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
        self.write_csv()
        # write pdf output
        self.write_opt_pdf()
        x_opt, R_opt, _, _ = self.res
        return x_opt, R_opt

    def updated_intermediate_output(self):
        if len(self.x) > 0:
            self.write_csv()
        if len(self.x) >= 3:
            self.write_opt_pdf()

    def _annotate_opt_pdf(self, fig, ax):
        # do not call super in this case
        ax.plot(self.x, self.R, 'o', c='darkslategray', ls = '-')

        # plot R + var(R) line
        var_R = self._calc_var_r()
        x_range = (min(self.x), max(self.x))
        logger.log(5, f"R + var(R): {min(self.R) + var_R:.4f}")
        plot_r_plus_var_r(ax, self.x, self.R, min(self.R), max(self.R),
                          (ax.get_xlim()[0], ax.get_xlim()[1]), var_R)

        # TODO: adapt error class to do this
        # # draw error estimates
        # x_sorted, R_sorted = self.x_R_sorted
        # min_id = np.argmin(self.x_R_sorted[1])
        # R_minus_R_min_var_R = R_sorted - np.min(self.R) - var_R
        # # lower bound
        # if get_n_zero_crossings(R_minus_R_min_var_R[:min_id+1]) == 1:
        #     draw_error(ax,
        #                get_zero_crossing(x_sorted[:min_id+1], R_minus_R_min_var_R[:min_id+1]),
        #                R_sorted[min_id], var_R, x_sorted[min_id],
        #                ax.get_ylim()[1] - ax.get_ylim()[0])
        # # upper bound
        # if get_n_zero_crossings(R_minus_R_min_var_R[min_id:]) == 1:
        #     draw_error(ax,
        #                get_zero_crossing(x_sorted[min_id:], R_minus_R_min_var_R[min_id:]),
        #                R_sorted[min_id], var_R, x_sorted[min_id],
        #                ax.get_ylim()[1] - ax.get_ylim()[0])


    def _calc_var_r(self):
        energy_range = self.rparams.total_energy_range()
        _, temp_rp = self._apply_params((self.x_R_min[0],))

        return (np.sqrt(8*np.abs(temp_rp.V0_IMAG) / energy_range)
                * self.x_R_min[1])


class SingleParameterMinimizer(OneDimensionalFDOptimizer):
    required_settings = {
        'method': str,
    }
    optional_settings = {
        'tol': float,
        'options': dict,
    }

    def __init__(self, fd_parameter, sl, rp, method, tol=None, options=None):
        super().__init__(fd_parameter, sl, rp)

        if method not in AVAILABLE_MINIMIZERS:
            raise ValueError(f"Minimizer method {method} not "
                             f"available. Available methods are: "
                             f"{AVAILABLE_MINIMIZERS}")
        self.minimizer_method = method
        self.tol = tol
        self.minimizer_options = options if options is not None else {}
        if "disp" not in self.minimizer_options:
            self.minimizer_options["disp"] = logger.level <= logging.DEBUG

    def optimize(self, x0):
        logger.info(f"Starting full dynamic calculation using "
                    f"{self.minimizer_method} minimizer method.")
        if len(x0) > 1:
            raise ValueError("Multiple start values x0 given for "
                             "SingleParameterMinimizer.")
        self.res = scipy.optimize.minimize(
            fun=self.evaluate,
            x0=float(x0[0]),
            method=self.minimizer_method,
            bounds=(self.bounds,),
            options=self.minimizer_options,
            tol=self.tol,
            callback=self.updated_intermediate_output,
        )

    def finalize(self):
        if not self.res.success:
            logger.warning("Minimization failed with message: "
                           f"{self.res.message}")
        logger.info(f"Minimizer {self.minimizer_method} finished successfully "
                    f"in {self.res.nit} iterations.")
        self.write_opt_pdf()

        x_opt, R_opt, _ = self.x_R_min
        return x_opt, R_opt

    def _annotate_opt_pdf(self, fig, ax):
        super()._annotate_opt_pdf(fig, ax)

    def evaluate(self, x_val):
        result = super().evaluate(x_val)
        logger.info(
            f"{self.param_name} = {x_val[0]:.4f}: "
            f"R = {result:.4f}")
        return result

    def updated_intermediate_output(self,intermediate_result):
        logger.debug(f"Complete iteration of minimizer {self.minimizer_method}."
                     " Updating CSV output.")
        self.write_csv()
        self.write_opt_pdf()


class ParameterMinimizer(FDOptimizer):

    required_settings = {
        'method': str,
    }
    optional_settings = {
        'tol': float,
        'options': dict,
    }

    def __init__(self, fd_parameters, sl, rp, method, tol=None, options=None):
        super().__init__(fd_parameters, sl, rp)

        if method not in AVAILABLE_MINIMIZERS:
            raise ValueError(f"Minimizer method {method} not "
                             f"available. Available methods are: "
                             f"{AVAILABLE_MINIMIZERS}")
        self.minimizer_method = method
        self.tol = tol  # tolerance for minimizer, passed to scipy
        self.minimizer_options = options  # kwargs for minimizer, passed to scipy
        self.scipy_bounds = scipy.optimize.Bounds(
            lb=[param.bounds[0] for param in self.fd_parameters],
            ub=[param.bounds[1] for param in self.fd_parameters],
        )
        self.minimizer_options = options if options is not None else {}
        if "disp" not in self.minimizer_options:
            self.minimizer_options["disp"] = logger.level <= logging.DEBUG

    def optimize(self, x0):
        logger.info(f"Starting full dynamic calculation using "
                    f"{self.minimizer_method} minimizer method.")
        self.res = scipy.optimize.minimize(
            fun=self.evaluate,
            x0=x0,
            method=self.minimizer_method,
            bounds=self.scipy_bounds,
            options=self.minimizer_options,
            tol=self.tol,
            callback=self.updated_intermediate_output,
        )

    def finalize(self):
        if not self.res.success:
            logger.warning("Minimization failed with message: "
                           f"{self.res.message}")
        logger.info(f"Minimizer {self.minimizer_method} finished successfully "
                    f"in {self.res.nit} iterations.")

        x_opt, R_opt = self.res.x, self.res.fun
        return x_opt, R_opt


    def evaluate(self, x_val):
        result = super().evaluate(x_val)
        logger.info(
            f"{self._eval_point_string(x_val)}: R = {result:.4f}")
        return result

    def _eval_point_string(self, x_val):
        return ", ".join(
            f"{param.name} = {x:.4f}"
            for param, x in zip(self.fd_parameters, x_val)
        )


    def updated_intermediate_output(self,intermediate_result):
        logger.debug(f"Complete iteration of minimizer {self.minimizer_method}."
                     " Updating CSV output.")
        self.write_csv()

    def write_fd_params_to_file(self):
        for par_id, param in enumerate(self.fd_parameters):
            store_fd_param_to_file(
                self.slab, self.rparams, param.name, self.res.x[par_id]
            )


def _update_parabola_settings(fd_parabola_settings, param_name):

    # check whether step is set; if not, choose default value
    if fd_parabola_settings["step"] == 0.:
        if param_name in ["v0i", "theta"]:
            step = 0.5
        elif param_name == "phi":
            step = 2.
        else:   # unit cell size
            step = 0.01
        logger.debug("Initial step size undefined, defaulting to {}"
                     .format(fd_parabola_settings["step"]))
    else:
        step = fd_parabola_settings["step"]

    if fd_parabola_settings["maxstep"] == 0.:
        max_step = 3 * fd_parabola_settings["step"]
    else :
        max_step = fd_parabola_settings["maxstep"]

    if fd_parabola_settings["convergence"] == 0.:
        convergence = 0.1 * fd_parabola_settings["step"]
    else:
        convergence = fd_parabola_settings["convergence"]

    min_points = fd_parabola_settings["minpoints"]
    if min_points <= 1:
        raise ValueError("Minimum number of points for parabola fit must be 1")

    max_points = fd_parabola_settings["maxpoints"]

    max_step = abs(max_step)  # always positive
    convergence = abs(convergence)

    logger.log(10,
               f"Parabola fit settings: step={step}, max_step={max_step}, "
               f"convergence={convergence}, min_points={min_points}, "
               f"max_points={max_points}."
               )

    return step, max_step, convergence, min_points, max_points


def store_fd_param_to_file(sl, rp, which, new_min):
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
        poscar.write(sl, filename=f"POSCAR_OUT_{rp.timestamp}", comments="all")

    # TODO: create parameter to toggle this

    # clean up tmpdirs
    # for path in tmpdirs:
    #     try:
    #         shutil.rmtree(path)
    #     except Exception as e:
    #         logger.warning("Failed to delete temporary directory {}: "
    #                        .format(path) + str(e))
