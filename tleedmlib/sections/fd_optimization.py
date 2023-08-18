# -*- coding: utf-8 -*-

"""
Created on Oct 22 2021

@author: Alexander M. Imre
@author: Florian Kraushofer

Tensor LEED Manager section Full-dynamic Optimization
"""
from abc import ABC, abstractmethod
from collections.abc import Iterable
from functools import lru_cache
from pathlib import Path
import copy
import csv
import logging
import os
import shutil

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    plt.style.use('viperleed.tleedm')
    from matplotlib import cm
except Exception:
    _CAN_PLOT = False
else:
    _CAN_PLOT = True

import numpy as np
from numpy.polynomial import Polynomial
import scipy.optimize

from viperleed.tleedmlib import psgen
from viperleed.tleedmlib.classes.r_error import get_zero_crossing, get_n_zero_crossings
from viperleed.tleedmlib.files import iofdopt as tl_io
from viperleed.tleedmlib.files.iorfactor import read_rfactor_columns
from viperleed.tleedmlib.files.ivplot import plot_iv
from viperleed.tleedmlib.files.ioerrorcalc import plot_r_plus_var_r, draw_error
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
            rfaclist = section_rfactor(sl, rp, 11)                              # TODO: disable plotting
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

    if not rp.FD_PARAMS:
        rp.setHaltingLevel(3)
        raise RuntimeError(
            "No parameters defined for optimization in the full-dynamic "
            "calculation. Set FD in the PARAMETERS file to choose a parameter!"
        )

    logger.info(f"FD optimization of parameters: {rp.FD_PARAMS}")

    # make sure there's a compiler ready, and we know the number of cores:      # TODO: this is repeated in multiple locations; refactor into base
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise
    rp.updateCores()  # if number of cores is not defined, try to find it

    if len(rp.FD_PARAMS) == 1:
        # single parameter optimization
        param_name = rp.FD_PARAMS[0]
        if param_name == "parabola":
            _update_parabola_settings(rp.FD_PARABOLA, param_name)
            # TODO: continue here; there is some confusion about min/max points and steps yet
            optimizer = None
        elif param_name == "error":
            pass #TODO
    else:
        # multiple parameters
        pass #TODO
        optimizer = None


    # perform optimization
    optimizer.optimize(x0=fd_param.start_value)
    # finalize by plotting and printing results
    opt_x, opt_R = optimizer.finalize()


    logger.info(
        f"Optimization finished. Best {which} = {opt_x:.4f} with "
        f"R = {opt_R:.4f}")

    # update PARAMETERS and POSCAR_OUT
    store_fd_param_to_file(sl, rp, fd_param.name, opt_x)

    # fetch I(V) from all and plot together
    known_points = np.array([optimizer.x, optimizer.R]).T
    try:
        optimizer.write_beams_pdf()
    except Exception as exc:
        logger.warning(f"Failed to plot I(V) curves from optimization: {exc}")
        raise exc # TODO: remove this raise

    # get FDParameter object
    fd_param = FDParameter(
        which,
        FD_PARAMETERS[which]["x0"](rp),
        FD_PARAMETERS[which]["bounds"],
        FD_PARAMETERS[which]["eval"],
    )

    optimizer = SingleParameterParabolaFit(
        fd_parameter=fd_param,
        sl=sl,
        rp=rp,
        convergence=convergence,
        step=step,
        max_points=max_points,
        min_points=min_points,
    )
    optimizer = SingleParameterBruteForceOptimizer(
        fd_parameter=fd_param,
        sl=sl,
        rp=rp,
        min_val=1,
        max_val=7,
        steps=3,
    )


def _update_parabola_settings(fd_parabola_settings, param_name):

    # check whether step is set; if not, choose default value
    if fd_parabola_settings.OPTIMIZE["step"] == 0.:
        if param_name in ["v0i", "theta"]:
            step = 0.5
        elif param_name == "phi":
            step = 2.
        else:   # unit cell size
            step = 0.01
        logger.debug("Initial step size undefined, defaulting to {}"
                     .format(fd_parabola_settings.OPTIMIZE["step"]))
    else:
        step = fd_parabola_settings.OPTIMIZE["step"]

    if fd_parabola_settings.OPTIMIZE["maxstep"] == 0.:
        max_step = 3 * fd_parabola_settings.OPTIMIZE["step"]
    else :
        max_step = fd_parabola_settings.OPTIMIZE["maxstep"]

    if fd_parabola_settings.OPTIMIZE["convergence"] == 0.:
        convergence = 0.1 * fd_parabola_settings.OPTIMIZE["step"]
    else:
        convergence = fd_parabola_settings.OPTIMIZE["convergence"]

    min_points = fd_parabola_settings.OPTIMIZE["min_points"]
    if min_points <= 1:
        raise ValueError("Minimum number of points for parabola fit must be 1")

    max_step = abs(max_step)  # always positive
    convergence = abs(convergence)

    return step, max_step, convergence, min_points




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
            return self.points[x]

        tmp_slab, tmp_rparams = self._apply_params(x)

        r, r_per_beam = get_fd_r(tmp_slab,
                                 tmp_rparams,
                                 work_dir=self._eval_dir(x),
                                 home_dir=self.rparams.workdir)

        # read spectra
        try:
            theo_spec, exp_spec = read_rfactor_columns(self._eval_dir(x))
        except:
            logger.warning("Failed to read spectrum data.")  # TODO: longer error message
            theo_spec, exp_spec = None, None

        self.points[x] = r, r_per_beam, theo_spec, exp_spec                     # TODO: why do we need to store expspec for all params?

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
        by_best_r = np.argsort(self.R)
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

    def __init__(self, fd_parameter, sl, rp):
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


class SingleParameterParabolaFit(OneDimensionalFDOptimizer):
    def __init__(self, fd_parameter, sl, rp, convergence, step, min_points, max_points, max_step=None):
        super().__init__(fd_parameter, sl, rp)
        self.convergence = convergence
        self.min_points = min_points
        self.max_points = max_points
        self.step = step
        self.max_step = max_step if max_step is not None else np.inf
        self.curvature_fails = 0


    def optimize(self, x0):
        if len(x0) == 0:
            raise RuntimeError("Must provide at least one evaluation point to "
                               "start.")
        for val in x0:
            if not self._in_bounds(val):
                raise ValueError("Initial point must be within bounds")
            self.evaluate(val)

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
        if len(self.x) >= 2:
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

    def __init__(self, fd_parameter, sl, rp, min_val, max_val, steps):
        super().__init__(fd_parameter, sl, rp)
        self.min_val = min_val
        self.max_val = max_val
        self.steps = steps

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
        self.write_opt_pdf()
        x_opt, R_opt, _, _ = self.res
        return x_opt, R_opt

    def updated_intermediate_output(self):
        pass

    def _annotate_opt_pdf(self, fig, ax):
        # do not call super in this case
        ax.plot(self.x, self.R, 'o', c='darkslategray', ls = '-')

        # plot R + var(R) line
        var_R = self._calc_var_r()
        x_range = (min(self.x), max(self.x))
        logger.log(5, f"R + var(R): {min(self.R) + var_R:.4f}")
        plot_r_plus_var_r(ax, self.x, self.R, min(self.R), max(self.R),
                          (ax.get_xlim()[0], ax.get_xlim()[1]), var_R)

        # draw error estimates
        x_sorted, R_sorted = self.x_R_sorted
        min_id = np.argmin(self.x_R_sorted[1])
        R_minus_R_min_var_R = R_sorted - np.min(self.R) - var_R
        # lower bound
        if get_n_zero_crossings(R_minus_R_min_var_R[:min_id+1]) == 1:
            draw_error(ax,
                       get_zero_crossing(x_sorted[:min_id+1], R_minus_R_min_var_R[:min_id+1]),
                       R_sorted[min_id], var_R, x_sorted[min_id],
                       ax.get_ylim()[1] - ax.get_ylim()[0])
        # upper bound
        if get_n_zero_crossings(R_minus_R_min_var_R[min_id:]) == 1:
            draw_error(ax,
                       get_zero_crossing(x_sorted[min_id:], R_minus_R_min_var_R[min_id:]),
                       R_sorted[min_id], var_R, x_sorted[min_id],
                       ax.get_ylim()[1] - ax.get_ylim()[0])



    def _calc_var_r(self):
        energy_range = self.rparams.total_energy_range()
        _, temp_rp = self._apply_params((self.x_R_min[0],))

        return (np.sqrt(8*np.abs(temp_rp.V0_IMAG) / energy_range)
                * self.x_R_min[1])


class SingleParameterMinimizer(OneDimensionalFDOptimizer):                      # TODO: can implement a callback method in the minimizer
    def __init__(self, fd_parameter, sl, rp, minimizer_method, tol=None):
        super().__init__(fd_parameter, sl, rp)

        if minimizer_method not in AVAILABLE_MINIMIZERS:
            raise ValueError(f"Minimizer method {minimizer_method} not "
                             f"available. Available methods are: "
                             f"{AVAILABLE_MINIMIZERS}")
        self.minimizer_method = minimizer_method
        self.tol = tol

    def optimize(self, x0):
        self.res = scipy.optimize.minimize(
            fun=self.evaluate,
            x0=float(x0),
            method=self.minimizer_method,
            bounds=(self.bounds,),
            options={"disp": logger.level <= logging.DEBUG},
            tol=self.tol
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
