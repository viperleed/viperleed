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

from viperleed.tleedmlib.classes.rparams._rparams import FD_PARAMETERS, AVAILABLE_MINIMIZERS, apply_scaling
from viperleed.tleedmlib.files import parameters, poscar
from viperleed.tleedmlib.sections.refcalc import refcalc as section_refcalc
from viperleed.tleedmlib.sections.rfactor import rfactor as section_rfactor
from viperleed.tleedmlib.classes.fd_optimizer import *

logger = logging.getLogger("tleedm.fdopt")

_DEFAULT_MINIMIZER_TOL = 1e-4



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

    logger.info("FD optimization of parameter(s): "
                f"{', '.join(rp.FD_PARAMS)}")

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
        fd_param = FDParameter(name=param_name,                                 # TODO: move object creation to point of dict creation
                               start_value=FD_PARAMETERS[param_name]["x0"](rp),
                               bounds=FD_PARAMETERS[param_name]["bounds"],
                               transform_func=FD_PARAMETERS[param_name]["eval"],)
        start_value = fd_param.start_value
        if rp.FD_METHOD.lower() == "parabola":
            (step, max_step, convergence, min_points, max_points
             ) = _update_parabola_settings(rp.FD_PARABOLA, param_name)
            optimizer = SingleParameterParabolaFit(
                fd_param,
                sl,
                rp,
                convergence=convergence,
                step=step,
                min_points=min_points,
                max_points=max_points,
                max_step=max_step,
            )
        elif rp.FD_METHOD.lower() == "error":
            optimizer = SingleParameterBruteForceOptimizer(
                fd_param,
                sl,
                rp,
                min_val=rp.FD_BRUTE_FORCE["min"],
                max_val=rp.FD_BRUTE_FORCE["max"],
                steps=rp.FD_BRUTE_FORCE["steps"],
                error=rp.FD_BRUTE_FORCE["error"],
            )
        elif rp.FD_METHOD.lower() in AVAILABLE_MINIMIZERS:
            optimizer = SingleParameterMinimizer(
                fd_param,
                sl,
                rp,
                minimizer_method=rp.FD_METHOD.lower(),
                tol=rp.FD_MINIMIZER["tol"],
            )
        else:
            raise NotImplementedError(
                "Unknown optimization method for single parameter "
                f"optimization: {rp.FD_METHOD}"
            )
    else:
        # multiple parameters
        fd_parameters = [
            FDParameter(name=param_name,
                        start_value=FD_PARAMETERS[param_name]["x0"](rp),
                        bounds=FD_PARAMETERS[param_name]["bounds"],
                        transform_func=FD_PARAMETERS[param_name]["eval"],)
            for param_name in rp.FD_PARAMS
        ]
        start_value = tuple(param.start_value for param in fd_parameters)
        if rp.FD_METHOD.lower() in AVAILABLE_MINIMIZERS:
            optimizer = ParameterMinimizer(
                fd_parameters,
                sl,
                rp,
                minimizer_method=rp.FD_METHOD.lower(),
                tol=rp.FD_MINIMIZER["tol"],
            )
        else:
            raise NotImplementedError(
                "Unknown optimization method for multiple parameter "
                f"optimization: {rp.FD_METHOD}"
            )


    # perform optimization
    optimizer.optimize(x0=start_value)
    # finalize by plotting and printing results
    opt_x, opt_R = optimizer.finalize()

    result_msg = f"Optimization finished. Best result with R = {opt_R:.4f}:"
    if len(rp.FD_PARAMS) == 1:
        result_msg += f" {rp.FD_PARAMS[0]} = {opt_x:.4f}"
    else:
        for param_name, param_value in zip(rp.FD_PARAMS, list(opt_x)):
            result_msg += f"\n{param_name} = {param_value:.4f}"

    logger.info(result_msg)

    # update PARAMETERS and POSCAR_OUT
    optimizer.write_fd_params_to_file()

    # fetch I(V) from all and plot together
    try:
        optimizer.write_beams_pdf()
    except Exception as exc:
        logger.warning(f"Failed to plot I(V) curves from optimization:\n{exc}")




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

    # clean up tmpdirs
    # for path in tmpdirs:
    #     try:
    #         shutil.rmtree(path)
    #     except Exception as e:
    #         logger.warning("Failed to delete temporary directory {}: "
    #                        .format(path) + str(e))



