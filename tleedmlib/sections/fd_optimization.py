# -*- coding: utf-8 -*-

"""
Created on Oct 22 2021

@author: Alexander M. Imre
@author: Florian Kraushofer

Tensor LEED Manager section Full-dynamic Optimization
"""
import logging

from viperleed.tleedmlib.classes.rparams._rparams import FD_PARAMETERS, AVAILABLE_MINIMIZERS
from viperleed.tleedmlib.classes.fd_optimizer import *

logger = logging.getLogger("tleedm.fdopt")


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
