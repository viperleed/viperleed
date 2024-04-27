"""Section Full-dynamic Optimization."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2021-10-22'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path
import shutil

import numpy as np
from numpy.polynomial import Polynomial

from viperleed.calc import psgen
from viperleed.calc.files import iofdopt
from viperleed.calc.files import parameters
from viperleed.calc.files import poscar
from viperleed.calc.sections.refcalc import refcalc as section_refcalc
from viperleed.calc.sections.rfactor import rfactor as section_rfactor
from viperleed.calc.classes.rparams.special.full_dynamic import FD_PARAMETERS
from viperleed.calc.classes.fd_optimizer.fd_optimizers import AVAILABLE_MINIMIZERS
from viperleed.calc.classes.fd_optimizer import *

logger = logging.getLogger(__name__)


def fd_optimization(slab, rp):
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

    if not rp.FD:
        rp.setHaltingLevel(3)
        raise RuntimeError(
            "No parameters defined for optimization in the full-dynamic "
            "calculation. Set FD in the PARAMETERS file to choose a parameter!"
        )

    logger.info("FD optimization of parameter(s): "
                f"{', '.join(param.name for param in rp.FD.parameters)}")

    # make sure there's a compiler ready, and we know the number of cores:      # TODO: this is repeated in multiple locations; refactor into base
    if rp.FORTRAN_COMP[0] == "":
        try:
            rp.getFortranComp()
        except Exception:
            logger.error("No fortran compiler found, cancelling...")
            raise
    rp.updateCores()  # if number of cores is not defined, try to find it

    # Initialize optimizer - it was already set during parameter interpretation
    optimizer = rp.FD.optimizer_class(
        rp.FD.parameters,
        slab,
        rp,
        **rp.FD.settings,
    )

    # Perform optimization
    optimizer.optimize(x0=rp.FD.start_values)

    # Finalize by plotting and printing results
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
