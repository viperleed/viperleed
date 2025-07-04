"""Module viperleed.calc.vlj.search.
"""

__authors__ = (
    "Alexander M. Imre (@amimre)",
)
__copyright__ = "Copyright (c) 2019-2024 ViPErLEED developers"
__created__ = "2025-04-17"
__license__ = "GPLv3+"

import copy
import logging

from viperleed.calc.files import poscar
from viperleed.calc.files.vibrocc import writeVIBROCC
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.vlj import VLJ_AVAILABLE

import numpy as np

if VLJ_AVAILABLE:
    from viperleed_jax.from_objects import (
        setup_tl_parameter_space,
        setup_tl_calculator,
    )
    from viperleed_jax import optimization
    from viperleed_jax.utils import benchmark_calculator


logger = logging.getLogger(__name__)


def vlj_search(slab, rpars):
    """Perform a search using the ViPErLEED-JAX backend.

    Parameters
    ----------
    slab : object
        The slab object to search.
    rpars : object
        The parameters for the search.

    Returns
    -------
    object
        The result of the search.
    """

    if not VLJ_AVAILABLE:
        err_msg = (
            "ViPErLEED-JAX is not available. Please install it to use this "
            "function."
        )
        logger.error(err_msg)
        raise RuntimeError(err_msg)

    # select tensors based on the tensor index
    tensor_path = (rpars.paths.home / DEFAULT_TENSORS
                    / f'{DEFAULT_TENSORS}_{rpars.TENSOR_INDEX:03d}.zip')
    if not tensor_path.exists():
        raise FileNotFoundError(f'Tensor {tensor_path} not found.')

    # get the next block from the displacements file – or return if we are done
    if rpars.last_R is None:
        search_block = rpars.vlj_displacements.next(np.inf)
    else:
        search_block = rpars.vlj_displacements.next(rpars.last_R)

    # create the raw parameter space (without applying DISPLACEMENTS yet)
    parameter_space = setup_tl_parameter_space(slab, rpars)

    # apply offsets to the parameter space
    offsets = rpars.vlj_displacements.offsets
    if offsets is not None:
        parameter_space.apply_offsets(offsets)

    # apply the current search block to the parameter space
    # (generates constraints and bounds)
    parameter_space.apply_search_segment(search_block)

    logger.info(
        "\nParameter space created\n"
        "-----------------------\n"
        f"{parameter_space.info}"
    )

    # initialize the calculator
    calculator = setup_tl_calculator(
        slab,
        rpars,
        tensor_path=tensor_path,
        phaseshifts_path=rpars.paths.home / 'PHASESHIFTS',
    )

    # apply the parameter space to the calculator
    calculator.set_parameter_space(parameter_space)
    logger.info('Calculator initialized with parameter space.')

    benchmark_results = benchmark_calculator(calculator, use_grad=False)
    logger.info(
        "\nBenchmark results\n"
        "-----------------\n"
        f"{benchmark_results}"
)

    cmaes_optimizer = optimization.CMAESOptimizer(
        fun=calculator.R,
            n_generations=200,
            pop_size=30,
            ftol=1e-3,
    )

    x0 = np.array([0.50]*calculator.n_free_parameters)
    cmaes_result = cmaes_optimizer(x0)

    logger.info(cmaes_result)

    # write to a copy of the slab object
    tmp_slab = copy.deepcopy(slab)
    calculator.apply_to_slab(tmp_slab, rpars, cmaes_result.best_x)
    # write intermediate results to files
    poscar.write(slab, "POSCAR_TL_intermediate", comments="all")
    writeVIBROCC(slab, "VIBROCC_TL_intermediate")

    # write result to file
    cmaes_result.write_to_file('CMAES_result.npz')

    slsqp_opt = optimization.SLSQPOptimizer(
        fun=calculator.R,
        grad=calculator.grad_R,
    )
    slsqp_result = slsqp_opt(x0=cmaes_result.best_x)


    logger.info(slsqp_result)

    # set last R

    rpars.last_R = slsqp_result.best_R

    # write the result to the slab object
    calculator.apply_to_slab(slab, rpars, slsqp_result.best_x)


    # write the updated POSCAR and VIBROCC files to work
    poscar.write(slab, 'POSCAR', comments='all')
    writeVIBROCC(slab, 'VIBROCC')

    # write result to file
    slsqp_result.write_to_file('SLSQP_result.npz')
