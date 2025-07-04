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
    from viperleed_jax.optimization.iterator import OptimizerIterator
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
    optimizer_iterator = OptimizerIterator(
        calculator=calculator,
        rpars=rpars,)

    # initial parameter vector
    x = optimizer_iterator.suggested_starting_point
    logger.debug(
        f'Initial parameter vector:\n{x}')

    for optimizer, result in optimizer_iterator:

        logger.info(
            f'Optimizer {optimizer.name} finished with best '
            'R = {result.best_R:.4f}'
        )
        logger.debug(f'Optimizer result:\n{result}')

        # write intermediate results to files
        tmp_slab = copy.deepcopy(slab)
        calculator.apply_to_slab(tmp_slab, rpars, result.best_x)

        poscar.write(slab, f'POSCAR_TL_{optimizer.name}_intermediate',
                     comments='all')
        writeVIBROCC(slab, f'VIBROCC_TL_{optimizer.name}_intermediate')

        # write result to file
        result.write_to_file(f'{optimizer.name}_result.npz')

        # update the last_R parameter
        rpars.last_R = result.best_R


    # Finished optimization
    logger.info(
        f'Finished optimization with best R = {result.best_R:.4f}'
    )
    logger.debug(
        f'Best parameter vector:\n{result.best_x}'
    )

    # apply the final result to the slab object
    logger.debug("Applying final result to slab object.")
    calculator.apply_to_slab(slab, rpars, result.best_x)

    # write the updated POSCAR and VIBROCC files to work
    poscar.write(slab, "POSCAR", comments="all")
    writeVIBROCC(slab, "VIBROCC")
