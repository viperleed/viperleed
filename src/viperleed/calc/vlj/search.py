"""Module viperleed.calc.vlj.search.
"""

__authors__ = (
    "Alexander M. Imre (@amimre)",
)
__copyright__ = "Copyright (c) 2019-2024 ViPErLEED developers"
__created__ = "2025-04-17"
__license__ = "GPLv3+"


import logging
from pathlib import Path


from viperleed.calc.constants import DEFAULT_TENSORS

try:
    from viperleed_jax.from_objects import calculator_from_objects
    from viperleed_jax import optimization
except ModuleNotFoundError:
    VLJ_AVAILABLE = False
else:
    VLJ_AVAILABLE = True

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
        raise RuntimeError(f'Tensor {tensor_path} not found.')

    # path to DISPLACEMENTS file
    displacements_file = Path() / 'DISPLACEMENTS'

    calculator = calculator_from_objects(slab, rpars,
                                        tensor_path=tensor_path,
                                        displacements_path=displacements_file,)
    logger.info(
        "\nParameter space created\n"
        "-----------------------\n"
        f"{calculator.parameter_space.info}"
    )
