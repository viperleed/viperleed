
__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-04-15"
__license__ = "GPLv3+"


import logging


_VLJ_DEPENDENCIES = (
    'viperleed_jax',
    'anytree',
    'jax',
    'spbessax',
    'tqdm',
)

logger = logging.getLogger(__name__)


def check_vlj_dependencies():
    """Check if the required dependencies for viperleed_jax are installed."""
    for dependency in _VLJ_DEPENDENCIES:
        try:
            __import__(dependency)
        except ImportError as e:
            msg = (
                f'Dependency {dependency} is required for use of the '
                'viperleed_jax plugin. If you have installed it, please check '
                'your PYTHONPATH.'
            )
            logger.error(msg)
            # raise ImportError with the original exception
            raise ImportError(f‘Missing dependency: {dependency}.‘) from e
