"""Module tl_backend of viperleed.files.displacements."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2024-10-15"
__license__ = "GPLv3+"

# Note by @michele-riva: (Part of) this module can probably be moved higher up
# the calc tree so we will be ready to support more backends, not only
# search-related. (See my suggestion of having a more generic BACKEND parameter
# in #333.)

from abc import ABC, abstractmethod

class TensorLEEDBackend(ABC):
    """Base class for the tensor LEED backends.

    These classes have information on the capabilities of a given backend, e.g.
    if it can handle a certain type of search block.
    """

    def __init__(self, name, handle_search_block_func):
        self.name = name
        self._handle_search_block_func = handle_search_block_func

    @abstractmethod
    def replace_search_block(self, offsets_block, search_block):
        """Check if the backend can handle the given search block.

        Parameters
        ----------
        offsets_block (OffsetsBlock): The offsets block.
        search_block (SearchBlock): The search block.

        Returns
        -------
        list(SearchBlock): If the backend can't handle the search block, but
            can provide a replacement, return the replacement blocks. Otherwise,
            or if the backend can handle the search block, return None.

        Raises
        ------
        IncompatibleBackendError: If the backend can't handle the search block
            and can't provide a replacement.
        """
        return self._handle_search_block_func(offsets_block, search_block)


def tenserleed_search_block_handler_func(offsets_block, search_block):
    """Handle the search block with the TensorLEED backend."""
    # TODO: xy shorthands exits!
    # For TensErLEED backend, this duplicates a search block xy -> xy[1 0] & xy[0 1]
    raise NotImplementedError


def viplerleed_jax_search_block_handler_func(offsets_block, search_block):
    """Handle the search block with the VIPERLEED backend."""
    # no special handling needed, ViPErLEED jax can handle all search blocks
    return


VIPERLEED_JAX_BACKEND = TensorLEEDBackend(
    'ViPErLEED jax', viplerleed_jax_search_block_handler_func
)

TENSERLEED_BACKENDS = TensorLEEDBackend(
    'TensErLEED', tenserleed_search_block_handler_func
)
