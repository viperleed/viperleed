"""Module tl_backend of viperleed.files.displacements."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-15'
__license__ = 'GPLv3+'

# Note by @michele-riva: (Part of) this module can probably be moved higher up
# the calc tree so we will be ready to support more backends, not only
# search-related. (See my suggestion of having a more generic BACKEND parameter
# in #333.)

from abc import ABC, abstractmethod


class BackendError(Exception):
    """Base class for errors in backends."""


class IncompatibleBackendError(BackendError):
    """Error when the backend is incompatible with user requests."""


class TensorLEEDBackend(ABC):
    """Base class for the tensor LEED backends.

    These classes have information on the capabilities of a given backend, e.g.
    if it can handle a certain type of search block.
    """

    def __init__(self, name, handle_search_block_func):
        self.name = name
        self._handle_search_block_func = handle_search_block_func

    @abstractmethod
    def replace_search_block(self, search_block):
        """Check if the backend can handle the given search block.

        Parameters
        ----------
        search_block: SearchBlock
        The search block to check.

        Returns
        -------
        tuple(SearchBlock)
            Returns a tuple with a valid set of replacement search
            blocks. If the original search block is valid, return a
            tuple with just the original search block.

        Raises
        ------
        IncompatibleBackendError
            If the backend can't handle the search block and can't
            provide a replacement.
        """
        return self._handle_search_block_func(search_block)

class TensErLEEDBackend(TensorLEEDBackend):
    """TensErLEED backend for ViPErLEED."""

    name = 'TensErLEED'

    @classmethod
    def replace_search_block(cls, search_block):
        """Replace the search block with a TensErLEED-compatible version."""
        # TODO: xy shorthands exits!

        # For TensErLEED backend, this duplicates a search block
        # for example: xy -> xy[1 0] & xy[0 1]
        raise NotImplementedError


class JAXBackend(TensorLEEDBackend):
    """JAX backend for ViPErLEED."""

    name = 'ViPErLEED JAX'

    @classmethod
    def replace_search_block(cls, search_block):
        """Replace the search block with a JAX-compatible version.

        This does not need special handling as ViPErLEED JAX can handle all
        currently implemented search blocks.
        """
        return (search_block,)
