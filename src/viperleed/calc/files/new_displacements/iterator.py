"""Module iterator of viperleed.calc.files.displacements."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-06-18"
__license__ = "GPLv3+"


from collections.abc import Iterator

class DisplacementsIterator(Iterator):
    """Class to iterate over search blocks from the DISPLACEMENTS file."""

    def __init__(self, rpars, displacements_file):
        """Initialize the iterator with a reader."""
        self._displacements_file = displacements_file

    def __iter__(self):
        """Return the iterator itself."""
        return self

    def __next__(self):
        """Return the next search block."""
        #TODO
