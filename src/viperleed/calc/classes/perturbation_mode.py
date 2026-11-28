"""Module perturbation_mode of viperleed.calc.classes."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-05-13'
__license__ = 'GPLv3+'

from enum import Enum


class PerturbationModeError(ValueError):
    """Error class for unknown perturbation mode errors."""


class PerturbationMode(str, Enum):
    """Enum for perturbation modes."""

    GEO = 'geo'
    VIB = 'vib'
    OCC = 'occ'
    DOM = 'dom'

    def __str__(self):
        """Return the string representation of the enum."""
        return self.value

    @classmethod
    def from_string(cls, mode_string):
        """Convert a string to a PerturbationMode enum.

        Parameters
        ----------
        mode_string : str
            The string to convert. Case-insensitive.

        Returns
        -------
        PerturbationMode
            The mode whose value is `mode_string`.

        Raises
        ------
        PerturbationModeError
            If `mode_string` names no known mode.
        """
        try:
            return cls(mode_string.lower())
        except ValueError:
            msg = f'Unknown perturbation mode: {mode_string!r}. '
            raise PerturbationModeError(msg) from None
