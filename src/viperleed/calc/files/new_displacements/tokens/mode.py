"""Module for the <mode> token in the DISPLACEMENTS file."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-04-10"
__license__ = "GPLv3+"


from viperleed.calc.classes.perturbation_mode import (
    PerturbationMode,
    PerturbationModeError,
)

from .base import DisplacementsFileToken, TokenParserError


class ModeTokenParserError(TokenParserError):
    """Class for errors during the ModeToken parsing."""

class ModeToken(DisplacementsFileToken):
    """Class for the ModeToken."""

    def __init__(self, mode_str):
        self.mode_str = mode_str
        try:
            self.mode = PerturbationMode(mode_str)
        except PerturbationModeError as err:
            msg = f'Unable to parse perturbation mode {self.mode_str}.'
            raise ModeTokenParserError(msg) from err

    def __eq__(self, other):
        """Compare self to other ModeToken."""
        if not isinstance(other, ModeToken):
            return False
        return self.mode is other.mode

    def __str__(self):
        """Return representation of ModeToken."""
        return f'ModeToken({self.mode})'
