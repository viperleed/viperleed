"""Module mode of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-10'
__license__ = 'GPLv3+'


from viperleed.calc.classes.perturbation_mode import (
    PerturbationMode,
    PerturbationModeError,
)

from .base import DisplacementsFileToken, TokenParserError


class ModeTokenParserError(TokenParserError):
    """Class for errors during the ModeToken parsing."""


class ModeToken(DisplacementsFileToken):
    """Class for the ModeToken."""

    def __init__(self, type_str):
        self.mode_str = type_str
        try:
            self.mode = PerturbationMode(type_str)
        except PerturbationModeError as err:
            msg = f'Unable to parse perturbation type {self.mode_str}.'
            raise ModeTokenParserError(msg) from err

    def __eq__(self, other):
        """Compare self to other ModeToken."""
        if not isinstance(other, ModeToken):
            return NotImplemented
        return self.mode is other.mode

    def __str__(self):
        """Return representation of ModeToken."""
        return f'ModeToken({self.mode})'
