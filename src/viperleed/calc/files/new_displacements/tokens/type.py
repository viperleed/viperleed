"""Module for the <type> token in the DISPLACEMENTS file."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-04-10"
__license__ = "GPLv3+"


from viperleed.calc.classes.perturbation_type import (
    PerturbationMode,
    PerturbationModeError,
)

from .base import DisplacementsFileToken, TokenParserError


class TypeTokenParserError(TokenParserError):
    """Class for errors during the TypeToken parsing."""

class TypeToken(DisplacementsFileToken):
    """Class for the TypeToken."""

    def __init__(self, type_str):
        self.type_str = type_str
        try:
            self.type = PerturbationMode(type_str)
        except PerturbationModeError as err:
            msg = f'Unable to parse perturbation type {self.type_str}.'
            raise TypeTokenParserError(msg) from err

    def __eq__(self, other):
        """Compare self to other TypeToken."""
        if not isinstance(other, TypeToken):
            return False
        return self.type is other.type

    def __str__(self):
        """Return representation of TypeToken."""
        return f'TypeToken({self.type})'
