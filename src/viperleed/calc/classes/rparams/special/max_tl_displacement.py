"""Module max_tl_displacement of viperleed.calc.classes.rparams.special.

Defines the MaxTLDisplacement class, a float with optional vib value.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-15'
__license__ = 'GPLv3+'

import numpy as np

from dataclasses import dataclass

from ._base import SpecialParameter
from .._defaults import NO_VALUE


@dataclass
class MaxTLDisplacement(SpecialParameter, param='MAX_TL_DISPLACEMENT'):
    """Maximum geometric and vibrational displacements relative to refcalc."""

    geo: float
    _vib: float = NO_VALUE
    action: str = 'refcalc'
    max_duration: float = 1800   # max. refcalc duration for 'continue' in s

    def __post_init__(self):
        """Convert non-float inputs and check their range."""
        for attr_name in ('geo', '_vib'):
            attr = getattr(self, attr_name)
            if attr is NO_VALUE:
                continue
            setattr(self, attr_name,
                    self._check_float_value(attr, extra_msg=f'{attr_name} '))

    @staticmethod
    def _check_float_value(value, extra_msg=''):
        """Return a float version of value. Raise if not acceptable."""
        try:
            float_v = float(value)
        except (ValueError, TypeError):
            raise TypeError(f'MAX_TL_DISPLACEMENT {extra_msg}value '
                            'must be float') from None
        if float_v <= 0:
            raise ValueError(f'MAX_TL_DISPLACEMENT {extra_msg}value '
                             'must be positive')
        return float_v

    @property
    def vib(self):
        """Return the maximum vibrational displacement."""
        return self.geo if self._vib is NO_VALUE else self._vib

    @classmethod
    def from_value(cls, value):
        """Return a MaxTLDisplacement from a 2-item tuple."""
        return cls(*value)

    def assign_float_values(self, values):
        """Assign unlabelled tuple of 1 or 2 values to geo [and _vib]."""
        self.assign_single_value('geo', values[0])
        if len(values) > 1:
            self.assign_single_value('_vib', values[1])

    def assign_single_value(self, flag, value):
        """Assign values to geo or _vib."""
        attr = flag
        if attr == 'vib':
            attr = '_vib'
        setattr(self, attr,
                self._check_float_value(value, extra_msg=f'{flag} '))

    def assign_action(self, values):
        """Assign values to an action. Also interprets the requested time
        if the action is 'refcalc'."""
        self.action = values[0]
        if values[0] == 'refcalc' and len(values) == 1:
            # disable checking for reference calculation time
            self.max_duration = NO_VALUE
        elif values[0] == 'refcalc' and len(values) > 1:
            time_str = values[-1].lower()
            multiplier = 1
            if time_str[-1] in {'h', 'm', 's'}:
                # if none of these letters, assume pure float in seconds
                if time_str[-1] == 'h':
                    multiplier = 3600
                if time_str[-1] == 'm':
                    multiplier = 60
                time_str = time_str[:-1]
            float_v = self._check_float_value(time_str,
                                              extra_msg='(flag "continue") ')
            self.max_duration = float_v*multiplier

    def is_too_far(self, atom):
        """Return whether `atom` was displaced too much."""
        if atom.distance(atom.oriState) > self.geo:
            return True
        if any(np.abs(atom.site.vibamp[el] - atom.site.oriState.vibamp[el])
               > self.vib for el in atom.site.vibamp):
            return True
        return False
