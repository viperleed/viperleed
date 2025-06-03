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

from dataclasses import dataclass
from enum import Enum
import logging

import numpy as np

from ..defaults import NO_VALUE
from .base import SpecialParameter

_LOGGER = logging.getLogger(__name__)
_MAP_TIME_MULTIPLIER = {
    # Assume pure float in seconds if none of these
    'h': 3600,
    'm': 60,
    's': 1,  # Explicitly include this so char is dropped if given
    }


class MaxTLAction(Enum):
    """What to do if we exceed MAX_TL_DISPLACEMENT."""

    IGNORE = 'ignore'
    STOP = 'stop'
    REFCALC = 'refcalc'

    @property
    def has_options(self):
        return self is MaxTLAction.REFCALC


@dataclass
class MaxTLDisplacement(SpecialParameter, param='MAX_TL_DISPLACEMENT'):
    """Maximum geometric and vibrational displacements relative to refcalc."""

    geo: float
    _vib: float = NO_VALUE
    action: MaxTLAction = MaxTLAction.REFCALC
    max_duration: float = 1800   # max. refcalc duration for 'refcalc' in s

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
        """Return a float version of `value`. Raise if not acceptable."""
        try:
            float_v = float(value)
        except (ValueError, TypeError):
            raise TypeError(f'MAX_TL_DISPLACEMENT {extra_msg}value '
                            'must be float.') from None
        if float_v <= 0:
            raise ValueError(f'MAX_TL_DISPLACEMENT {extra_msg}value '
                             'must be positive.')
        return float_v

    @property
    def vib(self):
        """Return the maximum vibrational displacement."""
        return self.geo if self._vib is NO_VALUE else self._vib

    @classmethod
    def from_value(cls, value):
        """Return a MaxTLDisplacement from a 2-item tuple."""
        return cls(*value)

    def assign_action(self, action, *options):
        """Assign an `action` with additional `options`.

        Parameters
        ----------
        action : MaxTLAction or str
            The action to be assigned. Should be one of the known actions.
        *options : str, optional
            Additional information for `action`. Only used for action 'refcalc'
            to set the `max_duration` attribute from the last item in
            `options`.

        Returns
        -------
        None.
        """
        try:
            self.action = MaxTLAction(action)
        except ValueError:
            raise ValueError(f'Unkonwn action {action!r}.') from None
        if self.action is MaxTLAction.REFCALC and not options:
            # disable checking for reference calculation time
            self.max_duration = NO_VALUE
            return
        if self.action is MaxTLAction.REFCALC:
            self._assign_max_duration(options[-1].lower().strip())

    def assign_float_values(self, values):
        """Assign unlabelled tuple of 1 or 2 values to geo [and _vib]."""
        self.assign_single_value('geo', values[0])
        if len(values) > 1:
            self.assign_single_value('_vib', values[1])

    def assign_single_value(self, flag, value):
        """Assign values to geo or _vib."""
        attr = flag
        # pylint: disable-next=magic-value-comparison
        if attr == 'vib':
            attr = '_vib'
        setattr(self, attr,
                self._check_float_value(value, extra_msg=f'{flag} '))

    def _assign_max_duration(self, time_str):
        """Set a maximum refcalc duration from a `time_str` interval."""
        try:
            multiplier = _MAP_TIME_MULTIPLIER[time_str[-1]]
            time_str = time_str[:-1]
        except KeyError:
            multiplier = 1
        float_v = self._check_float_value(time_str,
                                          extra_msg='(flag "refcalc") ')
        self.max_duration = float_v*multiplier

    def is_too_far(self, atom):
        """Return whether `atom` was displaced too much."""
        dist = atom.distance(atom.oriState)
        if dist > self.geo:
            _LOGGER.debug(
                f'MAX_TL_DISPLACEMENT: geometry: {atom} displaced by '
                f'{dist:.3f} A (> {self.geo} A)'
                )
            return True
        vib_diff = {el: np.abs(atom.site.vibamp[el]
                               - atom.site.oriState.vibamp[el])
                    for el in atom.site.vibamp}
        if any(v > self.vib for v in vib_diff.values()):
            max_dist_el = max(vib_diff, key=vib_diff.get)
            _LOGGER.debug(
                f'MAX_TL_DISPLACEMENT: vibration: {atom} displaced by '
                f'{vib_diff[max_dist_el]} A for element {max_dist_el} '
                f'(> {self.vib} A)'
                )
            return True
        return False
