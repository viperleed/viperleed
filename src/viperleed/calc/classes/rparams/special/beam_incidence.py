"""Module beam_incidence of viperleed.calc.classes.rparams.special.

Defines the BeamIncidence class, an immutable container for the
polar (THETA) and azimuthal (PHI) angles of the incident beam, as
given by the BEAM_INCIDENCE parameter.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from numbers import Real

from .base import SpecialParameter


FULL_TURN = 360.0
MAX_THETA = 90.0
_NORMAL_INCIDENCE_EPS = 1e-4  # degrees


class BeamIncidence(SpecialParameter, param='BEAM_INCIDENCE'):
    """The direction of the incident beam, in degrees.

    Instances are immutable. Use the constructor to obtain a
    new incidence (e.g., with a different .theta or .phi).

    Attributes
    ----------
    theta : int or float
        Polar angle of the incident beam with respect to the
        surface normal. Always between 0 and 90 degrees.
    phi : int or float
        Azimuthal angle of the incident beam, measured from
        the x axis. Between 0 and 360 degrees (both included).
    is_normal : bool
        Whether the beam is (numerically) at normal incidence.
    """

    __slots__ = ('_phi', '_theta')

    def __init__(self, theta=0.0, phi=0.0):
        """Initialize instance from angles in degrees.

        Parameters
        ----------
        theta : float, optional
            Polar angle with respect to the surface normal. Must
            be between -90 and 90 degrees. A negative value is
            equivalent to its absolute value with `phi` rotated
            by 180 degrees. Default is 0.
        phi : float, optional
            Azimuthal angle, measured from the x axis. Values
            outside [0, 360] are brought into this range modulo
            360 degrees. Default is 0.

        Raises
        ------
        TypeError
            If `theta` or `phi` are not real numbers.
        ValueError
            If `theta` is not between -90 and 90 degrees.
        """
        theta = self._to_real(theta, 'theta')
        phi = self._to_real(phi, 'phi')
        if abs(theta) > MAX_THETA:
            raise ValueError(f'Invalid theta={theta} for BeamIncidence. '
                             f'Must be between {-MAX_THETA} and '
                             f'{MAX_THETA} degrees')
        if not 0 <= phi <= FULL_TURN:
            phi %= FULL_TURN
        if theta < 0:
            theta = abs(theta)
            phi = (phi + FULL_TURN / 2) % FULL_TURN
        object.__setattr__(self, '_theta', theta)
        object.__setattr__(self, '_phi', phi)

    @classmethod
    def from_value(cls, value):
        """Return a BeamIncidence from a (theta, phi) pair."""
        if isinstance(value, cls):
            return value
        return cls(*value)

    @property
    def is_normal(self):
        """Return whether this is (numerically) normal incidence."""
        return self.theta <= _NORMAL_INCIDENCE_EPS

    @property
    def phi(self):
        """Return the azimuthal angle (degrees from x)."""
        return self._phi  # pylint: disable=no-member  # Slot

    @property
    def theta(self):
        """Return the polar angle (degrees from the surface normal)."""
        return self._theta  # pylint: disable=no-member  # Slot

    def __delattr__(self, name):
        """Raise AttributeError, as BeamIncidence is immutable."""
        raise AttributeError(f'{type(self).__name__} is immutable')

    def __eq__(self, other):
        """Return self == other."""
        if not isinstance(other, BeamIncidence):
            return NotImplemented
        return (self.theta, self.phi) == (other.theta, other.phi)

    def __hash__(self):
        """Return hash(self)."""
        return hash((self.theta, self.phi))

    def __repr__(self):
        """Return a representation string for self."""
        return f'BeamIncidence(theta={self.theta}, phi={self.phi})'

    def __reduce__(self):
        """Return information for copying and pickling self."""
        return type(self), (self.theta, self.phi)

    def __setattr__(self, name, value):
        """Raise AttributeError, as BeamIncidence is immutable."""
        raise AttributeError(f'{type(self).__name__} is immutable')

    @staticmethod
    def _to_real(value, name):
        """Return value as a real number. Raise if not possible.

        Real numbers are returned unchanged, so that an incidence
        given as integers (e.g., the default) keeps being written
        as integers (e.g., in experiment_symmetry.ini).
        """
        if isinstance(value, Real):
            return value
        if isinstance(value, (str, bytes)):
            raise TypeError(f'BeamIncidence {name} must be a real number')
        try:
            return float(value)
        except (TypeError, ValueError):
            raise TypeError(f'BeamIncidence {name} must be '
                            'a real number') from None
