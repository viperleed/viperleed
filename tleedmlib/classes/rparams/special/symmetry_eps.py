"""Module symmetry_eps of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-12-11

@author: Alexander Imre (@amimre)
@author: Michele Riva (@michele-riva)

Defines the class SymmetryEps, which is a float with optional z value.
"""

from functools import total_ordering
from numbers import Real

from ._base import SpecialParameter


# TODO: arithmetic operations return a float. This may be problematic.
# Dunder methods that may be needed: '__add__', '__eq__', '__floordiv__',
# '__mod__', '__mul__', '__neg__', '__pos__', '__pow__', '__radd__',
# '__rdivmod__', '__rfloordiv__', '__rmod__', '__rmul__', '__round__',
# '__rpow__', '__rsub__', '__rtruediv__', '__sub__', '__truediv__'
# TODO: forbid negative values!

@total_ordering
class SymmetryEps(float, SpecialParameter, param='SYMMETRY_EPS'):
    """SymmetryEps acts like a float but has an optional .z value.

    Used for interpreting the parameter SYMMETRY_EPS. The user
    can specify a second float value to be used for symmetry
    comparisons in the z direction (i.e., orthogonal to the
    surface). If no second value is given, the first value is
    used by default.
    """

    def __new__(cls, value, z=None):
        """Initialize instance."""
        try:
            float(value)
        except (ValueError, TypeError):
            raise TypeError('SymmetryEps value must be float') from None
        if z is not None:
            try:
                z = float(z)
            except (ValueError, TypeError):
                raise TypeError('SymmetryEps z value must be float') from None
        instance = super().__new__(cls, value)
        setattr(instance, '_z', z)
        return instance

    @property
    def has_z(self):
        """Return whether self has a z value defined."""
        # About the disable: the member exists, it's created in __new__
        return self._z is not None  # pylint: disable=no-member

    @property
    def z(self):  # pylint: disable=invalid-name
        """Return the z value of this SymmetryEps."""
        # About the disable: the member exists, it's created in __new__
        if not self.has_z:
            return float(self)
        return self._z  # pylint: disable=no-member

    def __eq__(self, other):
        """Return self == other."""
        if not isinstance(other, (SymmetryEps, Real)):
            return NotImplemented
        was_real = not isinstance(other, SymmetryEps)
        if was_real:
            other = SymmetryEps(other)
        _eq = (float(self), self.z) == (float(other), other.z)
        if was_real:
            # Do not transform a False into a NotImplemented to
            # prevent python from falling back onto trying the
            # reverse operation other.__eq__(self) that returns
            # the incorrect result
            return _eq
        return _eq or NotImplemented

    def __lt__(self, other):
        """Return self < other."""
        if not isinstance(other, (SymmetryEps, Real)):
            return NotImplemented
        was_real = not isinstance(other, SymmetryEps)
        if was_real:
            other = SymmetryEps(other)
        _lt = float(self) < float(other) and self.z < other.z
        if was_real:
            # Do not transform a False into a NotImplemented to
            # prevent python from falling back onto trying the
            # reverse operation other.__gt__(self) that returns
            # the incorrect result
            return _lt
        return _lt or NotImplemented

    def __ne__(self, other):
        """Return self != other."""
        # Notice that, while normally python should take care
        # correctly of the default __ne__ implementation, we
        # want to make sure that we do not fall back onto
        # checking other.__ne__(self) when comparing to a
        # simple Real value. Otherwise we get the wrong result.
        # Reimplementing __ne__ is generally suggested when
        # inheriting from a built-in objects like here (float).
        # See github.com/python/cpython/issues/48645
        _eq = self.__eq__(other)
        if _eq is NotImplemented:
            return _eq
        return not _eq

    def __hash__(self):
        """Return hash(self)."""
        has_no_z = not self.has_z or self.z == float(self)
        if has_no_z:
            return super().__hash__()
        return hash((float(self), self.z))
