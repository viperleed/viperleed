"""Module symmetry_eps of viperleed.calc.classes.rparams.special.

Defines the SymmetryEps class, a float with optional z value.
"""

__authors__ = (
    'Alexander Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-12-11'
__license__ = 'GPLv3+'

from functools import total_ordering
from numbers import Real

from .base import SpecialParameter


def _make_arithmetic(operation, allow_eps):
    """Return a direct and reverse version of operation.

    Parameters
    ----------
    operation : callable
        Arithmetic operation to be used for the float value and the
        z value of a SymmetryEps object. The methods returned are
        (essentially) operation(eps, other) and operation(other, eps),
        applied to both the main and z values of eps.
    allow_eps : bool
        Whether the second operand may be a SymmetryEps.

    Returns
    -------
    _direct : callable
        Returns operation(eps, other)
    _reverse : callable
        Returns operation(other, eps)
    """
    def _direct(self, other):
        if not allow_eps and isinstance(other, SymmetryEps):
            return NotImplemented
        if not isinstance(other, SymmetryEps):
            try:
                other = SymmetryEps(other)
            except (TypeError, ValueError):
                return NotImplemented
        value = operation(float(self), float(other))
        z_result = (operation(self.z, other.z) if self.has_z or other.has_z
                    else None)
        try:
            return self.__class__(value, z=z_result)
        except ValueError as exc:
            raise ValueError(exc) from None

    def _reverse(self, other):
        return _direct(self, other)
    return _direct, _reverse


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
        cls._check_float_value(value)
        if z is not None:
            z = cls._check_float_value(z, 'z ')
        instance = super().__new__(cls, value)
        setattr(instance, '_z', z)
        return instance

    @staticmethod
    def _check_float_value(value, extra_msg=''):
        """Return a float version of value. Raise if not acceptable."""
        try:
            float_v = float(value)
        except (ValueError, TypeError):
            raise TypeError(f'SymmetryEps {extra_msg}value '
                            'must be float') from None
        if float_v <= 0:
            raise ValueError(f'SymmetryEps {extra_msg}value must be positive')
        return float_v

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

    def __repr__(self):
        """Return a representation string for self."""
        txt = f'SymmetryEps({float(self)}'
        if self.has_z:
            txt += f', z={self.z}'
        return txt + ')'

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

    def __neg__(self):
        """Raise TypeError, as negative SymmetryEps are not valid."""
        raise TypeError('Negative SymmetryEps values are unsupported')

    def __pos__(self):
        """Return a copy of this SymmetryEps."""
        cls = self.__class__
        if self.has_z:
            return cls(float(self), z=self.z)
        return cls(float(self))

    def __round__(self, ndigits=None):
        """Return a rounded version of this SymmetryEps."""
        cls = self.__class__
        value = round(float(self), ndigits)
        if not self.has_z:
            return cls(value)
        return cls(value, round(self.z, ndigits))

    def __rpow__(self, other, modulo=None):
        """Raise TypeError not to allow other**self."""
        raise TypeError('unsupported operand type(s) for pow: '
                        f'{type(other).__name__!r} and '
                        f'{type(self).__name__!r}')

    def __pow__(self, other, modulo=None):
        """Return self**other."""
        cls = self.__class__
        if modulo is not None:
            # pow() supports the modulo argument only
            # when all three arguments are integers
            raise TypeError('3rd argument of pow() unsupported for '
                            f'{cls.__name__} objects')
        value = float(self)**other
        if not self.has_z:
            return cls(value)
        return cls(value, z=self.z**other)

    __add__, __radd__ = _make_arithmetic(float.__add__, allow_eps=True)
    __divmod__, _ = _make_arithmetic(float.__divmod__, allow_eps=False)
    __rdivmod__, _ = _make_arithmetic(float.__rdivmod__, allow_eps=False)
    __floordiv__, _ = _make_arithmetic(float.__floordiv__, allow_eps=False)
    __rfloordiv__, _ = _make_arithmetic(float.__rfloordiv__, allow_eps=False)
    __mod__, _ = _make_arithmetic(float.__mod__, allow_eps=False)
    __rmod__, _ = _make_arithmetic(float.__rmod__, allow_eps=False)
    __mul__, __rmul__ = _make_arithmetic(float.__mul__, allow_eps=False)
    __truediv__, _ = _make_arithmetic(float.__truediv__, allow_eps=False)
    __rtruediv__, _ = _make_arithmetic(float.__rtruediv__, allow_eps=False)
    __sub__, _ = _make_arithmetic(float.__sub__, allow_eps=True)
    __rsub__, _ = _make_arithmetic(float.__rsub__, allow_eps=True)
