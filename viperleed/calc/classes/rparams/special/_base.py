"""Module _base of viperleed.calc.classes.rparams.special.

Defines SpecialParameter, the base class for all other parameter
classes defined in this package, and a few exceptions.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-27'
__license__ = 'GPLv3+'


class SpecialParameterError(Exception):
    """Base exception for all special parameter classes."""


class NotASpecialParameterError(SpecialParameterError):
    """Not a known SpecialParameter, or incorrectly subclassed."""

    def __init__(self, param):
        """Initialize exception."""
        super().__init__(f'{param} is not a special parameter, or '
                         'it is not a subclass of SpecialParameter. '
                         f'known_subclasses={SpecialParameter._subclasses}')


class SpecialParameter:
    """Base class for all the PARAMETERS classes defined in this package."""

    _subclasses = {}  # {PARAM_NAME: cls}

    @classmethod
    def __init_subclass__(cls, *, param=None, **kwargs):
        """Register `cls` as being associated with `param`."""
        super().__init_subclass__(**kwargs)
        if param is not None:
            cls._subclasses[param.upper()] = cls

    @classmethod
    def from_value(cls, value):
        """Return an instance of this class from another `value`.

        This method is provided so that subclasses that need to treat
        `value` can reimplement this method. Examples are classes that
        need to tuple-unpack `value`, or that have other convenience
        methods that needs to be called beforehand.
        The default implementation creates an instance of this class
        with value as its only argument.

        Parameters
        ----------
        value : object
            Value to be used to create an instance of this class.

        Returns
        -------
        instance : SpecialParameter
            The instance created.
        """
        return cls(value)

    @classmethod
    def get_subclass(cls, param):
        """Return the subclass for `param`, or complain."""
        param = param.upper()
        try:
            return next(c for p, c in cls._subclasses.items()
                        if p.startswith(param))
        except StopIteration:
            raise NotASpecialParameterError(param) from None
