"""Module errors of viperleed.calc.files.parameters.

This module used to be called files/parameter_errors.py. Refactored
in October 2023.

Defines exceptions that may be raised when reading/writing/interpreting
a PARAMETERS file.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-06-07'
__license__ = 'GPLv3+'


class ParameterError(Exception):
    """Base class for errors raised during PARAMETERS interpretation."""

    _default_message = ''

    def __init__(self, parameter, message='', **__kwargs):
        """Initialize instance."""
        _message = f'PARAMETERS file: parameter {str(parameter)}:\n'
        _message += message or self._default_message
        self.parameter = parameter
        self.message = message or self._default_message
        super().__init__(_message)


class InconsistentParameterError(ParameterError):
    """A user parameter conflicts with the one derived by us."""


class MissingEqualsError(ParameterError):
    """A known parameter is present on a line without an '='."""


class ParameterConflictError(ParameterError):
    """Raised when the value of two or more parameters clash."""

    _default_message = 'Parameters have conflicting values'

    def __init__(self, *parameters, message=''):
        """Initialize instance."""
        super().__init__(', and '.join(parameters), message=message)


# base class for conversion errors
class ParameterConversionError(ParameterError):
    """Raised when a conversion fails."""

    _type = None

    def __init__(self, parameter, given_value=None, message=''):
        """Initialize instance."""
        if message:
            super().__init__(parameter, message)
            return
        message = 'Failed to convert '
        if given_value is not None:
            message += f'{given_value!r} to {self._type}'
        else:
            message += f'input to {self._type}. Check parameter syntax'
        super().__init__(parameter, message)


class ParameterBooleanConversionError(ParameterConversionError):
    """Raised when a boolean conversion fails."""

    _type = 'boolean'


class ParameterIntConversionError(ParameterConversionError):
    """Raised when an int conversion fails."""

    _type = 'integer'


class ParameterFloatConversionError(ParameterConversionError):
    """Raised when a float conversion fails."""

    _type = 'float'


class ParameterHasNoValueError(ParameterError):
    """Raised when a parameter is not recognized."""

    _default_message = 'Parameter appears to have no value'


class ParameterNeedsFlagError(ParameterError):
    """Raised when a flag is needed but not given."""

    _default_message = 'Parameter requires a flag'


class ParameterNeedsSlabError(ParameterError):
    """A parameter requiring a Slab was requested without a Slab present."""

    _default_message = 'Cannot interpret parameter without a slab'


class ParameterNotRecognizedError(ParameterError):
    """Raised when a parameter is not recognized."""

    _default_message = 'Parameter not recognized'


class ParameterNumberOfInputsError(ParameterError):
    """Raised when the number of inputs is unexpected."""

    def __init__(self, parameter, found_and_expected=None, message=''):
        """Initialize instance."""
        if message:
            pass
        elif found_and_expected:
            message = (f'Expected {found_and_expected[1]} inputs, '
                       f'but found {found_and_expected[0]}')
        else:
            message = 'Unexpected number of inputs'
        super().__init__(parameter, message)


class ParameterParseError(ParameterError):
    """Raised when parsing fails."""

    def __init__(self, parameter, message='',
                 supp_message='Check parameter syntax'):
        """Initialize instance."""
        if not message:
            message = f'Could not parse input. {supp_message}'
        super().__init__(parameter, message)


class ParameterRangeError(ParameterError):
    """Raised when the value is not in the allowed range."""

    def __init__(self, parameter, given_value=None,
                 allowed_range=None, message=None):
        """Initialize instance."""
        if message:
            pass
        elif given_value is not None and allowed_range is not None:
            message = (f'Value {given_value} is outside allowed range '
                       f'({allowed_range[0]} <= {parameter} <= '
                       f'{allowed_range[1]})')
        super().__init__(parameter, message)


class ParameterUnexpectedInputError(ParameterError):
    """Raised when unexpected input is encountered."""

    _default_message = 'Encountered unexpected input. Check parameter syntax'


class ParameterUnknownFlagError(ParameterError):
    """Raised when an unknown flag is encountered."""

    def __init__(self, parameter, message='', flag=''):
        """Initialize instance."""
        if not message:
            message = f'Unknown flag {flag!r} encountered'
        super().__init__(parameter, message)


class ParameterValueError(ParameterError):
    """Raised when the value is not allowed."""

    def __init__(self, parameter, given_value=None, message=''):
        """Initialize instance."""
        if not message:
            message = 'Could not interpret '
            message += f'{given_value!r}' if given_value else 'given value'
        super().__init__(parameter, message)


class SuperfluousParameterError(ParameterError):
    """A useless parameter was given in the current PARAMETERS file."""

    _default_message = 'Parameter cannot be used in this PARAMETERS file'
