# -*- coding: utf-8 -*-

class ParameterError(Exception):
    """Base class for errors raised during PARAMETERS interpretation"""

    def __init__(self, parameter, message="", **__kwargs):
        _message = f"PARAMETERS file: parameter {str(parameter)}:\n"
        if message:
            _message += message + " "  # add space before "input will be ignored"
        super().__init__(_message)


class ParameterNotRecognizedError(ParameterError):
    """Raised when a parameter is not recognized"""

    def __init__(self, parameter, message="Parameter not recognized."):
        super().__init__(parameter, message)


class ParameterUnexpectedInputError(ParameterError):
    """Raised when unexpected input is encountered"""

    def __init__(self, parameter, message=""):
        if not message:
            message = "Encountered unexpected input. Check parameter syntax."
        super().__init__(parameter, message)


# base class for conversion errors
class ParameterConversionError(ParameterError):
    """Raised when a conversion fails"""

    _type = None

    def __init__(self, parameter, given_value=None, message=""):
        if message:
            super().__init__(parameter, message)
            return
        message = 'Failed to convert '
        if given_value:
            message += f'"{given_value}" to {self._type}.'
        else:
            message += f'input to {self._type}. Check parameter syntax.'
        super().__init__(parameter, message)


class ParameterBooleanConversionError(ParameterConversionError):
    """Raised when a boolean conversion fails"""
    _type = 'boolean'


class ParameterIntConversionError(ParameterConversionError):
    """Raised when an int conversion fails"""
    _type = 'integer'


class ParameterFloatConversionError(ParameterConversionError):
    """Raised when a float conversion fails"""
    _type = 'float'


class ParameterValueError(ParameterError):
    """Raised when the value is not allowed"""

    def __init__(self, parameter, given_value=None, message=""):
        if not message:
            message = 'Could not interpret '
            message += f'"{given_value}".' if given_value else 'given_value.'
        super().__init__(parameter, message)


class ParameterParseError(ParameterError):
    """Raised when parsing fails"""

    def __init__(self, parameter, message="",
                 supp_message="Check parameter syntax."):
        if message:
            super().__init__(parameter, message)
            return
        super().__init__(parameter,
                         f"Could not parse input. {supp_message}")



class ParameterNumberOfInputsError(ParameterError):
    """Raised when the number of inputs is unexpected"""

    def __init__(self, parameter, found_and_expected=None, message=""):
        if message:
            super().__init__(parameter, message)
            return
        if found_and_expected:
            super().__init__(parameter,
                             f"Expected {found_and_expected[1]} inputs, "
                             "but found {found_and_expected[0]}.")
        else:
            super().__init__(parameter,
                            "Unexpected number of inputs.")


class ParameterRangeError(ParameterError):
    """Raised when the value is not in the allowed range"""

    def __init__(self, parameter, given_value=None,
                 allowed_range=None, message=None):
        if message:
            super().__init__(parameter, message)
            return
        if given_value is not None and allowed_range is not None:
            message = (f"Value {given_value} is outside allowed range "
                       f"({allowed_range[0]} <= {parameter} <= "
                       f"{allowed_range[1]}).")
        super().__init__(parameter, message)


class ParameterUnknownFlagError(ParameterError):
    """Raised when an unknown flag is encountered"""

    def __init__(self, parameter, message="", flag=""):
        if not message:
            message = f"Unknown flag {flag!r} encountered."
        super().__init__(parameter, message)


class ParameterNeedsFlagError(ParameterError):
    """Raised when a flag is needed but not given"""

    def __init__(self, parameter, message="Parameter requires a flag."):
        super().__init__(parameter, message)
