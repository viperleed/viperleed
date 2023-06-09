# -*- coding: utf-8 -*-

class ParameterError(Exception):
    '''Base class for errors raised during PARAMETERS interpretation'''

    def __init__(self, parameter, message):
        message = f"PARAMETERS file: parameter {str(parameter)}:\n"
        if message:
            message += message + " "  # add space before "input will be ignored"
        super().__init__(message)

class ParameterNotRecognizedError(ParameterError):
    '''Raised when a parameter is not recognized'''

    def __init__(self, parameter):
        super().__init__(parameter, "Parameter not recognized.")

class ParameterUnexpectedInputError(ParameterError):
    '''Raised when unexpected input is encountered'''

    def __init__(self, parameter):
        super().__init__(parameter,
                         "Encountered unexpected input."
                         "Check parameter syntax.")


# base class for conversion errors
class ParameterConversionError(ParameterError):
    '''Raised when a conversion fails'''
    _type = None

    def __init__(self, parameter, given_value=None):
        msg = 'Failed to convert '
        if given_value:
            msg += f'"{given_value}" to {self._type}.'
        else:
            msg = f'input to {self._type}. Check parameter syntax.'
        super().__init__(parameter, msg)


class ParameterBooleanConversionError(ParameterConversionError):
    '''Raised when a boolean conversion fails'''
    _type = 'boolean'


class ParameterIntConversionError(ParameterConversionError):
    '''Raised when an int conversion fails'''
    _type = 'integer'


class ParameterFloatConversionError(ParameterConversionError):
    '''Raised when a float conversion fails'''
    _type = 'float'


class ParameterValueError(ParameterError):
    '''Raised when the value is not allowed'''

    def __init__(self, parameter, given_value=None):
        if given_value:
            super().__init__(parameter,
                             f'Could not interpret "{given_value}".')
        else:
            super().__init__(parameter,
                             "Could not interpret given value.")

class ParameterParseError(ParameterError):
    '''Raised when parsing fails'''

    def __init__(self, parameter, supp_message=None):
        if supp_message:
            super().__init__(parameter,
                             f"Could not parse input. {supp_message}")
        else:
            super().__init__(parameter,
                            "Could not parse input."
                            "Check parameter syntax.")

class ParameterNumberOfInputsError(ParameterError):
    '''Raised when the number of inputs is unexpected'''

    def __init__(self, parameter, found_and_expected=None):
        if found_and_expected:
            super().__init__(parameter,
                             f"Expected {found_and_expected[1]} inputs, "
                             "but found {found_and_expected[0]}.")
        else:
            super().__init__(parameter,
                            "Unexpected number of inputs.")

class ParameterRangeError(ParameterError):
    '''Raised when the value is not in the allowed range'''

    def __init__(self, parameter, given_value, allowed_range):
        message = (f"Value {given_value} is outside allowed range "
                   f"({allowed_range[0]} <= {parameter} <= "
                   f"{allowed_range[1]}).")
        super().__init__(parameter, message)

class ParameterUnknownFlagError(ParameterError):
    '''Raised when an unknown flag is encountered'''

    def __init__(self, parameter, flag):
        super().__init__(parameter,
                         f"Unknown flag '{flag}' encountered.")

class ParameterNeedsFlagError(ParameterError):
    '''Raised when a flag is needed but not given'''

    def __init__(self, parameter):
        super().__init__(parameter,
                         "Parameter requires a flag.")
