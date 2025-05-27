"""Module winerrors of viperleed.gui.measure.camera.drivers.imagingsource.

Defines functionality to catch errors as return values of the
tisgrabber C code.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-21'
__license__ = 'GPLv3+'


from enum import Enum

# pylint: disable=missing-param-doc,missing-type-doc
# Unreported bug. Multiple arguments listed on the left
# of the docstring appear not to be correctly identified,
# although the numpy style allows this.
def check_dll_return(success='int', include_errors=tuple(),
                     exclude_errors=tuple(), valid_returns=()):
    """Return an appropriate callable for checking DLL errors.

    Parameters
    ----------
    success : str, optional
        Condition for being a successful return. Allowed values
        are:
        * 'int': check integer return values
        * 'pointer': check that the return is a valid pointer
        * '>NUMBER': check that the result is larger than NUMBER
        * '>=NUMBER': check that the result is larger than, or
            equal to NUMBER
        Default is 'int'.
    include_errors, exclude_errors : Sequence, optional
        Which error return values should be included/excluded while
        checking the return. Only one of the two should be given.
        Each element is (the name of) a DLLReturns. This is meaningful
        only for success == 'int', as there are multiple returns with
        the same value. Default is an empty tuple for both, i.e., take
        all the returns.
    valid_returns : Sequence, optional
        Which return values are to also be considered 'valid'. This
        argument is used only for 'int' checking, and makes sense to
        be given only if the function is a checker, returning SUCCESS
        when something is found, and one of these values when it is
        not. SUCCESS is alway considered a valid return. Default is
        an empty tuple.

    Returns
    -------
    checker : callable
        Error checker suitable for using with ctypes
        _FuncPtr.errcheck attribute.

    Raises
    ------
    ValueError
        If both include_errors and exclude_errors are given, or
        if 'success' is not one of the acceptable checker specifiers.
    """
    if include_errors and exclude_errors:
        raise ValueError("Can only give include_errors "
                         "or exclude_errors, not both.")
    if include_errors:
        include = [getattr(DLLReturns, e) for e in include_errors]
        errors = {e.value: e for e in DLLReturns if e in include}
    elif exclude_errors:
        exclude = [getattr(DLLReturns, e) for e in exclude_errors]
        errors = {e.value: e
                  for e in DLLReturns
                  if e not in exclude}
    else:
        errors = DLLReturns.as_dict()

    if valid_returns:
        valid_returns = [getattr(DLLReturns, e).value for e in valid_returns]

    limit = 0.0
    success = success.replace(' ', '')
    if success.startswith('>='):
        limit = float(success[2:])
    elif success.startswith('>'):
        limit = float(success[1:])

    def int_checker(result, func, args):
        """Check validity of the return of a int-returning function."""
        # print(f"{func.__name__}{args} returned {result}")
        if result in (DLLReturns.SUCCESS.value, *valid_returns):
            # All good
            return result

        error = errors.get(result, None)
        if error is None:
            raise ImagingSourceError(
                f"Unkonwn or excluded error code {result}"
                )
        raise ImagingSourceError(
            f"{func.__name__}{args} returned error "
            f"{error.name}: {error.message}", err_code=error
            )

    def ge_checker(result, func, args):
        """Check that the return is larger or equal than a limit."""
        err_txt = f"{func.__name__}{args} returned {result} < {limit}."
        error = errors.get(result, None)
        if error is not None:
            err_txt += f" This is error {error.name}: {error.message}."
        if result < limit:
            raise ImagingSourceError(err_txt, err_code=error)
        return result

    def gt_checker(result, func, args):
        """Check that the return is larger than a limit."""
        err_txt = f"{func.__name__}{args} returned {result} <= {limit}."
        error = errors.get(result, None)
        if error is not None:
            err_txt += f" This is error {error.name}: {error.message}"
        if result <= limit:
            raise ImagingSourceError(err_txt, err_code=error)
        return result

    def pointer_checker(result, func, args):
        """Check that the return of the function is a valid pointer."""
        if not result:
            # pointer to NULL
            raise ImagingSourceError(
                f"{func.__name__}{args} returned a pointer to NULL",
                err_code = DLLReturns.NULL_POINTER
                )
        return result

    if success == 'int':
        return int_checker
    if success == 'pointer':
        return pointer_checker
    if success.startswith('>='):
        return ge_checker
    if success.startswith('>'):
        return gt_checker
    raise ValueError(f"Invalid success={success!r}. Should be 'int', "
                     "'pointer', '>NUMBER', or '>=NUMBER'")


class DLLReturns(tuple, Enum):
    """Enum of return values from the firmware DLL."""

    SUCCESS = (1, "No error")
    ERROR = (0, "A generic error occurred")
    NO_HANDLE = (-1, "Grabber handle is invalid. Call create_grabber().")
    NO_DEVICE = (-2,
                 "Method requires an open device, but no device is open. "
                 "Call open(dev_name).")
    CAMERA_PROPERTY_NOT_AVAILABLE = (-2, "Property not available.")
    NOT_AVAILABLE = (-3, "Device does not support a {!r} property.")
    NO_PROPERTYSET = (-3,
                      "The property set was not queried for the device. "
                      "Call IC_QueryPropertySet(), then retry.")
    DEFAULT_WINDOW_SIZE_SET = (-3,
                               "Failed setting size of live-view window. "
                               "Call IC_SetDefaultWindowPosition(false), "
                               "then retry.")
    NOT_IN_LIVE_MODE = (-3,
                        "Property {!r} can only/cannot be set while "
                        "in live mode. Switch mode, then retry.")
    PROPERTY_ITEM_NOT_AVAILABLE = (
        -4, "VCD property item {!r} is not supported by the device."
        )
    PROPERTY_ELEMENT_NOT_AVAILABLE = (
        -5, "VCD property {!r} has no element {!r}."
        )
    PROPERTY_ELEMENT_WRONG_INTERFACE = (
        -6, "VCD property/element {!r}/{!r} cannot be accessed this way."
        )
    INDEX_OUT_OF_RANGE = (-7, "Index {} is out of range.")
    WRONG_XML_FORMAT = (-1,
                        "File {} does not contain data or has invalid format.")
    WRONG_INCOMPATIBLE_XML = (-3, "File {} contains incompatible XML data.")
    NOT_ALL_PROPERTIES_RESTORED = (
        -4, "Some properties in file {} could not be restored."
        )
    DEVICE_NOT_FOUND = (
        -5, "Could not open device while loading settings from file."
    )
    FILE_NOT_FOUND = (35, "File {} does not exist.")
    UNKOWN = (-2000, "Unknown or unspecified error.")
    NULL_POINTER = (-2001, "Returned a pointer to NULL")
    INVALID_SINK_FORMAT = (-2002, "Invalid sink/video format")
    REQUIRES_REBOOT = (-2003, ".close(), reboot the camera, and retry")
    FAILED_TO_ENABLE = (-2004, "Failed to enable/disable properties.")

    @classmethod
    def as_dict(cls):
        """Return a dict of {value: message}.

        Notice that, since multiple errors have the same return
        value, this will return an incorrect dictionary most of
        the times.
        """
        return {e.value: e for e in cls}

    # pylint: disable=invalid-overridden-method
    # Bug? .value is a @types.DynamicClassAttribute, i.e., somewhat
    # similar to a @property, not a method/callable as the error
    # message suggests. Seemes solved from Issue #2306, but looks
    # like it actually is not.
    @property
    def value(self):
        """Reimplement .value to return only the numeric code."""
        return self[0]
    # pylint: enable=invalid-overridden-method

    @property
    def message(self):
        """Return the error message of the enum."""
        return self[1]


class ImagingSourceError(Exception):
    """An error specific to a Imaging Source camera."""

    def __init__(self, msg, *args, err_code=DLLReturns.UNKOWN, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        msg : str
            Error message
        *args : object
            Other arguments passed to Exception
        err_code : DLLReturns
            Error code. Accessible with the .err_code attribute
        **kwargs : object
            Other keyword arguments passed to Exception

        Returns
        -------
        None.
        """
        self.err_code = err_code
        super().__init__(msg, *args, **kwargs)
