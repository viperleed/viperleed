"""Module ...

Created: 2021-07-21
Author: Michele Riva
"""
from enum import Enum

def check_dll_return(success='int', include_errors=tuple(),
                     exclude_errors=tuple()):
    """Return an appropriate callable for checking DLL errors.

    Parameters
    ----------
    success : str, optional
        Condition for being a successful return. Allowed values
        are:
        * 'int': check integer return values
        * 'pointer': check that the return is a valid pointer
        * '>0': check that the result is positive
        Default is 'int'.
    include_errors, exclude_errors : Sequence, optional
        Which error return values should be included/excluded while
        checking the return. Only one of the two should be given.
        Each element is (the name of) a DLLReturns. This is meaningful
        only for success == 'int', as there are multiple returns with
        the same value. Default is tuple() for both, i.e., take all
        the returns.

    Returns
    -------
    checker : callable
        Error checker suitable for using with ctypes
        _FuncPtr.errcheck attribute.
    """
    if include_errors and exclude_errors:
        raise ValueError("Can only give include_errors "
                         "or exclude_errors, not both.")
    if include_errors:
        include = (getattr(DLLReturns, e) for e in include_errors)
        errors = {e.value: e for e in DLLReturns if e in include}
    elif exclude_errors:
        exclude = (getattr(DLLReturns, e) for e in exclude_errors)
        errors = {e.value: e
                  for e in DLLReturns
                  if e not in exclude}
    else:
        errors = DLLReturns.as_dict()

    limit = 0
    success = success.replace(' ', '')
    if success.startswith('>='):
            limit = float(success[2:])
    elif success.startswith('>'):
            limit = float(success[1:])

    def int_checker(result, func, args):
        """Check validity of the return of a int-returning function."""
        print(f"{func.__name__}{args} returned {result}")
        if result.value == DLLReturns.SUCCESS.value:
            # All good
            return

        error = errors.get(result.value, None)
        if error is None:
            raise ImagingSourceError(
                f"Unkonwn or excluded error code {result.value}"
                )
        raise ImagingSourceError(
            f"{func.__name__}{args} returned error "
            f"{error.name}: {error.message}"
            )

    def pointer_checker(result, func, args):
        """Check that the return of the function is a valid pointer."""
        print(f"{func.__name__}{args} returned {result}")
        if result.value is None:
            raise ImagingSourceError(
                f"{func.__name__}{args} returned a pointer to NULL"
                )
        return result

    def gt_checker(result, func, args):
        """Check that the return is larger than a limit."""
        print(f"{func.__name__}{args} returned {result}")
        err_txt = f"{func.__name__}{args} returned {result} <= {limit}."
        err = errors.get(result.value, None)
        if err is not None:
            err_txt += f" This is error {err.name}: {err.message}."
        if result.value <= limit:
            raise ImagingSourceError(err_txt)
        return result

    def ge_checker(result, func, args):
        """Check that the return is larger or equal than a limit."""
        print(f"{func.__name__}{args} returned {result}")
        err_txt = f"{func.__name__}{args} returned {result} < {limit}"
        err = errors.get(result.value, None)
        if err is not None:
            err_txt += f" This is error {err.name}: {err.message}."
        if result.value < limit:
            raise ImagingSourceError(err_txt)
        return result

    if success == 'int':
        return int_checker
    if success == 'pointer'
        return pointer_checker
    if success.startswith('>='):
        return ge_checker
    if success.startswith('>'):
        return gt_checker
    raise ValueError(f"Invalid {success=!r}. Should be 'int', 'pointer', "
                     "'>NUMBER', or '>=NUMBER'")


class DLLReturns(tuple, Enum):
    """Enum of return values from the firmware DLL."""

    SUCCESS = (1, "No error")
    ERROR = (0, "A generic error occurred")
    NO_HANDLE = (-1, "Grabber handle is invalid. Call create_grabber().")
    NO_DEVICE = (-2,
                 "Method requires an open device, but no device is open. "
                 "Call OpenVideoCaptureDevice().")
    NOT_AVAILABLE = (-3, "Device does not support a {!r} property.")
    NO_PROPERTYSET = (-3,
                      "The property set was not queried for the device. "
                      "Call IC_QueryPropertySet(), then retry.")
    DEFAULT_WINDOW_SIZE_SET = (-3,
                               "Failed setting size of live-view window. "
                               "Call IC_SetDefaultWindowPosition(false), "
                               "then retry.")
    NOT_IN_LIVEMODE = (-3, "Property {!r} can only be set while in live mode.")
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

    @classmethod
    def as_dict(cls):
        """Return a dict of {value: message}.

        Notice that, since multiple errors have the same return
        value, this will return an incorrect dictionary most of
        the times.
        """
        return {e.value: e for e in cls}

    @property
    def value(self):
        """Reimplement .value to return only the numeric code."""
        return self[0]

    @property
    def message(self):
        """Return the error message of the enum."""
        return self[1]


class ImagingSourceError(Exception):
    pass
