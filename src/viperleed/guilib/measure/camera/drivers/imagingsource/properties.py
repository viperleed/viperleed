"""Module cameraproperties of Imaging Source driver.

Collects utilities to access th VCD properties of a Imaging Source
camera.
"""

from enum import Enum
from ctypes import CFUNCTYPE, c_int, c_char_p, py_object

# Prototype of frame-ready callbacks:
PropertyDiscoveryCallbackType = CFUNCTYPE(
    c_int,     # return type
    c_char_p,  # name of property/element/interface
    py_object  # data (python) useful to store the name
    )

@PropertyDiscoveryCallbackType
def store_property(name, storage):
    """Store a property in a python dict.

    This function is used as the callback for listing properties.

    Parameters
    ----------
    name : bytes
        Name of property/element/interface
    storage : dict
        Property names will be stored as keys, values
        are empty dictionaries.

    Returns
    -------
    finished : int
        Whether property search is over or not. This function
        always returns 0, such that all properties can be discovered.
    """
    storage[name.decode()] = {}
    return 0


class CameraProperty(Enum):
    """Class holding possible camera properties.

    Does not include VCD properties, which are accessed by name.
    """

    PAN = 0
    TILT = 1
    ROLL = 2
    ZOOM = 3
    EXPOSURE = 4
    IRIS = 5
    FOCUS = 6
    MEGA = 65536  # Guard for max enum value, i.e., 2**sizeof(int)

    @classmethod
    def get(cls, cam_property):
        """Return an appropriate enum from value or name.

        Parameters
        ----------
        cam_property : str or int or CameraProperty
            Some reference to a camera property. If a string, lookup
            is done by name (not case sensitive), otherwise by value.

        Returns
        -------
        property : CameraProperty
            The attribute corresponding to cam_property

        Raises
        ------
        ValueError
            If no property could be found
        """
        if isinstance(cam_property, str):
            cam_property = getattr(cls, cam_property.upper(), None)
            if cam_property is None:
                raise ValueError(f"Unknown camera property {cam_property}. "
                                 "Perhaps you meant video property?")
        else:
            try:
                cam_property = cls(cam_property)
            except ValueError as err:
                raise ValueError(f"Unknown camera property {cam_property}. "
                                 "Perhaps you meant video property?") from err
        return cam_property


class VideoProperty(Enum):
    """Class holding possible video properties."""
    BRIGHTNESS = 0
    CONTRAST = 1
    HUE = 2
    SATURATION = 3
    SHARPNESS = 4
    GAMMA = 5
    COLORENABLE = 6
    WHITEBALANCE = 7
    BLACK_LIGHT_COMPENSATION = 8
    GAIN = 9
    MEGA = 65536  # Guard for max enum value, i.e., 2**sizeof(int)

class SwitchProperty(Enum):
    """On/Off VCD property enum."""

    ON = 1
    OFF = 0

    @classmethod
    def get(cls, on_off):
        """Return an appropriate enum from value or name.

        Parameters
        ----------
        on_off : str, SwitchProperty or object
            If a str, should be 'on'/'off' (not case sensitive).
            If an object, its truth value will be used.

        Returns
        -------
        switch_property : SwitchProperty
            The attribute corresponding to on_off.

        Raises
        ------
        ValueError
            If a string is passed that is neither 'on' nor 'off'
            (not case sensitive)
        """
        if isinstance(on_off, str):
            try:
                switch_property = getattr(cls, on_off.upper())
            except AttributeError as err:
                raise ValueError(
                    f"Invalid string {on_off} for SwitchProperty. "
                    "Should be 'on' or 'off' (not case sensitive)"
                    ) from err
        else:
            switch_property = cls(bool(on_off))
        return switch_property


class VCDPropertyInterface(Enum):
    """List of VCD property interfaces with their methods."""

    RANGE = int
    SWITCH = bool
    BUTTON = 'one_push'  # == Activate once
    MAPSTRINGS = str
    ABSOLUTEVALUE = float
    UNKNOWN = None

    @classmethod
    def get(cls, name_or_value):
        """Return attribute by name."""
        if isinstance(name_or_value, str):
            interface = getattr(cls, name_or_value.upper(), None)
            if interface is None:
                raise ValueError("Invalid property interface "
                                 f"name {name_or_value}")
            return interface
        try:
            interface = cls(name_or_value)
        except AttributeError as err:
            raise ValueError("Invalid property interface "
                             f"{name_or_value}") from err
        return interface
