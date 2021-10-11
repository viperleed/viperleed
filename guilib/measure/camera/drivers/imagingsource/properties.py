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
                    )
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
