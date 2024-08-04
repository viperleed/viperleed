"""Module dataclasses_utils of viperleed.calc.lib.

Collects useful functions that add up on top of the dataclasses
stdlib module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-25'
__license__ = 'GPLv3+'

from dataclasses import fields
from dataclasses import is_dataclass
from typing import Optional


def is_optional_field(field_, with_type=None):
    """Return whether a datclasses.Field is optional."""
    if with_type is None:
        with_type = field_.type
    # The following trick works because Optional removes 'repeated'
    # entries, so that Optional[Optional[t]] == Optional[t]
    return field_.type == Optional[with_type]


# TODO: this may be made stricter to ensure it is only used
# inside the source code of a dataclass and not from user code.
# One can use sys._get_frame (CPython only) as suggested in
# https://stackoverflow.com/questions/900392 or
# https://stackoverflow.com/questions/2654113
def set_frozen_attr(self, attr_name, attr_value):
    """Set self.attr_name to attr_value even if self is frozen.

    It is **VERY IMPORTANT** to realize that this function is a
    workaround that makes a frozen dataclass not actually frozen.
    It is thus critical that this function is used **only** in
    a context where it makes sense to actually modify an attribute
    of a dataclass that is supposed to be immutable. This usually
    means that this function should be used ONLY INTERNALLY in the
    implementation of a dataclass.

    Parameters
    ----------
    self : dataclass
        The dataclass whose attribute should be set.
    attr_name : str
        The name of the attribute to set.
    attr_value : object
        The value of the attribute.

    Raises
    ------
    AttributeError
        If self has no `attr_name` attribute.
    """
    if not is_dataclass(self):
        raise TypeError('Cannot use set_frozen_attr on non-dataclass objects')
    try:
        getattr(self, attr_name)
    except AttributeError:
        # It may still be that the attribute is a init=False field.
        if not any(f.name == attr_name for f in fields(self) if not f.init):
            raise
    object.__setattr__(self, attr_name, attr_value)
