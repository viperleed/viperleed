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

import functools
from dataclasses import dataclass
from dataclasses import fields as data_fields
from dataclasses import is_dataclass
import typing
import sys

try:  # Import stuff that was made available since py3.8
    from typing import get_args as ty_get_args
except ImportError:
    from typing_extensions import get_args as ty_get_args
    from typing_extensions import get_origin as ty_get_origin
else:
    from typing import get_origin as ty_get_origin

if sys.version_info < (3, 12):  # No frozen_default keyword in py3.11
    from typing_extensions import dataclass_transform
else:
    from typing import dataclass_transform


frozen = dataclass_transform(frozen_default=True)(
    functools.partial(dataclass, frozen=True)
    )


def check_types(self, init_only=False):
    """Raise TypeError if fields don't match their type hints."""
    for field in data_fields(self):
        if init_only and not field.init:
            continue
        attr = field.name
        value = getattr(self, attr)
        actual_types = tuple(_find_type_origin(field.type))
        if actual_types and not isinstance(value, actual_types):
            raise TypeError(
                    f'Expected type {field.type} for argument {attr!r} '
                    f'but received type {type(value).__name__!r} instead'
                    )


def is_optional_field(field, with_type=None):
    """Return whether a datclasses.Field is optional."""
    if with_type is None:
        with_type = field.type
    # The following trick works because Optional removes 'repeated'
    # entries, so that Optional[Optional[t]] == Optional[t]
    return field.type == typing.Optional[with_type]


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
        non_init = (f for f in data_fields(self) if not f.init)
        if not any(f.name == attr_name for f in non_init):
            raise
    object.__setattr__(self, attr_name, attr_value)


_SpecialType = tuple({  # These types are not pubic API
    type(typing.Optional),
    type(typing.Any),   # In Py3.12 they are different
    })


def _find_type_origin(type_hint):
    """Yield nested type origins from a `type_hint`."""
    # Adapted from https://stackoverflow.com/questions/50563546
    if isinstance(type_hint, _SpecialType):
        # Special types, without extra parameters
        return

    actual_type = ty_get_origin(type_hint) or type_hint
    if isinstance(actual_type, _SpecialType):
        # Case of typing.Union[...] or typing.ClassVar[...] or ...
        args = ty_get_args(type_hint)
        all_origins = tuple(tuple(_find_type_origin(arg)) for arg in args)
        for origins in all_origins:
            yield from origins
    else:
        yield actual_type
