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

from typing import Optional


def is_optional_field(field, with_type=None):
    """Return whether a datclasses.Field is optional."""
    if with_type is None:
        with_type = field.type
    # The following trick works because Optional removes 'repeated'
    # entries, so that Optional[Optional[t]] == Optional[t]
    return field.type == Optional[with_type]
