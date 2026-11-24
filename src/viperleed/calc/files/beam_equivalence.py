"""Functions for reading and writing the BEAM_EQUIVALENCE file.

BEAM_EQUIVALENCE assigns calculated beams to experimental beams
explicitly, one beam assignment per line, for example
    (1 0) = (1 0) + 0.5*(-1 0)
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import re

from quicktions import Fraction

from viperleed.calc.files.beams import _snap_to_fraction

# One beam index: an integer, a decimal, or a fraction like 1/2
_INDEX = r'[+-]?(?:\d+/\d+|\d+(?:\.\d*)?|\.\d+)'

# A beam label: two indices in parentheses, separated by
# whitespace, a comma, or a vertical bar, e.g., (1/2 | -1)
_LABEL_RE = re.compile(
    rf'\(\s*(?P<first>{_INDEX})(?:\s*[,|]\s*|\s+)(?P<second>{_INDEX})\s*\)'
    )


def parse_beam_label(label):
    """Return the (h, k) indices of a beam label as Fractions.

    Parameters
    ----------
    label : str
        A beam label such as '(1 0)', '(1/2, -1)', or '(0.5 | 0)'.
        Indices are integers, fractions, or decimals, separated by
        whitespace, a comma, or a vertical bar. Whitespace is allowed
        anywhere inside the parentheses and around them.

    Returns
    -------
    tuple
        The two indices (h, k) as Fractions. Decimals that are within
        the rounding to six decimals of a fraction with a small
        denominator are snapped to that fraction, as for IVBEAMS.

    Raises
    ------
    ValueError
        If `label` is not two indices in parentheses.
    """
    match = _LABEL_RE.fullmatch(label.strip())
    if not match:
        msg = (f'Invalid beam label {label.strip()!r}. Expected two '
               'indices in parentheses, e.g., (1/2 -1).')
        raise ValueError(msg)
    try:
        return tuple(_snap_to_fraction(Fraction(match[group]))
                     for group in ('first', 'second'))
    except ZeroDivisionError:
        msg = f'Invalid beam label {label.strip()!r}: division by zero.'
        raise ValueError(msg) from None
