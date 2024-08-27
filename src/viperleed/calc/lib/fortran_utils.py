"""Module fortran_utils of viperleed.calc.lib.

Collects functionality useful for handling FORTRAN code.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'


def fortranContLine(s):
    """Takes a sting that might be too long to fit on a single line of fortran
    code and splits it into continuation lines as necessary. Returns a string
    with ampersands and line breaks."""
    limit = 72  # fortran line length limit
    if len(s) <= limit:
        return s
    o = s[:6]
    s = s[6:]
    while len(s) > (limit-6):   # 6 characters for beginning of line
        o += s[:(limit-6)]
        o += "&\n     &"
        s = s[(limit-6):]
    o += s
    return o
