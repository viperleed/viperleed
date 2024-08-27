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


_FORTRAN_LINE_LENGTH = 72  # FORTRAN line-length limit
_F77_CONTINUATION_POS = 6  # Column of the continuation character


def wrap_fortran_line(string):
    """Wrap a FORTRAN string into continuation lines with ampersands."""
    if len(string) <= _FORTRAN_LINE_LENGTH:
        return string
    head, rest = string[:_F77_CONTINUATION_POS], string[_F77_CONTINUATION_POS:]
    chunk_size = _FORTRAN_LINE_LENGTH - _F77_CONTINUATION_POS
    continuation_lines = (rest[i:i + chunk_size]
                          for i in range(0, len(rest), chunk_size))
    sep = f'&\n{"&":>{_F77_CONTINUATION_POS}}'
    return head + sep.join(continuation_lines)
