"""Module fortran_utils of viperleed.calc.lib.

Collects functionality useful for handling FORTRAN code.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'

import re
import shutil
import subprocess

from viperleed.calc.lib.version import Version

_FORTRAN_LINE_LENGTH = 72  # FORTRAN line-length limit
_F77_CONTINUATION_POS = 6  # Column of the continuation character


class FortranCompilerError(Exception):
    """Base class for exceptions related to the Fortran compiler(s)."""


class CompilerNotFoundError(FortranCompilerError):
    """Raised when a Fortran compiler is not found."""


class NoCompilerVersionFoundError(FortranCompilerError):
    """The version of a Fortran compiler could not be determined."""


def get_mpifort_version():
    """Return the version of the mpifort compiler."""
    # Check if mpifort is installed
    if not shutil.which('mpifort'):
        raise CompilerNotFoundError

    # Get version number
    version_nr_call = 'mpifort', '--version'
    try:
        result = subprocess.run(version_nr_call,
                                capture_output=True,
                                check=True)
    except subprocess.CalledProcessError as exc:
        msg = str(exc)
        if exc.stdout:
            msg += f'\noutput:{exc.stdout.decode()}'
        if exc.stderr:
            msg += f'\nerrors:\n{exc.stderr.decode()}'
        raise NoCompilerVersionFoundError(msg) from None

    output = result.stdout.decode().strip()
    match_ = re.match(r'GNU Fortran.*\) (?P<version>\d+\.\d+\.\d+)',
                      output)
    if not match_:
        raise NoCompilerVersionFoundError('Could not determine version '
                                          f'from\n{output}')
    return Version(match_['version'])


def wrap_fortran_line(string):
    """Wrap a FORTRAN string into continuation lines with ampersands."""
    if len(string) <= _FORTRAN_LINE_LENGTH:
        return string
    # The first line is _FORTRAN_LINE_LENGTH characters, the others
    # need to be wrapped to chunk_size below in order to fit the
    # continuation character at _F77_CONTINUATION_POS. Pull out
    # the first _F77_CONTINUATION_POS characters from the beginning
    # so we can simply split the rest into chunk_size-long parts.
    head, rest = string[:_F77_CONTINUATION_POS], string[_F77_CONTINUATION_POS:]
    chunk_size = _FORTRAN_LINE_LENGTH - _F77_CONTINUATION_POS
    continuation_lines = (rest[i:i + chunk_size]
                          for i in range(0, len(rest), chunk_size))
    sep = f'&\n{"&":>{_F77_CONTINUATION_POS}}'
    return head + sep.join(continuation_lines)
