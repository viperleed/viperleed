"""Module fortran_utils of viperleed.calc.lib.

Collects functionality useful for handling FORTRAN code.
"""

__authors__ = (
    "Florian Kraushofer (@fkraushofer)",
    "Alexander M. Imre (@amimre)",
    "Michele Riva (@michele-riva)",
)
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'

import re
import subprocess
import shutil

from .version import Version

_FORTRAN_LINE_LENGTH = 72  # FORTRAN line-length limit
_F77_CONTINUATION_POS = 6  # Column of the continuation character


class FortranCompilerError(Exception):
    """Base class for exceptions related to the Fortran compiler(s)."""


class CompilerNotFoundError(FortranCompilerError):
    """Raised when the Fortran compiler is not found."""


class NoCompilerVersionFoundError(FortranCompilerError):
    """Raised when the Fortran version could not be determined."""


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


def get_mpifort_version():
    """Check the version of the mpifort compiler."""
    # use sed to extract the version number; standard GNU util
    version_nr_call = ['mpifort', '--version']
    # check if mpifort is installed
    if not shutil.which('mpifort'):
        raise CompilerNotFoundError

    # get version number
    try:
        result = subprocess.run(
            version_nr_call, shell=False, check=True, capture_output=True
        )
    except subprocess.CalledProcessError:
        raise NoCompilerVersionFoundError
    output = result.stdout.decode().strip()
    version_nr_str = re.search(r"GNU Fortran.*\) (\d+\.\d+\.\d+)", output)
    mpifort_version = Version(version_nr_str.strip())

    return mpifort_version
