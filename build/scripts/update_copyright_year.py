"""Module update_copyright_year for viperleed files.

Sets the end year of the __copyright__ module dunder of given files
to the current year. Meant to be used as a pre-commit hook: prints
the files it changed, and exits with a non-zero code if any file was
changed. Only dunders of the exact form

    __copyright__ = 'Copyright (c) <start>-<end> ViPErLEED developers'

(single or double quotes) are considered, and only those whose end
year is older than the current one are rewritten. Everything else in
the files, including line endings and encoding, is left untouched.

Usage: python update_copyright_year.py FILE [FILE ...]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-22'
__license__ = 'GPLv3+'

from datetime import date
from pathlib import Path
import re
import sys


COPYRIGHT_RE = re.compile(
    rb'^(?P<head>__copyright__ = (?P<quote>[\'"])'
    rb'Copyright \(c\) [0-9]{4}-)'
    rb'(?P<end>[0-9]{4})'
    rb'(?P<tail> ViPErLEED developers(?P=quote))',
    re.MULTILINE,
    )


def update_copyright_year(path, year):
    """Set a stale copyright end year in file `path` to `year`.

    Parameters
    ----------
    path : str or Path
        The file to update.
    year : int
        The year to use as end year in the __copyright__ dunder.

    Returns
    -------
    changed : bool
        Whether the file was modified.
    """
    path = Path(path)
    contents = path.read_bytes()
    new_year = str(year).encode('ascii')

    def _replace(match):
        if int(match['end']) >= year:
            return match[0]
        return match['head'] + new_year + match['tail']

    new_contents = COPYRIGHT_RE.sub(_replace, contents)
    if new_contents == contents:
        return False
    path.write_bytes(new_contents)
    return True


def main(argv=None):
    """Update the copyright end year of files in argv. Return exit code."""
    if argv is None:
        argv = sys.argv[1:]
    year = date.today().year
    exit_code = 0
    for filename in argv:
        if update_copyright_year(filename, year):
            print(f'Updated copyright year in {filename}')
            exit_code = 1
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
