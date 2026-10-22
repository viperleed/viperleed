"""Module update_copyright_year for viperleed files.

Sets the end year of the __copyright__ module dunder of given files
to the current year. Meant to be used as a pre-commit hook: prints
the files it changed, and exits with a non-zero code if any file was
changed. Only files that are staged in git as added, modified or
renamed are considered; other arguments are skipped silently, so
that `pre-commit run --all-files` leaves untouched files alone. Only
dunders of the exact form

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
import os
from pathlib import Path
import re
import subprocess
import sys


COPYRIGHT_RE = re.compile(
    rb'^(?P<head>__copyright__ = (?P<quote>[\'"])'
    rb'Copyright \(c\) [0-9]{4}-)'
    rb'(?P<end>[0-9]{4})'
    rb'(?P<tail> ViPErLEED developers(?P=quote))',
    re.MULTILINE,
    )


def get_staged_files():
    """Return the resolved paths of files staged as added/modified/renamed.

    Returns
    -------
    staged : set of Path
        Absolute paths of the files in the git index that differ
        from HEAD (or all of them before the first commit), except
        deleted ones.
    """
    def _git(*args):
        return subprocess.run(('git', *args), stdout=subprocess.PIPE,
                              check=True).stdout
    toplevel = Path(os.fsdecode(_git('rev-parse', '--show-toplevel').strip()))
    names = _git('diff', '--cached', '--name-only', '--diff-filter=AMR', '-z')
    return {(toplevel / os.fsdecode(name)).resolve()
            for name in names.split(b'\0') if name}


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


def main():
    """Update the copyright end year of staged files in sys.argv.

    Returns
    -------
    exit_code : int
        1 if any file was changed, 0 otherwise.
    """
    year = date.today().year
    staged = get_staged_files()
    exit_code = 0
    for filename in sys.argv[1:]:
        if Path(filename).resolve() not in staged:
            continue
        if update_copyright_year(filename, year):
            print(f'Updated copyright year in {filename}')
            exit_code = 1
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
