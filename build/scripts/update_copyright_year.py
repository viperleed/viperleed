"""Module update_copyright_year for viperleed files.

Sets the end year of the __copyright__ module dunder of given files
to the current year. Meant to be used as a pre-commit hook: prints
the files it changed, and exits with a non-zero code if any file was
changed. Only files that are staged in git (except deleted ones) are
considered; other arguments are skipped silently, so that
`pre-commit run --all-files` leaves untouched files alone. When run
via `pre-commit run --from-ref A --to-ref B` (which sets the
PRE_COMMIT_FROM_REF and PRE_COMMIT_TO_REF environment variables),
the files changed between the two refs (except deleted ones) are
considered instead of the staged ones. git is
required: if it is not found, or not run inside a repository, the
script exits with an error without touching any file. Only
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


class GitError(Exception):
    """git could not be run, or a git command failed."""


def _git(*args):
    """Return the standard output of git called with `args`.

    Raises
    ------
    GitError
        If git is not found or exits with a non-zero code.
    """
    try:
        result = subprocess.run(('git', *args), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, check=False)
    except OSError as exc:
        raise GitError(f'could not run git ({exc})') from exc
    if result.returncode:
        error = os.fsdecode(result.stderr).strip().splitlines()
        reason = error[-1] if error else f'exit code {result.returncode}'
        command = ' '.join(('git', *args))
        raise GitError(f'"{command}" failed: {reason}')
    return result.stdout


def _get_toplevel():
    """Return the root directory of the current git repository.

    Raises
    ------
    GitError
        If git is not found or fails, e.g., outside a repository.
    """
    return Path(os.fsdecode(_git('rev-parse', '--show-toplevel').strip()))


def _resolve_dirs(path):
    """Return `path` as absolute, with symlinks resolved except its name.

    The file name itself is not resolved, so that a symlink and its
    target remain different paths.
    """
    path = Path(path).absolute()
    return path.parent.resolve() / path.name


def _resolve_names(names, toplevel):
    """Return absolute paths for the NUL-separated git paths `names`."""
    return {_resolve_dirs(toplevel / os.fsdecode(name))
            for name in names.split(b'\0') if name}


def get_changed_files(from_ref, to_ref, toplevel):
    """Return the resolved paths of files changed between two refs.

    The files are those pre-commit picks for `pre-commit run
    --from-ref from_ref --to-ref to_ref`, i.e., those changed on
    `to_ref` since its merge base with `from_ref`, except deleted ones.

    Parameters
    ----------
    from_ref : str
        The git ref to compare against.
    to_ref : str
        The git ref whose changes are taken.
    toplevel : Path
        The root directory of the repository.

    Returns
    -------
    changed : set of Path
        Absolute paths of the changed files.

    Raises
    ------
    GitError
        If git is not found or fails, e.g., for an unknown ref.
    """
    diff = ('diff', '--name-only', '--no-ext-diff', '--diff-filter=d', '-z')
    try:
        names = _git(*diff, f'{from_ref}...{to_ref}')
    except GitError:
        # Like pre-commit: refs without a merge base need a full diff
        names = _git(*diff, f'{from_ref}..{to_ref}')
    return _resolve_names(names, toplevel)


def get_staged_files(toplevel):
    """Return the resolved paths of all staged files except deleted ones.

    Parameters
    ----------
    toplevel : Path
        The root directory of the repository.

    Returns
    -------
    staged : set of Path
        Absolute paths of the files in the git index that differ
        from HEAD (or all of them before the first commit), except
        deleted ones.

    Raises
    ------
    GitError
        If git is not found or fails.
    """
    names = _git('diff', '--cached', '--name-only', '--diff-filter=d', '-z')
    return _resolve_names(names, toplevel)


def get_files_to_update():
    """Return the resolved paths of the files that may be updated.

    These are the files changed between PRE_COMMIT_FROM_REF and
    PRE_COMMIT_TO_REF if pre-commit sets both, the staged files
    otherwise. Deleted files are never included.

    Raises
    ------
    GitError
        If git is not found or fails, e.g., outside a repository.
    """
    toplevel = _get_toplevel()
    from_ref = os.environ.get('PRE_COMMIT_FROM_REF')
    to_ref = os.environ.get('PRE_COMMIT_TO_REF')
    if from_ref and to_ref:
        return get_changed_files(from_ref, to_ref, toplevel)
    return get_staged_files(toplevel)


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
    """Update the copyright end year of changed files in sys.argv.

    Only files that are staged, or changed between the refs given
    to `pre-commit run --from-ref --to-ref`, are considered.

    Returns
    -------
    exit_code : int
        1 if any file was changed, 0 otherwise. 2 if git is
        unavailable or fails, in which case no file is touched.
    """
    year = date.today().year
    try:
        candidates = get_files_to_update()
    except GitError as exc:
        print(f'update_copyright_year: error: {exc}', file=sys.stderr)
        return 2
    exit_code = 0
    for filename in sys.argv[1:]:
        if _resolve_dirs(filename) not in candidates:
            continue
        if update_copyright_year(filename, year):
            print(f'Updated copyright year in {filename}')
            exit_code = 1
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
