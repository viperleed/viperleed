"""Script lock_released_branch.

Adds the branch of a released version to the branches locked by a
repository ruleset. Used by .github/workflows/create_version_branch.yml.

The ruleset locks all version branches up to the latest release, with
as few patterns as possible. For example, after v0.14.2 its patterns
that look like versions are

    refs/heads/0.[0-9].*
    refs/heads/0.1[0-3].*
    refs/heads/0.14.[0-2]

On release of VERSION, these patterns are replaced by those for
VERSION, but only if they are exactly those for one of the earlier
releases (the final vX.Y.Z tags below VERSION, as found by git in the
current directory; not necessarily the highest one, as releases made
from version branches may not have updated the ruleset). Other
patterns (e.g., ~DEFAULT_BRANCH) are always kept. If an admin edited
the patterns that look like versions, they are left alone: the branch
of VERSION is appended as an exact pattern (unless a pattern already
matches it), and a warning is printed.

The ruleset is found by name and written back whole, with only its
patterns changed. It is read again right before it is written, and
the update is refused if it changed in the meantime. GitHub offers no
conditional write for rulesets, so an edit made in the instant
between that read and the write would still be lost: the ruleset is
read once more after the write, and the script fails if it differs
from what was written.

Usage: python lock_released_branch.py --repo OWNER/REPO --ruleset NAME
                                      VERSION
VERSION is X.Y.Z (no leading 'v'). Needs the GitHub CLI, with GH_TOKEN
holding a token with Administration: write on the repository.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-19'
__license__ = 'GPLv3+'

import argparse
import copy
from fnmatch import fnmatchcase
import json
import re
import subprocess
import sys

BRANCH_PREFIX = 'refs/heads/'
_FINAL_TAG_RE = re.compile(r'v(\d+)\.(\d+)\.(\d+)')
# Fields of a ruleset that change whenever it is written
_CHANGED_ON_WRITE = ('updated_at',)
# Fields of a ruleset that are written back when updating it: those
# documented as parameters of "Update a repository ruleset". The
# others (e.g., id, _links) are read-only
_RULESET_FIELDS = (
    'bypass_actors', 'conditions', 'enforcement', 'name', 'rules', 'target',
    )
_VERSION_LIKE_RE = re.compile(r'refs/heads/[0-9.*?\[\]-]+')


class LockError(Exception):
    """The ruleset cannot be updated safely."""


def lock_patterns(version):
    """Return the fewest patterns that match all versions <= version.

    Parameters
    ----------
    version : tuple of int
        The (major, minor, patch) version.

    Returns
    -------
    patterns : list of str
        Ruleset patterns, starting with BRANCH_PREFIX, matching the
        branches of all versions up to and including `version` (e.g.,
        refs/heads/0.[0-9].* for all 0.0 to 0.9 versions).
    """
    major, minor, patch = version
    patterns = [f'{g}.*' for g in _globs_below(major)]
    patterns.extend(f'{major}.{g}.*' for g in _globs_below(minor))
    patterns.extend(f'{major}.{minor}.{g}' for g in _globs_below(patch + 1))
    return [BRANCH_PREFIX + p for p in patterns]


def updated_patterns(include, version, released):
    """Return the patterns that also lock the branch of version.

    Parameters
    ----------
    include : list of str
        The current patterns of the ruleset.
    version : tuple of int
        The (major, minor, patch) version just released.
    released : iterable of tuple of int
        The (major, minor, patch) versions released before `version`.

    Returns
    -------
    new_include : list of str
        The new patterns. Equal to `include` if the branch is already
        locked.
    warning : str or None
        Why the patterns were not compacted, or None if they were
        (or did not need to be).
    """
    version_like = [p for p in include if _VERSION_LIKE_RE.fullmatch(p)]
    new = lock_patterns(version)
    if set(version_like) == set(new):
        return list(include), None
    # No version patterns at all are expected only before any release
    expected = [set(lock_patterns(v)) for v in released] or [set()]
    if set(version_like) in expected:
        others = [p for p in include if p not in version_like]
        return new + others, None
    branch = BRANCH_PREFIX + '.'.join(str(v) for v in version)
    if any(fnmatchcase(branch, p) for p in include):
        return list(include), None
    warning = (f'Version patterns {version_like} are not those of any '
               f'earlier release, so they were left alone and {branch} '
               'was appended. Consider compacting them by hand.')
    return [*include, branch], warning


def update_ruleset(repo, ruleset_name, version):
    """Lock the branch of version in the ruleset named ruleset_name.

    Parameters
    ----------
    repo : str
        The OWNER/REPO repository.
    ruleset_name : str
        The name of the repository ruleset to update.
    version : tuple of int
        The (major, minor, patch) version just released.

    Raises
    ------
    LockError
        If the ruleset does not exist, if it changed between reading
        and writing it, or if it differs from what was written after
        writing it.
    """
    ruleset_id = _find_ruleset(repo, ruleset_name)
    url = f'repos/{repo}/rulesets/{ruleset_id}'
    ruleset = _gh_api(url)
    include = ruleset['conditions']['ref_name']['include']
    new_include, warning = updated_patterns(include, version,
                                            _released_before(version))
    if warning:
        print(f'::warning::{warning}')
    if new_include == include:
        print(f'Ruleset {ruleset_name!r} already locks {version}')
        return
    if _gh_api(url) != ruleset:
        raise LockError(f'Ruleset {ruleset_name!r} changed while being '
                        'updated. Run the workflow again.')
    expected = copy.deepcopy(ruleset)
    expected['conditions']['ref_name']['include'] = new_include
    _gh_api(url, method='PUT',
            body={key: expected[key] for key in _RULESET_FIELDS})
    # Compare all fields, not just those written, so that any field
    # that the write did not keep (e.g., one added to the API later)
    # is noticed
    written = _gh_api(url)
    differ = sorted(k for k in expected.keys() | written.keys()
                    if k not in _CHANGED_ON_WRITE
                    and written.get(k) != expected.get(k))
    if differ:
        raise LockError(f'Ruleset {ruleset_name!r} differs from what was '
                        f'written in {differ}. Check it by hand.')
    print(f'Ruleset {ruleset_name!r}: patterns {include} -> {new_include}')


def _find_ruleset(repo, name):
    """Return the ID of the repository ruleset called name."""
    rulesets = _gh_api(f'repos/{repo}/rulesets?per_page=100'
                       '&includes_parents=false')
    ids = [r['id'] for r in rulesets if r['name'] == name]
    if len(ids) != 1:
        raise LockError(f'Expected one ruleset named {name!r} in {repo}, '
                        f'found {len(ids)}')
    return ids[0]


def _gh_api(endpoint, method='GET', body=None):
    """Call the GitHub API via the gh CLI and return the JSON reply."""
    cmd = ['gh', 'api', '--method', method, endpoint]
    if body is not None:
        cmd.extend(('--input', '-'))
    result = subprocess.run(cmd,
                            input=None if body is None else json.dumps(body),
                            capture_output=True, text=True, check=False)
    if result.returncode:
        raise LockError(f'gh api {method} {endpoint} failed: '
                        f'{result.stderr.strip()}')
    return json.loads(result.stdout)


def _globs_below(number):
    """Return the fewest globs matching integers 0 to number - 1."""
    if number <= 0:
        return []
    globs = []
    largest = number - 1
    for n_digits in range(1, len(str(largest)) + 1):
        start = 0 if n_digits == 1 else 10**(n_digits - 1)
        stop = min(largest, 10**n_digits - 1)
        globs.extend(_globs_between(str(start), str(stop)))
    return globs


def _globs_between(low, high):
    """Return the fewest globs matching integers low to high.

    Parameters
    ----------
    low, high : str
        Decimal integers with as many digits, with low <= high.

    Returns
    -------
    globs : list of str
        Patterns matching exactly the integers from low to high.
    """
    if len(low) == 1:
        return [_digit_class(int(low), int(high))]
    if low[0] == high[0]:
        return [low[0] + g for g in _globs_between(low[1:], high[1:])]
    n_rest = len(low) - 1
    head, tail = [], []
    first, last = int(low[0]), int(high[0])
    if low[1:] != '0' * n_rest:
        head = [low[0] + g for g in _globs_between(low[1:], '9' * n_rest)]
        first += 1
    if high[1:] != '9' * n_rest:
        tail = [high[0] + g for g in _globs_between('0' * n_rest, high[1:])]
        last -= 1
    middle = []
    if first <= last:
        middle = [_digit_class(first, last) + '[0-9]' * n_rest]
    return head + middle + tail


def _digit_class(first, last):
    """Return a glob matching one digit from first to last."""
    return str(first) if first == last else f'[{first}-{last}]'


def _released_before(version):
    """Return the final releases below version, from the git tags."""
    tags = subprocess.run(['git', 'tag', '--list', 'v*'],
                          capture_output=True, text=True, check=True)
    releases = (_FINAL_TAG_RE.fullmatch(t) for t in tags.stdout.split())
    versions = (tuple(int(n) for n in m.groups()) for m in releases if m)
    return sorted(v for v in versions if v < version)


def main():
    """Lock the branch of a released version in a ruleset."""
    parser = argparse.ArgumentParser(
        description='Add the branch of a released version to the '
                    'branches locked by a repository ruleset.',
        )
    parser.add_argument('--repo', required=True, help='OWNER/REPO')
    parser.add_argument('--ruleset', required=True, help='ruleset name')
    parser.add_argument('version', help="released version, e.g., '0.15.0'")
    args = parser.parse_args()
    match = re.fullmatch(r'(\d+)\.(\d+)\.(\d+)', args.version)
    if not match:
        sys.exit(f'::error::Version {args.version!r} is not X.Y.Z')
    version = tuple(int(n) for n in match.groups())
    try:
        update_ruleset(args.repo, args.ruleset, version)
    except LockError as exc:
        sys.exit(f'::error::{exc}')


if __name__ == '__main__':
    main()
