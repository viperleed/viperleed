"""Script check_dist_version.

Checks that a release tag matches the version of the built sdist and
wheel of viperleed. Used by .github/workflows/distribute_built_dist.yml.

Usage: python check_dist_version.py TAG DIST_DIR
Exits with a GitHub Actions error annotation unless DIST_DIR contains
exactly one sdist and one wheel, and TAG == 'v' + the version of each.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-21'
__license__ = 'GPLv3+'

import argparse
from email.parser import Parser
from pathlib import Path
import sys
import tarfile
import zipfile


class DistVersionError(Exception):
    """The built distributions do not match the tag."""


def check_dist_version(tag, dist_dir):
    """Raise DistVersionError unless tag == 'v' + each dist's version.

    Parameters
    ----------
    tag : str
        The release tag, e.g., 'v0.15.0'.
    dist_dir : str or Path
        Directory containing exactly one sdist (*.tar.gz) and one
        wheel (*.whl).

    Raises
    ------
    DistVersionError
        If there is not exactly one sdist and one wheel, if one of
        them has no top-level metadata file, or if the tag does not
        match the version of either of them.
    """
    dist_dir = Path(dist_dir)
    sdists = sorted(dist_dir.glob('*.tar.gz'))
    wheels = sorted(dist_dir.glob('*.whl'))
    if len(sdists) != 1 or len(wheels) != 1:
        raise DistVersionError(
            'Expected exactly one sdist and one wheel, '
            f'found {[p.name for p in dist_dir.iterdir()]}'
            )
    with tarfile.open(sdists[0]) as sdist:
        sdist_version = _version_from(
            lambda n: sdist.extractfile(n).read().decode(),
            sdist.getnames(),
            '/PKG-INFO',
            )
    with zipfile.ZipFile(wheels[0]) as wheel:
        wheel_version = _version_from(
            lambda n: wheel.read(n).decode(),
            wheel.namelist(),
            '.dist-info/METADATA',
            )
    print(f'Tag: {tag}; sdist: {sdist_version}; wheel: {wheel_version}')
    for kind, version in (('sdist', sdist_version),
                          ('wheel', wheel_version)):
        if tag != f'v{version}':
            raise DistVersionError(
                f'Tag {tag} does not match {kind} version {version}'
                )


def _version_from(read_member, names, suffix):
    """Return the Version field of the top-level metadata file."""
    for name in names:
        # One level deep: <name-version>/PKG-INFO in the sdist,
        # <name-version>.dist-info/METADATA in the wheel, not
        # copies nested further down (e.g., in *.egg-info/).
        if name.count('/') == 1 and name.endswith(suffix):
            return Parser().parsestr(read_member(name))['Version']
    raise DistVersionError(f'No {suffix} found')


def main():
    """Check the distributions in a directory against a tag."""
    parser = argparse.ArgumentParser(
        description='Check that TAG matches the version of the built '
                    'sdist and wheel in DIST_DIR.',
        )
    parser.add_argument('tag', help="release tag, e.g., 'v0.15.0'")
    parser.add_argument('dist_dir', help='directory with sdist and wheel')
    args = parser.parse_args()
    try:
        check_dist_version(args.tag, args.dist_dir)
    except DistVersionError as exc:
        sys.exit(f'::error::{exc}')


if __name__ == '__main__':
    main()
