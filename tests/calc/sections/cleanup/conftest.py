"""Test configuration for tests/calc/sections/cleanup.

Fixtures
--------
manifest
    A sample ManifestFile.
rpars
    An empty Rparams object.
workdir
    A temporary work directory.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files.manifest import ManifestFile

_MODULE = 'viperleed.calc.sections.cleanup'


@fixture(name='manifest')
def fixture_manifest():
    """Return the contents of a sample manifest file."""
    return ManifestFile('file1.txt', 'file2.txt', 'dir1')


@fixture(name='rpars')
def fixture_rpars():
    """Return an empty Rparams."""
    return Rparams()


@fixture(name='workdir')
def fixture_workdir(tmp_path):
    """Create a temporary work directory."""
    workdir = tmp_path / 'work_directory'
    workdir.mkdir()
    return workdir
