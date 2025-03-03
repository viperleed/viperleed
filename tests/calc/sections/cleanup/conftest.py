"""Test configuration for tests/calc/sections/cleanup.

Fixtures
--------
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

from viperleed.calc.classes.rparams import Rparams

_MODULE = 'viperleed.calc.sections.cleanup'


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
