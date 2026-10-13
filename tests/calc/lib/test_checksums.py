"""Tests for module checksums of viperleed.calc.lib."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.calc.files.tenserleed import KNOWN_TL_VERSIONS
from viperleed.calc.lib.checksums import VALID_CHECKSUMS


@parametrize(version=KNOWN_TL_VERSIONS)
def test_known_versions_have_checksums(version):
    """Check that _checksums.dat has entries for each known version."""
    checksums = VALID_CHECKSUMS.get(version)
    assert checksums
    assert all(checksums.values())
