"""Tests for module viperleed.calc.files.vibrocc."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-06-18'
__license__ = 'GPLv3+'

import pytest

from viperleed.calc.files import tenserleed
from viperleed.calc.lib.version import Version

# check validity of old and new version name formats
@pytest.mark.parametrize('version_str', ('1.6.0', '1.7.3', '1.73'))
def test_detect_version_number_from_path(tmp_path, version_str):
    tensorleed_super_dir = tmp_path
    version_dir = tensorleed_super_dir / f'TensErLEED-v{version_str}'
    version_dir.mkdir()

    # check detection of directory
    detected_path = tenserleed.get_tensorleed_path(tmp_path)
    assert detected_path.is_dir()

    # check detection of version
    source = tenserleed.TensErLEEDSource(version_dir)
    assert source.version == Version(
        tenserleed.OLD_TL_VERSION_NAMES.get(version_str, version_str))
