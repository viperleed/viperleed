"""Tests for cleanup.get_rpars_from_manifest of viperleed.calc.section."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-16'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.sections.cleanup import get_rpars_from_manifest as get


def test_manifest(manifest):
    """Check the return value when a ManifestFile object is given."""
    rpars = get(manifest)
    assert rpars.manifest is manifest
    # pylint: disable-next=protected-access               # OK in tests
    assert rpars.timer._stopped


@parametrize(obj=(tuple(), [], '1', type))
def test_raises(obj):
    """Check complaints when the wrong object type is given."""
    with pytest.raises(TypeError):
        get(obj)


def test_rpars(rpars):
    """Check the return value when an Rparams object is given."""
    assert get(rpars) is rpars
