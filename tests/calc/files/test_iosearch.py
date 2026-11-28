"""Tests for module iosearch of viperleed.calc.files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-27'
__license__ = 'GPLv3+'

import pytest

from viperleed.calc.files.iosearch import generateSearchInput
from viperleed.calc.lib.context import execute_in_dir


class _StopError(Exception):
    """Raised by a fake generateSearchPars to stop the search input."""


class TestGenerateSearchInput:
    """Tests for the generateSearchInput function."""

    def test_search_pars_first(self, tmp_path, mocker):
        """Check that search parameters come before writing anything.

        generateSearchInput relies on generateSearchPars to merge the
        VIBROCC and DISPLACEMENTS offsets of every search atom (via
        checkDelta, see #602), so it must call it before using them.
        """
        rpars = mocker.MagicMock()
        rpars.generateSearchPars.side_effect = _StopError
        slab = mocker.MagicMock()
        with execute_in_dir(tmp_path), pytest.raises(_StopError):
            generateSearchInput(slab, rpars)
        rpars.generateSearchPars.assert_called_once_with(slab)
        assert not any(tmp_path.iterdir())
