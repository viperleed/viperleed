"""Test configuration for viperleed.calc.bookkeeper.history.

Fixtures
--------
mock_path
    A fake path-like object.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-13'
__license__ = 'GPLv3+'

from pathlib import Path
from unittest.mock import MagicMock

from pytest_cases import fixture


@fixture
def mock_path():
    """Return a fake pathlib.Path."""
    return MagicMock(spec=Path)
