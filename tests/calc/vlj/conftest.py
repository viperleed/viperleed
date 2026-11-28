"""Fixtures for tests of viperleed.calc.vlj."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

import pytest

from .dnl_scope import scoped_dnl_patch


@pytest.fixture(autouse=True)
def fixture_scoped_dnl_patch():
    """Keep the plugin's swap of the numerical backend within one test."""
    with scoped_dnl_patch():
        yield
