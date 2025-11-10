"""Tests for viperleed.calc section refcalc.

This module contains tests for an actual execution of a reference
calculation for a multi-domain setup.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-10-22'
__license__ = 'GPLv3+'

import pytest

from ..conftest import _NON_INIT_TL_VERSION
from ..conftest import DomainsCalcFilesSetup as Setup


def test_explicit_refcalc_executed(tensorleed_path, tmp_path):
    """Ensure that an explicitly requested refcalc is always executed."""
    refcalc = Setup(
        surface_dir='telluride_stacking',
        tmp_test_path=tmp_path,
        copy_dirs=['deltas'],
        )
    tensors_before = refcalc.find_most_recent_tensors()
    refcalc.run_calc_from_setup(
        source=tensorleed_path,
        preset_params={'RUN': [0, 1],  # init, and explicit refcalc
                       'TL_VERSION': _NON_INIT_TL_VERSION,}
        )
    tensors_after = refcalc.find_most_recent_tensors()
    assert not refcalc.failed
    for src, before in tensors_before.items():
        # All domains should now have a new tensor
        assert tensors_after[src] == before + 1
