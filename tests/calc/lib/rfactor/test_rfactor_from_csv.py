"""Tests for the R factor from CSV files, unmocked, on real data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import numpy as np

from viperleed.calc import from_ase


def test_rfactor_from_csv(fe2o3_012_1x1_dataset):
    """Check that the R factor of the Fe2O3(012) beam files is finite.

    This runs the Python R factor end to end, through CachedSpline,
    unlike the other tests in this directory.
    """
    best_r, best_v0r, r_beams = from_ase.rfactor_from_csv(
        (fe2o3_012_1x1_dataset['exp_file'],
         fe2o3_012_1x1_dataset['theo_file']),
        fe2o3_012_1x1_dataset['v0i'],
        beam_correspondence=fe2o3_012_1x1_dataset['beam_correspondence'],
        )
    assert np.isfinite(best_r)
    assert np.isfinite(best_v0r)
    assert r_beams.shape == (len(fe2o3_012_1x1_dataset['hk']),)
    assert np.all(np.isfinite(r_beams))
