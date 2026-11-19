"""Tests for the R factor from CSV files, unmocked, on real data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import numpy as np
import pytest

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


def _with_duplicated_beam(csv_file):
    """Return the CSV text of csv_file with its first beam repeated."""
    lines = csv_file.read_text().splitlines()
    return '\n'.join(
        f'{line},{line.split(",")[1]}' for line in lines
    )


@pytest.mark.parametrize('duplicated', (0, 1))
@pytest.mark.parametrize('with_correspondence', (False, True))
def test_rfactor_from_csv_duplicate_labels(
    fe2o3_012_1x1_dataset,
    duplicated,
    with_correspondence,
):
    """Check that repeated beam labels in either file raise ValueError."""
    files = [fe2o3_012_1x1_dataset['exp_file'],
             fe2o3_012_1x1_dataset['theo_file']]
    contents = [f.read_text() for f in files]
    contents[duplicated] = _with_duplicated_beam(files[duplicated])
    correspondence = None
    if with_correspondence:
        correspondence = fe2o3_012_1x1_dataset['beam_correspondence']
        if duplicated:
            correspondence = (*correspondence, 0)
    with pytest.raises(ValueError, match='same index'):
        from_ase.rfactor_from_csv(
            contents,
            fe2o3_012_1x1_dataset['v0i'],
            beams_file_is_content=(True, True),
            beam_correspondence=correspondence,
        )
