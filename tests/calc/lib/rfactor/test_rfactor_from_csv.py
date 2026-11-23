"""Tests for the R factor from CSV files, unmocked, on real data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import numpy as np
import pytest

from viperleed.calc import from_ase
from viperleed.calc.lib.rfactor import r_1


@pytest.mark.parametrize('with_correspondence', [False, True])
def test_rfactor_from_csv(fe2o3_012_1x1_dataset, with_correspondence):
    """Check that the R factor of the Fe2O3(012) beam files is finite.

    This runs the Python R factor end to end, through CachedSpline,
    unlike the other tests in this directory.
    """
    kwargs = {}
    if with_correspondence:
        kwargs['beam_correspondence'] = (
            fe2o3_012_1x1_dataset['beam_correspondence']
        )
    best_r, best_v0r, r_beams = from_ase.rfactor_from_csv(
        (fe2o3_012_1x1_dataset['exp_file'],
         fe2o3_012_1x1_dataset['theo_file']),
        fe2o3_012_1x1_dataset['v0i'],
        **kwargs,
    )
    assert np.isfinite(best_r)
    assert np.isfinite(best_v0r)
    assert isinstance(r_beams, np.ndarray)
    assert r_beams.ndim == 1
    if with_correspondence:
        assert r_beams.shape == (len(fe2o3_012_1x1_dataset['hk']),)
    assert np.all(np.isfinite(r_beams))


def test_rfactor_from_csv_invalid_correspondence_length(
    fe2o3_012_1x1_dataset,
):
    """Check that a too-short beam_correspondence raises ValueError."""
    correspondence = fe2o3_012_1x1_dataset['beam_correspondence']
    with pytest.raises(ValueError, match='beam_correspondence'):
        from_ase.rfactor_from_csv(
            (fe2o3_012_1x1_dataset['exp_file'],
             fe2o3_012_1x1_dataset['theo_file']),
            fe2o3_012_1x1_dataset['v0i'],
            beam_correspondence=correspondence[:-1],
        )


def test_rfactor_from_csv_which_r_r1(fe2o3_012_1x1_dataset, mocker):
    """Check that which_r='r1' computes the R1 R factor."""
    r_1_func = mocker.patch(
        'viperleed.calc.lib.rfactor.r_1.r_1',
        wraps=r_1.r_1,
    )
    from_ase.rfactor_from_csv(
        (fe2o3_012_1x1_dataset['exp_file'],
         fe2o3_012_1x1_dataset['theo_file']),
        fe2o3_012_1x1_dataset['v0i'],
        which_r='r1',
    )
    r_1_func.assert_called()


def test_rfactor_from_csv_which_r_code_raises(fe2o3_012_1x1_dataset):
    """Check that an integer code is not an R-factor name."""
    with pytest.raises(ValueError, match='Unknown R-factor name'):
        from_ase.rfactor_from_csv(
            (fe2o3_012_1x1_dataset['exp_file'],
             fe2o3_012_1x1_dataset['theo_file']),
            fe2o3_012_1x1_dataset['v0i'],
            which_r='1',
        )


def _with_duplicated_beam(csv_file):
    """Return the CSV text of csv_file with its first beam repeated."""
    lines = csv_file.read_text().splitlines()
    return '\n'.join(
        f'{line},{line.split(",")[1]}' for line in lines
    )


@pytest.mark.parametrize('duplicated', [0, 1])
@pytest.mark.parametrize('with_correspondence', [False, True])
def test_rfactor_from_csv_duplicate_labels(
    fe2o3_012_1x1_dataset,
    duplicated,
    with_correspondence,
):
    """Check that repeated beam labels in a file raise ValueError."""
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


def test_rfactor_from_csv_accepts_lists(fe2o3_012_1x1_dataset):
    """Check that lists work in place of tuples."""
    kwargs = {
        'beam_correspondence': fe2o3_012_1x1_dataset['beam_correspondence'],
        'v0r_shift_range': (-2.0, 2.0),
    }
    files = (fe2o3_012_1x1_dataset['exp_file'],
             fe2o3_012_1x1_dataset['theo_file'])
    v0i = fe2o3_012_1x1_dataset['v0i']
    from_tuples = from_ase.rfactor_from_csv(files, v0i, **kwargs)
    kwargs['v0r_shift_range'] = list(kwargs['v0r_shift_range'])
    from_lists = from_ase.rfactor_from_csv(
        list(files),
        v0i,
        beams_file_is_content=[False, False],
        **kwargs,
    )
    assert from_lists[0] == from_tuples[0]
    assert from_lists[1] == from_tuples[1]
    assert np.array_equal(from_lists[2], from_tuples[2])
