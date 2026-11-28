"""Regression tests for R-factor calculations on real beam data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-04-07'
__license__ = 'GPLv3+'

import json

import numpy as np
import pytest

from .conftest import EXPECTED_FILE


_R_TYPES = ('pendry', 'r1', 'r2', 'rs', 'zj')
_R_KEYS = ('overall',
           'per_beam',
           'overall_from_splines',
           'per_beam_from_splines')


def _assert_regression_block(actual, expected):
    """Check that the R factors of actual and expected blocks agree."""
    for key in _R_KEYS:
        assert actual[key] == pytest.approx(expected[key],
                                            rel=1e-11,
                                            abs=1e-13)


def test_realdata_regression_payload(request, regression_payload):
    """Check the R factors of real data against the stored values."""
    if request.config.getoption('--update-rfactor-expected'):
        text = json.dumps(regression_payload, indent=2, sort_keys=True)
        EXPECTED_FILE.write_text(text + '\n')
        pytest.skip(f'Updated expected regression data in {EXPECTED_FILE}')

    if not EXPECTED_FILE.is_file():
        pytest.fail(
            'Missing expected regression file. Generate it with:\n'
            'pytest tests/calc/lib/rfactor/test_rfactors_regression.py '
            '--update-rfactor-expected'
            )

    expected = json.loads(EXPECTED_FILE.read_text())
    actual = regression_payload
    assert actual['label'] == expected['label']
    assert actual['v0i'] == pytest.approx(expected['v0i'])
    assert actual['energy_step'] == pytest.approx(expected['energy_step'])
    assert actual['n_beams'] == expected['n_beams']
    assert actual['hk'] == expected['hk']
    assert actual['beam_correspondence'] == expected['beam_correspondence']
    for r_type in _R_TYPES:
        _assert_regression_block(actual[r_type], expected[r_type])


def test_realdata_spline_and_precomputed_paths_match(regression_payload):
    """Check that splines and precomputed arrays give the same R."""
    for r_type in _R_TYPES:
        block = regression_payload[r_type]
        assert block['overall_from_splines'] == pytest.approx(
            block['overall'],
            rel=1e-11,
            abs=1e-13,
            )
        assert block['per_beam_from_splines'] == pytest.approx(
            block['per_beam'],
            rel=1e-11,
            abs=1e-13,
            )


def test_realdata_outputs_are_finite(regression_payload):
    """Check that R is finite, with one value per beam."""
    for r_type in _R_TYPES:
        block = regression_payload[r_type]
        assert np.isfinite(block['overall'])
        assert np.all(np.isfinite(block['per_beam']))
        assert len(block['per_beam']) == regression_payload['n_beams']
