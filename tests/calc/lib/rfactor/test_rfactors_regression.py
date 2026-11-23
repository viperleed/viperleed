"""Regression tests for R-factor calculations on real beam data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-04-07'
__license__ = 'GPLv3+'

import json

import numpy as np
import pytest


def _assert_regression_block(actual, expected):
    """Check that the R factors of actual and expected blocks agree."""
    assert actual['overall'] == pytest.approx(
        expected['overall'],
        rel=1e-11,
        abs=1e-13,
    )
    assert actual['per_beam'] == pytest.approx(
        expected['per_beam'],
        rel=1e-11,
        abs=1e-13,
    )
    assert actual['overall_from_splines'] == pytest.approx(
        expected['overall_from_splines'],
        rel=1e-11,
        abs=1e-13,
    )
    assert actual['per_beam_from_splines'] == pytest.approx(
        expected['per_beam_from_splines'],
        rel=1e-11,
        abs=1e-13,
    )


def test_realdata_regression_payload(request,
                                     regression_payload,
                                     expected_regression_file):
    """Check the R factors of real data against the stored values."""
    actual = regression_payload

    if request.config.getoption('--update-rfactor-expected'):
        text = json.dumps(actual, indent=2, sort_keys=True)
        expected_regression_file.write_text(text + '\n')
        pytest.skip(
            f'Updated expected regression data in {expected_regression_file}'
        )

    if not expected_regression_file.is_file():
        pytest.fail(
            'Missing expected regression file. Generate it with:\n'
            'pytest tests/calc/lib/rfactor/test_rfactors_regression.py '
            '--update-rfactor-expected'
        )

    expected = json.loads(expected_regression_file.read_text())

    assert actual['label'] == expected['label']
    assert actual['v0i'] == pytest.approx(expected['v0i'])
    assert actual['energy_step'] == pytest.approx(expected['energy_step'])
    assert actual['n_beams'] == expected['n_beams']
    assert actual['hk'] == expected['hk']
    assert actual['beam_correspondence'] == expected['beam_correspondence']

    for key in ('pendry', 'r1', 'r2', 'rs', 'zj'):
        _assert_regression_block(actual[key], expected[key])


def test_realdata_spline_and_precomputed_paths_match(regression_payload):
    """Check that splines and precomputed arrays give the same R."""
    payload = regression_payload

    for key in ('pendry', 'r1', 'r2', 'rs', 'zj'):
        assert payload[key]['overall_from_splines'] == pytest.approx(
            payload[key]['overall'],
            rel=1e-11,
            abs=1e-13,
        )
        assert payload[key]['per_beam_from_splines'] == pytest.approx(
            payload[key]['per_beam'],
            rel=1e-11,
            abs=1e-13,
        )


def test_realdata_outputs_are_finite(regression_payload):
    """Check that R is finite, with one value per beam."""
    payload = regression_payload

    for key in ('pendry', 'r1', 'r2', 'rs', 'zj'):
        assert np.isfinite(payload[key]['overall'])
        assert np.all(np.isfinite(payload[key]['per_beam']))
        assert len(payload[key]['per_beam']) == payload['n_beams']
