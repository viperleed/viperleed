"""Regression tests for R-factor calculations on real beam data."""

import numpy as np
import pytest

from viperleed.calc.lib import rfactor

from .conftest import (
    EXPECTED_FILE,
    build_regression_payload,
    load_expected_payload,
    write_expected_payload,
)


def _assert_regression_block(actual, expected):
    assert actual["overall"] == pytest.approx(expected["overall"], rel=1e-11, abs=1e-13)
    assert actual["per_beam"] == pytest.approx(expected["per_beam"], rel=1e-11, abs=1e-13)
    assert actual["overall_from_splines"] == pytest.approx(
        expected["overall_from_splines"], rel=1e-11, abs=1e-13
    )
    assert actual["per_beam_from_splines"] == pytest.approx(
        expected["per_beam_from_splines"], rel=1e-11, abs=1e-13
    )


def test_realdata_regression_payload(request, fe2o3_012_1x1_dataset):
    actual = build_regression_payload(fe2o3_012_1x1_dataset, rfactor)

    if request.config.getoption("--update-rfactor-expected"):
        write_expected_payload(actual)
        pytest.skip(f"Updated expected regression data in {EXPECTED_FILE}")

    if not EXPECTED_FILE.is_file():
        pytest.fail(
            "Missing expected regression file. Generate it with:\n"
            "pytest tests/calc/lib/rfactor/test_rfactors_regression.py "
            "--update-rfactor-expected"
        )

    expected = load_expected_payload()

    assert actual["label"] == expected["label"]
    assert actual["v0i"] == pytest.approx(expected["v0i"])
    assert actual["energy_step"] == pytest.approx(expected["energy_step"])
    assert actual["n_beams"] == expected["n_beams"]
    assert actual["hk"] == expected["hk"]
    assert actual["beam_correspondence"] == expected["beam_correspondence"]

    for key in ("pendry", "r1", "r2", "rs", "zj"):
        _assert_regression_block(actual[key], expected[key])


def test_realdata_spline_and_precomputed_paths_match(fe2o3_012_1x1_dataset):
    payload = build_regression_payload(fe2o3_012_1x1_dataset, rfactor)

    for key in ("pendry", "r1", "r2", "rs", "zj"):
        assert payload[key]["overall_from_splines"] == pytest.approx(
            payload[key]["overall"], rel=1e-11, abs=1e-13
        )
        assert payload[key]["per_beam_from_splines"] == pytest.approx(
            payload[key]["per_beam"], rel=1e-11, abs=1e-13
        )


def test_realdata_outputs_are_finite(fe2o3_012_1x1_dataset):
    payload = build_regression_payload(fe2o3_012_1x1_dataset, rfactor)

    for key in ("pendry", "r1", "r2", "rs", "zj"):
        assert np.isfinite(payload[key]["overall"])
        assert np.all(np.isfinite(payload[key]["per_beam"]))
        assert len(payload[key]["per_beam"]) == payload["n_beams"]
