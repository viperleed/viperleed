"""Shared fixtures and helpers for rfactor tests."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files.iorfactor import beamlist_to_array
from viperleed.calc.lib.rfactor.utils import average_beam_array
from viperleed.calc.lib.spline_interpolation import interpolate_ragged_array


_MOCK_ROUGH_E_AXIS = np.arange(15.0, 300.0, 3.0)
_MOCK_FINE_E_STEP = 0.5
_MOCK_FINE_E_AXIS = np.arange(15.0, 300.0, _MOCK_FINE_E_STEP)

_MOCK_SIN_INTENSITY = np.sin(_MOCK_ROUGH_E_AXIS / (3.0 * np.pi)) ** 2
_MOCK_COS_INTENSITY = np.cos(_MOCK_ROUGH_E_AXIS / (3.0 * np.pi)) ** 2
_MOCK_CONST_INTENSITY = np.ones_like(_MOCK_ROUGH_E_AXIS)

TESTS_DIR = Path(__file__).resolve().parents[3]
BEAM_DATA_DIR = TESTS_DIR / "_test_data" / "beam-data"
EXP_BEAMS_FILE = BEAM_DATA_DIR / "EXPBEAMS-fe2o3-012-1x1.csv"
THEO_BEAMS_FILE = BEAM_DATA_DIR / "THEOBEAMS-fe2o3-012-1x1.csv"
BEAMCORR_FILE = BEAM_DATA_DIR / "beamcorr-fe2o3-012-1x1.json"
EXPECTED_FILE = BEAM_DATA_DIR / "rfactor-fe2o3-012-1x1.expected.json"


def pytest_addoption(parser):
    parser.addoption(
        "--update-rfactor-expected",
        action="store_true",
        default=False,
        help="Regenerate frozen expected values for real-data rfactor regression tests.",
    )


@pytest.fixture
def sin_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_SIN_INTENSITY)


@pytest.fixture
def cos_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_COS_INTENSITY)


@pytest.fixture
def const_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_CONST_INTENSITY)


@pytest.fixture
def synthetic_pair(sin_spline, cos_spline):
    return {
        "v0i": 3.0,
        "energy_step": _MOCK_FINE_E_STEP,
        "energies": _MOCK_FINE_E_AXIS,
        "data_spline_1": sin_spline,
        "data_spline_2": cos_spline,
    }


def _evaluate_spline_triplet(data_spline, energy_grid):
    deriv_1 = data_spline.derivative()
    deriv_2 = deriv_1.derivative()
    return (
        data_spline(energy_grid),
        deriv_1(energy_grid),
        deriv_2(energy_grid),
    )


def _jsonify_hk(hk_values):
    return [list(hk) for hk in hk_values]


def _load_beam_correspondence():
    data = json.loads(BEAMCORR_FILE.read_text())
    return tuple(data["beam_correspondence"])


@pytest.fixture(scope="session")
def fe2o3_012_1x1_dataset():
    """Real beam data with fixed beam correspondence applied to THEO beams."""
    exp_out_beams = readOUTBEAMS(EXP_BEAMS_FILE)
    theo_out_beams = readOUTBEAMS(THEO_BEAMS_FILE)

    exp_grid, _, _, exp_beams = beamlist_to_array(exp_out_beams)
    theo_grid, _, _, theo_beams = beamlist_to_array(theo_out_beams)

    exp_hk = [beam.hk for beam in exp_out_beams]
    theo_hk = [beam.hk for beam in theo_out_beams]

    beam_corr = _load_beam_correspondence()

    if len(beam_corr) != len(exp_hk):
        raise ValueError(
            f"Beam correspondence length {len(beam_corr)} does not match "
            f"number of experimental beams {len(exp_hk)}"
        )
    if len(beam_corr) != len(theo_hk):
        raise ValueError(
            f"Beam correspondence length {len(beam_corr)} does not match "
            f"number of theoretical beams {len(theo_hk)}"
        )

    theo_beams_effective = average_beam_array(theo_beams, beam_corr)
    assert theo_beams_effective.shape[1] == len(exp_hk)

    exp_spline = interpolate_ragged_array(exp_grid, exp_beams)
    theo_spline = interpolate_ragged_array(theo_grid, theo_beams_effective)

    emin = max(np.nanmin(exp_grid), np.nanmin(theo_grid))
    emax = min(np.nanmax(exp_grid), np.nanmax(theo_grid))
    energy_step = 0.5
    energies = np.arange(emin, emax + 0.25 * energy_step, energy_step)

    exp_triplet = _evaluate_spline_triplet(exp_spline, energies)
    theo_triplet = _evaluate_spline_triplet(theo_spline, energies)

    return {
        "label": "fe2o3-012-1x1",
        "v0i": 3.0,
        "energy_step": energy_step,
        "energies": energies,
        "hk": exp_hk,
        "beam_correspondence": beam_corr,
        "exp_spline": exp_spline,
        "theo_spline": theo_spline,
        "exp_triplet": exp_triplet,
        "theo_triplet": theo_triplet,
    }


def build_regression_payload(dataset, rfactor_module):
    funcs = {
        "pendry": rfactor_module.R_pendry,
        "r1": rfactor_module.R_1,
        "r2": rfactor_module.R_2,
        "rs": rfactor_module.R_s,
        "zj": rfactor_module.R_zj,
    }

    payload = {
        "label": dataset["label"],
        "v0i": dataset["v0i"],
        "energy_step": dataset["energy_step"],
        "n_beams": len(dataset["hk"]),
        "hk": _jsonify_hk(dataset["hk"]),
        "beam_correspondence": list(dataset["beam_correspondence"]),
    }

    base_kwargs_spline = {
        "v0_imag": dataset["v0i"],
        "energy_step": dataset["energy_step"],
        "energy_grid": dataset["energies"],
        "data_spline_1": dataset["exp_spline"],
        "data_spline_2": dataset["theo_spline"],
    }
    base_kwargs_arrays = {
        "v0_imag": dataset["v0i"],
        "energy_step": dataset["energy_step"],
        "energy_grid": dataset["energies"],
        "data_and_derivatives_1": dataset["exp_triplet"],
        "data_and_derivatives_2": dataset["theo_triplet"],
    }

    for name, func in funcs.items():
        overall = func(**base_kwargs_arrays)
        per_beam = func(**base_kwargs_arrays, groups="beam")
        overall_from_splines = func(**base_kwargs_spline)
        per_beam_from_splines = func(**base_kwargs_spline, groups="beam")

        payload[name] = {
            "overall": float(np.asarray(overall)),
            "per_beam": np.asarray(per_beam, dtype=float).tolist(),
            "overall_from_splines": float(np.asarray(overall_from_splines)),
            "per_beam_from_splines": np.asarray(
                per_beam_from_splines, dtype=float
            ).tolist(),
        }

    return payload


def load_expected_payload():
    return json.loads(EXPECTED_FILE.read_text())


def write_expected_payload(payload):
    EXPECTED_FILE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
