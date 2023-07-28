"""Test module beamgen of viperleed.tests.

Created on 2023-06-09

@author: Alexander M. Imre

"""

import pytest
from pathlib import Path
import os, sys
from copy import deepcopy
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

from viperleed.tleedmlib.beamgen import get_beam_scattering_subsets
from viperleed.tleedmlib.beamgen import make_beamlist_string
from viperleed.tleedmlib.beamgen import calc_and_write_beamlist

_FIXTURES_PATH = Path('tests/fixtures/')

class TestBeamScatteringSubsets:
    def test_get_beam_scattering_subsets_integer_only(self):
        beam_indices_raw = [(0, 0), (1, 0), (2, 1), (1, 5)]
        expected_subset_classes = [(0, 0),]
        expected_reduced_indices = [(0, 0), (0, 0), (0, 0), (0, 0)]

        subset_classes, reduced_indices = get_beam_scattering_subsets(beam_indices_raw)
        assert subset_classes == expected_subset_classes
        assert reduced_indices == expected_reduced_indices


    def test_get_beam_scattering_subsets_fractional(self):
        beam_indices_raw = [(0, 0), (1, 0), (1/2, 1/2), (1, 1/5)]
        expected_subset_classes = [(0, 0),(0, 1/5), (1/2, 1/2)] # sorted by |k| !
        expected_reduced_indices = [(0, 0), (0, 0), (1/2, 1/2), (0, 1/5)]

        subset_classes, reduced_indices = get_beam_scattering_subsets(beam_indices_raw)
        assert subset_classes == expected_subset_classes
        assert reduced_indices == expected_reduced_indices



def test_generate_beamlist(tmp_path, ag100_parameters_example):
    rp, sl = ag100_parameters_example
    beamlist_path = tmp_path / 'BEAMLIST'
    sl.createLayers(rp)
    calc_and_write_beamlist(sl, rp, beamlist_name = beamlist_path)
    assert beamlist_path.exists()