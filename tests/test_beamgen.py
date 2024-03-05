"""Tests for module viperleed.calc.files.beamgen.

Created on 2023-06-09

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

from pytest_cases import fixture, parametrize_with_cases

from viperleed.calc.files import beamgen
from viperleed.calc.lib import symmetry

from .helpers import CaseTag, exclude_tags
from .poscar_slabs import CasePOSCARSlabs


class TestBeamScatteringSubsets:
    """Collection of tests for sorting of beams into subsets."""

    def test_get_beam_scattering_subsets_integer_only(self, subtests):
        """Test correct identification for integer-order beams."""
        beam_indices_raw = [(0, 0), (1, 0), (2, 1), (1, 5)]
        expected_subset_classes = [(0, 0),]
        expected_reduced_indices = [(0, 0), (0, 0), (0, 0), (0, 0)]

        (classes,
         indices) = beamgen.get_beam_scattering_subsets(beam_indices_raw)
        with subtests.test('Integer beams - classes'):
            assert classes == expected_subset_classes
        with subtests.test('Integer beams - indices'):
            assert indices == expected_reduced_indices

    def test_get_beam_scattering_subsets_fractional(self, subtests):
        """Test correct identification for integer and fractional beams."""
        beam_indices_raw = [(0, 0), (1, 0), (1/2, 1/2), (1, 1/5)]
        expected_classes = [(0, 0), (0, 1/5), (1/2, 1/2)]  # By |k|
        expected_indices = [(0, 0), (0, 0), (1/2, 1/2), (0, 1/5)]

        (classes,
         indices) = beamgen.get_beam_scattering_subsets(beam_indices_raw)
        with subtests.test('Integer and fractional beams - classes'):
            assert classes == expected_classes
        with subtests.test('Integer and fractional beams - indices'):
            assert indices == expected_indices


_BEAMGEN_CASES = {
    'cases': CasePOSCARSlabs,
    'filter': exclude_tags(CaseTag.NO_INFO),                                     # TODO: not great to exclude these, but there's a few that really need PARAMETERS
    'scope': 'class',
    }


class TestGenerateBeamlist:
    """Collection of tests for the generation of beam lists."""                 # TODO: check energy sorting

    @fixture(name='make_beamlist', scope='class')
    @parametrize_with_cases('args', **_BEAMGEN_CASES)
    def fixture_make_beamlist(self, args, tmp_path_factory, tensorleed_path):
        """Return slab, parameters, info and the path to a 'BEAMLIST'."""
        slab, param, info = args
        slab.create_layers(param)
        slab.make_bulk_slab(param)
        symmetry.findSymmetry(slab, param)
        symmetry.findSymmetry(slab.bulkslab, param, bulk=True)
        symmetry.findBulkSymmetry(slab.bulkslab, param)

        param.initTheoEnergies()
        param.source_dir = tensorleed_path
        param.updateDerivedParams()  # for TL_VERSION

        tmp_dir = tmp_path_factory.mktemp(
            basename=f'{info.poscar.name}_beamlist',
            numbered=True
            )
        beamlist = tmp_dir / 'BEAMLIST'
        beamgen.calc_and_write_beamlist(slab, param, beamlist_name=beamlist)
        return slab, param, info, beamlist

    def test_generate_beamlist(self, make_beamlist):
        """Check successful generation of a 'BEAMLIST' file."""
        *_, beamlist = make_beamlist
        assert beamlist.exists()
        # check that file is not empty
        with open(beamlist, 'r') as file:
            assert file.read()

    def test_beamlist_is_ordered_consistently(self, make_beamlist):
        """Check that the beams are ordered consistently (See issue #184)"""
        beamlist_contents = []
        for _ in range(3):
            *_, beamlist = make_beamlist
            with open(beamlist, 'r') as file:
                beamlist_contents.append(file.read())
        assert all(contents == beamlist_contents[0]
                   for contents in beamlist_contents)
