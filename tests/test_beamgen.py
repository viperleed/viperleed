"""Tests for module viperleed.tleedmlib.beamgen.

Created on 2023-06-09

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""
from pathlib import Path
import sys

from pytest_cases import fixture, parametrize_with_cases

VPR_PATH = str(Path(__file__).resolve().parents[2])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable
from viperleed.tleedmlib import beamgen, symmetry
from viperleed.tleedmlib.classes.rparams import EnergyRange

from .helpers import CaseTag, exclude_tags
from .poscar_slabs import CasePOSCARSlabs
# pylint: enable=wrong-import-position


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
            assert list(classes) == expected_subset_classes
        with subtests.test('Integer beams - indices'):
            assert indices == expected_reduced_indices

    def test_get_beam_scattering_subsets_fractional(self, subtests):
        """Test correct identification for integer and fractional beams."""
        beam_indices_raw = [(0, 0), (1, 0), (1/2, 1/2), (1, 1/5)]
        expected_classes = [(0, 0), (1/2, 1/2), (0, 1/5)]  # Same order
        expected_indices = [(0, 0), (0, 0), (1/2, 1/2), (0, 1/5)]

        (classes,
         indices) = beamgen.get_beam_scattering_subsets(beam_indices_raw)
        with subtests.test('Integer and fractional beams - classes'):
            assert list(classes) == expected_classes
        with subtests.test('Integer and fractional beams - indices'):
            assert indices == expected_indices


_BEAMGEN_CASES = {
    'cases': CasePOSCARSlabs,
    'filter': exclude_tags(CaseTag.NO_INFO),                                    # TODO: not great to exclude these, but there's a few that really need PARAMETERS
    'scope': 'class',
    }


class TestGenerateBeamlist:
    """Collection of tests for the generation of beam lists."""                 # TODO: check energy sorting

    @fixture(name='prepare_for_beamlist', scope='class')
    @parametrize_with_cases('args', **_BEAMGEN_CASES)
    def fixture_prepare_for_beamlist(self, args, tensorleed_path):
        """Return slab, parameters, info and the path to a 'BEAMLIST'."""
        slab, param, info = args
        slab.create_layers(param)
        slab.make_bulk_slab(param)
        symmetry.findSymmetry(slab, param)
        symmetry.findSymmetry(slab.bulkslab, param, bulk=True)
        symmetry.findBulkSymmetry(slab.bulkslab, param)

        # Limit energy range for performance reasons. 300 eV is more
        # than enough to generate a sizeable number of beams for any
        # reasonable system
        param.THEO_ENERGIES = EnergyRange(stop=300)
        param.initTheoEnergies()
        param.source_dir = tensorleed_path
        param.updateDerivedParams()  # for TL_VERSION
        return slab, param, info

    @fixture(name='make_beamlist', scope='class')
    def fixture_make_beamlist(self, prepare_for_beamlist, tmp_path_factory):
        """Create a 'BEAMLIST' file and return slab, parameters, info and path."""
        slab, param, info = prepare_for_beamlist
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

    def test_beamlist_is_correct(self, make_beamlist, data_path):
        """Check that the BEAMLIST file matches the expected contents."""
        *_, info, beamlist = make_beamlist
        with open(beamlist, 'r') as file:
            contents = file.read()
        # get the expected contents
        beamlist_name = info.poscar.name.replace('POSCAR', 'BEAMLIST')
        with open(data_path / 'BEAMLISTs' / beamlist_name, 'r') as file:
            expected_contents = file.read()
        assert contents == expected_contents
