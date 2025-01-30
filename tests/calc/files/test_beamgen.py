"""Tests for module viperleed.calc.files.beamgen."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-06-09'
__license__ = 'GPLv3+'

from contextlib import nullcontext

import pytest
from pytest_cases import fixture, parametrize_with_cases

from viperleed.calc import symmetry
from viperleed.calc.classes.rparams import EnergyRange
from viperleed.calc.files import beamgen
from viperleed.calc.lib.version import Version

from ...helpers import exclude_tags
from ..poscar_slabs import CasePOSCARSlabs
from ..tags import CaseTag


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
        slab, rpars, info = args
        slab.create_layers(rpars)
        slab.make_bulk_slab(rpars)
        symmetry.findSymmetry(slab, rpars)
        symmetry.findSymmetry(slab.bulkslab, rpars, bulk=True)
        symmetry.findBulkSymmetry(slab.bulkslab, rpars)

        # Limit energy range for performance reasons. 300 eV is more
        # than enough to generate a sizeable number of beams for any
        # reasonable system
        rpars.THEO_ENERGIES = EnergyRange(stop=300)
        rpars.initTheoEnergies()
        rpars.paths.source = tensorleed_path
        rpars.updateDerivedParams()  # for TL_VERSION
        return slab, rpars, info

    @fixture(name='get_expected_beamlist', scope='class')
    def fixture_get_expected_beamlist(self, data_path):
        """Return a path to the test-data file for the right beam list."""
        # Find a mapping between subfolders and the
        # range of applicable TensErLEED versions
        folder_to_range = {}
        for folder in (data_path / 'BEAMLISTs').iterdir():
            version_min_max = [Version(v) for v in folder.name.split('-') if v]
            # pylint: disable-next=magic-value-comparison
            if len(version_min_max) < 2:
                version_min_max.append(None)
            folder_to_range[folder] = tuple(version_min_max)

        def _get_tl_version_path(rpars):
            """Return which of the 'BEAMLISTs' subfolders applies."""
            version = rpars.TL_VERSION
            for folder, (v_min, v_max) in folder_to_range.items():
                if version < v_min:
                    continue
                if v_max and version > v_max:
                    continue
                return folder
            raise ValueError('Found no subfolder of BEAMLISTs with '
                             f'data for TensErLEED v{version}')

        def _get_path_to_beamlist_file(rpars, info):
            beamlist_name = info.poscar.name.replace('POSCAR', 'BEAMLIST')
            folder = _get_tl_version_path(rpars)
            beamlist_file = folder / beamlist_name
            if not beamlist_file.is_file():
                raise FileNotFoundError(f'Found no {beamlist_name} '
                                        f'file in {folder}')
            return beamlist_file
        return _get_path_to_beamlist_file

    @fixture(name='make_beamlist', scope='class')
    def fixture_make_beamlist(self, prepare_for_beamlist,
                              get_expected_beamlist,
                              tmp_path_factory):
        """Create a 'BEAMLIST', and return slab, parameters, info, and path."""
        slab, rpars, info = prepare_for_beamlist
        try:
            expected = get_expected_beamlist(rpars, info)
        except FileNotFoundError:
            expected = None
            context = pytest.raises(ValueError, match='Too many beams')
        else:
            context = nullcontext()

        tmp_dir = tmp_path_factory.mktemp(
            basename=f'{info.poscar.name}_beamlist',
            numbered=True,
            )
        beamlist = tmp_dir / 'BEAMLIST'
        with context as exc_info:
            beamgen.calc_and_write_beamlist(slab,
                                            rpars,
                                            beamlist_name=beamlist)
        if expected is None:
            # Skip all tests that use this fixture if there was
            # an understandable exception. This typically means
            # that the structure is unsupported in the TensErLEED
            # version under test.
            pytest.skip(str(exc_info.value))
        return slab, rpars, info, expected, beamlist

    def test_generate_beamlist(self, make_beamlist):
        """Check successful generation of a 'BEAMLIST' file."""
        *_, beamlist = make_beamlist
        assert beamlist.exists()
        # Check that file is not empty
        assert beamlist.read_text(encoding='utf-8').strip()

    def test_beamlist_is_correct(self, make_beamlist):
        """Check that the BEAMLIST file matches the expected contents."""
        *_, path_to_expected, beamlist = make_beamlist
        contents = beamlist.read_text(encoding='utf-8')
        expected = path_to_expected.read_text(encoding='utf-8')
        assert contents == expected
