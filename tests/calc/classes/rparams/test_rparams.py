"""Tests for module rparams of viperleed.calc.classes.rparams."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-27'
__license__ = 'GPLv3+'

import logging

import numpy as np
import pytest
from pytest_cases import fixture

from viperleed.calc.classes.rparams.domain_params import DomainParameters
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files import displacements
from viperleed.calc.files import vibrocc
from viperleed.calc.files.iodeltas import checkDelta
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version


_MODULE = 'viperleed.calc.classes.rparams.rparams'
_DOMAIN_OFFSETS = {  # VIBROCC search offsets, more than OCC 1
    'first': '\nGEO 1 = Ag 0 0 0.05\nVIB 1 = Ag 0.04\n',
    'second': '\nGEO 2 = Ag 0.01 0 0\nVIB 2 = Ag 0.01\n',
    }
_DOMAIN_GEO_OFFSETS = {  # DISPLACEMENTS geo offsets
    'first': 'Ag 2 z offset = 0.03',
    }


class TestGenerateSearchPars:
    """Tests for the Rparams.generateSearchPars method."""

    def test_header_from_other_version(self, tmp_path, mocker):
        """Check that a search refuses a delta of another version."""
        rpars = Rparams()
        rpars.TL_VERSION = Version('1.6.0')
        delta = tmp_path/'DEL_1_Ag'
        delta.write_text('title\n    8   12    5\n', encoding='utf-8')
        atom = mocker.MagicMock(is_bulk=False, num=1,
                                disp_occ={'Ag': [1.0]})
        slab = mocker.MagicMock(deltas_initialized=False)
        slab.__iter__.return_value = [atom]
        with execute_in_dir(tmp_path), pytest.raises(RuntimeError):
            rpars.generateSearchPars(slab)

    @fixture(name='prepare_search')
    def factory_prepare_search(self, data_path):
        """Return a function that gets a slab ready for a search.

        The slab has VIBROCC offsets (Ag(100)/mergeDisp/VIBROCC plus
        more_offsets lines), DISPLACEMENTS (DISPLACEMENTS_mixed, plus
        a geo_offset line in its GEO_DELTA block) and fake delta
        files, as before a search-only run (RUN = 3).
        """
        inputs_path = data_path / 'Ag(100)/mergeDisp'

        def _prepare(slab, rpars, workdir, more_offsets='', geo_offset=''):
            vibrocc_path = workdir / 'VIBROCC'
            vibrocc_path.write_text(
                (inputs_path / 'VIBROCC').read_text(encoding='utf-8')
                + more_offsets,
                encoding='utf-8',
                )
            vibrocc.readVIBROCC(rpars, slab, str(vibrocc_path))
            rpars.fileLoaded['VIBROCC'] = True
            slab.full_update(rpars)
            geo_block = '= GEO_DELTA\n'
            displacements_path = workdir / 'DISPLACEMENTS'
            displacements_path.write_text(
                (inputs_path / 'DISPLACEMENTS_mixed')
                .read_text(encoding='utf-8')
                .replace(geo_block, f'{geo_block}        {geo_offset}\n'),
                encoding='utf-8',
                )
            displacements.readDISPLACEMENTS(rpars, str(displacements_path))
            displacements.readDISPLACEMENTS_block(rpars,
                                                  slab,
                                                  rpars.disp_blocks[0])
            rpars.TL_VERSION = Version('2.0.0')
            for atom in slab:
                for element in (atom.el, 'vac'):
                    delta = workdir/f'DEL_{atom.num}_{element}_1'
                    delta.write_text('', encoding='utf-8')
        return _prepare

    @fixture(name='accept_deltas')
    def fixture_accept_deltas(self, mocker):
        """Make checkDelta accept any delta file, after merging offsets.

        The real checkDelta merges offsets, then rejects the empty
        delta files. Accept them, to get to the search parameters.
        """
        mocker.patch(f'{_MODULE}.checkDelta',
                     side_effect=lambda *args: checkDelta(*args) or True)

    @pytest.mark.usefixtures('accept_deltas')
    def test_merges_offsets(self, ag100, prepare_search, tmp_path):
        """Check that searching for deltas merges VIBROCC offsets.

        generateSearchInput does not merge offsets itself: it relies
        on generateSearchPars, run before it also in search-only runs
        (RUN = 3), calling checkDelta, thus Atom.mergeDisp. See #602.
        """
        slab, rpars, *_ = ag100
        prepare_search(slab, rpars, tmp_path)
        atom = slab.atlist[0]
        offset = atom.offset_occ[atom.el]
        merged = [occ + offset for occ in atom.disp_occ[atom.el]]
        assert atom.oriState is None
        with execute_in_dir(tmp_path):
            rpars.generateSearchPars(slab)
        assert atom in rpars.search_atlist
        assert atom.disp_occ[atom.el] == pytest.approx(merged)
        for search_atom in rpars.search_atlist:
            assert search_atom.oriState is not None
            assert not search_atom.offset_occ

    @pytest.mark.usefixtures('accept_deltas')
    def test_merges_offsets_all_domains(self, ag100, make_poscar,
                                        prepare_search, tmp_path, caplog):
        """Check that a multi-domain search merges offsets everywhere.

        Like for a single domain, generateSearchInput relies on this
        for every search atom of every domain. See #602.
        """
        *_, info = ag100
        rpars = Rparams()
        for name, offsets in _DOMAIN_OFFSETS.items():
            workdir = tmp_path / name
            workdir.mkdir()
            domain = DomainParameters(workdir, name)
            domain.slab, domain.rpars, *_ = make_poscar(info)
            prepare_search(domain.slab, domain.rpars, workdir,
                           offsets, _DOMAIN_GEO_OFFSETS.get(name, ''))
            assert sum(bool(atom.offset_geo and atom.offset_vib)
                       for atom in domain.slab) == 1
            rpars.domainParams.append(domain)
        assert sum(np.any(atom.disp_geo_offset['all'][0])
                   for domain in rpars.domainParams
                   for atom in domain.slab) == 1
        rpars.generateSearchPars(None)
        assert not [r for r in caplog.records if r.levelno >= logging.ERROR]
        for domain in rpars.domainParams:
            assert domain.rpars.search_atlist
            for atom in domain.rpars.search_atlist:
                assert atom.oriState is not None
                assert not atom.offset_geo
                assert not atom.offset_vib
                assert not atom.offset_occ
                assert list(atom.disp_geo_offset) == ['all']
                assert np.allclose(atom.disp_geo_offset['all'], 0)
