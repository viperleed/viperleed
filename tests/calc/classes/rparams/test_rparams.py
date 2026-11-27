"""Tests for module rparams of viperleed.calc.classes.rparams."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-27'
__license__ = 'GPLv3+'

import pytest

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.files import displacements
from viperleed.calc.files import vibrocc
from viperleed.calc.files.iodeltas import checkDelta
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version


_MODULE = 'viperleed.calc.classes.rparams.rparams'


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

    def test_merges_offsets(self, ag100, data_path, tmp_path, mocker):
        """Check that searching for deltas merges VIBROCC offsets.

        generateSearchInput does not merge offsets itself: it relies
        on generateSearchPars, run before it also in search-only runs
        (RUN = 3), calling checkDelta, thus Atom.mergeDisp. See #602.
        """
        slab, rpars, *_ = ag100
        inputs_path = data_path / 'Ag(100)/mergeDisp'
        vibrocc.readVIBROCC(rpars, slab, str(inputs_path / 'VIBROCC'))
        rpars.fileLoaded['VIBROCC'] = True
        slab.full_update(rpars)
        displacements.readDISPLACEMENTS(
            rpars,
            str(inputs_path / 'DISPLACEMENTS_mixed'),
            )
        displacements.readDISPLACEMENTS_block(rpars,
                                              slab,
                                              rpars.disp_blocks[0])
        rpars.TL_VERSION = Version('2.0.0')
        atom = slab.atlist[0]
        offset = atom.offset_occ[atom.el]
        merged = [occ + offset for occ in atom.disp_occ[atom.el]]
        assert atom.oriState is None
        for delta_atom in slab:
            for element in (atom.el, 'vac'):
                delta = tmp_path/f'DEL_{delta_atom.num}_{element}_1'
                delta.write_text('', encoding='utf-8')
        # The real checkDelta merges offsets, then rejects the empty
        # files. Accept them, to get to the search parameters.
        mocker.patch(f'{_MODULE}.checkDelta',
                     side_effect=lambda *args: checkDelta(*args) or True)
        with execute_in_dir(tmp_path):
            rpars.generateSearchPars(slab)
        assert atom in rpars.search_atlist
        assert atom.disp_occ[atom.el] == pytest.approx(merged)
        for search_atom in rpars.search_atlist:
            assert search_atom.oriState is not None
            assert not search_atom.offset_occ
