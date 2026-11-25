"""Tests for BEAM_EQUIVALENCE in the initialization of domain runs."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import shutil

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files.iosearch import writeRfInfo
from viperleed.calc.files.tenserleed import CURRENT_TL_VERSION
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.leedbase import getBeamCorrespondence
from viperleed.calc.run import run_calc
from viperleed.calc.sections.calc_section import BEAM_EQUIVALENCE_NAME
from viperleed.calc.sections.run_sections import section_loop

pytestmark = pytest.mark.tensorleed


def _log_text(work):
    """Return the concatenated calc logs found in work."""
    return '\n'.join(p.read_text() for p in work.glob('viperleed-calc*.log'))


@fixture(name='run_domains_init')
def factory_run_domains_init(mocker, data_path, make_section_tempdir,
                             tensorleed_path):
    """Return a function running calc on silver_and_bismuth."""
    loop = mocker.patch('viperleed.calc.run.section_loop',
                        wraps=section_loop)
    source = data_path / 'domains' / 'silver_and_bismuth' / 'initialization'

    def _run(tag, beam_equivalence=None, run=(0,), expect_fail=False):
        home = make_section_tempdir('silver_and_bismuth', tag)
        for path in source.iterdir():
            if path.is_file():
                shutil.copy(path, home)
            elif path.name != '_expect':
                shutil.copytree(path, home / path.name)
        for sub in ('silver', 'bismuth'):
            shutil.copy(source / '_expect' / 'work' / sub / 'PHASESHIFTS',
                        home / sub)
        if beam_equivalence is not None:
            (home / BEAM_EQUIVALENCE_NAME).write_text(beam_equivalence)
        work = home / DEFAULT_WORK
        shutil.copytree(home, work)
        with execute_in_dir(work):
            failed, _ = run_calc(
                system_name='silver_and_bismuth',
                source=tensorleed_path,
                preset_params={'RUN': list(run),
                               'TL_VERSION': CURRENT_TL_VERSION},
                home=home,
                )
        assert bool(failed) is expect_fail
        rpars, slab = loop.call_args.args
        return slab or rpars.pseudoSlab, rpars, work
    return _run


class TestDomainsBeamEquivalence:
    """Tests for init_domains with a BEAM_EQUIVALENCE file."""

    def test_roundtrip(self, run_domains_init):
        """Check that the file written by calc reproduces the result."""
        slab, rpars, work = run_domains_init('auto')
        vector = getBeamCorrespondence(slab, rpars)
        with execute_in_dir(work):
            rf_info = writeRfInfo(slab, rpars, work / 'rf.info')
        text = (work / DEFAULT_SUPP / BEAM_EQUIVALENCE_NAME).read_text()

        slab, rpars, work = run_domains_init('file', text)
        with execute_in_dir(work):
            rf_info_file = writeRfInfo(slab, rpars, work / 'rf.info')
        assert rpars.fileLoaded['BEAM_EQUIVALENCE']
        assert getBeamCorrespondence(slab, rpars) == vector
        assert rf_info_file == rf_info

    def test_file_filters_beams(self, run_domains_init):
        """Check that the file restricts the experimental beams."""
        text = '(1 0) = (1 0)\n(1 1) = (1 1)\n(9 9) = (1 0)\n'
        _, rpars, _ = run_domains_init('differs', text)
        assert len(rpars.expbeams) == 2

    def test_rhs_not_in_beamlist(self, run_domains_init):
        """Check the error for a calculated beam missing in BEAMLIST."""
        text = '(1 0) = (1 0) + (7 7)\n'
        *_, work = run_domains_init('no_beamlist', text, expect_fail=True)
        msg = 'in BEAM_EQUIVALENCE are not in BEAMLIST: ( 7| 7)'
        assert msg in _log_text(work)

    @parametrize(run=((0, 4), (0, 1)))
    def test_weights_error(self, run_domains_init, run):
        """Check that weights are rejected for runs needing a search."""
        text = '(1 0) = 1*(1 0) + 2*(1 1)\n(2 0) = (2 0)\n'
        _, rpars, work = run_domains_init(
            'weighted', text, run=run, expect_fail=True)
        assert 'BEAM_EQUIVALENCE contains' in _log_text(work)
        assert rpars.runHistory == [0]

    def test_all_lines_ignored(self, run_domains_init):
        """Check the error if no left-hand side is in EXPBEAMS."""
        text = '(9 9) = (1 0)\n(8 8) = (1 1)\n'
        *_, work = run_domains_init('all_ignored', text,
                                    expect_fail=True)
        assert 'No left-hand side of BEAM_EQUIVALENCE' in _log_text(work)
