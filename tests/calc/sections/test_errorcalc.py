"""Tests for module errorcalc of viperleed.calc.sections."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.sections import errorcalc

_MODULE = 'viperleed.calc.sections.errorcalc'


class TestErrorsFigures:
    """Tests for the plots of the error calculation."""

    @staticmethod
    def _mock_rpars(mocker, r_factor_type):
        """Return fake Rparams with one atom varied once per mode."""
        atom = mocker.MagicMock()
        search_par = mocker.MagicMock(atom=atom, steps=2, el='vac')
        return mocker.MagicMock(
            domainParams=[],
            best_v0r=0.,
            search_index=0,
            disp_blocks=[None],
            STOP=False,
            halt=0,
            HALTING=2,
            search_atlist=[atom],
            searchpars=[search_par],
            R_FACTOR_TYPE=RFactorType(r_factor_type),
            )

    @parametrize(r_factor_type=('pendry', 'smooth'))
    def test_figures_get_r_factor_type(self, r_factor_type, mocker,
                                       tmp_path, monkeypatch):
        """Check that error plots are labelled with R_FACTOR_TYPE."""
        monkeypatch.chdir(tmp_path)
        rpars = self._mock_rpars(mocker, r_factor_type)
        for name in ('R_Error',
                     'section_deltas',
                     'section_superpos'):
            mocker.patch(f'{_MODULE}.{name}')
        mocker.patch(f'{_MODULE}.readDISPLACEMENTS_block', return_value=True)
        mocker.patch(f'{_MODULE}.section_rfactor', return_value=[0.1])
        io_mock = mocker.patch(f'{_MODULE}.ioerrorcalc')
        io_mock.generate_errors_csv.return_value = ('', {})

        errorcalc.errorcalc(mocker.MagicMock(), rpars)

        make_figs = io_mock.make_errors_figs
        make_figs.assert_called_once()
        assert (make_figs.call_args.kwargs['r_factor_type']
                is RFactorType(r_factor_type))
