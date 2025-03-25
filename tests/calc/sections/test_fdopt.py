"""Tests for module fd_optimization of viperleed.calc.sections."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-25'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.fd_optimization import get_fd_r

_MODULE = 'viperleed.calc.sections.fd_optimization'


class TestGetFdR:
    """Tests for the get_fd_r function."""

    @fixture(name='mock_args')
    def fixture_mock_args(self, mocker, tmp_path):
        """Return arguments for get_fd_r."""
        return (mocker.MagicMock(name='slab'),
                Rparams(),
                tmp_path/'fd_op_partial')

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        return {
            'refcalc': mocker.patch(f'{_MODULE}.section_refcalc'),
            'r_factor': mocker.patch(f'{_MODULE}.section_rfactor',
                                     return_value=mocker.MagicMock()),
            'copy': mocker.patch('shutil.copy2'),
            }

    def test_success(self, mock_args, mocks, tmp_path, caplog):
        """Check results of a successful call."""
        slab, rpars, workdir = mock_args
        with execute_in_dir(tmp_path):
            rfactor, rfactor_per_beam = get_fd_r(*mock_args)
        assert rfactor is rpars.last_R
        assert rfactor_per_beam is mocks['r_factor'].return_value
        assert not any(rpars.TENSOR_OUTPUT)
        assert workdir.is_dir()

        mocks['refcalc'].assert_called_once_with(slab,
                                                 rpars,
                                                 parent_dir=tmp_path)
        mocks['r_factor'].assert_called_once_with(slab, rpars, 11)
        mocks['copy'].assert_any_call(tmp_path/'BEAMLIST', workdir)
        mocks['copy'].assert_any_call(tmp_path/'PHASESHIFTS', workdir)
        assert not caplog.text

    def test_negative_theta(self, mock_args, mocks, tmp_path, caplog):
        """Check successful call with a negative THETA."""
        _, rpars, _ = mock_args
        rpars.THETA = -30
        rpars.PHI = 50
        self.test_success(mock_args, mocks, tmp_path, caplog)
        expect_beamincidence = (30, 230)
        assert (rpars.THETA, rpars.PHI) == expect_beamincidence

    _raises = {
        'copy': OSError,
        'refcalc': Exception,
        'r_factor': Exception,
        }

    @parametrize('callee,exc', _raises.items(), ids=_raises)
    # pylint: disable-next=too-many-arguments  # 4/7 are fixtures
    def test_raises(self, callee, exc, mock_args, mocks, tmp_path, caplog):
        """Check complaints when failing to execute a `callee`."""
        mocks[callee].side_effect = exc
        with pytest.raises(exc):
            with execute_in_dir(tmp_path):
                get_fd_r(*mock_args)
        assert caplog.text
