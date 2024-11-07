"""Tests for module fortran_utils of viperleed.calc.lib."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'

import shutil
import subprocess

import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.fortran_utils import wrap_fortran_line
from viperleed.calc.lib.fortran_utils import get_mpifort_version
from viperleed.calc.lib.fortran_utils import MpifortNotFoundError
from viperleed.calc.lib.fortran_utils import (
    CouldNotDeterminMpifortVersionError
    )

from ...helpers import not_raises


with_mpifort = pytest.mark.skipif(not shutil.which('mpifort'),
                                  reason='no mpifort installed')


class TestWrapFortranLine:  # pylint: disable=too-few-public-methods
    """Tests for the wrap_fortran_line function."""

    _valid = {
        'at limit': ('x' * 72, 'x' * 72),
        'empty': ('', ''),
        'non ASCII': (
            'é' * 156,
            'é' * 72 + '&\n     &' + 'é' * 66 + '&\n     &' + 'é' * 18,
            ),
        'short': ('12345', '12345'),
        'spaces': ('123456' + ' ' * 76,
                   '123456' + ' ' * 66 + '&\n     &' + ' ' * 10),
        'wrap once': ('x' * 78, 'x' * 72 + '&\n     &' + 'x' * 6),
        'wrap twice exactly': (
            'x' * 72 + 'y' * 66,
            'x' * 72 + '&\n     &' + 'y' * 66,
            ),
        'wrap multiple': (
            'x' * 72 + 'y' * 66 + 'z' * 12,
            'x' * 72 + '&\n     &' + 'y' * 66 + '&\n     &' + 'z' * 12,
            ),
        'code line': (
            '      DATA NLMS(1), NLMS(2), NLMS(3), NLMS(4), NLMS(5), NLMS(6), '
            'NLMS(7), NLMS(8), NLMS(9), NLMS(10), NLMS(11), NLMS(12), '
            'NLMS(13), NLMS(14), NLMS(15), NLMS(16), NLMS(17), NLMS(18)/ 0, '
            '76, 284, 809, 1925, 4032, 7680, 13593, 22693, 36124, 55276, '
            '81809, 117677, 165152, 226848, 305745, 405213, 529036/',
            '''\
      DATA NLMS(1), NLMS(2), NLMS(3), NLMS(4), NLMS(5), NLMS(6), NLMS(7)&
     &, NLMS(8), NLMS(9), NLMS(10), NLMS(11), NLMS(12), NLMS(13), NLMS(1&
     &4), NLMS(15), NLMS(16), NLMS(17), NLMS(18)/ 0, 76, 284, 809, 1925,&
     & 4032, 7680, 13593, 22693, 36124, 55276, 81809, 117677, 165152, 22&
     &6848, 305745, 405213, 529036/'''
            )
        }

    @parametrize('string,expect', _valid.values(), ids=_valid)
    def test_valid(self, string, expect):
        """Check expected outcome with valid arguments."""
        assert wrap_fortran_line(string) == expect


class TestGetMpifortVersion:
    """Tests for the get_mpifort_version function."""

    mock_version = '13.2.0'
    check_exists_cmd = 'mpifort --version'
    pull_version_cmd = 'mpifort --version | sed'

    @staticmethod
    def mock_run_fails(cmd, **_):
        """Simulate a failed subprocess.run."""
        raise subprocess.CalledProcessError(1, cmd)

    def mock_run_success(self, cmd, **_):
        """Simulate a successful subprocess.run."""
        cmd_str = ' '.join(cmd)
        if self.pull_version_cmd in cmd_str:
            # Simulate the expected output for version extraction
            result = f'{self.mock_version}\n'.encode()
            return subprocess.CompletedProcess(args=cmd,
                                               returncode=0,
                                               stdout=result)
        return self.mock_mpifort_exists_ok(cmd)

    def mock_mpifort_exists_ok(self, cmd, **_):
        """Simulate a successful check for the existence of mpifort."""
        if self.check_exists_cmd in cmd:
            # Simulate the initial check call
            return subprocess.CompletedProcess(args=cmd, returncode=0)
        raise ValueError(f'Unexpected command {cmd}')

    def test_success(self, monkeypatch):
        """Test successful retrieval of mpifort version."""
        monkeypatch.setattr(subprocess, 'run', self.mock_run_success)
        assert get_mpifort_version() == self.mock_version

    @with_mpifort
    def test_success_dont_mock(self):
        """Run a realistic version check when mpifort exists."""
        with not_raises(CouldNotDeterminMpifortVersionError):
            get_mpifort_version()

    def test_not_installed(self, monkeypatch):
        """Test behavior when mpifort is not installed."""
        monkeypatch.setattr(subprocess, 'run', self.mock_run_fails)
        with pytest.raises(MpifortNotFoundError):
            get_mpifort_version()

    def test_could_not_determine(self, monkeypatch):
        """Test behavior when mpifort version cannot be determined."""
        def mock_run_could_not_determine(cmd, **_):
            """Fail only on version check."""
            try:
                return self.mock_mpifort_exists_ok(cmd)
            except ValueError:
                pass
            return self.mock_run_fails(cmd)

        monkeypatch.setattr(subprocess, 'run', mock_run_could_not_determine)
        with pytest.raises(CouldNotDeterminMpifortVersionError):
            get_mpifort_version()
