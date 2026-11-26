"""Tests for module fortran_utils of viperleed.calc.lib."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexandra Mia Imre (@alexmiame)',
    'Michele Riva (@michele-riva)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'

import shutil
from subprocess import CalledProcessError
from subprocess import CompletedProcess
from subprocess import TimeoutExpired

import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.fortran_utils import CompilerNotFoundError
from viperleed.calc.lib.fortran_utils import NoCompilerVersionFoundError
from viperleed.calc.lib.fortran_utils import get_mpifort_version
from viperleed.calc.lib.fortran_utils import gfortran_finds_openblas
from viperleed.calc.lib.fortran_utils import wrap_fortran_line

from ...helpers import not_raises


class TestGetMpifortVersion:
    """Tests for the get_mpifort_version function."""

    mock_version = '13.2.0'
    check_version_cmd = 'mpifort', '--version'

    @staticmethod
    def mock_run_fails(cmd, **_):
        """Simulate a failed subprocess.run call."""
        raise CalledProcessError(1,
                                 cmd,
                                 'fake stdout'.encode(),
                                 'fake stderr'.encode())

    @staticmethod
    def mock_run_no_version(cmd, **_):
        """Simulate the failed outcome of a version check."""
        return CompletedProcess(args=cmd,
                                returncode=0,
                                stdout=b'Unknown output\n')

    def mock_run_success(self, cmd, **_):
        """Simulate a successful subprocess.run."""
        if cmd == self.check_version_cmd:
            # Simulate the expected output for mpifort --version
            result = (f'GNU Fortran (Homebrew GCC {self.mock_version}) '
                      f'{self.mock_version}\n')
            return CompletedProcess(args=cmd,
                                    returncode=0,
                                    stdout=result.encode())
        raise ValueError(f'Unexpected command {cmd}')

    def test_success(self, mocker):
        """Test successful retrieval of mpifort version."""
        mocker.patch('shutil.which', return_value=True)
        mocker.patch('subprocess.run', new=self.mock_run_success)
        version = get_mpifort_version()
        assert version == self.mock_version

    def test_success_dont_mock(self, missing_dependency):
        """Run a realistic version check when mpifort exists."""
        if not shutil.which('mpifort'):
            missing_dependency('mpifort')
        with not_raises(NoCompilerVersionFoundError):
            version = get_mpifort_version()
        assert any(part > 0 for part in version)

    def test_not_installed(self, mocker):
        """Test behavior when mpifort is not installed."""
        mocker.patch('shutil.which', return_value=None)
        with pytest.raises(CompilerNotFoundError):
            get_mpifort_version()

    @parametrize(fake_run=('mock_run_no_version', 'mock_run_fails'))
    def test_could_not_determine(self, fake_run, mocker):
        """Test behavior when mpifort version cannot be determined."""
        mocker.patch('shutil.which', return_value=True)
        mocker.patch('subprocess.run', new=getattr(self, fake_run))
        with pytest.raises(NoCompilerVersionFoundError):
            get_mpifort_version()


class TestGfortranFindsOpenblas:
    """Tests for the gfortran_finds_openblas function."""

    lib_name = 'libopenblas.dll.a'
    probe_cmd = 'gfortran', f'-print-file-name={lib_name}'
    on_windows = 'viperleed.calc.lib.fortran_utils._on_windows'

    @pytest.fixture(autouse=True)
    def fake_windows(self, mocker):
        """Pretend to be on Windows, where the probe is meant to run."""
        mocker.patch(self.on_windows, return_value=True)

    def mock_run_prints(self, mocker, output):
        """Make subprocess.run print output, and return the mock."""
        if not isinstance(output, bytes):
            output = f'{output}\n'.encode()
        def _run(cmd, **_):
            assert cmd == self.probe_cmd
            return CompletedProcess(args=cmd, returncode=0, stdout=output)
        return mocker.patch('subprocess.run', side_effect=_run)

    def test_found(self, tmp_path, mocker):
        """Check that an existing absolute path means found."""
        lib_path = tmp_path / self.lib_name
        lib_path.touch()
        # MinGW gfortran writes its stdout in text mode: CRLF
        run = self.mock_run_prints(mocker, f'{lib_path}\r\n'.encode())
        assert gfortran_finds_openblas()
        run.assert_called_once()

    def test_not_utf8(self, mocker):
        """Check that output that is not UTF-8 does not raise."""
        output = b'C:\\Users\\M\xfcller\\libopenblas.dll.a\r\n'  # cp1252
        self.mock_run_prints(mocker, output)
        with not_raises(UnicodeDecodeError):
            assert not gfortran_finds_openblas()

    def test_bare_name(self, mocker):
        """Check that an echoed bare name means not found."""
        self.mock_run_prints(mocker, self.lib_name)
        assert not gfortran_finds_openblas()

    def test_bare_name_in_cwd(self, tmp_path, mocker, monkeypatch):
        """Check that a file of that name in the cwd does not count."""
        monkeypatch.chdir(tmp_path)
        (tmp_path / self.lib_name).touch()
        self.mock_run_prints(mocker, self.lib_name)
        assert not gfortran_finds_openblas()

    def test_missing_path(self, tmp_path, mocker):
        """Check that an absolute path to no file means not found."""
        self.mock_run_prints(mocker, tmp_path / self.lib_name)
        assert not gfortran_finds_openblas()

    _fails = {
        'non-zero exit': CalledProcessError(1, probe_cmd),
        'no gfortran': FileNotFoundError('gfortran'),
        'timeout': TimeoutExpired(probe_cmd, 30),
        }

    @parametrize(exc=_fails.values(), ids=_fails)
    def test_probe_fails(self, exc, mocker):
        """Check that a failing probe returns False without raising."""
        mocker.patch('subprocess.run', side_effect=exc)
        assert not gfortran_finds_openblas()

    def test_probe_args(self, mocker):
        """Check that the probe raises on failure and has a timeout."""
        run = self.mock_run_prints(mocker, self.lib_name)
        gfortran_finds_openblas()
        _, kwargs = run.call_args
        assert kwargs['check']
        assert kwargs['timeout'] > 0

    def test_not_windows(self, mocker):
        """Check that the probe refuses to run off Windows."""
        mocker.patch(self.on_windows, return_value=False)
        run = mocker.patch('subprocess.run')
        with pytest.raises(RuntimeError, match='Windows'):
            gfortran_finds_openblas()
        run.assert_not_called()


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
