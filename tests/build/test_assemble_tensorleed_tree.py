"""Tests for the script build/scripts/assemble_tensorleed_tree.py."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-18'
__license__ = 'GPLv3+'

import importlib.util
from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

_SCRIPT = (Path(__file__).resolve().parents[2]
           / 'build/scripts/assemble_tensorleed_tree.py')


def _load_script():
    """Return the script, imported as a module."""
    spec = importlib.util.spec_from_file_location(_SCRIPT.stem, _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


script = _load_script()

_MAKE = {  # Command calling make, by whether on Windows
    True: ['cmd', '/c', 'make.bat'],
    False: ['make'],
    }
_RANDOM_POSIX = ['-O3', '-c', 'random_.c', '-o', 'random_.o']


@fixture(name='on_windows')
@parametrize(windows=(True, False), ids=('windows', 'posix'))
def fixture_on_windows(windows, mocker):
    """Pretend to run on Windows or not."""
    mocker.patch.object(script, '_on_windows', return_value=windows)
    return windows


@fixture(name='mock_assemble')
def fixture_mock_assemble(mocker):
    """Replace assemble_tree in the script with a mock."""
    return mocker.patch.object(script, 'assemble_tree')


class TestCommandLine:
    """Tests for parsing the command-line arguments in main."""

    def test_defaults(self, mock_assemble):
        """Check the default compilers."""
        script.main(['tree'])
        mock_assemble.assert_called_once_with('tree',
                                              compilers='gcc',
                                              c_compiler=None,
                                              fortran_compiler=None)

    def test_options(self, mock_assemble):
        """Check that compiler options are passed on."""
        script.main(['--compilers', 'intel', '--cc', 'icc',
                     '--fc', 'ifx', 'tree'])
        mock_assemble.assert_called_once_with('tree',
                                              compilers='intel',
                                              c_compiler='icc',
                                              fortran_compiler='ifx')

    def test_unknown_compilers(self, mock_assemble, capsys):
        """Check that an unknown compiler family is refused."""
        with pytest.raises(SystemExit):
            script.main(['--compilers', 'clang', 'tree'])
        mock_assemble.assert_not_called()
        assert 'invalid choice' in capsys.readouterr().err

    def test_error_exits(self, mocker):
        """Check that a TreeAssemblyError ends with a message."""
        mocker.patch.object(script,
                            'assemble_tree',
                            side_effect=script.TreeAssemblyError('bad'))
        with pytest.raises(SystemExit, match='Error: bad'):
            script.main(['tree'])


class TestCommands:
    """Tests for the compile commands built from the options."""

    @pytest.mark.usefixtures('on_windows')
    def test_default_random(self):
        """Check the default command compiling random_.c."""
        assert script._random_command('gcc', None) == ['gcc', *_RANDOM_POSIX]

    def test_default_make(self, on_windows):
        """Check the default command compiling eeasisss."""
        assert (script._make_command('gcc', None)
                == [*_MAKE[on_windows], 'gcc'])

    def test_intel_random(self, on_windows):
        """Check the command compiling random_.c with Intel."""
        expect = (['icx', '/c', 'random_.c', '/Fo:random_.o', '/O3']
                  if on_windows else ['icx', *_RANDOM_POSIX])
        assert script._random_command('intel', None) == expect

    def test_intel_make(self, on_windows):
        """Check the command compiling eeasisss with Intel."""
        assert (script._make_command('intel', None)
                == [*_MAKE[on_windows], 'intel'])

    @parametrize(compilers=('gcc', 'intel'))
    @pytest.mark.usefixtures('on_windows')
    def test_c_compiler(self, compilers):
        """Check that a C compiler given replaces the default."""
        command = script._random_command(compilers, 'my-cc')
        assert command[0] == 'my-cc'
        assert command[1:] == script._random_command(compilers, None)[1:]

    @parametrize(compilers=('gcc', 'intel'))
    def test_fortran_compiler_posix(self, compilers, mocker):
        """Check that a Fortran compiler given is passed to make."""
        mocker.patch.object(script, '_on_windows', return_value=False)
        assert script._make_command(compilers, 'my-fc') == [
            'make', compilers, f'{compilers}_fcomp=my-fc',
            ]

    def test_fortran_compiler_windows(self, mocker, tmp_path):
        """Check that a Fortran compiler is refused on Windows."""
        mocker.patch.object(script, '_on_windows', return_value=True)
        target = tmp_path / 'tree'
        with pytest.raises(script.TreeAssemblyError, match=r'make\.bat'):
            script.assemble_tree(target, fortran_compiler='ifx')
        assert not target.exists()

    def test_commands_passed_on(self, mocker, tmp_path):
        """Check that assemble_tree uses the commands built."""
        mocker.patch.object(script, '_on_windows', return_value=False)
        mock_assemble = mocker.patch.object(script, '_assemble_tree')
        script.assemble_tree(tmp_path,
                             compilers='intel',
                             fortran_compiler='ifx')
        mock_assemble.assert_called_once_with(
            tmp_path.resolve(),
            ['icx', *_RANDOM_POSIX],
            ['make', 'intel', 'intel_fcomp=ifx'],
            )
