"""Tests for the script build/scripts/assemble_tensorleed_tree.py."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-18'
__license__ = 'GPLv3+'

import http.client
import importlib.util
import io
from pathlib import Path
import urllib.error

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

    @parametrize(compilers=('gcc', 'intel'))
    def test_make(self, compilers, on_windows):
        """Check the command compiling eeasisss without --fc."""
        assert (script._make_command(compilers, None)
                == [*_MAKE[on_windows], compilers])

    def test_intel_random(self, on_windows):
        """Check the command compiling random_.c with Intel."""
        expect = (['icx', '/c', 'random_.c', '/Fo:random_.o', '/O3']
                  if on_windows else ['icx', *_RANDOM_POSIX])
        assert script._random_command('intel', None) == expect

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
        mocker.patch.object(script, '_warn_about_newer_tags')
        mock_assemble = mocker.patch.object(script, '_assemble_tree')
        script.assemble_tree(tmp_path,
                             compilers='intel',
                             fortran_compiler='ifx')
        mock_assemble.assert_called_once_with(
            tmp_path.resolve(),
            ['icx', *_RANDOM_POSIX],
            ['make', 'intel', 'intel_fcomp=ifx'],
            )


class TestNewerTags:
    """Tests for the warning about tags newer than those pinned."""

    _pinned = ('v2.1.0', 'v2.0.0', 'v1.7.6', 'v1.7.5', 'v1.6.1',
               'EEASiSSS-v1.0.2')

    @parametrize(new_tags=(
        ('v2.2.0',),
        ('v2.1.1', 'v3.0'),
        ('EEASiSSS-v1.1.0',),
        ))
    def test_newer(self, new_tags, mocker, capsys):
        """Check the warning when newer tags exist."""
        mocker.patch.object(script,
                            '_list_tags',
                            return_value=[*new_tags, *self._pinned])
        script._warn_about_newer_tags()
        err = capsys.readouterr().err
        assert err.startswith('Warning:')
        assert ', '.join(sorted(new_tags)) in err

    @parametrize(tags=(
        (),
        _pinned,
        ('v1.7.7', 'v2.0.1', 'EEASiSSS-v1.0.1', 'latest', 'Other-v9'),
        ))
    def test_none_newer(self, tags, mocker, capsys):
        """Check that tags not newer than pinned give no output.

        Older versions than the newest pinned one of their kind (e.g.,
        v1.7.7 with v2.1.0 pinned), tags without a version, and kinds
        that are not pinned do not count.
        """
        mocker.patch.object(script, '_list_tags', return_value=list(tags))
        script._warn_about_newer_tags()
        assert not any(capsys.readouterr())

    @parametrize(exc=(
        urllib.error.URLError('no network'),
        TimeoutError('timed out'),
        ValueError('not JSON'),
        KeyError('name'),
        http.client.IncompleteRead(b'[{"na'),
        ))
    def test_listing_fails(self, exc, mocker, capsys):
        """Check that a failure to list the tags only gives a note."""
        mocker.patch.object(script, '_list_tags', side_effect=exc)
        script._warn_about_newer_tags()
        assert capsys.readouterr().err.startswith('Note:')

    def test_list_tags(self, mocker):
        """Check that tag names are read from the JSON listing."""
        listing = b'[{"name": "v2.1.0", "commit": {}}, {"name": "v2.0.0"}]'
        urlopen = mocker.patch.object(script.urllib.request,
                                      'urlopen',
                                      return_value=io.BytesIO(listing))
        assert script._list_tags() == ['v2.1.0', 'v2.0.0']
        assert urlopen.call_args[0] == (script.TAGS_URL,)

    def test_warned_after_assembly(self, mocker, tmp_path):
        """Check that assemble_tree checks for newer tags."""
        mocker.patch.object(script, '_assemble_tree')
        warn = mocker.patch.object(script, '_warn_about_newer_tags')
        script.assemble_tree(tmp_path)
        warn.assert_called_once_with()
