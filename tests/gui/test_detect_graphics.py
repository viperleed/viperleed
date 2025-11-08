"""Tests for module detect_graphics of viperleed.gui."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-22'
__license__ = 'GPLv3+'

import os
from pathlib import Path
from subprocess import SubprocessError


import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import _QtUnixPlatformsGetter
from viperleed.gui.detect_graphics import _TIMEOUT
from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.detect_graphics import PyQtSanityChecker
from viperleed.gui.detect_graphics import Qt5DependencyFinder
from viperleed.gui.detect_graphics import _import_dummy_pyqt_and_run_app
from viperleed.gui.detect_graphics import _init_dummy_qapp
from viperleed.gui.detect_graphics import find_missing_qt_dependencies
from viperleed.gui.detect_graphics import has_graphics
from viperleed.gui.detect_graphics import has_pyqt
from viperleed.gui.detect_graphics import suppress_file_permission_warnings

_MODULE = 'viperleed.gui.detect_graphics'


# TODO: would be better to find a way to mock at least
# _import_dummy_pyqt_and_run_app instead of mocking
# multiprocessing.Process. My tests seem to reveal
# that, while _check_pyqt_sanity gets the correct
# target, that target is never actually called.


@fixture(name='mock_multiprocessing')
def fixture_mock_process(mocker):
    """Replace multiprocessing.Process with a mock.

    Notice that it is necessary to do so in order to correctly
    test the behavior of has_display. In fact, we can't mock
    detect_graphics.QApplication as the module is re-imported
    by multiprocessing. We could, in principle, mock the target
    of multiprocessing.Process, i.e., _init_dummy_qapp. However,
    I tried without success to do so.
    """
    def _mock(has_display):
        mock_process = mocker.MagicMock()
        mock_process_cls = mocker.MagicMock(return_value=mock_process)
        # Mimic the expected behavior when a display is not present
        mock_process.is_alive.return_value = not has_display
        mocker.patch(f'{_MODULE}.Process', mock_process_cls)
        return mock_process, mock_process_cls
    return _mock


def test_init_dummy_qapp(mocker):
    """Check calls in _init_dummy_qapp."""
    mock_qtwidgets = mocker.patch('importlib.import_module').return_value
    mock_qapp = mock_qtwidgets.QApplication
    mock_suppress = mocker.patch(
        f'{_MODULE}.suppress_file_permission_warnings'
        )
    with pytest.raises(SystemExit):
        _init_dummy_qapp()
    mock_suppress.assert_called_once_with()
    mock_qapp.assert_called_once_with([])


@parametrize(pyqt_present=(True, False))
@parametrize(has_display=(True, False))
def test_has_graphics(pyqt_present, has_display, mocker):
    """Ensure correct identification of graphics capabilities."""
    mocker.patch(f'{_MODULE}.has_pyqt', return_value=pyqt_present)
    check_sanity = mocker.patch(
        f'{_MODULE}.check_pyqt_sanity',
        return_value=PyQtSanity.OK if has_display else PyQtSanity.NO_DISPLAY
        )
    result = has_graphics()
    assert result == (pyqt_present and has_display)

    # Check also internal calls
    if pyqt_present:
        check_sanity.assert_called_once_with()
    else:
        check_sanity.assert_not_called()


@parametrize(pyqt_sanity=PyQtSanity)
def test_has_pyqt(pyqt_sanity, mocker):
    """Ensure correct detection of PyQt5 presence."""
    mocker.patch(f'{_MODULE}.check_pyqt_sanity', return_value=pyqt_sanity)
    pyqt_present = pyqt_sanity in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}
    assert has_pyqt() == pyqt_present


class TestImportDummyPyqtAndRunApp:
    """Tests for the _import_dummy_pyqt_and_run_app subprocess worker."""

    # The modules that we attempt to import
    _modules = (
        'PyQt5.QtCore',
        'PyQt5.QtWidgets',
        )

    @fixture(name='mock_dummy_qapp', autouse=True)
    def fixture_mock_dummy_qapp(self, mocker):
        """Replace _init_dummy_qapp with a mock."""
        # This is to ensure tests don't hang forever on a system
        # without graphics.
        return mocker.patch(f'{_MODULE}._init_dummy_qapp')

    def test_cant_import_pyqt(self, mocker):
        """Check the expected exit code when importing PyQt5 fails."""
        mocker.patch('importlib.import_module',
                     side_effect=ImportError)
        with pytest.raises(SystemExit) as excinfo:
            _import_dummy_pyqt_and_run_app()
        assert excinfo.value.code == PyQtSanity.IMPORT_ERROR.value

    @parametrize(module=_modules)
    def test_cant_import_pyqt_submodule(self, module, mocker):
        """Check the expected exit code when importing PyQt5 fails."""
        def _fail_module_import(mod_name):
            if mod_name == module:
                raise ImportError
        mocker.patch('importlib.import_module',
                     side_effect=_fail_module_import)
        with pytest.raises(SystemExit) as excinfo:
            _import_dummy_pyqt_and_run_app()
        assert excinfo.value.code == PyQtSanity.IMPORT_ERROR.value

    def test_pyqt_not_found(self, mocker):
        """Check the expected exit code when PyQt5 is not found."""
        mocker.patch('importlib.import_module',
                     side_effect=ModuleNotFoundError)
        with pytest.raises(SystemExit) as excinfo:
            _import_dummy_pyqt_and_run_app()
        assert excinfo.value.code == PyQtSanity.NOT_FOUND.value

    def test_success(self, mocker, mock_dummy_qapp):
        """Test an execution that does not cause the interpreter to exit."""
        # NB: normally, the interpreter always exits when
        # _import_dummy_pyqt_and_run_app is called. The only
        # reason why we don't exit here is the patching of
        # _init_dummy_qapp via mock_dummy_qapp.
        mock_import = mocker.patch('importlib.import_module')
        _import_dummy_pyqt_and_run_app()

        expect_calls = [
            mocker.call('PyQt5'),
            *(mocker.call(m) for m in self._modules)
            ]
        mock_import.assert_has_calls(expect_calls, any_order=True)
        mock_dummy_qapp.assert_called_once_with()


class TestFindMissingDependencies:
    """Tests for the find_missing_qt_dependencies function."""

    def test_raises(self, mocker):
        """Check complaints when calling the function without Qt."""
        mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.NOT_FOUND)
        with pytest.raises(NotImplementedError):
            find_missing_qt_dependencies()

    def test_windows(self, mocker):
        """Check the result of finding dependencies on Windows."""
        mocker.patch('sys.platform', 'win32')
        missing = find_missing_qt_dependencies()
        assert not missing

    def test_not_implemented(self, mocker):
        """Check the result when dependency finding is not implemented."""
        mocker.patch('sys.platform', 'some-unsupported-platform')
        mock_finder = mocker.MagicMock()
        mocker.patch(f'{_MODULE}.Qt5DependencyFinder',
                     return_value=mock_finder)
        mock_finder.find_missing_dependencies.side_effect = NotImplementedError
        missing = find_missing_qt_dependencies()
        assert not missing

    def test_missing_deps(self, mocker):
        """Check the expected result when missing dependencies are found."""
        mocker.patch('sys.platform', 'some-supported-platform')
        mock_finder = mocker.MagicMock()
        mocker.patch(f'{_MODULE}.Qt5DependencyFinder',
                     return_value=mock_finder)
        missing = find_missing_qt_dependencies()
        assert missing is mock_finder.find_missing_dependencies.return_value


class TestSuppressFilePermissionWarnings:
    """Tests for the suppress_file_permission_warnings function."""

    def test_no_xdg_runtime_dir(self, mocker):
        """Check no changes if XDG_RUNTIME_DIR is not set (e.g., not UNIX)."""
        mocker.patch(f'{_MODULE}._UNIX_RUNTIME_DIR', '')
        mock_mkdir = mocker.patch('pathlib.Path.mkdir')
        mocker.patch.dict('os.environ', clear=True)

        suppress_file_permission_warnings()
        mock_mkdir.assert_not_called()
        # pylint: disable-next=magic-value-comparison
        assert 'XDG_RUNTIME_DIR' not in os.environ

    def test_xdg_set(self, tmp_path, mocker):
        """Check implementation when XDG_RUNTIME_DIR is set."""
        expect_kwargs = {'mode': 0o0700, 'parents': True, 'exist_ok': True}
        mkdir = Path.mkdir
        def mock_mkdir(*args, **kwargs):
            assert kwargs == expect_kwargs
            mkdir(*args, **kwargs)

        mocker.patch(f'{_MODULE}._UNIX_RUNTIME_DIR', str(tmp_path))
        mocker.patch('pathlib.Path.mkdir', mock_mkdir)
        expect_path = (tmp_path / 'viperleed-gui').resolve()
        assert not expect_path.exists()

        suppress_file_permission_warnings()
        assert expect_path.is_dir()
        assert os.environ['XDG_RUNTIME_DIR'] == str(expect_path)

    def test_mkdir_fails(self, tmp_path, mocker):
        """Check setting of fallback directory if the default fails."""
        fails = tmp_path/'fails'
        fallback = tmp_path/'fallback'
        fallback.mkdir()  # Otherwise kwargs['parents'] is mismatched
        expect_kwargs = {'mode': 0o0700, 'parents': True, 'exist_ok': True}
        not_created = fails/'viperleed-gui'
        expect_path = fallback/'viperleed-gui'

        mocker.patch(f'{_MODULE}._UNIX_RUNTIME_DIR', str(fails))
        mock_qtc = mocker.patch('importlib.import_module').return_value
        mock_qtc.QStandardPaths.writableLocation.return_value = str(fallback)

        mkdir = Path.mkdir
        def mkdir_raises(path, *args, **kwargs):
            check_kwargs = kwargs.copy()
            if path in (not_created, expect_path):
                if args:
                    check_kwargs['mode'] = args[0]
                assert check_kwargs == expect_kwargs
            if path == not_created:
                raise OSError
            mkdir(path, *args, **kwargs)

        mocker.patch('pathlib.Path.mkdir', mkdir_raises)

        suppress_file_permission_warnings()
        assert not not_created.exists()
        assert expect_path.is_dir()
        assert os.environ['XDG_RUNTIME_DIR'] == str(expect_path)


class TestPyQtSanityChecker:
    """Tests for the PyQtSanityChecker class."""

    def test_checks_only_once(self, mocker):
        """Ensure that checking is skipped when already done."""
        mock_sanity = mocker.MagicMock()
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity',
                                            return_value=mock_sanity)
        # Reset status manually, just in case it was checked already
        mocker.patch.object(PyQtSanityChecker, '_status', None)
        checker = PyQtSanityChecker()

        # Call it once: performs the check.
        assert checker() is mock_sanity
        # pylint: disable-next=protected-access
        assert PyQtSanityChecker._status is mock_sanity

        # Call it again: should skip the check
        assert checker() is mock_sanity
        mock_do_check.assert_called_once_with()

    def test_forget_status(self, mocker):
        """Ensure resetting re-checks."""
        # Reset status manually, just in case it was checked already
        mocker.patch.object(PyQtSanityChecker, '_status', None)
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity')
        checker = PyQtSanityChecker()
        # Check once to ensure forget_status does the right thing
        checker()

        checker.forget_status()
        checker()

        # pylint: disable-next=magic-value-comparison
        assert mock_do_check.call_count == 2

    _exit_codes = {
        0: PyQtSanity.OK,
        1: PyQtSanity.NOT_FOUND,
        2: PyQtSanity.IMPORT_ERROR,
        3: PyQtSanity.RUNTIME_CRASH,
        4: PyQtSanity.NO_DISPLAY,
        -6: PyQtSanity.RUNTIME_CRASH,   # SIGABRT
        134: PyQtSanity.RUNTIME_CRASH,  # SIGABRT on Unix
        'any code': PyQtSanity.RUNTIME_CRASH
        }

    @parametrize(has_display=(True, False))
    @parametrize('exit_code,expect', _exit_codes.items())
    # pylint: disable-next=too-many-arguments   # 2/6 fixtures
    def test_check_sanity(self, has_display, exit_code, expect,
                          mock_multiprocessing, mocker):
        """Check expected return value when no display is available."""
        proc, proc_cls = mock_multiprocessing(has_display)
        proc.exitcode = exit_code
        expect_calls = {
            'start': mocker.call(),
            'join': mocker.call(_TIMEOUT),
            'is_alive': mocker.call(),
            }
        if not has_display:
            expect_calls['terminate'] = mocker.call()
            expect = PyQtSanity.NO_DISPLAY
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result is expect
        proc_cls.assert_called_once_with(target=_import_dummy_pyqt_and_run_app)
        for method_name, call in expect_calls.items():
            method = getattr(proc, method_name)
            assert method.mock_calls == [call]


@fixture(name='mock_pyqt_root')
def fixture_mock_pyqt_root(mocker):
    """Fake the existence of PyQt5's root path."""
    mock_root = '/path/to/pyqt'
    mock_pyqt = mocker.patch('importlib.import_module').return_value
    mock_pyqt.__file__ = f'{mock_root}/__init__.py'
    mocker.patch('pathlib.Path.exists', resturn_value=True)
    return Path(mock_root).resolve()


class TestQt5DependencyFinder:
    """Tests for the Qt5DependencyFinder helper class."""

    def test_init(self, mock_pyqt_root):
        """Check initialization of a finder."""
        expect_root = mock_pyqt_root / 'Qt5/plugins/platforms'
        finder = Qt5DependencyFinder()
        # pylint: disable-next=protected-access
        assert finder._root == expect_root

    _install_linux = {
        'found': (('lib_one', 'lib_two'), 'sudo apt install lib_one lib_two'),
        'not found': ((), ''),
        }

    @parametrize('libs,expect', _install_linux.values(), ids=_install_linux)
    def test_install_libs_linux(self, libs, expect, mocker):
        """Check result of _find_install_for_libs_linux."""
        mock_list = mocker.patch(
            f'{_MODULE}.Qt5DependencyFinder._list_install_for_libs_linux',
            return_value=libs,
            )
        missing = mocker.MagicMock()
        # pylint: disable-next=protected-access
        result = Qt5DependencyFinder._find_install_for_libs_linux(missing)
        mock_list.assert_called_once_with(missing)
        assert result == expect
        mocker.patch('sys.platform', 'linux')
        assert Qt5DependencyFinder.find_install_for_libs(missing) == result

    _unsupported = (
        'aix',     # Python >= 3.8
        'aix5',    # Python < 3.8
        'aix7',    # Python < 3.8
        'atheos',
        'cygwin',
        'darwin',  # macOS
        'win32',
        'msys',
        'os2',
        'os2emx',
        'riscos',
        'freebsd7',
        'freebsd8',
        'freebsdN',
        'openbsd6',
        'invalidplatform',
        )
    _supported = (
        'linux',
        'validplatform',
        )

    @parametrize(platform=_supported)
    def test_find_install_for_libs_supported(self, platform, mocker):
        """Check result of find_install_for_libs on supported platforms."""
        mocker.patch('sys.platform', platform)
        mock_impl = mocker.patch(
            # Patch the class, as find_install_for_libs is @classmethod
            f'{_MODULE}.Qt5DependencyFinder._find_install_for_libs_{platform}',
            create=True,
            )
        missing = mocker.MagicMock()
        finder = Qt5DependencyFinder()
        result = finder.find_install_for_libs(missing)
        assert result is mock_impl.return_value
        mock_impl.assert_called_once_with(missing)

    @parametrize(platform=_unsupported)
    def test_find_install_for_libs_unsupported(self, platform, mocker):
        """Check result of find_install_for_libs on unsupported platforms."""
        mocker.patch('sys.platform', platform)
        result = Qt5DependencyFinder.find_install_for_libs(mocker.MagicMock())
        no_install_suggestion = ''
        assert result == no_install_suggestion

    @parametrize(platform=(*_supported, *_unsupported))
    def test_find_install_for_libs_no_missing(self, platform, mocker):
        """Check shortcut if nothing is missing."""
        mocker.patch('sys.platform', platform)
        mock_delegate = mocker.patch.object(
            Qt5DependencyFinder,
            f'_find_install_for_libs_{platform}',
            create=True,
            )
        none_missing = {}
        assert not Qt5DependencyFinder.find_install_for_libs(none_missing)
        mock_delegate.assert_not_called()

    @pytest.mark.usefixtures('mock_pyqt_root')
    @parametrize(platform=_supported)
    def test_find_missing_dependencies_supported(self, platform, mocker):
        """Check result of find_missing_dependencies on supported platforms."""
        mocker.patch('sys.platform', platform)
        finder = Qt5DependencyFinder()
        mock_impl = mocker.patch.object(
            finder,
            f'_find_missing_deps_{platform}',
            create=True,
            )
        result = finder.find_missing_dependencies()
        assert result is mock_impl.return_value
        mock_impl.assert_called_once_with()

    @pytest.mark.usefixtures('mock_pyqt_root')
    @parametrize(platform=_unsupported)
    def test_find_missing_dependencies_unsupported(self, platform, mocker):
        """Check find_missing_dependencies raises on unsupported platforms."""
        mocker.patch('sys.platform', platform)
        finder = Qt5DependencyFinder()
        with pytest.raises(NotImplementedError):
            finder.find_missing_dependencies()

    @parametrize(exc=(FileNotFoundError, SubprocessError))
    def test_find_missing_dependencies_implementation_fails(self, exc, mocker):
        """Check find_missing_dependencies when implementation raises."""
        mocker.patch('sys.platform', 'linux')
        finder = Qt5DependencyFinder()
        mocker.patch.object(finder,
                            '_find_missing_deps_linux',
                            side_effect=exc)
        assert not finder.find_missing_dependencies()


class TestQt5DependencyFinderLinux:
    """Tests for finding missing dependencies on Linux."""

    @pytest.mark.usefixtures('mock_pyqt_root')
    def test_find_missing_deps_linux(self, mocker):
        """Check the result of _find_missing_deps_linux."""
        finder = Qt5DependencyFinder()
        mocker.patch(f'{_MODULE}._QtUnixPlatformsGetter.get',
                     return_value=('missing_deps', 'no_missing_deps'))
        root = finder._root  # pylint: disable=protected-access
        so_files = (
            root/'libmissing_deps.so',
            root/'anothermissing_deps.so.999',
            root/'nested/missing_deps.so.12',
            root/'no_missing_deps.so.9',
            )
        dependencies = {
            str(root/'libmissing_deps.so'): '''\
dep_one.so => not found
dep_one.so => not found, it's a duplicate
dep_two.so.3 => dep was found, here the info\n''',
            str(root/'anothermissing_deps.so.999'): '''\
dep_three.so.9 (0x00007ffd6f7d7000) => not found, hex code at left\n''',
            str(root/'nested/missing_deps.so.12'): '''\
nested/deeper/dep_four => not found\n''',
            str(root/'no_missing_deps.so.9'): '''\
dep_five => was found
dep_six.so => this one is also already installed''',
            }
        expect = {
            'libmissing_deps.so': ['dep_one'],
            'anothermissing_deps.so.999': ['dep_three'],
            'nested/missing_deps.so.12': ['dep_four'],
            }

        def mock_glob(_, pattern):
            pattern = pattern.replace('*', '')
            for file in so_files:
                if pattern in file.name:
                    yield file

        def mock_ldd(cmd, **_):
            """Fake the result of subprocess.run(['ldd' ...])."""
            _, so_name = cmd
            mock_process = mocker.MagicMock()
            mock_process.stdout = dependencies[so_name].encode()
            return mock_process

        mocker.patch('pathlib.Path.rglob', mock_glob)
        mocker.patch('subprocess.run', mock_ldd)

        # pylint: disable-next=protected-access
        missing = finder._find_missing_deps_linux()
        assert missing == expect

    def test_list_install_libs_success(self, mocker):
        """Check the successful listing of apt modules."""
        missing = ['multiple', 'devonly', 'another']
        expect_rgx = 'multiple|devonly|another'
        def mock_apt(cmd, **__):
            """Fake the result of calling 'apt-cache search'."""
            *_, rgx = cmd
            assert rgx == expect_rgx
            mock_proc = mocker.MagicMock()
            mock_proc.stdout = '''\
multiplelib-one - description -- some more text
multiplelib-one4 - more recent than the previous - should be selected
multiplelib-one5-dev - this one looks more recent but it is a development one
devonly-three-dev - this comes only with developer options
another-lib - this has no number'''.encode()
            return mock_proc

        expect = ['multiplelib-one4',
                  'devonly-three-dev',
                  'another-lib']
        mocker.patch('subprocess.run', mock_apt)
        # pylint: disable-next=protected-access
        libs = Qt5DependencyFinder._list_install_for_libs_linux(missing)
        assert libs == expect

    def test_list_install_libs_one_missing(self, mocker):
        """Check that no info is returned if any install is missing."""
        missing = ['one', 'cant-find-install']
        def mock_apt(_, **__):
            """Fake the result of calling 'apt-cache search'."""
            mock_proc = mocker.MagicMock()
            mock_proc.stdout = 'onelib - this is found'.encode()
            return mock_proc

        mocker.patch('subprocess.run', mock_apt)
        # pylint: disable-next=protected-access
        assert not Qt5DependencyFinder._list_install_for_libs_linux(missing)

    @parametrize(exc=(FileNotFoundError, SubprocessError))
    def test_list_install_libs_fails(self, exc, mocker):
        """Check result when subprocess raises."""
        missing = 'a', 'b'
        mocker.patch('subprocess.run', side_effect=exc)
        # pylint: disable-next=protected-access
        assert not Qt5DependencyFinder._list_install_for_libs_linux(missing)


class TestQtUnixPlatformGetter:
    """Tests for the _QtUnixPlatformsGetter helper class."""

    _supported = (
        '5.15',
        *(f'5.15.{v}' for v in range(18)),
        )
    _unsupported = (
        '4.12.99',
        '5.15.18',
        '5.16',
        '6.0.0',
        )

    @fixture(name='mock_qtc')
    def fixture_mock_qtc(self, mocker):
        """Return a fake QtCore module."""
        return mocker.patch('importlib.import_module').return_value

    @fixture
    @parametrize(version=_supported)
    def mock_supported_version(self, version, mock_qtc):
        """Fake a supported PyQt5 version."""
        mock_qtc.QT_VERSION_STR = version

    def test_cant_import_qtcore(self, mocker):
        """Check outcome when QtCore fails to import."""
        mocker.patch('importlib.import_module', side_effect=ImportError)
        assert not _QtUnixPlatformsGetter.get()

    @parametrize(version=_unsupported)
    def test_unsupported_version(self, version, mock_qtc):
        """Check complaints for unsupported Qt versions."""
        mock_qtc.QT_VERSION_STR = version
        with pytest.raises(NotImplementedError):
            _QtUnixPlatformsGetter.get()

    @pytest.mark.usefixtures('mock_supported_version')
    def test_qt_qpa_platform(self, mocker):
        """Check result when QT_QPA_PLATFORM is set."""
        mock_qpa = str(mocker.MagicMock())
        mocker.patch.dict('os.environ',
                          QT_QPA_PLATFORM=mock_qpa,
                          clear=True)
        result = _QtUnixPlatformsGetter.get()
        assert result == {mock_qpa}

    @pytest.mark.usefixtures('mock_supported_version')
    @parametrize(default=(None, 'some-default'))
    def test_no_xdg_session(self, default, mocker):
        """Check result when no XDG_SESSION_TYPE is set."""
        kwargs = {'clear': True}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get()
        assert result == {default}

    @parametrize(default=(None, 'some-default'))
    @parametrize(gnome=({},
                        {'XDG_CURRENT_DESKTOP': 'some-gnome-value'},
                        {'XDG_SESSION_DESKTOP': 'some-gnome-value'}))
    @pytest.mark.usefixtures('mock_supported_version')
    def test_wayland(self, default, gnome, mocker):
        """Check result when XDG_SESSION_TYPE is set to wayland."""
        kwargs = {'clear': True,
                  'XDG_SESSION_TYPE': 'wayland',
                  **gnome}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get()
        expect = {default} if gnome else {default, 'wayland'}
        assert result == expect

    @parametrize(default=(None, 'some-default'))
    @parametrize(xdg=('valid', 'invalid', 'xcb'))
    @pytest.mark.usefixtures('mock_supported_version')
    def test_xdg_other(self, default, xdg, mocker):
        """Check result when XDG_SESSION_TYPE is set to non-wayland/x11."""
        kwargs = {'clear': True,
                  'XDG_SESSION_TYPE': xdg}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get()
        assert result == {xdg, default}

    @parametrize(default=(None, 'some-default'))
    @pytest.mark.usefixtures('mock_supported_version')
    def test_x_eleven(self, default, mocker):
        """Check result when XDG_SESSION_TYPE is set to x11."""
        kwargs = {'clear': True,
                  'XDG_SESSION_TYPE': 'x11'}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get()
        assert result == {default, 'xcb'}
