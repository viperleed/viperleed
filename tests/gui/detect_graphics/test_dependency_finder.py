"""Tests for the QtDependencyFinder class of detect_graphics."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

from pathlib import Path
from subprocess import SubprocessError
from types import SimpleNamespace

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import QtDependencyFinder
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError
from viperleed.gui.qt_compat.binding import get_binding

from .conftest import BOTH
from .conftest import MODULE


class TestQtDependencyFinder:
    """Tests for the QtDependencyFinder helper class."""

    @parametrize(binding=BOTH)
    def test_init(self, binding, mocker):
        """Check initialization of a finder for a given binding."""
        find = mocker.patch.object(QtDependencyFinder,
                                   '_find_installed_binding')
        finder = QtDependencyFinder(binding)
        assert finder._binding == binding  # pylint: disable=W0212
        find.assert_not_called()

    def test_init_no_binding(self, mocker):
        """Check that the installed binding is used if none is given."""
        find = mocker.patch.object(QtDependencyFinder,
                                   '_find_installed_binding')
        finder = QtDependencyFinder()
        assert finder._binding is find.return_value  # pylint: disable=W0212

    @parametrize(installed=BOTH)
    def test_find_installed_binding(self, installed, mocker):
        """Check that the first installed candidate is found."""
        mocker.patch(f'{MODULE}.get_candidates', return_value=BOTH)
        mocker.patch('importlib.util.find_spec',
                     side_effect=lambda name: name == installed or None)
        # pylint: disable-next=protected-access
        assert QtDependencyFinder._find_installed_binding() == installed

    def test_find_installed_binding_none(self, mocker):
        """Check complaints when no candidate is installed."""
        mocker.patch(f'{MODULE}.get_candidates', return_value=BOTH)
        mocker.patch('importlib.util.find_spec', return_value=None)
        with pytest.raises(QtBindingNotFoundError):
            # pylint: disable-next=protected-access
            QtDependencyFinder._find_installed_binding()

    # {name: (binding, folders in the package, expected platforms)}
    _platforms_folders = {
        'PyQt5': (PYQT5, ('Qt5/plugins/platforms',), 'Qt5/plugins/platforms'),
        'PyQt6': (PYQT6, ('Qt6/plugins/platforms',), 'Qt6/plugins/platforms'),
        'PyQt5<=5.15.2': (PYQT5, ('Qt/plugins/platforms',),
                          'Qt/plugins/platforms'),
        'both': (PYQT6, ('Qt/plugins/platforms', 'Qt6/plugins/platforms'),
                 'Qt6/plugins/platforms'),
        'other major': (PYQT6, ('Qt5/plugins/platforms',), None),
        'system Qt': (PYQT5, (), None),
        }

    @parametrize('binding,folders,expect', _platforms_folders.values(),
                 ids=_platforms_folders)
    def test_get_platforms_path(self, binding, folders, expect,
                                tmp_path, mocker):
        """Check that plugins are found where the wheels put them."""
        for folder in folders:
            (tmp_path / folder).mkdir(parents=True)
        spec = SimpleNamespace(submodule_search_locations=[str(tmp_path)])
        find_spec = mocker.patch('importlib.util.find_spec',
                                 return_value=spec)
        import_module = mocker.patch('importlib.import_module')
        # pylint: disable-next=protected-access
        path = QtDependencyFinder(binding)._get_platforms_path()
        assert path == (tmp_path / expect if expect else None)
        find_spec.assert_called_once_with(binding)
        import_module.assert_not_called()

    @parametrize(spec=(None, SimpleNamespace(submodule_search_locations=None)))
    def test_get_platforms_path_no_package(self, spec, mocker):
        """Check that there are no plugins without a binding package."""
        mocker.patch('importlib.util.find_spec', return_value=spec)
        # pylint: disable-next=protected-access
        assert QtDependencyFinder(PYQT5)._get_platforms_path() is None

    @pytest.mark.usefixtures('qtc')  # Skips without Qt
    def test_get_platforms_path_real(self):
        """Check the plugins path of the installed binding."""
        # pylint: disable-next=protected-access
        path = QtDependencyFinder(get_binding())._get_platforms_path()
        # pylint: disable-next=magic-value-comparison
        assert path.name == 'platforms'
        assert path.is_dir()

    _install_linux = {
        'found': (('lib_one', 'lib_two'), 'sudo apt install lib_one lib_two'),
        'not found': ((), ''),
        }

    @parametrize('libs,expect', _install_linux.values(), ids=_install_linux)
    def test_install_libs_linux(self, libs, expect, mocker):
        """Check result of _find_install_for_libs_linux."""
        mock_list = mocker.patch(
            f'{MODULE}.QtDependencyFinder._list_install_for_libs_linux',
            return_value=libs,
            )
        missing = mocker.MagicMock()
        # pylint: disable-next=protected-access
        result = QtDependencyFinder._find_install_for_libs_linux(missing)
        mock_list.assert_called_once_with(missing)
        assert result == expect
        mocker.patch('sys.platform', 'linux')
        assert QtDependencyFinder.find_install_for_libs(missing) == result

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
            f'{MODULE}.QtDependencyFinder._find_install_for_libs_{platform}',
            create=True,
            )
        missing = mocker.MagicMock()
        finder = QtDependencyFinder(PYQT5)
        result = finder.find_install_for_libs(missing)
        assert result is mock_impl.return_value
        mock_impl.assert_called_once_with(missing)

    @parametrize(platform=_unsupported)
    def test_find_install_for_libs_unsupported(self, platform, mocker):
        """Check result of find_install_for_libs on unsupported platforms."""
        mocker.patch('sys.platform', platform)
        result = QtDependencyFinder.find_install_for_libs(mocker.MagicMock())
        no_install_suggestion = ''
        assert result == no_install_suggestion

    @parametrize(platform=(*_supported, *_unsupported))
    def test_find_install_for_libs_no_missing(self, platform, mocker):
        """Check shortcut if nothing is missing."""
        mocker.patch('sys.platform', platform)
        mock_delegate = mocker.patch.object(
            QtDependencyFinder,
            f'_find_install_for_libs_{platform}',
            create=True,
            )
        none_missing = {}
        assert not QtDependencyFinder.find_install_for_libs(none_missing)
        mock_delegate.assert_not_called()

    @parametrize(platform=_supported)
    def test_find_missing_dependencies_supported(self, platform, mocker):
        """Check result of find_missing_dependencies on supported platforms."""
        mocker.patch('sys.platform', platform)
        finder = QtDependencyFinder(PYQT5)
        mock_impl = mocker.patch.object(
            finder,
            f'_find_missing_deps_{platform}',
            create=True,
            )
        result = finder.find_missing_dependencies()
        assert result is mock_impl.return_value
        mock_impl.assert_called_once_with()

    @parametrize(platform=_unsupported)
    def test_find_missing_dependencies_unsupported(self, platform, mocker):
        """Check find_missing_dependencies raises on unsupported platforms."""
        mocker.patch('sys.platform', platform)
        finder = QtDependencyFinder(PYQT5)
        with pytest.raises(NotImplementedError):
            finder.find_missing_dependencies()

    @parametrize(exc=(FileNotFoundError, SubprocessError))
    def test_find_missing_dependencies_implementation_fails(self, exc, mocker):
        """Check find_missing_dependencies when implementation raises."""
        mocker.patch('sys.platform', 'linux')
        finder = QtDependencyFinder(PYQT5)
        mocker.patch.object(finder,
                            '_find_missing_deps_linux',
                            side_effect=exc)
        assert not finder.find_missing_dependencies()


class TestQtDependencyFinderLinux:
    """Tests for finding missing dependencies on Linux."""

    _modules_ldd = {  # {module file: output of ldd}
        'QtCore.abi3.so': 'libQt6Core.so.6 => /some/path (0x00007f)\n',
        'QtGui.abi3.so': (
            'libGL.so.1 => not found\n'
            'libEGL.so.1 => not found\n'
            ),
        'QtWidgets.cpython-37m-x86_64-linux-gnu.so': (  # Not abi3
            'libxkbcommon.so.0 => not found\n'
            'libGL.so.1 => not found\n'
            ),
        # Not needed by the GUI, or not compiled modules
        'QtNetwork.abi3.so': 'libgssapi_krb5.so.2 => not found\n',
        'sip.cpython-313-x86_64-linux-gnu.so': 'libfoo.so => not found\n',
        'QtGui.pyi': 'not a library, not found\n',
        }
    _no_package = {
        'not installed': None,
        'not a package': SimpleNamespace(submodule_search_locations=None),
        }
    _platforms_ldd = {  # {plugin file: output of ldd}
        'libqxcb.so': (
            'dep_one.so => not found\n'
            'dep_one.so => not found, again (a duplicate)\n'
            'dep_two.so.3 => dep was found, here the info\n'
            ),
        'libqwayland-generic.so.5': (  # Qt < 6.10
            'dep_three.so.9 (0x00007ffd6f7d7000) => not found, hex code\n'
            ),
        'libqwayland.so': (  # Qt >= 6.10
            'dep_five => was found\n'
            'dep_six.so => this one is also already installed\n'
            ),
        'nested/libqxcb.so.12': 'nested/deeper/dep_four => not found\n',
        # Not plugins of the platforms looked for
        'libqwayland-egl.so': 'dep_seven => not found\n',
        'libqoffscreen.so': 'dep_eight => not found\n',
        'libqxcbsomething.so': 'dep_nine => not found\n',
        }

    @fixture(name='fake_libs')
    def fixture_fake_libs(self, tmp_path, mocker):
        """Return a function that makes fake libraries for ldd.

        The function takes a dict {relative path: output of ldd},
        makes the files in a temporary folder, and returns it. ldd
        gives the output for each file.
        """
        outputs = {}

        def _ldd(cmd, **_):
            """Fake the result of subprocess.run(['ldd' ...])."""
            _, so_name = cmd
            name = Path(so_name).relative_to(tmp_path).as_posix()
            return mocker.MagicMock(stdout=outputs[name].encode())

        mocker.patch('subprocess.run', _ldd)

        def _make(ldd_output):
            outputs.update(ldd_output)
            for name in ldd_output:
                lib = tmp_path / name
                lib.parent.mkdir(parents=True, exist_ok=True)
                lib.touch()
            return tmp_path
        return _make

    def test_find_missing_deps_linux(self, mocker):
        """Check that modules and platform plugins are both checked."""
        finder = QtDependencyFinder(PYQT5)
        mocker.patch.object(finder, '_find_missing_deps_of_modules_linux',
                            return_value={'QtGui.abi3.so': ['libGL']})
        mocker.patch.object(finder, '_find_missing_deps_of_platforms_linux',
                            return_value={'libqxcb.so': ['libxcb']})
        # pylint: disable-next=protected-access
        missing = finder._find_missing_deps_linux()
        assert missing == {'QtGui.abi3.so': ['libGL'],
                           'libqxcb.so': ['libxcb']}

    def test_find_missing_deps_linux_unknown_platforms(self, mocker):
        """Check that missing module dependencies survive unknown Qt."""
        finder = QtDependencyFinder('PyQt7')
        mocker.patch.object(finder, '_find_missing_deps_of_modules_linux',
                            return_value={'QtGui.abi3.so': ['libGL']})
        mocker.patch.object(finder, '_get_platforms_path',
                            return_value=Path('platforms'))
        # pylint: disable-next=protected-access
        missing = finder._find_missing_deps_linux()
        assert missing == {'QtGui.abi3.so': ['libGL']}

    def test_find_missing_deps_of_modules_linux(self, fake_libs, mocker):
        """Check the missing dependencies of the Qt modules."""
        root = fake_libs(self._modules_ldd)
        spec = SimpleNamespace(submodule_search_locations=[str(root)])
        find_spec = mocker.patch('importlib.util.find_spec',
                                 return_value=spec)
        finder = QtDependencyFinder(PYQT5)
        # pylint: disable-next=protected-access
        missing = finder._find_missing_deps_of_modules_linux()
        assert missing == {
            'QtGui.abi3.so': ['libEGL', 'libGL'],
            'QtWidgets.cpython-37m-x86_64-linux-gnu.so': ['libGL',
                                                          'libxkbcommon'],
            }
        find_spec.assert_called_once_with(PYQT5)

    @parametrize(spec=_no_package.values(), ids=_no_package)
    def test_find_missing_deps_of_modules_linux_none(self, spec, mocker):
        """Check that nothing is missing without a binding package."""
        mocker.patch('importlib.util.find_spec', return_value=spec)
        run = mocker.patch('subprocess.run')
        finder = QtDependencyFinder(PYQT5)
        # pylint: disable-next=protected-access
        assert not finder._find_missing_deps_of_modules_linux()
        run.assert_not_called()

    def test_find_missing_deps_of_platforms_linux(self, fake_libs, mocker):
        """Check the missing dependencies of the platform plugins."""
        root = fake_libs(self._platforms_ldd)
        finder = QtDependencyFinder(PYQT5)
        mock_get = mocker.patch(f'{MODULE}._QtUnixPlatformsGetter.get',
                                return_value={'xcb', 'wayland'})
        mocker.patch.object(finder, '_get_platforms_path', return_value=root)
        # pylint: disable-next=protected-access
        missing = finder._find_missing_deps_of_platforms_linux()
        assert missing == {
            'libqxcb.so': ['dep_one'],
            'libqwayland-generic.so.5': ['dep_three'],
            'nested/libqxcb.so.12': ['dep_four'],
            }
        mock_get.assert_called_once_with(PYQT5)

    def test_find_missing_deps_of_platforms_linux_none(self, mocker):
        """Check that nothing is looked for without platform plugins."""
        finder = QtDependencyFinder(PYQT5)
        get_platforms = mocker.patch(f'{MODULE}._QtUnixPlatformsGetter.get')
        mocker.patch.object(finder, '_get_platforms_path',
                            return_value=None)
        run = mocker.patch('subprocess.run')
        # pylint: disable-next=protected-access
        assert not finder._find_missing_deps_of_platforms_linux()
        get_platforms.assert_not_called()
        run.assert_not_called()

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
        libs = QtDependencyFinder._list_install_for_libs_linux(missing)
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
        assert not QtDependencyFinder._list_install_for_libs_linux(missing)

    @parametrize(exc=(FileNotFoundError, SubprocessError))
    def test_list_install_libs_fails(self, exc, mocker):
        """Check result when subprocess raises."""
        missing = 'a', 'b'
        mocker.patch('subprocess.run', side_effect=exc)
        # pylint: disable-next=protected-access
        assert not QtDependencyFinder._list_install_for_libs_linux(missing)
