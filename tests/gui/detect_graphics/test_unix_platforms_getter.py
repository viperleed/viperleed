"""Tests for the _QtUnixPlatformsGetter class of detect_graphics."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import _QtUnixPlatformsGetter
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6


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
        '7.0.0',
        )
    _qt6_env = {  # {environment: platforms}
        'no session': ({}, {'xcb'}),
        'x11': ({'XDG_SESSION_TYPE': 'x11', 'DISPLAY': ':0'}, {'xcb'}),
        'wayland session': ({'XDG_SESSION_TYPE': 'wayland'},
                            {'xcb', 'wayland'}),
        'wayland display': ({'WAYLAND_DISPLAY': ''}, {'xcb', 'wayland'}),
        'wayland gnome': ({'XDG_SESSION_TYPE': 'wayland',
                           'XDG_CURRENT_DESKTOP': 'GNOME'},
                          {'xcb', 'wayland'}),
        'requested': ({'QT_QPA_PLATFORM': 'offscreen',
                       'XDG_SESSION_TYPE': 'wayland'},
                      {'offscreen'}),
        'requested list': ({'QT_QPA_PLATFORM': 'wayland;xcb:arg=1'},
                           {'wayland', 'xcb'}),
        }

    @fixture(name='mock_qtc')
    def fixture_mock_qtc(self, mocker, mock_import_module):
        """Return a fake QtCore module."""
        mock_qtc = mocker.MagicMock()
        fake_modules = {'PyQt5.QtCore': mock_qtc}
        mock_import_module(**fake_modules)
        return mock_qtc

    @fixture
    @parametrize(version=_supported)
    def mock_supported_version(self, version, mock_qtc):
        """Fake a supported Qt version."""
        mock_qtc.QT_VERSION_STR = version

    def test_cant_import_qtcore(self, mocker):
        """Check outcome when QtCore fails to import."""
        mocker.patch('importlib.import_module', side_effect=ImportError)
        assert not _QtUnixPlatformsGetter.get(PYQT5)

    @parametrize(version=_unsupported)
    def test_unsupported_version(self, version, mock_qtc):
        """Check complaints for unsupported Qt versions."""
        mock_qtc.QT_VERSION_STR = version
        with pytest.raises(NotImplementedError):
            _QtUnixPlatformsGetter.get(PYQT5)

    @parametrize('env,expect', _qt6_env.values(), ids=_qt6_env)
    @parametrize(version=('6.0.0', '6.2.4', '6.5.3', '6.10.0', '6.11.0'))
    def test_qt6(self, env, expect, version, mocker, mock_import_module):
        """Check the platforms of Qt 6, via PyQt6.QtCore."""
        mock_qtc = mocker.MagicMock(QT_VERSION_STR=version)
        mock_import_module(**{'PyQt6.QtCore': mock_qtc})
        mocker.patch.dict('os.environ', clear=True, **env)
        assert _QtUnixPlatformsGetter.get(PYQT6) == expect

    def test_qt7(self, mocker, mock_import_module):
        """Check that Qt 7 is not supported yet, via PyQt7.QtCore."""
        mock_qtc = mocker.MagicMock(QT_VERSION_STR='7.0.0')
        mock_import_module(**{'PyQt7.QtCore': mock_qtc})
        with pytest.raises(NotImplementedError):
            _QtUnixPlatformsGetter.get('PyQt7')

    @pytest.mark.usefixtures('mock_supported_version')
    def test_qt_qpa_platform(self, mocker):
        """Check result when QT_QPA_PLATFORM is set."""
        mock_qpa = str(mocker.MagicMock())
        mocker.patch.dict('os.environ',
                          QT_QPA_PLATFORM=mock_qpa,
                          clear=True)
        result = _QtUnixPlatformsGetter.get(PYQT5)
        assert result == {mock_qpa}

    @pytest.mark.usefixtures('mock_supported_version')
    def test_qt_qpa_platform_list(self, mocker):
        """Check result when QT_QPA_PLATFORM lists several platforms."""
        mocker.patch.dict('os.environ', clear=True,
                          QT_QPA_PLATFORM='wayland:arg=1;xcb')
        result = _QtUnixPlatformsGetter.get(PYQT5)
        assert result == {'wayland', 'xcb'}

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
        result = _QtUnixPlatformsGetter.get(PYQT5)
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
        result = _QtUnixPlatformsGetter.get(PYQT5)
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
        result = _QtUnixPlatformsGetter.get(PYQT5)
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
        result = _QtUnixPlatformsGetter.get(PYQT5)
        assert result == {default, 'xcb'}
