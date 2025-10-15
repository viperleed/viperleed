"""Tests for module hardwarebase of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-10-13'
__license__ = 'GPLv3+'

import abc
import sys
import types

import pytest

from viperleed.gui.measure.hardwarebase import class_from_name
from viperleed.gui.measure.hardwarebase import get_devices
from viperleed.gui.measure.hardwarebase import import_with_sub_modules


def in_module(module, forced_name=''):
    """Register a class to belong to a module."""
    def _register_cls(cls):
        if forced_name:
            cls.__name__ = forced_name
        cls.__module__ = module.__name__
        setattr(module, cls.__name__, cls)
        return cls
    return _register_cls


@pytest.fixture(name='fake_pkg')
def fixture_fake_pkg(mocker):
    """
    Create a fake package structure like:
    fakepkg/
        sub_a.py
    and mock importlib.import_module and pkgutil.iter_modules so that
    import_with_sub_modules and related functions can discover them.
    """
    base = 'fakepkg'
    full_base = 'viperleed.gui.measure.fakepkg'

    pkg = types.ModuleType(base)
    pkg.__path__ = ['dummy_path']
    pkg.__package__ = full_base

    sys.modules[base] = pkg
    sys.modules[full_base] = pkg

    submodules = {}

    def add_submodule(name, module):
        sys.modules[f'{base}.{name}'] = module
        sys.modules[f'{full_base}.{name}'] = module
        submodules[name] = module

    # Mock importlib.import_module.
    def fake_import(name, package=None):
        return sys.modules[name]

    mocker.patch('importlib.import_module', side_effect=fake_import)

    # Mock pkgutil.iter_modules to reflect whatever submodules exist.
    def fake_iter_modules(path):
        return [types.SimpleNamespace(name=n) for n in submodules]

    mocker.patch('pkgutil.iter_modules', side_effect=fake_iter_modules)

    pkg.add_submodule = add_submodule
    pkg.submodules = submodules

    # Define submodule containing a class A.
    sub_a = types.ModuleType('fakepkg.sub_a')
    @in_module(sub_a)
    class A:
        pass
    pkg.add_submodule('sub_a', sub_a)
    return pkg


class TestImportWithSubModules:
    """Tests for the import_with_submodules iterator."""

    def test_invalid_import(self, mocker):
        mocker.patch('importlib.import_module', side_effect=ImportError)
        with pytest.raises(AttributeError):
            _ = list(import_with_sub_modules('nonexistent'))

    def test_yield_main_and_submodule(self, fake_pkg):
        names = [m.__name__ for m in import_with_sub_modules('fakepkg')]
        expect_names = [
            'fakepkg',        # First the top-level module
            'fakepkg.sub_a',  # Then all its submodules
            ]
        assert names == expect_names

    def test_yield_module_only(self, fake_pkg):
        delattr(fake_pkg, '__path__')
        names = [m.__name__ for m in import_with_sub_modules('fakepkg')]
        expect_names = ['fakepkg',]
        assert names == expect_names


class TestClassFromName:
    """Tests for the class_from_name function."""

    def test_find_class(self, fake_pkg):
        cls = class_from_name('fakepkg', 'A')
        assert cls.__name__ == 'A'

    def test_no_class_found(self, fake_pkg):
        with pytest.raises(ValueError):
            class_from_name('fakepkg', 'NotExist')

    def test_duplicate_classes(self, fake_pkg):
        # Add submodule that defines another A
        sub_b = types.ModuleType('fakepkg.sub_b')
        @in_module(sub_b, forced_name='A')
        class A2:
            pass
        fake_pkg.add_submodule('sub_b', sub_b)

        with pytest.raises(ImportError):
            class_from_name('fakepkg', 'A')


class TestGetDevices:
    """Tests for the get_devices function."""

    def test_device_getting(self, fake_pkg):
        driver_a = types.ModuleType('fakepkg.driver_a')

        class DummyDeviceInfo:
            def __init__(self, name):
                self.unique_name = name

        @in_module(driver_a)
        class DummyDevice:
            def list_devices(self):
                return [DummyDeviceInfo('device_a')]
        fake_pkg.add_submodule('driver_a', driver_a)

        devices = get_devices('fakepkg')
        assert len(devices) == 1
        assert 'device_a' in devices
        cls, dev = devices['device_a']
        assert cls is DummyDevice
        assert dev.unique_name == 'device_a'

    def test_ignore_non_device_classes(self, fake_pkg):
        sub = types.ModuleType('fakepkg.drivers')

        # Add abstract class with list_devices. 'A' already is a
        # non-abstract class without 'list_devices', so we do not
        # need to add a non-device class.
        @in_module(sub)
        class AbstractDevice(metaclass=abc.ABCMeta):
            @abc.abstractmethod
            def list_devices(self):
                pass
        fake_pkg.add_submodule('drivers', sub)

        devices = get_devices('fakepkg')
        assert not devices
