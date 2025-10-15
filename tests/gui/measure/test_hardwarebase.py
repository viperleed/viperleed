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


@pytest.fixture
def fake_pkg(mocker):
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
    class A:
        pass
    A.__module__ = sub_a.__name__
    sub_a.A = A
    pkg.add_submodule('sub_a', sub_a)

    return pkg


# Tests for import_with_sub_modules
def test_import_with_sub_modules_yields_main_and_submodules(fake_pkg):
    result = list(import_with_sub_modules('fakepkg'))
    names = [m.__name__ for m in result]
    assert 'fakepkg' in names
    assert 'fakepkg.sub_a' in names


def test_import_with_sub_modules_invalid_import(mocker):
    mocker.patch('importlib.import_module', side_effect=ImportError)
    with pytest.raises(AttributeError):
        tmp = list(import_with_sub_modules('nonexistent'))


# Tests for class_from_name
def test_class_from_name_finds_class(fake_pkg):
    cls = class_from_name('fakepkg', 'A')
    assert cls.__name__ == 'A'


def test_class_from_name_no_class_found(fake_pkg):
    with pytest.raises(ValueError):
        class_from_name('fakepkg', 'NotExist')


def test_class_from_name_duplicate_classes(fake_pkg):
    # Add submodule that defines another A
    sub_b = types.ModuleType('fakepkg.sub_b')
    class A2:
        pass
    A2.__name__ = 'A'
    A2.__module__ = sub_b.__name__
    sub_b.A = A2
    fake_pkg.add_submodule('sub_b', sub_b)

    with pytest.raises(ImportError):
        class_from_name('fakepkg', 'A')


# Tests for get_devices
class DummyDevice:
    def __init__(self, name):
        self.unique_name = name


def test_get_devices_collects_all_devices(fake_pkg):
    driver_a = types.ModuleType('fakepkg.driver_a')

    class DummyDevice:
        def __init__(self, name):
            self.unique_name = name

    class DummyDeviceClass:
        def list_devices(self):
            return [DummyDevice('device_a')]

    DummyDeviceClass.__module__ = driver_a.__name__
    driver_a.DummyDeviceClass = DummyDeviceClass
    fake_pkg.add_submodule('driver_a', driver_a)

    devices = get_devices('fakepkg')
    assert 'device_a' in devices
    cls, dev = devices['device_a']
    assert cls is DummyDeviceClass
    assert dev.unique_name == 'device_a'


def test_get_devices_ignores_non_device_classes(fake_pkg):
    sub = types.ModuleType('fakepkg.drivers')

    # Add abstract class with list_devices. 'A' already is a
    # non-abstract class without 'list_devices', so we do not
    # need to add a non-device class.
    class AbstractDevice(metaclass=abc.ABCMeta):
        @abc.abstractmethod
        def list_devices(self):
            pass

    AbstractDevice.__module__ = sub.__name__
    sub.AbstractDevice = AbstractDevice
    fake_pkg.add_submodule('drivers', sub)

    devices = get_devices('fakepkg')
    assert devices == {}
