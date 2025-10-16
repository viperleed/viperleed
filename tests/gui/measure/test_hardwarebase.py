"""Tests for module hardwarebase of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-10-13'
__license__ = 'GPLv3+'

import abc
import types

import pytest

from viperleed.gui.measure.hardwarebase import class_from_name
from viperleed.gui.measure.hardwarebase import get_devices
from viperleed.gui.measure.hardwarebase import import_with_sub_modules


# pylint: disable=magic-value-comparison, unused-argument
def in_module(module):
    """Register a class to belong to a module."""
    def _register_cls(cls):
        cls.__module__ = module.__name__
        setattr(module, cls.__name__, cls)
        return cls
    return _register_cls


@pytest.fixture(name='fake_pkg')
def fixture_fake_pkg(mocker):
    """Create a fake package structure like:
    fakepkg/
        sub_a.py
    and mock sys.modules and pkgutil.iter_modules so that
    import_with_sub_modules and related functions can discover them.
    """
    base = 'fakepkg'
    full_base = 'viperleed.gui.measure.fakepkg'

    pkg = types.ModuleType(base)
    pkg.__path__ = ['dummy_path']
    pkg.__package__ = full_base

    mocker.patch.dict('sys.modules', {base: pkg, full_base: pkg})
    submodules = {}

    def add_submodule(name, module):
        """Add submodule to fake package."""
        mocker.patch.dict('sys.modules', {f'{base}.{name}': module,
                                          f'{full_base}.{name}': module})
        submodules[name] = module

    # Mock pkgutil.iter_modules to reflect whatever submodules exist.
    def fake_iter_modules(path):
        """Iterate through fake modules."""
        return [types.SimpleNamespace(name=n) for n in submodules]

    mocker.patch('pkgutil.iter_modules', side_effect=fake_iter_modules)

    pkg.add_submodule = add_submodule
    pkg.submodules = submodules

    # Define submodule containing a class DummyA.
    sub_a = types.ModuleType('fakepkg.sub_a')
    @in_module(sub_a)
    class DummyA:   # pylint: disable=too-few-public-methods, unused-variable
        """Dummy class."""
    pkg.add_submodule('sub_a', sub_a)   # pylint: disable=no-member
    return pkg


class TestImportWithSubModules:
    """Tests for the import_with_submodules iterator."""

    def test_invalid_import(self, mocker):
        """Test import error."""
        mocker.patch('importlib.import_module', side_effect=ImportError)
        with pytest.raises(AttributeError):
            _ = list(import_with_sub_modules('nonexistent'))

    def test_yield_main_and_submodule(self, fake_pkg):
        """Test import of package and submodule."""
        names = [m.__name__ for m in import_with_sub_modules('fakepkg')]
        expect_names = [
            'fakepkg',        # First the top-level module
            'fakepkg.sub_a',  # Then all its submodules
            ]
        assert names == expect_names

    def test_yield_module_only(self, fake_pkg):
        """Test import of single module."""
        delattr(fake_pkg, '__path__')
        names = [m.__name__ for m in import_with_sub_modules('fakepkg')]
        expect_names = ['fakepkg',]
        assert names == expect_names

    def test_yield_recursive(self, fake_pkg, mocker):
        """Test recursive import of modules."""
        sub_a = fake_pkg.submodules['sub_a']
        sub_a.__path__ = ['dummy_sub_path']
        sub_a.__package__ = 'fakepkg.sub_a'
        nested = types.ModuleType('fakepkg.sub_a.nested_sub')
        fake_pkg.add_submodule('sub_a.nested_sub', nested)

        mocker.patch(
            'pkgutil.iter_modules',
            side_effect=[
                [types.SimpleNamespace(name='sub_a')],
                [types.SimpleNamespace(name='nested_sub')]
                ],
            )

        names = [m.__name__ for m in import_with_sub_modules('fakepkg',
                                                             recursive=True)]
        expect_names = [
            'fakepkg',
            'fakepkg.sub_a',
            'fakepkg.sub_a.nested_sub',   # Recursively imported module
            ]
        assert names == expect_names


class TestClassFromName:
    """Tests for the class_from_name function."""

    def test_find_class(self, fake_pkg):
        """Test finding class."""
        cls = class_from_name('fakepkg', 'DummyA')
        assert cls.__name__ == 'DummyA'

    def test_no_class_found(self, fake_pkg):
        """Test finding class failure."""
        with pytest.raises(ValueError):
            class_from_name('fakepkg', 'NotExist')

    def test_duplicate_classes(self, fake_pkg):
        """Test duplicate class declaration."""
        # Add submodule that defines another DummyA
        sub_b = types.ModuleType('fakepkg.sub_b')
        @in_module(sub_b)
        class DummyA: # pylint: disable=too-few-public-methods, unused-variable
            """Dummy class."""
        fake_pkg.add_submodule('sub_b', sub_b)

        with pytest.raises(RuntimeError):
            class_from_name('fakepkg', 'DummyA')

    def test_imported_class_does_not_count_as_duplicate(self, fake_pkg):
        """Test reimport of class."""
        # Import DummyA from sub_a into sub_b
        sub_a = fake_pkg.submodules['sub_a']
        sub_b = types.ModuleType('fakepkg.sub_b')

        # Re-export DummyA (imported reference)
        sub_b.DummyA = sub_a.DummyA
        fake_pkg.add_submodule('sub_b', sub_b)

        cls = class_from_name('fakepkg', 'DummyA')
        assert cls is sub_a.DummyA
        assert cls.__module__ == sub_a.__name__


class TestGetDevices:
    """Tests for the get_devices function."""

    def test_device_getting(self, fake_pkg):
        """Test getting devices."""
        driver_a = types.ModuleType('fakepkg.driver_a')

        class DummyDeviceInfo:  # pylint: disable=too-few-public-methods
            """Dummy device info."""
            def __init__(self, name):
                """Init dummy device info."""
                self.unique_name = name

        @in_module(driver_a)
        class DummyDevice:      # pylint: disable=too-few-public-methods
            """Dummy dwevice."""
            def list_devices(self):
                """Dummy list devices."""
                return [DummyDeviceInfo('device_a')]
        fake_pkg.add_submodule('driver_a', driver_a)

        devices = get_devices('fakepkg')
        assert len(devices) == 1
        assert 'device_a' in devices
        cls, dev = devices['device_a']
        assert cls is DummyDevice
        assert dev.unique_name == 'device_a'

    def test_ignore_non_device_classes(self, fake_pkg):
        """Test ignoring non-devices."""
        sub = types.ModuleType('fakepkg.drivers')

        # Add abstract class with list_devices. 'DummyA' already is
        # a non-abstract class without 'list_devices', so we do not
        # need to add a non-device class.
        @in_module(sub)
        # pylint: disable=too-few-public-methods, unused-variable
        class AbstractDevice(metaclass=abc.ABCMeta):
            """Dummy device ABC."""
            @abc.abstractmethod
            def list_devices(self):
                """Dummy list devices."""
        fake_pkg.add_submodule('drivers', sub)

        devices = get_devices('fakepkg')
        assert not devices
