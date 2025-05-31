"""Tests for module cli_base of viperleed."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-19'
__license__ = 'GPLv3+'

from argparse import ArgumentParser
from argparse import ArgumentTypeError
from argparse import Namespace
import logging
from pathlib import Path
import re
import sys

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.cli_base import NotAViPErLEEDCLIError
from viperleed.cli_base import StreamArgument
from viperleed.cli_base import ViPErLEEDCLI
from viperleed.cli_base import ViPErLEEDCLIWithAutoChildren
from viperleed.cli_base import float_in_zero_one
from viperleed.cli_base import length_choices
from viperleed.cli_base import maximum_length
from viperleed.cli_base import minimum_length
from viperleed.cli_base import positive_float
from viperleed.cli_base import required_length
from viperleed.cli_base import strip_cli_module
from viperleed.cli_base import _to_float


class TestToFloat:
    """Tests for the _to_float helper."""

    def test_valid_float(self):
        """Check conversion of a valid float string."""
        # pylint: disable-next=magic-value-comparison
        assert _to_float('3.14') == 3.14

    def test_invalid_float(self):
        """Check complaints when converting an invalid float string."""
        with pytest.raises(ArgumentTypeError):
            _to_float('not_a_float')


class TestFloatInZeroOne:
    """Tests for the float_in_zero_one function."""

    def test_valid_range(self):
        """Check conversion of a valid float string."""
        # pylint: disable-next=magic-value-comparison
        assert float_in_zero_one('0.5') == 0.5

    _invalid = {
        'too small': '0',
        'too large': '1',
        'outside range': '1.1',
        }

    @parametrize(str_value=_invalid.values(), ids=_invalid)
    def test_raises(self, str_value):
        """Check complaints for an out-of-range `str_value`."""
        with pytest.raises(ArgumentTypeError):
            float_in_zero_one(str_value)


class TestLengthChoices:
    """Tests for the length_choices action factory."""

    def test_valid_choices(self):
        """Check parsing a valid number of arguments."""
        parser = ArgumentParser()
        parser.add_argument('--numbers',
                            nargs='+',
                            action=length_choices(1, 2))
        args = parser.parse_args(['--numbers', '10', '20'])
        assert args.numbers == ['10', '20']

        args = parser.parse_args(['--numbers', '10'])
        assert args.numbers == ['10']

    def test_invalid_choice(self):
        """Test complaints when parsing the wrong number of arguments."""
        parser = ArgumentParser()
        parser.add_argument('--numbers',
                            nargs='+',
                            action=length_choices(2, 3))
        with pytest.raises(SystemExit):
            parser.parse_args(['--numbers', '10', '20', '30', '40'])
        with pytest.raises(SystemExit):
            parser.parse_args(['--numbers', '10'])


class TestMaximumLength:
    """Tests for the maximum_length action factory."""

    action = maximum_length

    def test_in_range(self):
        """Check the successful parsing of an acceptable number of values."""
        n_max = 3
        parser = ArgumentParser()
        action = type(self).action  # Use class: no self when calling
        parser.add_argument('--values', nargs='+', action=action(n_max=n_max))
        for n_args in range(1, n_max):
            values = ['1']*n_args
            args = parser.parse_args(['--values', *values])
            assert args.values == values

    def test_too_many(self):
        """Check complaints when too many arguments are given."""
        n_max = 3
        parser = ArgumentParser()
        action = type(self).action  # Use class: no self when calling
        parser.add_argument('--values', nargs='+', action=action(n_max=n_max))
        values = ['1']*4
        with pytest.raises(SystemExit):
            parser.parse_args(['--values', *values])


class TestMinimumLength:
    """Tests for the minimum_length action factory."""

    action = minimum_length

    def test_in_range(self):
        """Check the successful parsing of an acceptable number of values."""
        n_min = 3
        parser = ArgumentParser()
        action = type(self).action  # Use class: no self when calling
        parser.add_argument('--values', nargs='+', action=action(n_min=n_min))
        for n_args in range(n_min, 5):
            values = ['1']*n_args
            args = parser.parse_args(['--values', *values])
            assert args.values == values

    def test_too_few(self):
        """Check complaints when too few arguments are given."""
        n_min = 3
        parser = ArgumentParser()
        action = type(self).action  # Use class: no self when calling
        parser.add_argument('--values', nargs='+', action=action(n_min=n_min))
        values = ['1']*2
        with pytest.raises(SystemExit):
            parser.parse_args(['--values', *values])


class TestPositiveFloat:
    """Tests for the positive_float function."""

    def test_positive_value(self):
        """Check the conversion of a positive float string."""
        # pylint: disable-next=magic-value-comparison
        assert positive_float('10.5') == 10.5

    _invalid = {
        'negative': '-5',
        'zero': '0',
        }

    @parametrize(str_value=_invalid.values(), ids=_invalid)
    def test_raises(self, str_value):
        """Check complaints for an out-of-range `str_value`."""
        with pytest.raises(ArgumentTypeError):
            positive_float(str_value)


# pylint: disable-next=too-few-public-methods   # Bug?
class TestRequiredLength:
    """Tests for the required_length action factory."""

    action = required_length

    test_above_max_length = TestMaximumLength.test_too_many
    test_below_min_length = TestMinimumLength.test_too_few
    test_valid_min_length = TestMinimumLength.test_in_range
    test_valid_max_length = TestMaximumLength.test_in_range


class TestStreamArgument:
    """Tests for the StreamArgument argument type."""

    _terminal = ('stdin', 'stdout', 'stderr')

    def test_file_not_found(self):
        """Check complaints when given an invalid file path."""
        stream = StreamArgument('r')('does-not-exist.txt')
        with pytest.raises(FileNotFoundError):
            with stream:
                pytest.fail('Should not reach this point')

    def test_invalid_path(self):
        """Check complaints when given an invalid file path."""
        stream = StreamArgument('r')
        with pytest.raises(ArgumentTypeError, match='Cannot open'):
            stream(12345)  # Invalid type

    @parametrize(interactive=(True, False))
    @parametrize(terminal=_terminal)
    def test_is_interactive(self, interactive, terminal, mocker):
        """Check the is_interactive property."""
        mock_resource = mocker.MagicMock(closed=False)
        mock_resource.isatty.return_value = interactive
        mocker.patch(f'sys.{terminal}', mock_resource)
        stream = StreamArgument('r')(mock_resource)
        assert stream.is_interactive == interactive

    def test_is_stream_detection(self):
        """Check correct detection of the stream nature of __call__ arg."""
        # pylint: disable=protected-access                # OK in tests
        assert StreamArgument._is_stream(sys.stdout)
        assert not StreamArgument._is_stream('not_a_stream')

    @parametrize(terminal=_terminal)
    def test_terminal_stream(self, terminal, mocker):
        """Check correct (not) opening/closing of a stream to the terminal."""
        sys_terminal = mocker.patch(f'sys.{terminal}')
        stream = StreamArgument('r')(sys_terminal)
        assert stream.is_terminal
        with stream as open_stream:
            assert open_stream is sys_terminal
        sys_terminal.open.assert_not_called()

    def test_open_read_context(self, tmp_path):
        """Check reading from a file path."""
        expect_read = 'content'
        file_path = tmp_path / 'test.txt'
        file_path.write_text(expect_read)
        open_ = StreamArgument('r')
        with open_(file_path) as file:
            assert file.read() == expect_read

    def test_open_write_context(self, tmp_path):
        """Check writing to a file path."""
        expect_write = 'new content'
        file_path = tmp_path / 'test_write.txt'
        open_ = StreamArgument('w')
        with open_(file_path) as file:
            file.write(expect_write)
        assert file_path.read_text() == expect_write


class TestStripCLIModule:
    """Tests for the strip_cli_module function."""

    def test_strip_cli(self):
        """Check removal of a 'cli' suffix."""
        # pylint: disable-next=magic-value-comparison
        assert strip_cli_module('viperleed.cli') == 'viperleed'

    def test_no_cli_suffix(self):
        """Check nothing is removed if module does not end with 'cli'."""
        # pylint: disable-next=magic-value-comparison
        assert strip_cli_module('viperleed.module') == 'viperleed.module'


def _make_cli_cls(**cls_args):
    """Return two subclasses of ViPErLEEDCLI as parent and child."""
    parent = type('ParentCLI', (ViPErLEEDCLI,), {})
    child = type('ChildCLI', (ViPErLEEDCLI,), {}, **cls_args)
    parent.register_child(child)
    return parent, child


@fixture(name='make_cli_cls')
def factory_cli_cls():
    """Return two subclasses of ViPErLEEDCLI."""
    return _make_cli_cls


class TestViPErLEEDCLI:
    """Tests for the ViPErLEEDCLI class."""

    def test_add_alias(self, make_cli_cls, mocker):
        """Check addition of aliases of a CLI."""
        parent_cls, child_cls = make_cli_cls(cli_name='test')
        child_call = mocker.patch.object(child_cls, '__call__')
        cli = parent_cls()
        cli.add_child_aliases('test', 'alias')
        cli(['alias'])
        child_call.assert_called_once()

    def test_call_without_arguments(self, make_cli_cls, capsys, mocker):
        """Check result of calling a (child) CLI without arguments."""
        mocker.patch('sys.argv', ['test'])
        cli_cls, _ = make_cli_cls(cli_name='test')
        cli = cli_cls()
        with pytest.raises(SystemExit):
            cli()
        captured = capsys.readouterr()
        assert re.match(r'usage:.*\[-h\] \[--version\].*', captured.out)
        assert not captured.err

    _mock_cli, _ = _make_cli_cls()
    _members = {  # ((name_in_module, obj_in_module), exc_or_result)
        'not a class': (('not_a_class', '123'), NotAViPErLEEDCLIError),
        'cli base': (('ViPErLEEDCLI', ViPErLEEDCLI), NotAViPErLEEDCLIError),
        'cli base, other': (
            ('ViPErLEEDCLIWithAutoChildren', ViPErLEEDCLIWithAutoChildren),
            NotAViPErLEEDCLIError,
            ),
        'not a cli class': (('OtherClass', tuple), NotAViPErLEEDCLIError),
        'private cli': (('_Private', _mock_cli), NotAViPErLEEDCLIError),
        'valid cli': (('ValidCLI', _mock_cli), _mock_cli),
        }

    @parametrize('module_member,expect', _members.values(), ids=_members)
    def test_child_cls_from_module(self, module_member, expect, mocker):
        """Check fetching a CLI class from a module name."""
        mocker.patch('importlib.import_module')
        mocker.patch('inspect.getmembers', return_value=(module_member,))
        # pylint: disable-next=protected-access           # OK in tests
        get = ViPErLEEDCLI._child_class_from_module_name
        if expect is NotAViPErLEEDCLIError:
            with pytest.raises(expect):
                get('some_module')
        else:
            assert get('some_module') == expect

    def test_cli_name_default(self, make_cli_cls):
        """Check default assignment of a CLI name."""
        _, cli_cls = make_cli_cls()
        assert cli_cls().cli_name is not None  # From module

    def test_cli_name_custom(self, make_cli_cls):
        """Check custom assignment of a CLI name."""
        name = 'custom_name'
        _, cli_cls = make_cli_cls(cli_name=name)
        assert cli_cls().cli_name is name

    def test_dont_parse_twice(self, make_cli_cls):
        """Check that arguments are not re-parsed."""
        cli_cls, _ = make_cli_cls()
        already_parsed = Namespace()
        parsed = cli_cls().parse_cli_args(already_parsed)
        assert parsed is already_parsed

    def test_get_logger(self, make_cli_cls):
        """Check retrieval of a CLI's logger."""
        cli_cls, _ = make_cli_cls()
        logger = cli_cls.get_logger()
        assert isinstance(logger, logging.Logger)

    def test_keyboard_interrupt(self, make_cli_cls, capsys, mocker):
        """Check graceful handling of KeyboardInterrupt."""
        mocker.patch('sys.exit')
        mocker.patch('sys.argv', [])
        cli_cls, _ = make_cli_cls()
        mocker.patch.object(cli_cls, '__call__', side_effect=KeyboardInterrupt)
        cli_cls.run_as_script()
        captured = capsys.readouterr()
        sys.exit.assert_called_once()
        # pylint: disable-next=magic-value-comparison
        assert 'keyboard' in captured.err

    def test_parse_cli_args(self, make_cli_cls):
        """Check parsing of valid CLI arguments."""
        parent_cli, _ = make_cli_cls()
        class _ChildCLI(ViPErLEEDCLI, cli_name='test'):
            # We can't use self as the first argument, as it
            # would conflict with the 'self' of this test.
            # pylint: disable-next=no-self-argument,arguments-renamed
            def add_parser_arguments(self_, parser):
                super().add_parser_arguments(parser)
                self_.add_verbose_option(parser)
        parent_cli.register_child(_ChildCLI)
        cli = parent_cli()
        parsed_args = cli.parse_cli_args(['test', '-v'])
        assert isinstance(parsed_args, Namespace)
        assert parsed_args.command == _ChildCLI().cli_name
        assert parsed_args.verbose
        func = parsed_args.func
        # NB: func == _ChildCLI().__call__ would fail, as the instance
        # whose __call__ we register as func is created freshly at each
        # add_parser_arguments call.
        # pylint: disable-next=unidiomatic-typecheck
        assert type(func.__self__) is _ChildCLI
        # pylint: disable-next=magic-value-comparison
        assert func.__name__ == '__call__'

    def test_print_version(self, make_cli_cls, capsys, mocker):
        """Check printing of version to the terminal."""
        mocker.patch('sys.argv', ['test', '--version'])
        cli_cls, _ = make_cli_cls(cli_name='test')
        cli = cli_cls()
        with pytest.raises(SystemExit):
            cli()
        captured = capsys.readouterr()
        assert re.match(r'ViPErLEED.*v[\d.]+', captured.out)
        assert not captured.err

    @staticmethod
    def _check_has_child(parent, child):
        """Check that a `parent` CLI class has a `child` class registered."""
        for child_name, child_cls in parent().children.items():
            if child().cli_name in child_name:
                assert child_cls is child
                return
        pytest.fail('child not found')

    def test_register_child_cls(self, make_cli_cls):
        """Check registering a child using a ViPErLEEDCLI subclass."""
        self._check_has_child(*make_cli_cls())

    def test_register_child_instance(self, make_cli_cls):
        """Check registering a child using a ViPErLEEDCLI instance."""
        parent, _ = make_cli_cls()
        _, child = make_cli_cls()
        parent.register_child(child())
        self._check_has_child(parent, child)

    _invalid_child = {
        'invalid_module': ValueError,
        1: NotAViPErLEEDCLIError,
        }

    @parametrize('invalid_child,exc', _invalid_child.items())
    def test_register_child_invalid(self, invalid_child, exc, make_cli_cls):
        """Check complaints when registering an invalid child."""
        parent, _ = make_cli_cls()
        with pytest.raises(exc):
            parent.register_child(invalid_child)

    def test_run_as_script(self, make_cli_cls, mocker):
        """Check a successful execution of a CLI as a main."""
        mocker.patch('sys.exit')
        mocker.patch('sys.argv', [])
        _, cli_cls = make_cli_cls()
        cli_cls.run_as_script()
        sys.exit.assert_called()


class TestFindSubUtilities:
    """Tests for the find_sub_utilities class method of ViPErLEEDCLI."""

    @fixture(name='fake_package')
    def fixture_mock_package(self, mocker):
        """Return a fake package."""
        fake_package = mocker.Mock(ispkg=True)
        fake_package.name = 'a_package'
        fake_package.module_finder.path = 'path/to/fake/package'
        return fake_package

    @fixture(name='fake_modules')
    def fixture_mock_modules(self, mocker):
        """Return a list of fake modules."""
        fake_modules = []
        module_names = 'a_module', '__main__', 'cli'
        for module_name in module_names:
            fake_module = mocker.Mock(ispkg=False)
            fake_module.name = module_name
            fake_modules.append(fake_module)
        return fake_modules

    @fixture(autouse=True)
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        mocker.patch('importlib.util.find_spec',
                     return_value=mocker.Mock(origin='/fake/path/cli.py'))
        mocker.patch('pathlib.Path.exists', return_value=True)

    @fixture(name='make_cli_with_module')
    def fixture_make_cli_with_module(self, make_cli_cls):
        """Return a CLI subclass with its __module__ set."""
        def _make(module_name):
            cli_cls, _ = make_cli_cls()
            cli_cls.__module__ = module_name
            return cli_cls
        return _make

    def test_module_and_package(self,
                                fake_package,
                                fake_modules,
                                make_cli_with_module,
                                mocker):
        """Check correct finding of modules and packages."""
        def _mock_iter_modules(modules):
            module, *_ = modules
            _fake_pkg_path = Path(fake_package.module_finder.path)
            _fake_pkg_path /= fake_package.name
            if Path(module) == _fake_pkg_path:
                submodule = mocker.Mock()
                submodule.name = 'cli'
                return (submodule, *fake_modules)
            return (*fake_modules, fake_package)

        mocker.patch('pkgutil.iter_modules', _mock_iter_modules)
        test_cls = make_cli_with_module('root.cli')
        result = test_cls.find_sub_utilities()
        assert result == ('root.a_package.cli', 'root.a_module')

    def test_no_cli_in_package(self,
                               fake_package,
                               fake_modules,
                               make_cli_with_module,
                               mocker):
        """Check no modules are found if no cli exists in a package."""
        def _mock_iter_modules(modules):
            module, *_ = modules
            _fake_pkg_path = Path(fake_package.module_finder.path)
            _fake_pkg_path /= fake_package.name
            if Path(module) == _fake_pkg_path:
                return fake_modules[:-1]
            return (fake_package,)

        mocker.patch('pkgutil.iter_modules', _mock_iter_modules)
        test_cls = make_cli_with_module('root.cli')
        result = test_cls.find_sub_utilities()
        assert not result

    def test_subclass_not_in_cli_module(self, make_cli_with_module, mocker):
        """Check no search is done if a class is not in a cli.py module."""
        iter_modules = mocker.patch('pkgutil.iter_modules')
        test_cls = make_cli_with_module('package.module.that.is.not.a_cli')
        result = test_cls.find_sub_utilities()
        assert not result
        iter_modules.assert_not_called()


def test_auto_register_children(mocker):
    """Check automatic recognition of children CLIs."""
    child_module = 'child_module'
    modules = [child_module, 'invalid_module']
    ChildCLI = type('ChildCLI', (ViPErLEEDCLI,), {})
    ChildCLI.__module__ = child_module
    def _mock_child_from_module(module):
        if module == child_module:
            return ChildCLI
        return 'not a cli subclass'
    mocker.patch.object(ViPErLEEDCLIWithAutoChildren,
                        'find_sub_utilities',
                        return_value=modules)
    mocker.patch.object(ViPErLEEDCLIWithAutoChildren,
                        '_child_class_from_module_name',
                        _mock_child_from_module)
    ParentCLI = type('ParentCLI', (ViPErLEEDCLIWithAutoChildren,), {})
    assert child_module in ParentCLI().children
