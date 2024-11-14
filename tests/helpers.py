"""Module helpers of viperleed.tests.

Contains some useful general definitions that can be used when creating
or running tests.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-28'
__license__ = 'GPLv3+'

from contextlib import contextmanager
import copy
from dataclasses import dataclass, fields
import functools
import inspect
import os
from pathlib import Path

import pytest
from pytest_cases import case as case_decorator
from pytest_cases import fixture
from pytest_cases.case_funcs import CASE_FIELD
from pytest_cases.case_funcs import get_case_tags
from pytest_cases.case_funcs import is_case_class
from pytest_cases.case_funcs import is_case_function
from pytest_cases.filters import CaseFilter


# Think about a decorator for injecting fixtures.
# Some ideas at
# https://github.com/pytest-dev/pytest/issues/2424
# https://github.com/pytest-dev/pytest/issues/6322
# https://github.com/nteract/testbook/issues/4


TEST_DATA = Path(__file__).parent / '_test_data'
POSCAR_PATH = TEST_DATA / 'POSCARs'

# ##############################   EXCEPTIONS   ###############################

class CustomTestException(BaseException):
    """A custom exception for checking try...except blocks."""


# ##############################   DECORATORS   ###############################

def fixture_factory(*args, **fixture_kwargs):                                   # TODO: does not work for a callable with arguments
    """Return a decorator for fixtures, turning functions into factories.

    Decorating a function with this decorator will turn its body
    into a fixture factory. This means that the decorated functions
    can then be used in tests like so:

    @fixture_factory       # or @fixture(scope=..., name=...)
    def fixture(arguments):
        # ... do stuff ...
        return values

    def test(fixture):
        return_value = fixture(arguments)

    Parameters
    ----------
    *args : object
        Positional arguments. Only one positional argument is
        acceptable, i.e., the function to be turned into a
        fixture factory. If a single argument is given, there
        must be no fixture_kwargs.
    **fixture_kwargs : object
        The keyword arguments to give to pytest.fixture.
        The default is scope='session'.

    Returns
    -------
    decorator
    """

    def _fixture_wrapper(func):
        # We need an extra layer of wrapping because pytest will call
        # the fixture and pass its return value(s) to tests using it:
        # The one that will be called by pytest is _fixture.
        if 'scope' not in fixture_kwargs:
            fixture_kwargs['scope'] = 'session'
        if 'name' not in fixture_kwargs:
            fixture_kwargs['name'] = f'make_{func.__name__}'

        @fixture(**fixture_kwargs)
        @functools.wraps(func)
        def _fixture(*__args, **__kwargs):
            def _factory(*func_args, **func_kwargs):
                return func(*func_args, **func_kwargs)
            return _factory
        return _fixture

    if len(args) > 1:
        raise ValueError
    if len(args) == 1 and (fixture_kwargs or not callable(args[0])):
        raise ValueError
    if len(args) == 1:
        # Decorated as @fixture_factory
        return _fixture_wrapper(args[0])
    # Decorated as @fixture_factory(**kwargs)
    return _fixture_wrapper


def with_case_tags(*tags):
    """Attach `tags` to all cases defined in the decorated class."""
    def _decorator(cls):
        if is_case_function(cls):
            raise ValueError(
                'Cannot use with_case_tags on a case '
                'function. Use the @case decorator instead'
                )
        if not is_case_class(cls):
            raise ValueError('with_case_tags can only be applied to classes '
                             'defining a collection of cases')
        for case_name in dir(cls):
            case_ = getattr(cls, case_name)
            if not is_case_function(case_):  # Not a case
                continue
            try:
                case_info = getattr(case_, CASE_FIELD)
            except AttributeError:
                # Not explicitly decorated with @case. Do so now.
                case_ = case_decorator(case_)
                case_info = getattr(case_, CASE_FIELD)
            tags_to_add = tuple(t for t in tags if t not in case_info.tags)
            case_info.add_tags(tags_to_add)
        return cls
    return _decorator


# ##############################   FUNCTIONS   ################################

def duplicate_all(*sequence):
    """Return a deepcopy of all arguments passed."""
    return [copy.deepcopy(item) for item in sequence]


def inject_fixture(factory, *args, **kwargs):
    """Inject the return of factory in the caller's namespace.

    Parameters
    ----------
    factory : callable
        Should return a string and a fixture. The former
        is used as the name of the fixture to be injected
    *args : object
        Positional arguments for fixture_factory
    **kwargs : object
        Keyword arguments for fixture_factory

    Returns
    -------
    name : str
        The named of the injected fixture
    fixture : pytest.fixture
        The injected fixture
    """
    caller_globals = inspect.stack()[1][0].f_globals
    name, fixture_func = factory(*args, **kwargs)
    caller_globals[name] = fixture_func
    return name, fixture_func


def filesystem_from_dict(as_dict, root):
    """Create files and directories from a dictionary version of a tree.

    Parameters
    ----------
    as_dict : dict
        A dictionary representation of a file-system tree. Directories
        are dictionaries, whose keys are (relative) paths. File
        contents are dictionary values: they may be strings or
        None. In the latter case empty files are created.
    root : Path
        The path to the real file-system root directory where `as_dict`
        contents should be created.

    Returns
    -------
    created : set of Path
        Paths to all files/folders created.
    """
    created = set()
    for entry, contents in as_dict.items():
        path = root / entry
        created.add(path)
        # Directory:
        if isinstance(contents, dict):
            path.mkdir(parents=True, exist_ok=True)
            sub_paths = filesystem_from_dict(contents, path)
            created.update(sub_paths)
            continue
        # File:
        if contents is None:
            path.touch()
        else:
            path.write_text(contents, encoding='utf-8')
    return created


def flat_fixture(func, **fixture_args):                                         # TODO: better name. Only usable for parameter-less functions
    """Turn a function into a fixture with fixture_args.

    For use without decoration, like so:
    >>> fixture_from_func = flat_fixture(func, name='some_name')

    Parameters
    ----------
    func : function
        The function to be turned into a fixture.
    **fixture_args : object
        Argument to pass to the @fixture decorator.

    Returns
    -------
    fixture
    """
    _decorator = fixture(**fixture_args)
    return _decorator(func)


# ######################   CONTEXT MANAGERS   #########################

@contextmanager
def execute_in_dir(path):
    """Safely execute code in a specific directory."""
    home = Path().resolve()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(home)


@contextmanager
def not_raises(exc):
    """Fail a test if a specific exception is raised."""
    # Exclude this function when reporting the exception trace
    __tracebackhide__ = True  # pylint: disable=unused-variable
    try:
        yield
    except exc:
        pytest.fail(f'DID RAISE {exc.__name__}')


@contextmanager
def make_obj_raise(obj_or_dotted_name, exc, attr=None):
    """Temporarily make obj_or_dotted_name.attr raise exc.

    This is a wrapper around the pytest.monkeypatch fixture used
    as a context.

    Parameters
    ----------
    obj_or_dotted_name : object or str
        The object whose attribute should be set. May also be
        the dotted name of a module[.submodules][.class].attribute.
    exc : BaseException
        The exception to be raised.
    attr : str, optional
        The name of the attribute to be modified. May be omitted if
        a full dotted name is given as the first argument.

    Yields
    ------
    patch : fixture
        The pytest monkeypatch context.
    """
    # Exclude this function when reporting the exception trace
    __tracebackhide__ = True  # pylint: disable=unused-variable
    if inspect.ismodule(obj_or_dotted_name):
        obj_name = obj_or_dotted_name.__name__
    elif isinstance(obj_or_dotted_name, str):  # Likely a dotted name
        obj_name = obj_or_dotted_name
    elif inspect.isclass(obj_or_dotted_name):
        obj_name = obj_or_dotted_name.__name__
    else:
        obj_name = type(obj_or_dotted_name).__name__

    _err = ('When using the one-argument form, the argument must be a string '
            'corresponding to the dotted name of the attribute to set.')
    if attr is None and not isinstance(obj_or_dotted_name, str):
        raise TypeError(_err)
    if attr is not None and not isinstance(attr, str):
        raise TypeError('attribute name to be patched must be a string')
    # pylint: disable-next=magic-value-comparison  # OK for the dot
    if isinstance(obj_or_dotted_name, str) and not '.' in obj_or_dotted_name:
        raise ValueError(_err)
    if not issubclass(exc, BaseException):
        raise TypeError('exc must be a BaseException subclass')

    def _replaced(*args, **kwargs):
        patched = f'Replaced {obj_name}'
        if attr is not None:
            patched += f'.{attr}'
        raise exc(f'{patched} was called with args={args}, kwargs={kwargs}')

    monkeypatch = pytest.MonkeyPatch
    with monkeypatch.context() as patch:
        if attr is None:
            patch.setattr(obj_or_dotted_name, _replaced)
        else:
            patch.setattr(obj_or_dotted_name, attr, _replaced)
        yield patch


@contextmanager
def raises_exception(obj_or_dotted_name, exc, attr=None):
    """Make obj_or_dotted_name.attr raise exc, and ensure code raises it.

    This context manager is commonly used to test that a certain
    piece of code correctly catches specific exceptions.

    Parameters
    ----------
    obj_or_dotted_name : object or str
        The object whose attribute should be set. May also be
        the dotted name of a module[.submodules][.class].attribute.
    exc : BaseException
        The exception to be raised.
    attr : str, optional
        The name of the attribute to be modified. May be omitted if
        a full dotted name is given as the first argument.

    Yields
    ------
    patch : fixture
        The pytest monkeypatch context.
    raises : fixture
        The pytest.raises context.
    """
    # Exclude this function when reporting the exception trace
    __tracebackhide__ = True  # pylint: disable=unused-variable
    obj = obj_or_dotted_name  # Just to have 'with' fit one line
    with make_obj_raise(obj, exc, attr) as patch, pytest.raises(exc) as raises:
        yield patch, raises


@contextmanager
def raises_test_exception(obj_or_dotted_name, attr=None):
    """Temporarily make obj.attr raise CustomTestException.

    This context manager is commonly used to test that a certain
    piece of code does not use catch-all try...except blocks in
    either try...except or try...except Exception forms.

    Parameters
    ----------
    obj_or_dotted_name : object or str
        The object whose attribute should be set. May also be
        the dotted name of a module[.submodules][.class].attribute.
    attr : str, optional
        The name of the attribute to be modified. May be omitted if
        a full dotted name is given as the first argument.

    Yields
    ------
    patch : fixture
        The pytest monkeypatch context.
    raises : fixture
        The pytest.raises context.
    """
    # Exclude this function when reporting the exception trace
    __tracebackhide__ = True  # pylint: disable=unused-variable
    exc = CustomTestException
    with raises_exception(obj_or_dotted_name, exc, attr) as context:
        yield context


# ###########################   FILTERS   #############################

def exclude_tags(*tags):
    """Return a filter that excludes cases with given tags."""
    def _filter(case):
        """Return False if case has any of the tags."""
        case_tags = get_case_tags(case)
        return not any(tag in case_tags for tag in tags)
    return CaseFilter(_filter)


def has_any_tag(*tags):
    """Return a filter selecting cases with at least one of these tags."""
    def _filter(case):
        case_tags = set(get_case_tags(case))
        return any(t in case_tags for t in tags)

    if not tags:
        raise ValueError('Must select at least one tag')
    return CaseFilter(_filter)


# ###########################   CLASSES   #############################

@dataclass
class InfoBase:
    """Base class for all test-information classes."""

    def __iter__(self):
        """Yield values of the fields of this InfoBase."""
        for field_ in fields(self):
            yield getattr(self, field_.name)

    def __repr__(self):
        """Return a string representation of this InfoBase."""
        non_empty = []
        for field_ in fields(self):
            value = getattr(self, field_.name)
            if value is None:
                continue
            if (not isinstance(value, bool)
                    and not isinstance(value, InfoBase) and not value):
                continue
            if isinstance(value, InfoBase) and value.empty():
                continue
            non_empty.append((field_.name, value))
        return (f'{self.__class__.__name__}('
                + ', '.join(f'{n}={v}' for n, v in non_empty)
                + ')')

    def set(self, *args, **kwargs):
        """Set a bunch of values without recreating self."""
        for value, field_ in zip(args, fields(self)):
            setattr(self, field_.name, value)
        for name, value in kwargs.items():
            if not hasattr(self, name):
                raise AttributeError(name)
            setattr(self, name, value)
        try:
            self.__post_init__()
        except AttributeError:
            pass

    def empty(self):
        """Return whether any of the fields have truethy values."""
        return not any(self)
