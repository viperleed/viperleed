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
from pytest_cases import fixture
from pytest_cases.filters import get_case_tags


# Think about a decorator for injecting fixtures.
# Some ideas at
# https://github.com/pytest-dev/pytest/issues/2424
# https://github.com/pytest-dev/pytest/issues/6322
# https://github.com/nteract/testbook/issues/4


TEST_DATA = Path(__file__).parent / '_test_data'
POSCAR_PATH = TEST_DATA / 'POSCARs'

# ##############################   EXCEPTIONS   ###############################

class CustomTestException(Exception):
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


def exclude_tags(*tags):
    """Return a filter that excludes cases with given tags."""
    def _filter(case):
        """Return False if case has any of the tags."""
        case_tags = get_case_tags(case)
        return not any(tag in case_tags for tag in tags)
    return _filter


# #######################   CONTEX MANAGERS   #########################

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
def not_raises(exception):
    """Fail a test if a specific exception is raised."""
    # Exclude this function when reporting the exception trace
    __tracebackhide__ = True  # pylint: disable=unused-variable
    try:
        yield
    except exception:
        pytest.fail(f'DID RAISE {exception.__name__}')


@contextmanager
def raises_test_exception(obj, attr):
    """Temporarily make obj.attr raise CustomTestException."""
    # Exclude this function when reporting the exception trace
    __tracebackhide__ = True  # pylint: disable=unused-variable
    obj_name = (obj.__name__ if inspect.ismodule(obj)
                else type(obj).__name__)

    def _replaced(*args, **kwargs):
        raise CustomTestException(f'Replaced {obj_name}.{attr} was called '
                                  f'with args={args}, kwargs={kwargs}')

    monkey_patch = pytest.MonkeyPatch
    with monkey_patch.context() as patch, pytest.raises(CustomTestException):
        patch.setattr(obj, attr, _replaced)
        yield


# ###############################   CLASSES   #################################

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
