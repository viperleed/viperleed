"""Module helpers of viperleed.tests.

Created on 2023-02-28

@author: Michele Riva (@michele-riva)

Contains some useful general definitions that can be used when creating
or running tests.
"""

from contextlib import contextmanager
import copy
from dataclasses import dataclass, field, fields
from enum import IntEnum, auto
import functools
import inspect
import os
from pathlib import Path
from typing import Dict, List, Set, Tuple, Mapping

import numpy as np
import pytest
from pytest_cases import fixture
from pytest_cases.filters import get_case_tags

from viperleed.tleedmlib.classes.rparams import LayerCuts


# Think about a decorator for injecting fixtures.
# Some ideas at
# https://github.com/pytest-dev/pytest/issues/2424
# https://github.com/pytest-dev/pytest/issues/6322
# https://github.com/nteract/testbook/issues/4


TEST_DATA = Path(__file__).parent / '_test_data'
POSCAR_PATH = TEST_DATA / 'POSCARs'


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

# ###############################   CLASSES   #################################

class CaseTag(IntEnum):
    """Enumeration of tags to use for cases."""
    BULK = auto()
    BULK_PROPERTIES = auto()
    LAYER_INFO = auto()
    NEED_ROTATION = auto()
    NO_INFO = auto()
    NON_MINIMAL_CELL = auto()
    RAISES = auto()
    THICK_BULK = auto()


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


@dataclass(repr=False)
class POSCARInfo(InfoBase):
    """Container for information about POSCAR files."""
    name: str = ''
    n_atoms: int = None
    n_atoms_by_elem: dict = field(init=False, default_factory=dict)
    n_cells: int = 1  # How many 2D cells are in the POSCAR?

    def __post_init__(self):
        """Post-process initialization values."""
        if isinstance(self.n_atoms, Mapping):
            self.n_atoms_by_elem = self.n_atoms
            self.n_atoms = sum(self.n_atoms_by_elem.values())


# Avoid using atoms that are on planes: Easier to test.
@dataclass(repr=False)
class DisplacementInfo(InfoBase):
    """Information about an atom to be displaced."""
    atom_nr: int
    atom_on_axis: bool


LinkGroup = Dict[int, Set[int]]


@dataclass(repr=False)
class SymmetryInfo(InfoBase):
    """Symmetry information pertaining to atoms in a slab."""

    on_planes: Tuple[int] = ()        # .num of atoms on plane
    on_axes: Tuple[int] = ()          # .num of atoms on rotation axis
    link_groups: LinkGroup = field(   # {num: {linked_nums}}
        default_factory=dict
        )
    hermann: str = ''          # expected plane group (Hermann-Maugin)
    invalid_direction: str = ''       # For subgroup reduction


@dataclass(repr=False)
class BulkInfo(InfoBase):
    """Information exclusively valid for bulk slabs."""

    screw_orders: Set[int] = None    # Rotation orders of screw axes
    n_glide_planes: int = None       # Number of 3D glide planes
    repeat: (float,)*3 = None        # Repeat vector
    periods: list = None             # Candidate periods


@dataclass(repr=False)
class ParametersInfo(InfoBase):
    """A container of information for PARAMETERS tests."""
    param_path : str = '' # Relative to data_path
    expected: dict = field(default_factory=dict)


@dataclass(repr=False)
class TestInfo(InfoBase):
    """Container of various pieces of information useful while testing."""
    poscar: ... = field(default_factory=POSCARInfo)
    symmetry: ... = field(default_factory=SymmetryInfo)
    bulk: ... = field(default_factory=BulkInfo)
    displacements: List[DisplacementInfo] = field(default_factory=list)
    parameters: ... = field(default_factory=ParametersInfo)
    param_presets: dict = field(default_factory=dict)  # For Rparams
    debug: dict = field(default_factory=dict)

    def update_params(self, param):
        """Update an Rparam with the presets stored."""
        for attr_name, attr_value in self.param_presets.items():
            setattr(param, attr_name, attr_value)


@dataclass(repr=False)
class BulkSlabAndRepeatInfo(InfoBase):
    """Container for information about bulk atoms and repeat vector."""
    bulk_like_below: float
    # Here the expected values:
    bulk_repeat: np.ndarray
    n_bulk_atoms: int
    bulk_cuts: List[float]
    bulk_dist: float


@dataclass(repr=False)
class LayerInfo(InfoBase):
    """Container for information about expected layer properties."""
    layer_cuts: LayerCuts
    n_bulk_layers: int
    # Here the expected values:
    cuts: List[float]
    n_layers: int
    n_sublayers: int
    n_atoms_per_layer: List[int]
    n_atoms_per_sublayer: List[int]
    smallest_interlayer_spacing: float
