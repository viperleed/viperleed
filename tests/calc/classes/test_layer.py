"""Tests for module viperleed.calc.classes.layer."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__created__ = '2023-10-04'

from contextlib import contextmanager
from dataclasses import dataclass

import numpy as np
import pytest
import pytest_cases

from viperleed.calc.classes.atom import Atom
from viperleed.calc.classes.layer import Layer, SubLayer
from viperleed.calc.classes.layer import LayerHasNoAtomsError
from viperleed.calc.classes.slab import Slab

from ...helpers import InfoBase, duplicate_all


_NO_VALUE = object()

@dataclass(repr=False)
class LayerInfo(InfoBase):
    """Container of information useful for checking tests."""
    n_atoms: int = _NO_VALUE
    thickness: float = _NO_VALUE
    cartbotz: float = _NO_VALUE
    cartori: tuple = _NO_VALUE
    element: str = _NO_VALUE
    cartpos: tuple = _NO_VALUE
    pos: tuple = _NO_VALUE


@contextmanager
def skip_if_no_info(info, attr_name):
    """Skip test if info.attr_name is not available."""
    value = getattr(info, attr_name)
    if value is _NO_VALUE:
        pytest.skip(reason=f'No expected value available for {attr_name}')
    yield


@pytest_cases.fixture(name='make_empty_layer', scope='session')
def factory_empty_layer():
    """Return an empty layer from a Layer subclass."""
    def _make(cls, **kwargs):
        mock_slab = Slab()
        mock_slab.ucell = np.array([[1.1, 0.0, 0.0],
                                    [0.0, 3.9, 0.0],
                                    [0.3, 0.5, 4.8]]).T
        return cls(mock_slab, 0, **kwargs), LayerInfo()
    return _make


def add_atom(layer):
    """Add one atom to layer."""
    atom = Atom('C', (0.1, 0.2, 0.3), 1, layer.slab)
    atom.cartpos = atom.pos.dot(atom.slab.ucell.T)
    layer.atlist.append(atom)
    atom.layer = layer
    if atom.el not in atom.slab.n_per_elem:
        atom.slab.n_per_elem[atom.el] = 0
    atom.slab.n_per_elem[atom.el] += 1
    layer.update_position()
    return atom


def move_atom(atom, delta_cartpos):
    """Move an atom by delta_cartpos."""
    u_inv = np.linalg.inv(atom.slab.ucell.T)
    atom.cartpos += delta_cartpos
    atom.pos = atom.cartpos.dot(u_inv)


@pytest_cases.parametrize(cls=(Layer, SubLayer))
def case_one_atom_layer(cls, make_empty_layer):
    """Return layer containing one atom, and some info."""
    layer, info = make_empty_layer(cls)
    add_atom(layer)
    info.n_atoms = 1
    info.thickness = 0
    info.cartbotz = 1.44  # 0.3*4.8
    info.cartori = (.09, .15, 1.44)
    if isinstance(layer, SubLayer):
        info.element = 'C'
        info.pos = (0.1, 0.2, 0.3)
        info.cartpos = (0.2, 0.93, 1.44)
    return layer, info


def case_more_atoms_layer(make_empty_layer):
    """Return layer containing a few atoms, and some info."""
    layer, info = make_empty_layer(Layer)
    atom = add_atom(layer)

    displacements = ((0.1, 0.2, -0.3),
                     (0.0, -0.1, 2.2),  # topmost atom
                     (-.08, 1.5, 1.4))
    for displacement in displacements:
        new_atom = atom.duplicate()
        move_atom(new_atom, displacement)
    layer.update_position()
    info.n_atoms = 1 + len(displacements)
    info.thickness = 2.5  # 2.2-(-0.3)
    info.cartbotz = 1.14  # 1.44 -0.3
    info.cartori = (.2275, .3791667, 3.64)
    return layer, info


class TestLayerRaises:
    """Collection of tests which check that exceptions are raised."""

    _properties = ('element', 'pos', 'cartpos')
    _methods = ('update_position',)

    @pytest_cases.parametrize(attribute=_properties)
    def test_raises_no_atoms_property(self, attribute, make_empty_layer):
        """Assert that accessing attribute raises without atoms."""
        empty_sublayer, *_ = make_empty_layer(SubLayer)
        with pytest.raises(LayerHasNoAtomsError):
            getattr(empty_sublayer, attribute)

    @pytest_cases.parametrize(method_name=_methods)
    def test_raises_no_atoms_method(self, method_name, make_empty_layer):
        """Assert that calling a method raises without atoms."""
        empty_sublayer, *_ = make_empty_layer(SubLayer)
        method = getattr(empty_sublayer, method_name)
        with pytest.raises(LayerHasNoAtomsError):
            method()


class TestLayer:
    """Collection of various tests for a Layer object."""

    _attributes = ('n_atoms', 'thickness', 'cartbotz', 'cartori',
                   'element', 'pos', 'cartpos')

    @pytest_cases.parametrize(attribute=_attributes)
    @pytest_cases.parametrize_with_cases('layer_and_info', cases='.')
    def test_attribute(self, attribute, layer_and_info):
        """Check correctness of the attribute of a layer."""
        layer, info = layer_and_info
        with skip_if_no_info(info, attribute):
            value = getattr(layer, attribute)
            expected = getattr(info, attribute)
            assert value == pytest.approx(expected)

    @pytest_cases.parametrize_with_cases('layer_and_info',
                                         cases=case_one_atom_layer)
    def test_atom_displaced_no_update(self, layer_and_info):
        """Check that layer positions do not change if atoms are moved."""
        layer, _ = layer_and_info
        ori_before, z_before = duplicate_all(layer.cartori, layer.cartbotz)
        move_atom(layer.atlist[0], (0.1, 0.1, 0.1))
        assert all(layer.cartori == ori_before)
        assert layer.cartbotz == z_before
