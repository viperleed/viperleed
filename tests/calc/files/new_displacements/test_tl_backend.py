"""Tests for module viperleed.calc.files.new_displacements.tl_backend."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

import pytest

from viperleed.calc.files.new_displacements.tl_backend import (
    TensErLEEDBackend,
    TensorLEEDBackend,
)


def test_replace_search_block_is_abstract():
    """Check that the base class declares an abstract class method."""
    assert 'replace_search_block' in TensorLEEDBackend.__abstractmethods__
    assert isinstance(
        TensorLEEDBackend.__dict__['replace_search_block'],
        classmethod,
    )


def test_base_class_not_instantiable():
    """Check that the abstract backend cannot be instantiated."""
    with pytest.raises(TypeError):
        TensorLEEDBackend()


def test_subclass_needs_name():
    """Check that subclasses must define a backend name."""
    with pytest.raises(TypeError, match='name'):
        type('Nameless', (TensorLEEDBackend,), {})


def test_tensorleed_backend_replace_not_implemented():
    """Check the TensErLEED backend's placeholder implementation."""
    with pytest.raises(NotImplementedError):
        TensErLEEDBackend.replace_search_block(None)
