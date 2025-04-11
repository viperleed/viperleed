"""Tests for module files/displacements/lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'

import pytest

from viperleed_jax.files.displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    OccDeltaLine,
    OffsetsLine,
    VibDeltaLine,
)


def test_geo_delta_line_repr():
    line = GeoDeltaLine('Si', '1', 'z', 0.0, 1.0, 0.1)
    assert 'Si 1' in repr(line)
    assert '0.0' in repr(line)
    assert '1.0' in repr(line)


def test_vib_delta_line_equality():
    line1 = VibDeltaLine('L', None, 0.0, 0.5, 0.1)
    line2 = VibDeltaLine('L', None, 0.0, 0.5, 0.1)
    assert line1 == line2


def test_constraint_line_without_direction_for_occ():
    line = ConstraintLine('occ', 'L', None, 0.5)
    assert line.direction is None


def test_constraint_line_invalid_direction_for_occ():
    with pytest.raises(ValueError, match='not allowed'):
        ConstraintLine('occ', 'L', 'z', 0.5)


def test_offsets_line_occ_with_invalid_direction():
    with pytest.raises(ValueError, match='not allowed'):
        OffsetsLine('occ', 'L', 'z', 0.3)
