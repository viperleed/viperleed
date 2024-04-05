"""Tests for the viperleed poscar sort by z utility."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-05'
__license__ = 'GPLv3+'

import random

import numpy as np

import pytest_cases

from viperleed.utilities.poscar.sort_by_z import SortByZCLI

from ... import poscar_slabs

@pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs)
def test_sort_by_z(test_slab):
    slab, *_ = test_slab
    parser = SortByZCLI().parser
    args = parser.parse_args([])
    # shuffle the slab
    slab.atlist = random.sample(slab.atlist, len(slab.atlist))

    # Sort the slab by z-coordinate
    sorted_slab = SortByZCLI().process_slab(slab, args)

    # Check if the slab is sorted by z-coordinate
    assert np.all(np.diff([at.pos[2] for at in sorted_slab]) <= 0)

@pytest_cases.parametrize_with_cases('test_slab', cases=poscar_slabs)
def test_sort_by_z_reversed(test_slab):
    slab, *_ = test_slab
    parser = SortByZCLI().parser
    args = parser.parse_args(['-r'])
    # shuffle the slab
    slab.atlist = random.sample(slab.atlist, len(slab.atlist))

    # Sort the slab by z-coordinate
    sorted_slab = SortByZCLI().process_slab(slab, args)

    # Check if the slab is sorted by z-coordinate
    assert np.all(np.diff([at.pos[2] for at in sorted_slab]) >= 0)
