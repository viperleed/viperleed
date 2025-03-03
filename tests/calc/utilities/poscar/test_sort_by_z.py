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
from viperleed.calc.files import poscar

from ....helpers import POSCAR_PATH
from ... import poscar_slabs
from ...tags import CaseTag as Tag

infoless = pytest_cases.parametrize_with_cases('test_slab',
                                               cases=poscar_slabs,
                                               has_tag=Tag.NO_INFO)

@infoless
def test_sort_by_z(test_slab):
    slab, *_ = test_slab
    parser = SortByZCLI().parser
    args = parser.parse_args([])
    # shuffle the slab
    slab.atlist = random.sample(slab.atlist, len(slab.atlist))

    # Sort the slab by z-coordinate
    sorted_slab = SortByZCLI().process_slab(slab, args)

    # Check if the slab is sorted by z-coordinate
    assert slab_is_sorted(sorted_slab)

@infoless
def test_sort_by_z_reversed(test_slab):
    slab, *_ = test_slab
    parser = SortByZCLI().parser
    args = parser.parse_args(['-r'])
    # shuffle the slab
    slab.atlist = random.sample(slab.atlist, len(slab.atlist))

    # Sort the slab by z-coordinate
    sorted_slab = SortByZCLI().process_slab(slab, args)

    # Check if the slab is sorted by z-coordinate
    assert slab_is_sorted(sorted_slab, reversed=True)

@infoless
def test_overwrite_file(test_slab, tmp_path):
    *_, info = test_slab
    poscar_path = POSCAR_PATH/info.poscar.name

    # copy to a temporary file because we will overwrite the file
    tmp_poscar_path = tmp_path / info.poscar.name
    tmp_poscar_path.write_text(poscar_path.read_text())

    sort_by_z_cli = SortByZCLI()
    # give the temporary file as input and output
    sort_by_z_cli([
        '-i', str(tmp_poscar_path),
        '-o', str(tmp_poscar_path),
        ])

    # read the modified POSCAR file
    modified_slab = poscar.read(tmp_poscar_path)

    # Check if the slab is sorted by z-coordinate
    assert slab_is_sorted(modified_slab)

def slab_is_sorted(slab, reversed=False):
    """Check if the slab is sorted by z-coordinate."""
    if reversed:
        z_order = np.diff([at.pos[2] for at in slab]) >= 0
    else:
        z_order = np.diff([at.pos[2] for at in slab]) <= 0

    not_same_el = np.array(at1.el != at2.el
                      for at1, at2 in zip(slab.atlist[:-1], slab.atlist[1:]))
    print(z_order, not_same_el)
    return np.all(np.logical_or(z_order, not_same_el))
