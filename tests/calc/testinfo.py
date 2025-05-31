"""Module testinfo of viperleed.tests.calc.

Defines various InfoBase subclasses useful to define expected
information in tests for the calc package of viperleed. This
module used to be part of viperleed.tests.helpers.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-03-19'
__license__ = 'GPLv3+'

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Set, Tuple

import numpy as np

from viperleed.calc.classes.rparams.special.layer_cuts import LayerCuts
from viperleed.calc.lib.dataclass_utils import non_init_field

from ..helpers import InfoBase


@dataclass(repr=False)
class BulkInfo(InfoBase):
    """Information exclusively valid for bulk slabs."""

    screw_orders: Set[int] = None    # Rotation orders of screw axes
    n_glide_planes: int = None       # Number of 3D glide planes
    repeat: (float,)*3 = None        # Repeat vector
    periods: list = None             # Candidate periods


@dataclass(repr=False)
class BulkSlabAndRepeatInfo(InfoBase):
    """Container for information about bulk atoms and repeat vector."""
    bulk_like_below: float
    # Here the expected values:
    bulk_repeat: np.ndarray   # From bulk to surface
    n_bulk_atoms: int
    bulk_cuts: List[float]
    bulk_dist: float
    bulk_ucell: np.ndarray


# Avoid using atoms that are on planes: Easier to test.
@dataclass(repr=False)
class DisplacementInfo(InfoBase):
    """Information about an atom to be displaced."""
    atom_nr: int
    atom_on_axis: bool


@dataclass(repr=False)
class LayersInfo(InfoBase):
    """Container for information about expected layer properties."""
    layer_cuts: LayerCuts
    n_bulk_layers: int
    # Here the expected values:
    cuts: List[float]
    n_layers: int
    n_sublayers: int
    n_atoms_per_layer: List[int]
    n_atoms_per_sublayer: List[int]
    smallest_interlayer_gap: float


@dataclass(repr=False)
class NearestNeighborInfo(InfoBase):
    """Container for information about nearest neighbors."""
    # Keys are atom number
    nearest_neighbor_distances: Dict[int, float]


@dataclass(repr=False)
class ParametersInfo(InfoBase):
    """A container of information for PARAMETERS tests."""
    param_path : str = '' # Relative to data_path
    expected: dict = field(default_factory=dict)


@dataclass(repr=False)
class POSCARInfo(InfoBase):
    """Container for information about POSCAR files."""
    name: str = ''
    n_atoms: int = None
    n_atoms_by_elem: dict = non_init_field(default_factory=dict, repr=True)
    n_cells: int = 1  # How many 2D cells are in the POSCAR?

    def __post_init__(self):
        """Post-process initialization values."""
        if isinstance(self.n_atoms, Mapping):
            self.n_atoms_by_elem = self.n_atoms
            self.n_atoms = sum(self.n_atoms_by_elem.values())


@dataclass(repr=False)
class SurfaceAtomInfo(InfoBase):
    """Container for information about which atoms are at the surface."""
    # Keys are atom number
    surface_atom_nums: Tuple[int]


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
