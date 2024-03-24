"""Module poscar_slabs of viperleed.tests.calc.

Contains definition of pytest cases generated from POSCAR files.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-09-05'
__license__ = 'GPLv3+'

import copy
import inspect

import numpy as np
from pytest_cases import parametrize, lazy_value, case

from viperleed.calc.classes.rparams import Rparams, LayerCuts, SymmetryEps
from viperleed.calc.files import poscar

from ..helpers import POSCAR_PATH, duplicate_all
from .tags import CaseTag as Tag
from .testinfo import BulkSlabAndRepeatInfo
from .testinfo import DisplacementInfo
from .testinfo import LayersInfo
from .testinfo import NearestNeighborInfo
from .testinfo import SurfaceAtomInfo
from .testinfo import TestInfo


def _get_poscar_info(*args):
    """Return a TestInfo object appropriately filled with args.

    Parameters
    ----------
    *args : object
        Values to be stored in the TestInfo object. The following
        values are expected:
        name : str
            File name of the POSCAR file.
        n_atoms : int or dict, optional
            Nr. of atoms that the structure is composed of. If a
            dict, it should have the form {elem: count} where elem
            is the element name in the POSCAR file, and count the
            number of atoms of that element.
        group : str, optional
            Hermann-Maugin name of the plane group of the
            structure
        displaced_atoms_info : Sequence, optional
            Each item has two elements. The first is the index
            (1-based!) of the atom to be displaced. The second
            a bool specifying whether the atom is free to move
            or not. If only one atom is to be displaced, it is
            acceptable to give only a single 2-tuple.
        param_presets : dict, optional
            Special values of PARAMETERS to be used when handling
            this slab. Keys are attributes of Rparams objects,
            values their values. A deep copy of `param_presets` is
            stored in `test_info`.

    Returns
    -------
    test_info : TestInfo
        The information provided, packed into an object that
        tests recognize.
    """
    test_info = TestInfo()
    test_info.poscar.set(*args[:2])  # name and n_atoms
    try:
        test_info.symmetry.hermann = args[2]
    except IndexError:  # No symmetry or other info available
        return test_info

    displace_info = args[3:4]
    if len(displace_info) == 2 and isinstance(displace_info[0], int):
        # Given as a single tuple. Only one atom to displace
        displace_info = (displace_info,)
    for atom_nr, on_axis in displace_info:
        test_info.displacements.append(DisplacementInfo(atom_nr, on_axis))

    params = args[4:5]
    if params:
        test_info.param_presets = copy.deepcopy(params[0])
    return test_info


def make_poscar_ids(suffix=None):
    """Return a callable that generates case IDs with a suffix."""

    def _clean_up_name(name):
        """Return a cleaned-up name."""
        name = name.replace('POSCAR_', '')
        if suffix is None:
            return name

        name = name.replace('case_', '')
        if 'thick' in suffix:
            return name.replace('_bulk', '')
        return f'{name}_{suffix}'

    def _make_poscar_ids(*args, **argnames_and_values):
        """Return case-IDs for the case parametrizations given."""
        # There should only be a single relevant argument.
        args += tuple(argnames_and_values.values())
        assert len(args) == 1

        # Its value can be a bunch of TestInfo objects...
        info = args[0]
        if isinstance(info, TestInfo):
            return _clean_up_name(info.poscar.name)

        # ...or a case
        return _clean_up_name(info.__name__)
    return _make_poscar_ids


# PARAMETERS presets for slabs
_PRESETS = {
    'Fe3O4': {'LAYER_CUTS': LayerCuts.from_string('0.1 0.2 <dz(1.0)'),
              'N_BULK_LAYERS': 2,
              'SYMMETRY_EPS': SymmetryEps(0.3),
              'BULK_REPEAT': np.array([0.0, -4.19199991, 4.19199991]),
              'SUPERLATTICE': np.array(((1, -1), (1, 1))),
              'superlattice_defined': True},
    }

POSCARS_WITH_LITTLE_SYMMETRY_INFO = (
    _get_poscar_info('POSCAR_TiO2_supercell', 540, 'pmm', (540, False)),
    _get_poscar_info('POSCAR_Fe3O4_SCV', 83, 'cmm', (51, False),
                     _PRESETS['Fe3O4']),
    _get_poscar_info('POSCAR_Al2O3_NiAl(111)_cHole_20061025',
                     {'Ni': 402, 'Al': 134+132+188, 'O': 188+213},
                     'p3'),
    _get_poscar_info('POSCAR_Cu2O(111)_1x1_surplus_oxygen', 66, 'p3m1'),
    )


POSCARS_WITHOUT_INFO = [
    _get_poscar_info(f.name) for f in POSCAR_PATH.glob('POSCAR*')
    if 'duplicate' not in f.name
    ]

WITH_DUPLICATE_ATOMS = [
    _get_poscar_info(f.name) for f in POSCAR_PATH.glob('POSCAR*duplicate*')
    ]


def get_info_by_name(name, container=POSCARS_WITH_LITTLE_SYMMETRY_INFO):
    """Return a TestInfo object by name."""
    return next(i for i in container if name in i.poscar.name)


def _add_known_bulk_properties(info, bulk_info):
    """Add bulk properties to a TestInfo object."""
    info.bulk_properties = bulk_info
    return info


def _add_known_layer_properties(info, layer_info):
    """Add layer properties to a TestInfo object."""
    info.layer_properties = layer_info
    return info


def _add_nearest_neighbor_info(info, nn_info):
    """Add nearest-neighbour information to a TestInfo object."""
    info.nearest_neighbors = nn_info
    return info


def _add_surface_atom_info(info, surface_atom_info):
    """Add info about surface atoms to a TestInfo object."""
    info.surface_atoms = surface_atom_info
    return info


POSCAR_WITH_KNOWN_BULK_REPEAT = (
    _add_known_bulk_properties(
        _get_poscar_info('POSCAR_Cu2O_111', 22+43),
        BulkSlabAndRepeatInfo(
            bulk_like_below=0.55,
            bulk_repeat=np.array([0, -3.4766, 2.4584]),
            n_bulk_atoms=6,
            bulk_cuts=[0.08333],
            bulk_dist=0.61459,  # Closest z distance between Cu and O
            bulk_ucell=np.array([[ 6.02172136, -3.01086068,  0.        ],
                                 [ 0.        ,  5.21496367, -3.4766    ],
                                 [ 0.        ,  0.        ,  2.45835787]]),
            ),
        ),
    )
POSCAR_WITH_LAYER_INFO = ()
POSCAR_WITH_NEAREST_NEIGHBOR_INFO = ()
POSCAR_WITH_SURFACE_ATOM_INFO = ()

SLAB_Cu2O_111 = get_info_by_name('Cu2O(111)')


class CasePOSCARSlabs:
    """Collection of test setups various structures read from POSCAR files."""

    @staticmethod
    def _read(name):
        """Return a Slab from a POSCAR file name."""
        return poscar.read(POSCAR_PATH / name)

    @staticmethod
    def _updated_param(info):
        """Return a fresh Rparams objects updated from info."""
        param = Rparams()
        info.update_params(param)
        return param

    @parametrize(info=POSCARS_WITH_LITTLE_SYMMETRY_INFO,
                 idgen=make_poscar_ids())
    def case_poscar(self, info):
        """Return a slab, an Rparams and information for tests."""
        slab = self._read(info.poscar.name)
        param = self._updated_param(info)
        slab.full_update(param)
        return slab, param, info

    @parametrize(info=POSCARS_WITHOUT_INFO, idgen=make_poscar_ids())
    @case(tags=Tag.NO_INFO)
    def case_infoless_poscar(self, info):
        """Return a slab, an Rparams and an (essentially) empty info."""
        return self.case_poscar(info)

    @parametrize(info=POSCAR_WITH_KNOWN_BULK_REPEAT, idgen=make_poscar_ids())
    @case(tags=Tag.BULK_PROPERTIES)
    def case_bulk_repeat_poscar(self, info):
        """Return a slab, an Rparams and info on expected bulk properties."""
        slab, rpars, info = self.case_poscar(info)
        rpars.BULK_REPEAT = info.bulk_properties.bulk_repeat
        return slab, rpars, info

    @parametrize(info=POSCAR_WITH_LAYER_INFO, idgen=make_poscar_ids())
    @case(tags=Tag.LAYER_INFO)
    def case_layer_info_poscar(self, info):
        """Return a slab, an Rparams and info on expected layers."""
        slab, rpars, info = self.case_poscar(info)
        rpars.LAYER_CUTS = info.layer_properties.layer_cuts
        rpars.N_BULK_LAYERS = info.layer_properties.n_bulk_layers
        slab.create_layers(rpars)
        return slab, rpars, info

    @parametrize(info=POSCAR_WITH_NEAREST_NEIGHBOR_INFO,
                 idgen=make_poscar_ids())
    @case(tags=Tag.NEAREST_NEIGHBOURS)
    def case_nearest_neighbors_poscar(self, info):
        """Return a slab, an Rparams and info on expected nearest neighbors."""
        return self.case_poscar(info)

    @parametrize(info=POSCAR_WITH_SURFACE_ATOM_INFO, idgen=make_poscar_ids())
    @case(tags=Tag.SURFACE_ATOMS)
    def case_surface_atom_poscar(self, info):
        """Return a slab, an Rparams and info on expected surface atoms."""
        return self.case_poscar(info)

    @case(tags=(Tag.BULK_PROPERTIES,
                Tag.LAYER_INFO,
                Tag.NEAREST_NEIGHBOURS,
                Tag.SURFACE_ATOMS,))
    def case_poscar_ag100(self):
        """Return a Ag(100) slab."""
        info = _get_poscar_info('POSCAR_Ag(100)', 6, 'p4m', (1, True))
        info.param_presets = {'N_BULK_LAYERS': 1,
                              'BULK_REPEAT': np.array([-1.44, -1.44, 2.03647])}
        info.bulk_properties = BulkSlabAndRepeatInfo(
            bulk_like_below=0.65,
            bulk_repeat=info.param_presets['BULK_REPEAT'],
            n_bulk_atoms=1,
            bulk_cuts=[0.35],
            bulk_dist=0.0,
            bulk_ucell=np.array([[ 2.88   ,  0.     , -1.44   ],
                                 [ 0.     ,  2.88   , -1.44   ],
                                 [ 0.     ,  0.     ,  2.03647]]),
            )
        info.layer_properties = LayersInfo(
            layer_cuts=LayerCuts.from_string('dz(1.2)'),
            n_bulk_layers=info.param_presets['N_BULK_LAYERS'],
            cuts=[0.35, 0.45, 0.55, 0.65, 0.75],
            n_layers=6,
            n_atoms_per_layer=[1, 1, 1, 1, 1, 1],
            n_atoms_per_sublayer=[1, 1, 1, 1, 1, 1],
            n_sublayers=6,
            smallest_interlayer_gap=2.03646,
            )
        info.nearest_neighbors = NearestNeighborInfo(
            nearest_neighbor_distances={1: 2.88, 2: 2.88, 3: 2.88,
                                        4: 2.88, 5: 2.88, 6: 2.88,},
            )
        info.surface_atoms = SurfaceAtomInfo(surface_atom_nums=(1, 2),)
        return self.case_poscar(info)

    @case(tags=Tag.VACUUM_GAP_SMALL)
    def case_poscar_36carbon_atoms_p6m(self):
        """Return a 6-fold symmetric slab with 36 carbon atoms."""
        info = _get_poscar_info('POSCAR_36C_p6m', 36, 'p6m', (1, False))
        return self.case_poscar(info)

    @case(tags=Tag.VACUUM_GAP_SMALL)
    def case_poscar_36carbon_atoms_cm(self):
        """Return a slanted slab with 36 carbon atoms and cm group."""
        info = _get_poscar_info('POSCAR_36C_cm', 36, 'cm', (1, False))
        return self.case_poscar(info)

    @case(tags=Tag.NON_MINIMAL_CELL)
    def case_poscar_diamond(self):
        """Return a non-minimal diamond(111) slab."""
        info = _get_poscar_info('POSCAR_diamond', 96, 'rcm', (90, False))
        info.poscar.n_cells = 8
        info.param_presets = {'SUPERLATTICE': np.array([[2, -2], [2, 2]])}
        return self.case_poscar(info)

    def case_poscar_fe2o3_012_1x1(self):
        """Return a non-optimized Fe2O3(012) slab."""
        info = _get_poscar_info('POSCAR_Fe2O3_012_1x1_not_optimized', 40, 'pg')
        return self.case_poscar(info)

    @case(tags=(Tag.NON_MINIMAL_CELL,
                Tag.BULK_PROPERTIES,
                Tag.VACUUM_GAP_SMALL))
    def case_poscar_fe3o4_001_cod(self):
        """Return a slab from a bulk-truncated Fe3O4(001) POSCAR."""
        info = TestInfo()
        info.poscar.set('POSCAR_Fe3O4_(001)_cod1010369', 112)
        info.poscar.n_cells = 2
        info.param_presets = copy.deepcopy(_PRESETS['Fe3O4'])
        info.symmetry.hermann = 'cmm'
        info.symmetry.on_planes = (
            1, 2, 5, 6, 13, 14, 17, 18, 21, 22, 23, 24, 27, 28, 31, 32, 41, 42,
            43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 59, 60, 61, 62, 69, 70, 75,
            76, 79, 80, 83, 84, 89, 90, 91, 92, 93, 94, 95, 96, 99, 100, 103,
            104, 109, 110, 111, 112
            )
        info.symmetry.on_axes = (33, 34, 35, 36, 37, 38, 39, 40)
        info.symmetry.link_groups = {
            1: {1, 17}, 2: {2, 18}, 3: {3, 9, 15, 25}, 4: {4, 10, 16, 26},
            5: {5, 13}, 6: {6, 14}, 7: {7, 11, 19, 29}, 8: {8, 12, 20, 30},
            21: {21, 31}, 22: {22, 32}, 23: {23, 27}, 24: {24, 28}, 33: {33},
            34: {34}, 35: {35, 39}, 36: {36, 40}, 37: {37}, 38: {38},
            41: {41, 45}, 42: {42, 46}, 43: {43, 47}, 44: {44, 48},
            49: {49, 79}, 50: {50, 80}, 51: {51, 71, 77, 97},
            52: {52, 72, 78, 98}, 53: {53, 75}, 54: {54, 76},
            55: {55, 73, 81, 107}, 56: {56, 74, 82, 108},
            57: {57, 65, 85, 101}, 58: {58, 66, 86, 102}, 59: {59, 83},
            60: {60, 84}, 61: {61, 69}, 62: {62, 70}, 63: {63, 67, 87, 105},
            64: {64, 68, 88, 106}, 89: {89, 111}, 90: {90, 112}, 91: {91, 103},
            92: {92, 104}, 93: {93, 109}, 94: {94, 110}, 95: {95, 99},
            96: {96, 100}
            }
        info.displacements.extend((DisplacementInfo(8, False),
                                   DisplacementInfo(36, True)))
        info.bulk.repeat = tuple(info.param_presets['BULK_REPEAT'])
        info.bulk_properties = BulkSlabAndRepeatInfo(
            bulk_like_below=0.7,
            bulk_repeat=np.array([0, 4.192, 4.192]),
            n_bulk_atoms=14,
            bulk_cuts=[0.21885],
            bulk_dist=1.048,  # Z distance between Fe_tet and Fe_oct
            bulk_ucell=np.array([[ 4.192, -4.192,  0.000],
                                 [ 4.192,  4.192, -4.192],
                                 [ 0.000,  0.000,  4.192]]),
            )
        return self.case_poscar(info)

    @case(tags=(Tag.NON_MINIMAL_CELL, Tag.SURFACE_ATOMS))
    def case_poscar_lsmo_001_rt2(self):
        """Return a sqrt(2) x sqrt(2) LSMO (001) A-site bulk-truncated slab."""
        info = _get_poscar_info('POSCAR_LSMO_001_A_site', 60)
        info.poscar.n_cells = 2
        info.surface_atoms = SurfaceAtomInfo((2, 6, 9, 11,      # A
                                              19, 26, 38, 47))  # O
        slab, rpars, info = self.case_poscar(info)

        # Assign mixed sites. Notice that, due to a bug in initSites,
        # we have to clear the current atomic sites before making new
        # ones. Otherwise we'd get an empty sitelist.
        rpars.ELEMENT_MIX = {'A': ['La', 'Sr']}
        for atom in slab:
            atom.site = None
        slab.initSites(rpars)
        a_site = next(s for s in slab.sitelist if s.el == 'A')
        a_site.occ = {'La': 0.8, 'Sr': 0.2}
        return slab, rpars, info

    @case(tags=(Tag.NON_MINIMAL_CELL, Tag.VACUUM_GAP_SMALL))
    def case_poscar_mgo(self):
        """Return a non-minimal Mg(001) slab."""
        info = _get_poscar_info('POSCAR_MgO_cod_9006456',
                                {'Mg': 14, 'O': 14}, 'p4m')
        info.poscar.n_cells = 2
        return self.case_poscar(info)

    @case(tags=(Tag.NON_MINIMAL_CELL, Tag.VACUUM_GAP_SMALL))
    def case_poscar_sb_si_111(self):
        """Return a non-minimal, rectangular slab of Sb/Si(111)."""
        info = _get_poscar_info('POSCAR_Sb_Si(111)_rect', 72+12, 'pm')
        # While it looks like there are 4 cells (rt3 of Si), some
        # of the Sb atoms are not quite in the right position
        info.poscar.n_cells = 2
        return self.case_poscar(info)

    @case(tags=Tag.VACUUM_GAP_MIDDLE)
    def case_poscar_sto110_4x1(self):
        """Return a symmetric slab of SrTiO3(110)-4x1."""
        info = _get_poscar_info('POSCAR_STO(110)-4x1', 136, 'pm', (1, False))
        return self.case_poscar(info)

    @case(tags=Tag.BULK_PROPERTIES)
    def case_poscar_tio2_small(self):
        """Return a rutile TiO2(110) with a few bulk layers."""
        info = _get_poscar_info('POSCAR_TiO2_small', 29)
        _add_known_bulk_properties(
            info,
            BulkSlabAndRepeatInfo(
                bulk_like_below=0.43,
                bulk_repeat=np.array([-3.24845, 0.00000, 3.24876]),
                n_bulk_atoms=6,
                bulk_cuts=[0.2281, 0.1628],
                bulk_dist=1.2682,
                bulk_ucell=np.array([[6.4969, 0.000, -3.24845],
                                     [0.0000, 2.959,  0.00000],
                                     [0.0000, 0.000,  3.24876]]),
                ),
            )
        slab, rpars, info = self.case_bulk_repeat_poscar(info)
        rpars.N_BULK_LAYERS = 2
        slab.full_update(rpars)
        return slab, rpars, info


class CaseBulkSlabs:
    """Collection of bulk slabs."""

    @staticmethod
    def _make_bulk(slab, param):
        """Return a bulk slab and parameters from a surface slab."""
        param = copy.deepcopy(param)
        bulk_slab = slab.make_bulk_slab(param)
        return bulk_slab, param

    @classmethod
    def all_slabs(cls):
        """Return an iterable of bulk slabs as lazy values."""
        instance = cls()
        return [lazy_value(method)
                for name, method in inspect.getmembers(instance)
                if name.startswith('case_')]

    @case(tags=Tag.BULK)
    def case_fe3o4_bulk(self):
        """Return a bulk Fe3O4(001) slab."""
        *surf, _ = CasePOSCARSlabs().case_poscar_fe3o4_001_cod()
        slab, param = self._make_bulk(*surf)
        info = TestInfo()
        info.symmetry.set(
            (3, 5, 13, 15, 49, 55, 57, 61, 69, 73, 79, 85),
            (33, 43),
            {3: {3, 15}, 5: {5, 13}, 33: {33}, 43: {43}, 49: {49, 79},
             55: {55, 73}, 57: {57, 85}, 61: {61, 69}},
            'pmm'
            )
        info.bulk.screw_orders = {4}  # {2, 4} would make more sense!
        info.bulk.n_glide_planes = 2
        info.bulk.periods = [3]
        return slab, param, info


@parametrize(bulk=CaseBulkSlabs.all_slabs(), ids=make_poscar_ids('thick_bulk'))
@case(tags=(Tag.BULK, Tag.THICK_BULK))
def case_double_bulk(bulk):
    """Return a bulk slab with twice the thickness."""
    slab, param, info = bulk
    thick_slab = slab.with_double_thickness()

    param, info = duplicate_all(param, info)
    param.BULK_REPEAT *= 2
    info.bulk.repeat = tuple(param.BULK_REPEAT)

    # Propagate sublayer periods information
    info.bulk.periods += [2*p for p in info.bulk.periods]

    # Propagate symmetry information by using duplicate atoms
    sym_info = info.symmetry
    _map = {at.duplicate_of.num: at.num
            for at in thick_slab if at.duplicate_of}
    on_planes = sym_info.on_planes
    on_planes += tuple(new_at for base_at, new_at in _map.items()
                       if base_at in sym_info.on_planes)
    on_axes = sym_info.on_axes
    on_axes += tuple(new_at for base_at, new_at in _map.items()
                     if base_at in sym_info.on_axes)
    links = sym_info.link_groups.copy()
    links.update({new_at: {_map[linked] for linked in links[base_at]}
                  for base_at, new_at in _map.items() if base_at in links})
    sym_info.set(on_planes, on_axes, links)

    return thick_slab, param, info
