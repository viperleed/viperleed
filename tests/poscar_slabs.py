"""Module poscar_slabs of viperleed.tests.

Created on 2023-09-05

@author: Michele Riva (@michele-riva)

Contains definition of pytest cases generated from POSCAR files.
"""

from pathlib import Path
import copy
import inspect
import io
import sys

import numpy as np
from pytest_cases import parametrize, lazy_value, case

VPR_PATH = str(Path(__file__).resolve().parents[2])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Cannot do anything about it until we make viperleed installable
from viperleed.calc.classes.rparams import Rparams, LayerCuts, SymmetryEps
from viperleed.calc.files import poscar

from .helpers import POSCAR_PATH, TestInfo, DisplacementInfo, CaseTag as Tag
from .helpers import duplicate_all
# pylint: enable=wrong-import-position

POSCAR_FILES = tuple(POSCAR_PATH.glob('POSCAR*'))

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
              'SUPERLATTICE': np.array(((1, 1), (1, -1))),
              'superlattice_defined': True},
    'Ag': {'N_BULK_LAYERS': 1,
           'BULK_REPEAT': np.array([1.44, 1.44, -2.03647])},
    }

POSCARS_WITH_LITTLE_SYMMETRY_INFO = (
    _get_poscar_info('POSCAR_Ag(100)', 6, 'p4m', (1, True), _PRESETS['Ag']),
    _get_poscar_info('POSCAR_STO(110)-4x1', 136, 'pm', (1, False)),
    _get_poscar_info('POSCAR_TiO2', 540, 'pmm', (540, False)),
    _get_poscar_info('POSCAR_diamond', 96, 'rcm', (90, False)),
    _get_poscar_info('POSCAR_36C_p6m', 36, 'p6m', (1, False)),
    _get_poscar_info('POSCAR_36C_cm', 36, 'cm', (1, False)),
    _get_poscar_info('POSCAR_Fe3O4_SCV', 83, 'cmm', (51, False),
                     _PRESETS['Fe3O4']),
    _get_poscar_info('POSCAR_Al2O3_NiAl(111)_cHole_20061025',
                     {'Ni': 402, 'Al': 134+132+188, 'O': 188+213},
                     'p3'),
    _get_poscar_info('POSCAR_Cu2O(111)_1x1_surplus_oxygen', 66, 'p3m1'),
    )

def _get_info_by_name(name):
    """Return a TestInfo object by name."""
    return next(i for i in POSCARS_WITH_LITTLE_SYMMETRY_INFO
                if name in i.poscar.name)

AG_100 = _get_info_by_name('Ag(100)')
SLAB_36C_cm = _get_info_by_name('36C_cm')
SLAB_Cu2O_111 = _get_info_by_name('Cu2O(111)')

POSCARS_WITHOUT_INFO = [
    _get_poscar_info(f.name) for f in POSCAR_PATH.glob('POSCAR*')
    ]


class CasePOSCARFiles:
    @staticmethod
    def _string(file):
        """Return the string of a POSCAR file name."""
        with open(POSCAR_PATH / file, 'r', encoding='utf-8') as file:
            contents = file.read()
        return contents

    @parametrize(file=POSCAR_FILES)
    def case_poscar_file(self, file):
        """Return a POSCAR file name."""
        content = self._string(file)
        return file, content


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
        slab.fullUpdate(param)
        return slab, param, info

    @parametrize(info=POSCARS_WITHOUT_INFO, idgen=make_poscar_ids())
    @case(tags=Tag.NO_INFO)
    def case_infoless_poscar(self, info):
        """Return a slab, an Rparams and an (essentially) empty info."""
        return self.case_poscar(info)

    def case_poscar_fe3o4_001_cod(self):
        """Return a slab from a bulk-truncated Fe3O4(001) POSCAR."""
        info = TestInfo()
        info.poscar.set('POSCAR_Fe3O4_(001)_cod1010369', 112)
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
        return self.case_poscar(info)


class CaseBulkSlabs:
    """Collection of bulk slabs."""

    @staticmethod
    def _make_bulk(slab, param):
        """Return a bulk slab and parameters from a surface slab."""
        param = copy.deepcopy(param)
        bulk_slab = slab.makeBulkSlab(param)
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
        return slab, param, info


@parametrize(bulk=CaseBulkSlabs.all_slabs(), ids=make_poscar_ids('thick_bulk'))
@case(tags=(Tag.BULK, Tag.THICK_BULK))
def case_double_bulk(bulk):
    """Return a bulk slab with twice the thickness."""
    slab, param, info = bulk
    thick_slab = slab.doubleBulkSlab()

    param, info = duplicate_all(param, info)
    param.BULK_REPEAT *= 2
    info.bulk.repeat = tuple(param.BULK_REPEAT)

    # Propagate symmetry information by using duplicate atoms
    sym_info = info.symmetry
    _map = {at.duplicate_of.num: at.num
            for at in thick_slab.atlist if at.duplicate_of}
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
