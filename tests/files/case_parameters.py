"""Cases for test_parameters.

Created on 2023-09-07

@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

from pytest_cases import parametrize

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.files import parameters

from ..poscar_slabs import POSCARS_WITHOUT_INFO, AG_100
# pylint: enable=wrong-import-position


_POSCAR_INFO = {
    'Ag': AG_100,
    'Ir': next(info for info in POSCARS_WITHOUT_INFO
               if info.poscar.name.endswith('Ir(100)-(2x1)-O')),
    }
_READ = {
    'Ag': {'V0_IMAG': 5.0, 'THEO_ENERGIES': [50, 350, 3]},
    'Ir': {'RUN': [0, 1, 2, 3],
           'THEO_ENERGIES': [49, 700, 3], 'LMAX': [8, 14],
           'BULK_LIKE_BELOW': 0.35, 'T_DEBYE': 420, 'T_EXPERIMENT': 100,
           'SITE_DEF': {'Ir': {'surf': {3, 2}}, 'O': {'ads': {1}}},
           'VIBR_AMP_SCALE': ['*surf 1.3',], 'R_FACTOR_TYPE': 1},
    }
_PATHS = {
    'Ag': 'Ag(100)/initialization/PARAMETERS',
    'Ir': 'parameters/PARAMETERS_Ir(100)-(2x1)-O',
    }

for label, test_info in _POSCAR_INFO.items():
    test_info.parameters.param_path = _PATHS[label]
    test_info.parameters.expected = _READ[label]


@parametrize(info=_POSCAR_INFO.values(), ids=_POSCAR_INFO)
def case_parameters_slab(info, make_poscar, data_path):
    """Return a slab, an Rparam (read from file), and TestInfo."""
    slab, _, info = make_poscar(info)
    rpars = parameters.read(data_path / info.parameters.param_path)
    return slab, rpars, info
