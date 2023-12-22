"""Cases for test_parameters.

Created on 2023-09-07

@author: Michele Riva (@michele-riva)
"""

from pathlib import Path
import sys

import numpy as np
from pytest_cases import case, parametrize

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.classes.rparams import LMax
from viperleed.tleedmlib.files import parameters

from ...helpers import TestInfo, CaseTag
from ...poscar_slabs import POSCARS_WITHOUT_INFO, AG_100
# pylint: enable=wrong-import-position


_POSCAR_INFO = {
    'Ag': AG_100,
    'Ir': next(info for info in POSCARS_WITHOUT_INFO
               if info.poscar.name.endswith('Ir(100)-(2x1)-O')),
    }
_READ = {
    'Ag': {'V0_IMAG': 5.0, 'THEO_ENERGIES': [50, 350, 3],
           'RUN': [0], 'LOG_LEVEL': 10, 'N_BULK_LAYERS': 1,
           'BULK_REPEAT': np.array([1.44, 1.44, -2.03646753]),
           'SITE_DEF' : {'Ag': {'surf': {1}}}, 'LMAX': LMax(8, 12),
           'T_DEBYE': 330, 'T_EXPERIMENT': 100,
           'VIBR_AMP_SCALE': ['*surf 1.3',],
           'SUPERLATTICE': np.identity(2)},
    'Ir': {'RUN': [0, 1, 2, 3],
           'THEO_ENERGIES': [49, 700, 3], 'LMAX': LMax(8, 14),
           'BULK_LIKE_BELOW': 0.35, 'T_DEBYE': 420, 'T_EXPERIMENT': 100,
           'SITE_DEF': {'Ir': {'surf': {3, 2}}, 'O': {'ads': {1}}},
           'VIBR_AMP_SCALE': ['*surf 1.3',], 'V0_IMAG': 5.0},
    'empty': {}, 'stop': {}, 'no_stop': {}, 'missing_equals': {},
    'left empty': {},
    }
for key in ('stop', 'no_stop', 'left empty'):
    _READ[key].update(**_READ['Ir'])  # Copies of the Ir file
_READ['missing_equals'].update(**_READ['Ag'])
del _READ['missing_equals']['SUPERLATTICE']

# FORTRAN_COMP added only for the 'Ag'
_READ['Ag']['FORTRAN_COMP'] = ['gfortran', '-llapack -lpthread -lblas']

_PATHS = {
    'Ag': 'Ag(100)/initialization/PARAMETERS',
    'Ir': 'parameters/PARAMETERS_Ir(100)-(2x1)-O',
    'empty': 'parameters/PARAMETERS_empty',
    'stop': 'parameters/PARAMETERS_stop',
    'no_stop': 'parameters/PARAMETERS_stop_false',
    'missing_equals': 'parameters/PARAMETERS_missing_equals',
    'unknown_param': 'parameters/PARAMETERS_unknown',
    'left empty': 'parameters/PARAMETERS_left_side_empty',
    'no value': 'parameters/PARAMETERS_no_value',
    'typo': 'parameters/PARAMETERS_typo',
    }


def _fill_test_info(label_, info=None, needs_expected=True):
    """Fill parameters test info into a given, or new, TestInfo."""
    if info is None:
        info = TestInfo()
    try:
        info.parameters.expected = _READ[label_]
    except KeyError:
        if needs_expected:
            raise
    info.parameters.param_path = _PATHS[label_]
    return info


for label, test_info in _POSCAR_INFO.items():
    _fill_test_info(label, test_info)


def _get_all_parameters_infos():
    """Yield all TestInfo objects of known PARAMETERS files."""
    yield from _POSCAR_INFO.items()

    for label_ in _PATHS:
        if label_ in _POSCAR_INFO:
            continue
        try:
            info = _fill_test_info(label_)
        except KeyError:
            continue
        yield label_, info


_ALL_INFOS = dict(_get_all_parameters_infos())


class CasesParametersFile:
    """Collection of simple PARAMETERS file cases."""

    @parametrize(info=_ALL_INFOS.values(), ids=_ALL_INFOS)
    def case_parameters_file(self, info, data_path):
        """Return information regarding a PARAMETERS file."""
        return info, data_path/info.parameters.param_path

    def case_stop(self, data_path):
        """Return one PARAMETERS file with a line without an '=' sign."""
        info = _fill_test_info('stop')
        return self.case_parameters_file(info, data_path)

    def case_misses_equals(self, data_path):
        """Return one PARAMETERS file with a line without an '=' sign."""
        info = _fill_test_info('missing_equals')
        return self.case_parameters_file(info, data_path)

    @case(tags=CaseTag.RAISES)
    def case_unknown_parameter(self, data_path):
        """Return one PARAMETERS file with an unknown parameter."""
        info = _fill_test_info('unknown_param', needs_expected=False)
        return self.case_parameters_file(info, data_path)

    @case(tags=CaseTag.RAISES)
    def case_no_value_parameter(self, data_path):
        """Return one PARAMETERS file with a value-less parameter."""
        info = _fill_test_info('no value', needs_expected=False)
        return self.case_parameters_file(info, data_path)

    @case(tags=CaseTag.RAISES)
    def case_typo(self, data_path):
        """Return one PARAMETERS file with a typo in it."""
        info = _fill_test_info('typo', needs_expected=False)
        return self.case_parameters_file(info, data_path)


@parametrize(info=_POSCAR_INFO.values(), ids=_POSCAR_INFO)
def case_parameters_slab(info, make_poscar, data_path):
    """Return a slab, an Rparam (read from file), and TestInfo."""
    slab, _, info = make_poscar(info)
    rpars = parameters.read(data_path / info.parameters.param_path)
    return slab, rpars, info
