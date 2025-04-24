"""Cases for test_parameters."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-09-07'
__license__ = 'GPLv3+'

from copy import deepcopy

import numpy as np
from pytest_cases import case, parametrize

from viperleed.calc.classes.rparams.special.l_max import LMax
from viperleed.calc.files import parameters

from ...poscar_slabs import CasePOSCARSlabs
from ...poscar_slabs import POSCARS_WITHOUT_INFO
from ...poscar_slabs import get_info_by_name
from ...tags import CaseTag
from ...testinfo import TestInfo


def _get_poscar_info(name):
    """Return a POSCAR info from a multitude of possible sources."""
    try:
        return get_info_by_name(name)
    except StopIteration:  # Not one of those in the default container
        pass

    # See if there's an explicit case with that name
    try:
        case_func = getattr(CasePOSCARSlabs(), name)
    except AttributeError:
        pass
    else:
        *_, info = case_func()
        return info

    try:
        return get_info_by_name(name, container=POSCARS_WITHOUT_INFO)
    except StopIteration:  # Also not one without info
        pass
    raise ValueError(f'{name!r} not found in poscar_slabs')


_POSCAR_INFO = {
    'Ag': _get_poscar_info('case_poscar_ag100'),
    'Ir': _get_poscar_info('Ir(100)-(2x1)-O'),
    'PtRh': _get_poscar_info('POSCAR_Pt25Rh75(100)-p(3x1)-O'),
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
    'PtRh': {'RUN': [0],
           'THEO_ENERGIES': [49, 550, 3], 'ELEMENT_MIX': {'Me': ['Rh', 'Pt']},
           'BULK_LIKE_BELOW': 0.45, 'T_DEBYE': 250, 'T_EXPERIMENT': 300,
           'SITE_DEF': {
               'O': {'surf': {1, 2}},
               'Me': {'layer_1': {3, 4, 5}, 'layer_2': {6, 7, 8},
                      'layer_3': {9, 10, 11}, 'layer_4': {12, 13, 14}},
               },
           'VIBR_AMP_SCALE': ['*surf 1.3',],
           },
    'domains': {'RUN': [4], 'LOG_LEVEL': 10, 'THEO_ENERGIES': [50, 152, 3],
                'V0_IMAG': 5.0, 'LMAX': LMax(8, 12),
                'DOMAIN': [('', 'silver'), ('Bi', 'bismuth')]},
    'domains_identical': {},
    'empty': {}, 'stop': {}, 'no_stop': {}, 'missing_equals': {},
    'left empty': {},
    }
for key in ('stop', 'no_stop', 'left empty'):
    _READ[key].update(**_READ['Ir'])  # Copies of the Ir file
_READ['missing_equals'].update(**_READ['Ag'])
del _READ['missing_equals']['SUPERLATTICE']
_READ['domains_identical'].update(**_READ['domains'])
_READ['domains_identical']['DOMAIN'] = [('', 'silver'), ('', 'silver')]

# FORTRAN_COMP added on a copy of the Ag(100) file
_READ['fortran comp'] = deepcopy(_READ['Ag'])
_READ['fortran comp']['FORTRAN_COMP'] = ['gfortran',
                                         '-llapack -lpthread -lblas']

_PATHS = {
    'Ag': 'Ag(100)/initialization/PARAMETERS',
    'Ir': 'parameters/PARAMETERS_Ir(100)-(2x1)-O',
    'PtRh': 'parameters/PARAMETERS_Pt25Rh75(100)-p(3x1)-O',
    'domains': 'parameters/PARAMETERS_domains',
    'domains_identical': 'parameters/PARAMETERS_domains_identical',
    'empty': 'parameters/PARAMETERS_empty',
    'fortran comp': 'parameters/PARAMETERS_fortran_comp',
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

    def case_domains(self, data_path):
        """Return the main PARAMETERS file of a multi-domain calculation."""
        info = _fill_test_info('domains')
        return self.case_parameters_file(info, data_path)

    def case_domains_identical(self, data_path):
        """Return the main PARAMETERS file of a two-same-domain calculation."""
        info = _fill_test_info('domains_identical')
        return self.case_parameters_file(info, data_path)

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
