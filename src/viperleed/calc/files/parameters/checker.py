"""Module checker of viperleed.calc.files.parameters.

Defines the ParametersChecker class, useful for checking that parameters
read from a PARAMETERS file do not clash with one another.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-25'
__license__ = 'GPLv3+'

import logging
import operator

from viperleed.calc.lib.base import available_cpu_count   # For N_CORES
from viperleed.calc.lib.string_utils import parent_name

from .errors import ParameterConflictError


_LOGGER = logging.getLogger(parent_name(__name__))
_PYTHON_RFACTOR_INTPOL_DEG = 3


# About the disable: ParametersChecker has only one single job:
# calling its own '_check_xxx' methods. Hence it makes sense to
# have a single entry point (check_parameter_conflicts).
# pylint: disable-next=too-few-public-methods
class ParametersChecker:
    """A container of methods for checking clashes between PARAMETERS.

    To add a new check, it is sufficient to implement a method
    whose name starts with '_check_'. The methods can use the
    self._rpars object to access parameter values.
    """

    def __init__(self):
        """Initialize instance."""
        self._rpars = None

    def check_parameter_conflicts(self, rpars):
        """Make sure the parameters in `rpars` are consistent."""
        self._rpars = rpars
        checker_names = (method_name for method_name in dir(self)
                         if method_name.startswith('_check_'))
        for checker_name in checker_names:
            checker = getattr(self, checker_name)
            checker()

    def _get_attribute_values(self, *attr_names):
        """Yield values of rpars attributes with given names."""
        yield from (getattr(self._rpars, attr) for attr in attr_names)

    def _check_element_name_collisions(self):
        """Make sure ELEMENT_MIX and ELEMENT_RENAME are not in conflict."""
        params = 'ELEMENT_RENAME', 'ELEMENT_MIX'
        values = self._get_attribute_values(*params)
        conflict = operator.and_(*(v.keys() for v in values))
        if conflict:
            message = ('The same POSCAR element cannot appear in both. '
                       'Conflicting elements: ' + ', '.join(conflict))
            raise ParameterConflictError(*params, message=message)

    def _check_and_update_fortran_comp(self):                                   # TODO: it may be necessary to add a .exe on Windows
        """Update FORTRAN_COMP parameter with info from pre and post."""
        # First the 'standard' compiler
        pre, post = self._rpars.FORTRAN_COMP
        if not post and pre in {'ifort', 'gfortran'}:
            self._rpars.getFortranComp(comp=pre)

        # Then MPI
        pre, post = self._rpars.FORTRAN_COMP_MPI
        if not post and pre in {'mpifort', 'mpiifort'}:
            self._rpars.getFortranMpiComp(comp=pre)

    def _check_intpol_deg_python_rfactor(self):
        """Make sure INTPOL_DEG is cubic if using the Python R factor."""
        intpol_deg, legacy = self._get_attribute_values('INTPOL_DEG',
                                                        'R_FACTOR_LEGACY')
        if legacy or intpol_deg == _PYTHON_RFACTOR_INTPOL_DEG:
            return
        message = (
            'The Python R factor supports only cubic interpolation '
            f'(INTPOL_DEG = {_PYTHON_RFACTOR_INTPOL_DEG}). Set INTPOL_DEG '
            f'= {_PYTHON_RFACTOR_INTPOL_DEG}, or set R_FACTOR_LEGACY = True '
            'to use the TensErLEED R factor.'
            )
        raise ParameterConflictError('INTPOL_DEG', 'R_FACTOR_LEGACY',
                                     message=message)

    def _check_smoothing_ignored_python_rfactor(self):
        """Warn if R_FACTOR_SMOOTHING is set but the R factor is Python."""
        smoothing, legacy = self._get_attribute_values('R_FACTOR_SMOOTHING',
                                                       'R_FACTOR_LEGACY')
        if smoothing and not legacy:
            _LOGGER.warning(
                f'R_FACTOR_SMOOTHING is set to {smoothing}, but it is '
                'ignored by the Python R factor. Set R_FACTOR_LEGACY = '
                'True to smooth the experimental curves with the '
                'TensErLEED R factor.'
                )

    def _check_n_cores(self):
        """Warn if N_CORES is larger than the number of available CPUs."""
        n_cores = self._rpars.N_CORES
        available = available_cpu_count()
        if not n_cores or available <= 0:
            # Will auto-detect, or detection failed
            return
        if n_cores > available:
            _LOGGER.warning(
                f'The N_CORES parameter is set to {n_cores}, which is larger '
                f'than the number of available processors ({available}). This '
                'may not be ideal for performance, or may render execution '
                'impossible. Consider reducing N_CORES.'
                )
