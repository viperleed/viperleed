# -*- coding: utf-8 -*-
"""Module _checker of viperleed.tleedmlib.files.parameters.

Created on 2023-10-25

@author: Michele Riva (@michele-riva)

Defines the ParametersChecker class, useful for checking that parameters
read from a PARAMETERS file do not clash with one another.
"""

import operator

from .errors import ParameterConflictError


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
