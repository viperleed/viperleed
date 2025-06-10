"""Package special of viperleed.calc.classes.rparams.

This package defines classes for not-so-easy parameters to be used
as attributes for an Rparams object. These are, for example:
- parameters that would otherwise need to have variable types and
  require a bunch of isinstance checks when used
- parameters that have peculiar structure
- parameters that have multiple values that are cleaner when stored
  as attributes
- parameters that have a non-so-simple conversion to/from string
  or Assignment, i.e., those for which a ParameterInterpreter
  interpret_<...> method would be very complex.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-21'
__license__ = 'GPLv3+'
