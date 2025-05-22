"""Package lib of viperleed.calc.

Defines generic functionality used throughout the viperleed.calc
package.

Modules
-------
base
    Non-LEED-related classes and functions.
checksums
    Functions used to prevent user tinkering with FORTRAN source files.
    The actual checksums are contained in _checksums.dat and should
    only be edited by expert contributors.
coordinates
    Functions to handle atomic coordinates.
dataclass_utils
    Extension of the datclasses standard-library module,
    with additions for backward compatibility to Python 3.7.
fortran_utils
    Functions for interacting and manipulating FORTRAN code.
itertools_utils
    Extension of the itertools standard-library module,
    with additions for backward compatibility to Python 3.7.
leedbase
    LEED- and TensErLEED-related base functions used in calc.
log_utils
    Functions useful to handle the logging standard-library
    module and its `Logger`s.
math_utils
    Mathematical functions.
matplotlib_utils
    Extension of the matplotlib package, with additions for
    backward compatibility to Python 3.7.
matrix
    Functions and exceptions related to matrices, as well as
    matrix transforms.
parallelization
    Functions to help parallelize code execution.
periodic_table
    Collection of basic atom data.
sequence_utils
    Functions specifically intended for generic sequences.
string_utils
    Functions for reading, writing, and manipulating strings.
time_utils
    Time formats used in viperleed.calc and related functions.
    Timers for code execution monitoring.
"""
# NB: woods_notation is undocumented on purpose, as it will be
# replaced with functionality in guilib.

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-21'
__license__ = 'GPLv3+'
