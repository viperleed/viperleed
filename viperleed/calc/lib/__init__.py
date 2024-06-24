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
leedbase
    LEED- and TensErLEED-related base functions used in calc.
periodic_table
    Collection of basic atom data.
"""
# NB: woods_notation is undocumented on purpose, as it will be
# replaced with functionality in guilib.

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-21'
__license__ = 'GPLv3+'
