"""Module constants of viperleed.calc.

Defines constants used throughout the viperleed.calc package. Typical
examples are default names of files and folders. Part of the contents
of this module used to be in calc.__init__ and calc.sections.cleanup.
They have been moved here to prevent cyclic imports.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-03'
__license__ = 'GPLv3+'

COMPILE_LOGS_DIRNAME = 'compile_logs'    # Save compilation logs here
DEFAULT_DELTAS = 'Deltas'    # Where delta-amplitude files are stored
DEFAULT_DOMAIN_FOLDER_PREFIX = 'Domain_'  # Auto-generated subfolders
DEFAULT_HISTORY = 'history'
DEFAULT_OUT = 'OUT'
DEFAULT_SUPP = 'SUPP'
DEFAULT_TENSORS = 'Tensors'  # Where Tensor files are stored
DEFAULT_WORK = 'work'
DEFAULT_WORK_HISTORY = 'workhistory'
LOG_PREFIX = 'viperleed-calc'
ORIGINAL_INPUTS_DIR_NAME = 'original_inputs'
SKIP_IN_DOMAIN_MAIN = (
    # Input files that should not be present in the
    # main work directory in a multi-domain calculation.
    # They are all those that are dependent on the structure.
    'PHASESHIFTS',
    'POSCAR',
    'VIBROCC',
    )
