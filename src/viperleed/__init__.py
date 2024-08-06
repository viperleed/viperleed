"""=================
    ViPErLEED
=================
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-03'
__license__ = 'GPLv3+'
__version__ = '0.12.2'


GLOBALS = {
    'USE_GUI': None,
    'version': __version__,
    'StateRecorder': None, # Recorder for calculation states
    'version_message': ('ViPErLEED (Vienna Package for Erlangen LEED) '
                        f'v{__version__}')
    }
# Name of environment variable specifying the path to the tensor-LEED
# source code. Typically the viperleed-tensorleed repository.
VIPERLEED_TENSORLEED_ENV = 'VIPERLEED_TENSORLEED'
