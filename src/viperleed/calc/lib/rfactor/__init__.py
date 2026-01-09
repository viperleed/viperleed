"""Module R-factor."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-02-21'
__license__ = 'GPLv3+'


from .pendry import R_pendry
from .r_1 import R_1
from .r_2 import R_2
from .smooth import R_s
from .zannazi_jona import R_zj

_R_FACTOR_SYNONYMS = {
    R_pendry: ('pendry', 'r_p', 'r_pendry', 'rp', 'pendry r-factor', 'p'),
    R_1: ('r1', 'r_1', 'r1 factor'),
    R_2: ('r2', 'r_2', 'r2 factor'),
    R_s: ('s', 'rs', 'r_s', 'smooth', 'schmid'),
    R_zj: (
        'zj',
        'zj factor',
        'zannazi',
        'zannazi jona',
        'zannazi-jona',
    ),
}

def select_rfactor(name):
    """Select R-factor function by name using synonyms."""
    _rfactor_name = name.lower().strip()
    for func, synonyms in _R_FACTOR_SYNONYMS.items():
        if _rfactor_name in synonyms:
            return func
    err_msg = f'Unknown R-factor name: {name}'
    raise ValueError(err_msg)
