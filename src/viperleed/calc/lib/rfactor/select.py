"""Module select of viperleed.calc.lib.rfactor.

Defines functions that pick the function computing an R factor.
"""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType

from . import pendry
from . import r_1
from . import r_2
from . import smooth
from . import zannazi_jona

_R_FACTOR_FUNCTIONS = {
    RFactorType.PENDRY: pendry.r_pendry,
    RFactorType.R2: r_2.r_2,
    RFactorType.ZJ: zannazi_jona.r_zjj,
    RFactorType.SMOOTH: smooth.r_s,
    }

# The R1 R factor is kept only for legacy reasons. R_FACTOR_TYPE
# cannot select it (RFactorType has no R1 member): only these names,
# passed to select_rfactor_by_name (i.e., which_r of from_ase), can.
_R1_NAMES = ('r1', 'r_1', 'r1 factor')


def select_rfactor(r_factor_type):
    """Return the function that computes an R factor of a given type.

    Parameters
    ----------
    r_factor_type : RFactorType or str or int
        The type of R factor. Any value accepted by RFactorType.

    Returns
    -------
    r_func : callable
        The function computing the R factor.

    Raises
    ------
    ValueError
        If `r_factor_type` is not a valid RFactorType.
    """
    return _R_FACTOR_FUNCTIONS[RFactorType(r_factor_type)]


def select_rfactor_by_name(name):
    """Return the function that computes an R factor called `name`.

    Contrary to select_rfactor, this function can also select the
    R1 R factor. It is meant for user-facing functions that take a
    name rather than the R_FACTOR_TYPE parameter.

    Parameters
    ----------
    name : str
        The name of the R factor. Any value accepted by RFactorType,
        or one of 'r1', 'r_1', 'r1 factor' (case-insensitive).

    Returns
    -------
    r_func : callable
        The function computing the R factor.

    Raises
    ------
    ValueError
        If `name` is not a known R-factor name.
    """
    if name.strip().lower() in _R1_NAMES:
        return r_1.r_1
    try:
        return select_rfactor(name)
    except ValueError:
        raise ValueError(f'Unknown R-factor name: {name}') from None
