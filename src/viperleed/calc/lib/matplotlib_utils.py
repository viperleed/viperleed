"""Functions useful for handling matplotlib.

In particular, stuff that makes new matplotlib features available in
older python versions.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-10'
__license__ = 'GPLv+3'

from functools import wraps
import logging
from pathlib import Path
import sys
from traceback import format_exc

if sys.version_info >= (3, 10):
    import importlib.resources as importlib_resources
else:
    # According to matplotlib.style.core, even though Py3.9 has
    # importlib.resources, it doesn't properly handle modules
    # added in sys.path.
    import importlib_resources

try:
    from matplotlib import __version__ as mpl_version
except ModuleNotFoundError:
    CAN_PLOT = False
else:
    import matplotlib
    from matplotlib import _rc_params_in_file
    from matplotlib import style as mpl_style
    from matplotlib.style.core import STYLE_EXTENSION
    CAN_PLOT = True

from viperleed.calc.lib.version import Version


_LOGGER = logging.getLogger(__name__)


class HasNoMatplotlibError(Exception):
    """No matplotlib installation was found."""

    _default_msg = ('Cannot proceed without matplotlib installed. Run '
                    'python -m pip install matplotlib, then try again')

    def __init__(self, message=None):
        super().__init__(message or self._default_msg)


# ###########################  DECORATORS  ############################

def log_without_matplotlib(logger, level='debug', msg=''):
    """Emit a message to logger if matplotlib is not available, then return.

    Parameters
    ----------
    logger : Logger
        The logger that should emit the message.
    level : str, optional
        Which of the methods of `logger` to use. Default is 'debug'.
    msg : str, optional
        Extra information to be logged. Default is an empty string.

    Returns
    -------
    decorator : callable
        A decorator ready to be applied to any function. If
        no matplotlib is found, the decorated function logs
        the desired message, then returns None. Otherwise,
        it returns the same as the decorated function.
    """
    log_emit = getattr(logger, level)
    log_msg = 'Necessary modules for plotting not found.'
    if msg:
        log_msg += f' {msg}'

    def _decorator(func):
        @wraps(func)
        def _wrapper(*args, **kwargs):
            if not CAN_PLOT:
                log_emit(log_msg)
                return None
            return func(*args, **kwargs)
        return _wrapper
    return _decorator


def skip_without_matplotlib(func):
    """Decorate func, and return early if matplotlib is not available."""
    @wraps(func)
    def _wrapper(*args, **kwargs):
        if not CAN_PLOT:
            return None
        return func(*args, **kwargs)
    return _wrapper


def raise_without_matplotlib(func):
    """Decorate func, and complain if matplotlib is not available."""
    @wraps(func)
    def _wrapper(*args, **kwargs):
        if not CAN_PLOT:
            raise HasNoMatplotlibError
        return func(*args, **kwargs)
    return _wrapper


# ############################  FUNCTIONS  ############################

@raise_without_matplotlib
def close_figures(pyplot, *figures):
    """Close pyplot figures safely."""
    _log_msg = ('Failed to close figure. Please report this to the '
                'ViPErLEED developers: %s')
    for figure in figures:
        try:
            pyplot.close(figure)
        except Exception:                                                       # TODO: should catch correct ones, but let's have people tell us.
            _LOGGER.warning(_log_msg, format_exc())


@skip_without_matplotlib
def prepare_matplotlib_for_calc():
    """Prepare matplotlib for use in viperleed.calc."""
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')                                                       # TODO: check with Michele if this causes conflicts
    use_calc_style()


@skip_without_matplotlib
def use_calc_style():
    """Use viperleed.calc matplotlib style for plotting."""
    if Version(mpl_version) >= '3.7':
        mpl_style.use('viperleed.calc')
        return
    # In the early versions of matplotlib, dotted names cannot be
    # used to locate styles. Use here a simplified version of the
    # implementation of more recent matplotlib (v3.7.5)
    vpr_path = Path(importlib_resources.files('viperleed'))
    path = vpr_path / f'calc.{STYLE_EXTENSION}'

    # Implementation of _rc_params_in_file differs slightly in older
    # matplotlib versions. Hopefully it works anyway.
    mpl_style.use(_rc_params_in_file(path))
