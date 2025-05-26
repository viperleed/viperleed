"""Module mpl_graphics of viperleed.gui.

Defines functionality useful for the interaction between the ViPErLEED
Graphical User Interface and matplotlib. Notice that this module
behaves somewhat specially concerning imports: all imports are
dynamic. This is meant such that Qt-related modules are only
imported if appropriate.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-26'  # Parts of this were in basewidgets
__license__ = 'GPLv3+'

from enum import Enum
import importlib

from viperleed.gui.detect_graphics import has_pyqt


NO_CAIRO = (0, 0, 0)
MIN_CAIRO = (1, 18, 0)  # See Issue #9. Surely 1.17.2 fails


class MatplotLibBackend(Enum):
    """Supported back-ends for matplotlib."""

    # module-to-import, string-for-matplotlib.use
    AGG = (None, 'qt5agg')
    CAIRO = ('mplcairo', 'module://mplcairo.qt')
    DEFAULT = CAIRO

    @property
    def module(self):
        """Return the name of the module to be imported."""
        return self.value[0]

    @property
    def mpl_use(self):
        """Return the string for matplotlib.use."""
        return self.value[1]


def find_and_import_backend():
    """Determine a suitable backend for matplotlib."""
    # Notes on mplcairo: It used to be necessary to import mplcairo
    # before matplotlib. However, as of 2025, I see no more references
    # in the documentation concerning this need. If this turns out to
    # be a problem, we can test for matplotlib being imported using
    # sys.modules, and importlib.reload() it after importing mplcairo.
    target = MatplotLibBackend.DEFAULT
    if target.module:
        try:
            importlib.import_module(target.module)
        except ImportError:
            # Back-end not found. Fall back to AGG
            return MatplotLibBackend.AGG
    if target is MatplotLibBackend.CAIRO:
        cairo_version = get_cairo_version()
        if cairo_version < MIN_CAIRO:
            return MatplotLibBackend.AGG
    return target


def get_cairo_version():
    """Return the cairo version used by mplcairo."""
    try:
        cairo = importlib.import_module('mplcairo')
    except ImportError:
        return NO_CAIRO
    try:
        version_str = cairo.get_versions()['cairo']
    except KeyError:
        return NO_CAIRO
    version_str, *_ = version_str.split('@')
    version_tuple = tuple(version_str.strip().split('.'))
    try:
        return tuple(int(p) for p in version_tuple)
    except ValueError:
        return NO_CAIRO


def import_figure_canvas():
    """Return an appropriate FigureCanvas from the imported matplotlib."""
    if not has_pyqt():
        raise ImportError('PyQt5 not found.')
    backend = import_matplotlib()
    if backend is MatplotLibBackend.AGG:
        module = importlib.import_module('matplotlib.backends.backend_qt5agg')
        return module.FigureCanvas
    if backend is MatplotLibBackend.CAIRO:
        module = importlib.import_module('mplcairo.qt')
        return module.FigureCanvasQTCairo
    raise ImportError(f'Unsupported backend {backend}')


def import_matplotlib():
    """If PyQt5 is available, import matplotlib and set its backend.

    Returns
    -------
    backend : MatplotLibBackend or None
        The backend that is currently in use.
        None if PyQt5 is not present.
    """
    if not has_pyqt():
        return None
    backend = find_and_import_backend()
    matplotlib = importlib.import_module('matplotlib')
    matplotlib.use(backend.mpl_use)
    return backend
