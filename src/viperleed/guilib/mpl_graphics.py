"""Module mpl_graphics of viperleed.guilib.

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
import sys

from viperleed.guilib.detect_graphics import has_pyqt


PY38 = (3, 8)


class MatplotLibBackend(Enum):
    """Supported back-ends for matplotlib."""

    # module-to-import, string-for-matplotlib.use
    AGG = (None, 'qt5agg')
    CAIRO = ('mplcairo', 'module://mplcairo.qt')

    @property
    def module(self):
        """Return the name of the module to be imported."""
        return self.value[0]

    @property
    def mpl_use(self):
        """Return the string for matplotlib.use."""
        return self.value[1]

    @classmethod
    def get_default(cls):
        """Return the default backend."""
        return cls.CAIRO if sys.version_info < PY38 else cls.AGG


def find_and_import_backend():
    """Determine a suitable backend for matplotlib."""
    # Notes on mplcairo: It used to be necessary to import mplcairo
    # before matplotlib. However, as of 2025, I see no more references
    # in the documentation concerning this need. If this turns out to
    # be a problem, we can test for matplotlib being imported using
    # sys.modules, and importlib.reload() it after importing mplcairo.

    # TODO: there seems to be a bug of some sort in the way we're
    # using mplcairo. The UI crashes without any information when
    # setting mplcairo as a backend upon, e.g., loading an input
    # file. Tested on:
    #   - py3.7  (mpl 3.5.3,  cairo 0.5)
    #   - py3.8  (mpl 3.7.5,  cairo 0.6.1)
    #   - py3.9  (mpl 3.9.0,  cairo 0.6.1)
    #   - py3.12 (mpl 3.10.3, cairo 0.6.1)
    # For now, use the AGG backend for everything. When fixed, replace
    # the next line with
    #     target = MatplotLibBackend.get_default()
    target = MatplotLibBackend.AGG
    if target.module:
        try:
            importlib.import_module(target.module)
        except ImportError:
            # Back-end not found. Fall back to AGG
            return MatplotLibBackend.AGG
    return target


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
