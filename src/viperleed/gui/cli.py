"""Module cli of viperleed.gui.

Defines the ViPErLEEDGUICLI class, the main entry point for starting
the ViPErLEED Graphical User Interface.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2020-01-11'
__license__ = 'GPLv3+'

from importlib import import_module
from itertools import chain
from pathlib import Path
import sys

from viperleed.cli_base import ViPErLEEDCLI
from viperleed.gui.base import catch_gui_crash
from viperleed.gui.constants import LOGO
from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.detect_graphics import QtDependencyFinder
from viperleed.gui.detect_graphics import check_pyqt_sanity
from viperleed.gui.detect_graphics import find_missing_qt_dependencies
from viperleed.gui.detect_graphics import has_pyqt
from viperleed.gui.detect_graphics import suppress_file_permission_warnings
from viperleed.gui.helpers import resources_path
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import SUPPORTED_BINDINGS
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError
from viperleed.gui.qt_compat.binding import get_binding
from viperleed.gui.qt_compat.binding import get_candidates

_QT_COMPAT = 'viperleed.gui.qt_compat'

# Part of error message emitted if a likely ABI mismatch is detected
# on the current Qt installation (e.g., Anaconda vs. pip)
_ABI_MISMATCH_MSG = '''
If you are executing viperleed in a conda environment, try one of the
following:
    1. Creating a new, clean environment without Qt by calling
       conda create with the --no-default-packages flag, then
            pip install "viperleed[GUI]"
       there.
    2. Deactivating the current environment first, if you have
       installed viperleed globally.
If none of the above works, or you're not in a conda environment, please
open an issue under https://github.com/viperleed/viperleed/issues.'''

# Part of error message emitted if the result of the sanity check
# may be wrong (e.g., on a slow machine with a display)
_SKIP_CHECK_MSG = '''
If your system has a display, you can start the GUI anyway with
    viperleed gui --skip-graphics-check'''

# Error messages emitted if the sanity check fails. {qt} is replaced
# by the Qt bindings that were looked for (e.g., 'PyQt6 or PyQt5').
_SANITY_TO_ERR_MSG = {
    PyQtSanity.NOT_FOUND: 'no Qt binding is installed (looked for {qt}).',
    # All of the following may still be caused by missing dependencies
    PyQtSanity.OK: '',
    PyQtSanity.IMPORT_ERROR: (
        'a Qt binding ({qt}) was found, but it could not be loaded.'
        + _ABI_MISMATCH_MSG
        ),
    PyQtSanity.RUNTIME_CRASH: (
        'a Qt binding ({qt}) was found, but the GUI crashed at startup.'
        + _ABI_MISMATCH_MSG
        ),
    PyQtSanity.NO_DISPLAY: (
        'the system appears to have no graphics capability (i.e., no '
        'monitor was detected).' + _SKIP_CHECK_MSG
        ),
    PyQtSanity.UNRESPONSIVE: (
        'loading Qt ({qt}) did not complete in time. '
        'Try once more, or check whether your system is very busy.'
        + _SKIP_CHECK_MSG
        ),
    }


class ViPErLEEDGUICLI(ViPErLEEDCLI, cli_name='gui'):
    """Main entry point for GUI. Both command-line and graphical."""

    def __call__(self, args=None):
        """Call either the CLI or graphical versions of the GUI."""
        args = self.parse_cli_args(args)
        if args.nogui:
            return commandline_main()
        if args.skip_graphics_check:
            # Also prevent later checks, e.g., in has_pyqt()
            check_pyqt_sanity.assume_ok()
        else:
            self.check_can_run_gui()
        return gui_main()

    def add_parser_arguments(self, parser):
        """Add CLI arguments for viperleed.gui to `parser`."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            '--nogui',
            help=('run the ViPErLEED graphical user interface in '
                  'command-line mode, i.e., without any windows'),
            action='store_true',
            )
        parser.add_argument(
            '--skip-graphics-check',
            help=('do not check whether Qt (PyQt5 or PyQt6) works and '
                  'the system has a display before starting the '
                  'graphical user interface. Use this if the check '
                  'wrongly reports that no display is available. If Qt '
                  'is broken, Python may crash without further '
                  'explanation'),
            action='store_true',
            )

    def check_can_run_gui(self):
        """Raise SystemExit if the GUI cannot execute.

        Raises
        ------
        SystemExit
            If the GUI cannot execute because
            - neither PyQt5 nor PyQt6 is installed,
            - Qt misses any system dependencies,
            - loading Qt or opening the GUI would crash Python,
            - the current machine has no graphics capability,
            - the QT_API environment variable has an invalid value.
        """
        err_msg = (
            'Cannot execute the ViPErLEED graphical user interface because {}'
            )
        # The check may take a few seconds, e.g., on systems without
        # a display. Tell users, and end the line before any errors.
        print('Checking graphics capabilities...', flush=True, end='')
        try:
            sanity = check_pyqt_sanity()
        except QtBindingNotFoundError as exc:  # Invalid QT_API
            print('Done', flush=True)
            self.parser.error(err_msg.format(exc))
        print('Done', flush=True)
        details = _SANITY_TO_ERR_MSG[sanity]  # Let KeyError propagate
        details = details.format(qt=_describe_candidates())
        if sanity is PyQtSanity.NOT_FOUND:
            # This is the only case in which it makes no sense to
            # check for missing dependencies. All other non-sane
            # configurations may be because of missing dependencies.
            self.parser.error(err_msg.format(details))
        missing = tuple(chain.from_iterable(
            find_missing_qt_dependencies().values()
            ))
        if missing:
            sep = '\n    '
            fmt_missing = sep.join(missing)
            fmt_missing = sep + fmt_missing
            how_to_install = QtDependencyFinder.find_install_for_libs(missing)
            install_msg = 'Try again after installing them'
            install_msg += ('.' if not how_to_install
                            else f' via{sep}{how_to_install}')
            self.parser.error(err_msg.format(
                'the following Qt dependencies are missing on your system:'
                f'{fmt_missing}\n{install_msg}'
                ))
        if sanity is not PyQtSanity.OK:
            self.parser.error(err_msg.format(details))


def _describe_candidates():
    """Return the names of the Qt bindings looked for, for messages."""
    try:
        candidates = get_candidates()
    except QtBindingNotFoundError:  # Invalid QT_API
        candidates = SUPPORTED_BINDINGS
    return ' or '.join(candidates)


def commandline_main():
    """Body of the command-line version of the GUI."""
    print('Running in command line version...', flush=True, end='')
    print('Not implemented yet.', flush=True)
    return 0


def gui_main():
    """Body of the GUI when running with graphics capability.

    Body of the functionality that invokes the ViPErLEED
    Graphical User Interface.
    """
    (qtc,
     qtg,
     qtw,
     plugin_selector) = import_graphics_modules()

    catch_gui_crash()

    print('Loading GUI...', flush=True, end='')
    if get_binding() == PYQT5:  # Always on (and not available) in Qt 6
        app_attribute = qtc.Qt.ApplicationAttribute
        qtg.QGuiApplication.setAttribute(app_attribute.AA_EnableHighDpiScaling)
        qtg.QGuiApplication.setAttribute(app_attribute.AA_UseHighDpiPixmaps)
    suppress_file_permission_warnings()  # Next line emits the warnings
    app = qtw.QApplication(sys.argv)
    app.setWindowIcon(qtg.QIcon(LOGO))

    # Import some fonts from ./fonts folder
    font_path = resources_path('gui/fonts')
    # * Text: family =  'DejaVu Sans'
    qtg.QFontDatabase.addApplicationFont(
        str(Path(font_path, 'DejaVuSans.ttf').resolve())
        )
    # * Math: family =  'CMU Serif'
    qtg.QFontDatabase.addApplicationFont(
        str(Path(font_path, 'cmunrm.otf').resolve())
        )

    plugin_selector_window = plugin_selector()
    plugin_selector_window.show()
    print('Done', flush=True)
    return app.exec()


def import_graphics_modules():
    """Dynamically import modules for graphical capability."""
    if not has_pyqt():
        raise ImportError('No Qt binding (PyQt5 or PyQt6) is available.')
    gui_select = import_module('viperleed.gui.selectplugin')
    return (
        import_module(f'{_QT_COMPAT}.QtCore'),  # Not graphics per-se
        import_module(f'{_QT_COMPAT}.QtGui'),
        import_module(f'{_QT_COMPAT}.QtWidgets'),
        gui_select.ViPErLEEDSelectPlugin,
        )
