"""Module lib of viperleed.gui.widgets.

Library of functions that are common to several Qt objects.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-12'
__license__ = 'GPLv3+'

from datetime import datetime
import logging
import re
import sys
import traceback as _m_traceback
import warnings

import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed.gui.decorators import print_thread


###############################################################################
#                                   FUNCTIONS                                 #
###############################################################################


# See also https://fman.io/blog/pyqt-excepthook/
def catch_gui_crash(base_log_path):
    """Show unhandled exceptions to the user."""

    sys._excepthook = sys.excepthook

    class _HOOK(qtc.QObject):
        """Exception hook that runs in the correct thread."""

        # Prepare a QMessageBox
        _msg = qtw.QMessageBox()
        _msg.setWindowTitle("Unhandled exception")
        _msg.setIcon(_msg.Critical)
        btn = _msg.addButton(_msg.Ok)
        btn.setText("Close Application")

        # Prepare the info text
        __repo = 'viperleed'
        __issue = ('<a href=\"https://github.com/viperleed/'
                   f'{__repo}/issues\">GitHub Issue</a>')
        __email = ('<a href="mailto:riva@iap.tuwien.ac.at">'
                   'riva@iap.tuwien.ac.at</a>')
        __base_text = ("An unhandled exception occurred. This is most likely "
                       " a bug.<br><br>Please provide us with an appropriate "
                       f"bug report as a {__issue} or via email to {__email}."
                       "<br><br> The application will now close.")


        # Signal to call GUI-related methods in the right thread
        __show_msg_requested = qtc.pyqtSignal(str, str)

        def __init__(self):
            """Initialize instance."""
            super().__init__()
            self._msg.finished.connect(self.__quit_app)

            self.excinfo = ()

            # Logger is set up and created only when an exception occurs
            self.loginfo = {'log': None, 'fname': None}

        @print_thread
        def __call__(self, exctype, value, traceback):
            """Show the message box to the user, then exit normally."""
            print("called exception hook", flush=True)
            self.excinfo = (exctype, value, traceback)
            info = "".join(
                   _m_traceback.format_exception(exctype, value, traceback)
                   )
            logger = self.loginfo['log']
            if not logger:
                _now = datetime.now().strftime("%Y%m%d_%H%M%S")
                _logname = base_log_path / f"exceptions_{_now}.log"
                self.loginfo['fname'] = _logname
                logging.basicConfig(filename=_logname, level=logging.ERROR)      # TODO: probably also some more complicated form
                self.loginfo['log'] = logger = logging.getLogger()
            logger.error(info)
            text = (self.__base_text
                    + "<br><br>Please include file "
                    f"{self.loginfo['fname']} in the report.")
            try:
                self.__call_impl(text, info)
            except RuntimeError:
                # Probably the C++ QObject underlying self is gone
                pass

        def __call_impl(self, text, info):
            """Finish calling self."""
            try:
                self.__show_msg_requested.disconnect()                          # TODO: use safe_disconnect
            except TypeError:
                # Not connected
                pass
            self.__show_msg_requested.connect(self.__show_message)
            self.__show_msg_requested.emit(text, info)

        @qtc.pyqtSlot(str, str)
        def __show_message(self, text, detailed_text):
            """Show a message box to the user."""
            self._msg.setText(text)
            self._msg.setDetailedText(detailed_text)
            self._msg.exec_()

        @qtc.pyqtSlot(int)
        def __quit_app(self, _):
            """Quit application after exception."""
            sys._excepthook(*self.excinfo)
            qtw.qApp.quit()

    sys.excepthook = _HOOK()


def change_control_text_color(ctrl, color):
    """Change colour of the text in a control.

    Typically useful for signalling an error
    on user input on some of the controls.

    Parameters
    ----------
    ctrl_name : QWidget
        The control whose colour has to be changed.
    color : QColor or str
        Color to be used

    Returns
    -------
    None.
    """
    color = qtg.QColor(color)
    palette = ctrl.palette()
    if isinstance(ctrl, qtw.QPushButton):
        palette.setColor(palette.ButtonText, color)
    elif isinstance(ctrl, qtw.QLabel):
        palette.setColor(palette.WindowText, color)
    else:
        palette.setColor(palette.Text, color)
    ctrl.setPalette(palette)


def editStyleSheet(qwidget, new_entries):
    """
    Append or replace new_entries in the styleSheet of qwidget

    Parameters
    ----------
    qwidget: QWidget
             QWidget whose style sheet will be edited
    new_entries: str
                 list of properties and values that are to be inserted in the
                 style sheet.
                 Format: "property1: value1; property2: value2; ..."
    """
    if not isinstance(qwidget, qtw.QWidget):
        raise ValueError("Argument 0 of editStyleSheet must be a QWidget "
                         f"subclass. Found {type(qwidget)} instead")

    if not isinstance(new_entries, str):
        raise TypeError("Argument 1 of editStyleSheet must be a string. "
                        f"Found {type(new_entries)} instead")

    entry_re = r"""
        (?P<property>[-a-z ]+)              # property:
                                            #    lowercase letters, '-' and ' '
        :                                   # colon
        (?P<value>[\w\%\#\(\),:.-_\\\/ ]+)  # value:
                                            #    alphanumeric, %, #,
                                            #    parentheses, ' ' and
                                            #    .,:-_ \/
        """
    checkFormat = re.match(r"".join((r'^(', entry_re, r';)+$')),
                           new_entries, re.VERBOSE)
    if checkFormat is None:
        raise ValueError("Argument 1 of editStyleSheet must be a string "
                         "of the form 'property1: value1; "
                         "property2: value2; ...'")

    if len(qwidget.styleSheet()) == 0: # no style sheet present
        qwidget.setStyleSheet(new_entries)
        return

    # prepare a compiled re that will be used on all entries
    entry_format = re.compile(entry_re, re.VERBOSE)

    new_entries = (entry
                   for entry in re.split(';', new_entries)
                   if len(entry) > 0)

    for entry in new_entries:
        e_match = entry_format.match(entry)
        property = e_match.group('property')
        value = e_match.group('value')
        # check whether the entry is already present
        if property in qwidget.styleSheet():
            # update what's after "property:" till ";" (included) with value
            start = qwidget.styleSheet().index(property) + len(property) + 1
            # +1 to skip the ':'
            end = qwidget.styleSheet().index(';', start) + 1
            old = qwidget.styleSheet()
            qwidget.setStyleSheet("".join((old[:start], value, old[end:])))
        else: # otherwise simply append the new property
            qwidget.setStyleSheet(qwidget.styleSheet() + entry)


def drawText(painter, text, transform=None, combine=False):
    # reimplementation of text painting function that does not give awful
    # output. type(painter)==QPainter, type(p)==QPointF, type(text)==QString

    rawFont = qtg.QRawFont.fromFont(painter.font())
    indexes = rawFont.glyphIndexesForString(text)

    painter.save()
    paths = [rawFont.pathForGlyph(index) for index in indexes]
    advances = rawFont.advancesForGlyphIndexes(indexes,
                                               qtg.QRawFont.UseDesignMetrics
                                               | qtg.QRawFont.KernedAdvances)
    if transform is not None:
        painter.setWorldTransform(transform, combine=combine)
    for (path,advance) in zip(paths,advances):
        painter.fillPath(path, painter.pen().brush())
        painter.translate(advance)
    painter.restore()


def get_all_children_widgets(parent, exclude_parent=False, recursive=True):
    """
    This is an extension of the QObject.children() method that finds
    all the children QWidgets of a QObject.

    Parameters
    ----------
    parent: QWidget or QLayout (or subclasses)
    exclude_parent: bool, default = False
                    If True parent is not contained in the returned set
    recursive: bool, default = True
               If False, only first-generation widgets are returned
               If True, all-generations widgets are included

    Returns
    -------
    list of QWidget
        - if parent is a QWidget:
          All QWidgets that are children of any depth of parent
        - if parent is a QLayout:
          All the QWidgets managed by parent and all their children of any depth
    """
    if not isinstance(parent, (qtw.QWidget, qtw.QLayout)):
        return set()

    children = set(parent.children())

    # find all widgets that are directly children of parent
    childrenWidgs = {child for child in [*children, parent]
                      if isinstance(child, qtw.QWidget)}
    # find all layouts that are children of parent
    childrenLays = {child for child in [*children, parent]
                      if isinstance(child, qtw.QLayout)}

    # add the widgets that are managed by the layouts in the list of children
    for lay in childrenLays:
        to_add = (lay.itemAt(idx).widget()
                  for idx in range(lay.count())
                  if isinstance(lay.itemAt(idx), qtw.QWidgetItem))
        childrenWidgs.update(widg
                             for widg in to_add
                             if isinstance(widg, qtw.QWidget))

    if recursive:
        # and run recursively to find all the nested children
        for child in childrenWidgs.copy():
            if child != parent:
                childrenWidgs.update(get_all_children_widgets(child))

    if exclude_parent:
        childrenWidgs.discard(parent)

    return childrenWidgs


def move_to_front(window):  # TODO: move to a nicer place
    """Move a window to the front."""
    window.show()
    # Move window to front
    window.raise_()
    # Show as a window, also in case it is minimized
    window.setWindowState(window.windowState()
                          & ~qtc.Qt.WindowMinimized
                          | qtc.Qt.WindowActive)


def screen_fraction(obj, size):
    """Return the fraction of the screen of obj occupied by size."""
    try:
        scr_size = obj.window().windowHandle().screen().availableSize()
    except AttributeError:
        return -1
    return max(size.width()/scr_size.width(), size.height()/scr_size.height())


def raise_on_qt_messages():
    """Produce warnings and Exceptions instead of Qt messages.

    Warns
    -----
    QtDebug
        For QtCore.QtDebugMsg
    QtWarning
        For QtCore.QtWarningMsg
    QtInfo
        For QtCore.QtInfoMsg

    Raises
    ------
    QtCritical
        For QtCore.QtCriticalMsg or QtCore.QtSystemMsg
    QtFatal
        For QtCore.QtFatalMsg
    """

    _map = {qtc.QtDebugMsg: QtDebug,
            qtc.QtWarningMsg: QtWarning,
            qtc.QtInfoMsg: QtInfo,
            qtc.QtCriticalMsg: QtCritical,
            qtc.QtSystemMsg: QtCritical,
            qtc.QtFatalMsg: QtFatal}

    def __handler(severity, context, message):
        """Handler function for Qt messages.

        Changes the default behaviour to merely printing to stderr
        and rather raise appropriate exceptions/warnings.

        Parameters
        ----------
        severity : QtCore.QtMsgType
            Severity level of the message. Can be QtDebugMsg, QtInfoMsg,
            QtWarningMsg, QtCriticalMsg, QtFatalMsg, QtSystemMsg.
        context : QtCore.QMessageLogContext
            Information on the context in which the message was
            generated. It normally does not contain any information,
            except when using debug builds of Qt (and PyQt, and python,
            and all modules).
        message : str
            The original Qt message.

        Returns
        -------
        None.
        """
        if context.file:
            # There's sensible context information
            message += (f"line: {context.line}, func: {context.function}"
                        f"file: {context.file}")
        msg = _map[severity](message)

        if isinstance(msg, (QtDebug, QtWarning, QtInfo)):
            warnings.warn(msg, stacklevel=2)
            return
        raise msg

    qtc.qInstallMessageHandler(__handler)


def remove_spacing_and_margins(layout):
    """Remove spacing and margins from a layout."""
    layout.setSpacing(0)
    layout.setContentsMargins(0, 0, 0, 0)


def retain_size_when_hidden(widget):
    """Retain widget size when widget is not visible."""
    policy = widget.sizePolicy()
    policy.setRetainSizeWhenHidden(True)
    widget.setSizePolicy(policy)


################################################################################
#                                   CLASSES                                    #
################################################################################


class QtDebug(RuntimeWarning):
    pass


class QtWarning(RuntimeWarning):
    pass


class QtInfo(RuntimeWarning):
    pass


class QtCritical(Exception):
    pass


class QtFatal(Exception):
    pass


class QDoubleValidatorNoDot(qtg.QDoubleValidator):
    def validate(self, text:str, cursor_pos:int):
        text = text.replace(',', '.')
        return super().validate(text, cursor_pos)


class AllGUIFonts():  ## > Will handle in a different way!
    # May be worth removing the __init__ completely since all instances
    # have anyway the same 'value' for all fonts
    def __init__(self):
        self.buttonFont = qtg.QFont()
        self.smallTextFont = qtg.QFont()
        self.labelFont = qtg.QFont()
        self.smallButtonFont = qtg.QFont()
        self.largeTextFont = qtg.QFont()
        self.plotTitleFont = qtg.QFont()

        allfonts = [self.buttonFont,
                    self.smallTextFont,
                    self.labelFont,
                    self.smallButtonFont,
                    self.largeTextFont,
                    self.plotTitleFont]

        [x.setFamily("DejaVu Sans") for x in allfonts]

        self.buttonFont.setPointSize(10)
        self.smallTextFont.setPointSize(9)
        self.labelFont.setPointSize(10)
        self.smallButtonFont.setPointSize(8)
        self.largeTextFont.setPointSize(12)
        self.plotTitleFont.setPointSize(14)

        self.mathFont = qtg.QFont()
        self.mathFont.setFamily("CMU Serif")
        self.mathFont.setPointSize(15)
        self.mathFont.setStyleStrategy(qtg.QFont.PreferAntialias)
