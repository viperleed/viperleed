"""Module pathselector of viperleed.gui.measure.widgets.

Defines the PathSelector class: a QWidget for picking the
path to a file or to a directory.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-09-07'
__license__ = 'GPLv3+'

import errno
import os
from pathlib import Path
import sys

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.widgets.lib import change_control_text_color


_ELLIPSIS = '\u2026'
_IS_WINDOWS = sys.platform.startswith('win')
_ILLEGAL_CHARS = '<>":|?*' if _IS_WINDOWS else ''
_ILLEGAL_NAMES = (
    'CON', 'PRN', 'AUX', 'NUL', *(f'COM{i+1}' for i in range(9)),
    *(f'LPT{i+1}' for i in range(9))
    ) if _IS_WINDOWS else ()
_ROOT = Path(os.environ.get('SYSTEMDRIVE', 'C:' if _IS_WINDOWS else '') + '/')
assert _ROOT.is_dir()


class PathSelector(qtw.QWidget):
    """A widget for selecting a path to file or directory."""

    path_changed = qtc.pyqtSignal(Path)  # New full path

    def __init__(self, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        select_file : bool, optional
            Whether this path selector is intended for picking the
            path to a file (if True) or the one to a directory (if
            False). Default is True.
        existing_file : bool, optional
            Whether this selector is intended for picking only files
            that already exist on the file system (if True), or whether
            it is used for saving files (if False). Default is True.
        read_only : bool, optional
            Whether the selected path is meant to be only read and
            not edited. Default is False.
        max_chars : int, optional
            Maximum number of characters that are to be displayed.
            Paths longer than this will be shown as elided in the
            middle. Values below 10 will be ignored. Default is 35.
        parent : QWidget, or None, optional
            Parent of this widget, as well as of the QFileDialog.
            Default is None.
        **kwargs : object
            Other optional keyword arguments passed on to QFileDialog

        Returns
        -------
        None.
        """
        self.__glob = {
            'select_file': kwargs.pop('select_file', True),
            'existing_file': kwargs.pop('existing_file', True),
            'read_only': kwargs.pop('read_only', False),
            'max_chars': max(10, kwargs.pop('max_chars', 35)),
            'kwargs': kwargs,
            'full_path': None,
            'elided_path': '',
            }
        super().__init__(kwargs.get('parent', None))

        _fdialog = qtw.QFileDialog
        selector = _fdialog.getExistingDirectory
        if self.__glob['select_file'] and self.__glob['existing_file']:
            selector = _fdialog.getOpenFileName
        elif self.__glob['select_file'] and not self.__glob['existing_file']:
            selector = _fdialog.getSaveFileName
        self.__selector = selector

        self.__lineedit = qtw.QLineEdit()
        self.__browse = qtw.QToolButton()
        browse_action = qtw.QAction('\u00b7'*3)
        browse_action.setToolTip('')
        self.__browse.setDefaultAction(browse_action)

        self.__compose()
        self.__connect()

    @property
    def parent_directory(self):
        """Return the directory containing the current selection."""
        if self.path is None:
            return Path().resolve()
        return self.path.parent.resolve()

    @property
    def path(self):
        """Return the currently selected path."""
        return self.__glob['full_path']

    @path.setter
    def path(self, new_path):
        """Set a new_path if acceptable."""
        if not self.is_acceptable_path(new_path):
            self.__show_path()  # restore old value
            return

        # If used for opening existing files, the
        # path should point to an existing file
        new_path = Path(new_path)
        if (self.__selector is qtw.QFileDialog.getOpenFileName
                and not new_path.is_file()):
            self.__show_path()  # restore old value
            return

        old_path = self.path
        self.set_path(new_path)
        self.__glob['kwargs']['directory'] = str(self.parent_directory)
        if old_path is not None and self.path != old_path:
            self.path_changed.emit(self.path)

    @property
    def read_only(self):
        """Return whether this selector is used only for display."""
        return self.__glob['read_only']

    @read_only.setter
    def read_only(self, read_only):
        """Set this widget as read-only or as edit-enabled."""
        was_read_only = self.read_only
        if was_read_only == bool(read_only):
            # No change
            return
        self.__glob['read_only'] = bool(read_only)
        self.__browse.setVisible(not self.read_only)
        self.__lineedit.setReadOnly(self.read_only)

        # Add/remove the spacer item between line edit and browse
        # It is always at position 1, 0 is self.__lineedit.
        layout = self.layout()
        if self.read_only:
            layout.removeItem(layout.itemAt(1))
        else:
            layout.insetSpacing(1, 2)

    def set_path(self, new_path):
        """Set the contents of this path selector."""
        self.__glob['full_path'] = Path(new_path)
        self.__lineedit.setToolTip(str(new_path))
        full_path = self.__glob['full_path']
        selects_file = self.__glob['select_file']
        if not selects_file:
            valid = full_path.is_dir()
        elif selects_file and self.__glob['existing_file']:
            valid = full_path.is_file()
        else:  # A file to be selected for saving
            valid = True
        color = 'black' if valid else 'red'
        change_control_text_color(self.__lineedit, color)
        self.__show_path()

    def get_posix_path(self):
        """Return the selected path as a string."""
        return str(self.path.as_posix())

    def __compose(self):
        """Place children widgets."""
        # Set up the line edit: min size from max chars
        max_ = self.__glob['max_chars']
        line = self.__lineedit
        edit_width = line.fontMetrics().boundingRect('m'*max_).width()
        line.setMinimumWidth(edit_width)

        # Set up the browse button: style, fix size to match
        # the neighbour lineEdit
        btn = self.__browse
        btn.setToolButtonStyle(qtc.Qt.ToolButtonTextOnly)
        btn_width = btn.fontMetrics().boundingRect(btn.text()).width() + 10
        # Modify the height to closely match the one of the lineEdit.
        # On Windows 10 we need 2 more pixels.  TODO: test other OSes
        btn_height = line.sizeHint().height() + 2
        btn.setFixedSize(btn_width, btn_height)

        # Place all in a layout.
        layout = qtw.QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(line)
        layout.addWidget(btn)
        self.setLayout(layout)

        # Update visibility and read/write property of children
        # This also updates the space between the children
        self.read_only = self.__glob['read_only']

    def __connect(self):
        """Connect appropriate signals."""
        self.__browse.clicked.connect(self.__on_browse_pressed)
        self.__lineedit.editingFinished.connect(self.__on_path_edit_finished)

    @qtc.pyqtSlot(bool)
    def __on_browse_pressed(self, _):
        """Open a modal file/directory browser, and store the selection."""
        kwargs = self.__glob['kwargs']
        if 'directory' not in kwargs:
            kwargs['directory'] = str(self.parent_directory)
        selection = self.__selector(**kwargs)
        if not selection:
            return
        if isinstance(selection, tuple):
            pathname, *_ = selection
        else:
            pathname = selection
        if not pathname:
            return
        self.path = pathname

    def __show_path(self):
        """Show selected path, possibly elided if too long."""
        parts = list(self.path.parts)
        if not parts:
            # Empty path
            self.__glob['elided_path'] = ''
            self.__lineedit.setText('')
            return

        path_left = Path(parts.pop(0))
        try:
            path_right = Path(parts.pop(-1))
        except IndexError:
            # There was only one element
            path_right = Path()

        max_ = self.__glob['max_chars']
        while parts:
            if len(str(path_left / path_right)) > max_:
                break
            *parts, right = parts
            path_right = right / path_right
            if not parts:
                break
            left, *parts = parts
            path_left /= left

        elided = ''
        if parts:
            # Still some parts missing. Add ellipsis.
            path_left /= _ELLIPSIS
            elided = str(Path(*parts))
        self.__glob['elided_path'] = elided
        self.__lineedit.setText(str(path_left / path_right))

    def __on_path_edit_finished(self):
        """Check that the new path typed makes sense."""
        _line = self.__lineedit
        _elided = self.__glob['elided_path']
        self.path = Path(_line.text().replace(_ELLIPSIS, _elided))

    # The next method is adapted from here, and the linked StackOverflow:
    # https://gist.github.com/mo-han/240b3ef008d96215e352203b88be40db
    def is_acceptable_path(self, pathname):
        """Check whether pathname is acceptable for the current filesystem."""
        try:
            pathname = Path(pathname)
        except (TypeError, ValueError):
            return False

        drive = pathname.anchor
        if not Path(drive).is_dir():
            # Invalid drive
            return False

        parts = pathname.relative_to(drive).parts
        if parts and parts[-1] in _ILLEGAL_NAMES:
            return False

        for part in parts:
            if not self.__is_acceptable_path_part(part):
                return False
        return True

    @staticmethod
    def __is_acceptable_path_part(part):
        """Return whether a path part is acceptable for this filesystem."""
        if any(c in part for c in _ILLEGAL_CHARS):
            return False

        try:
            (_ROOT / part).lstat()
        except ValueError:  # Null character
            return False
        except OSError as err:
            # Does not exist (OK) or had some fuck-up in it:
            # * Windows error 123, ERROR_INVALID_NAME;
            # * POSIX ENAMETOOLONG, ERANGE on other OSes
            if hasattr(err, 'winerror') and err.winerror == 123:
                return False
            if err.errno in (errno.ENAMETOOLONG, errno.ERANGE):
                return False
        return True
