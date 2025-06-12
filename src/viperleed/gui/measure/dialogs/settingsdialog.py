"""Module settingsdialog of viperleed.gui.measure.dialogs.

This module contains the definition of the SettingsDialog class
that proposes to the user viewing and editing common and advanced
options for devices and measurements. The Synchronization between
the ViPErLEEDSettings displayed and the widgets in the dialog is
performed via an instance of the SettingsHandler class defined here.

The widgets managed by SettingsHandler instances are grouped in
instances of SettingsDialogSectionBase, SettingsDialogSection and
SettingsDialogOption. They are viewed in insertion order.

SettingsDialogSectionBase is commonly subclassed, and used only in
case a single SettingsDialogSection containing multiple options is
not sufficient. This is typically the case when multiple entries in
the settings should be edited at once.

SettingsDialogSection instances contain options as instances of
SettingsDialogOption. All such instances are displayed as grouped.

SettingsDialogOption instances that do not belong to an explicit
section appear on their own.

Part of the code here is inspired by https://github.com/pythonguis/pyqtconfig
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-09-02'
__license__ = 'GPLv3+'

import ast
import collections
import copy
from dataclasses import dataclass
import enum
from pathlib import Path
from types import MethodType as _bound_method
import warnings

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure.widgets.fieldinfo import FieldInfo
from viperleed.gui.measure.widgets.pathselector import PathSelector
from viperleed.gui.measure.widgets.spinboxes import CoercingDoubleSpinBox
from viperleed.gui.widgets.buttons import QNoDefaultDialogButtonBox
from viperleed.gui.widgets.buttons import QNoDefaultPushButton

# TODO: find a proper mechanism to make "invalid" values disable
# "Apply" or "Ok" (but leave "Cancel" enabled). Probably equip
# widgets with some universally named, validity-checking method.

# TODO: context menu to "reset" each entry separately.

# pylint: disable=too-many-lines
# We can probably live with 1011 instead of 1000

_MSGBOX = qtw.QMessageBox


def __get_qbuttongroup(_self):
    """Return a list of (index, checked) for buttons in the group."""
    return [(i, b.isChecked()) for i, b in enumerate(_self.buttons())]


def __set_qbuttongroup(_self, values):
    """Set state of buttons from a list of (index, checked) values."""
    btns = _self.buttons()
    for idx, checked in values:
        btns[idx].setChecked(bool(checked))


def __notify_qbuttongroup(_self):
    """Return the .buttonClicked signal."""
    return _self.buttonClicked

# Hooks are: (getter, setter, notifier signal, converter from string)
# The converter is called in the setter to correctly convert the string
# values in a ViPErLEEDSettings file to values that can be managed by
# the widget. It is left to None if the setter is a callable, or if
# the widget is capable of accepting strings.
_DEFAULT_HOOKS = {
    qtw.QLabel: ('text', 'setText', None, None),
    qtw.QLineEdit: ('text', 'setText', 'textChanged', None),
    qtw.QCheckBox: ('isChecked', 'setChecked', 'stateChanged',
                    ast.literal_eval),
    CoercingDoubleSpinBox: ('cleanText', 'setValue', 'value_coerced', float),
    qtw.QSpinBox: ('cleanText', 'setValue', 'valueChanged', int),
    qtw.QDoubleSpinBox: ('cleanText', 'setValue', 'valueChanged', float),
    qtw.QButtonGroup: (__get_qbuttongroup, __set_qbuttongroup,
                       __notify_qbuttongroup, None),
    qtw.QSlider: ('value', 'setValue', 'valueChanged', int),
    qtw.QPushButton: ('isChecked', 'setChecked', 'toggled', bool),              # BUG: converter(value) == bool(value) == True whatever non-empty string value! -> perhaps literal_eval?
    qtw.QAction: ('isChecked', 'setChecked', 'toggled', bool),                  # BUG: converter(value) == bool(value) == True whatever non-empty string value! -> perhaps literal_eval?
    PathSelector: ('get_posix_path', 'set_path', 'path_changed', None),
    qtw.QTextEdit: ('toPlainText', 'setPlainText', 'textChanged', None),
    }


def _bind_default_methods(obj, **kwargs):
    """Bind default methods to an object given its type."""
    # Prefer exact type first, then accept also base classes
    hook = _get_hook(obj)
    if not hook:
        return False

    names = ('get_', 'set_', 'notify_')
    for name, method in zip(names, hook):
        flag = kwargs.get(name, False)
        if not flag:    # Don't add this method
            continue
        if not method:  # Requested, but not defined
            return False
        if isinstance(method, str):  # A method name
            bound = getattr(obj, method)
        else:                        # A callable
            bound = _bound_method(method, obj)
        setattr(obj, name, bound)

    # Fix the getter and setter in case we need to convert values
    getter, setter, _, converter = hook
    set_ = kwargs.get('set_', False)
    if set_ and isinstance(setter, str) and converter is not None:
        def __setter(_self, value):
            getattr(_self, setter)(converter(value))
        def __getter(_self):
            return str(getattr(_self, getter)())
        obj.set_ = _bound_method(__setter, obj)
        obj.get_ = _bound_method(__getter, obj)
    return True


def _get_hook(obj):
    """Return a hook for an object, inferring it from its type."""
    key = type(obj)
    try:
        hook = _DEFAULT_HOOKS[key]
    except KeyError:
        hook = next(
            (v for k, v in _DEFAULT_HOOKS.items() if isinstance(obj, k)),
            None)
    return hook


class SettingsTag(enum.Flag):
    """Flags that decide the display behaviour of a settings widget."""
    NONE = 0
    REGULAR = enum.auto()
    R = REGULAR
    ADVANCED = enum.auto()
    A = ADVANCED
    MEASUREMENT = enum.auto()
    READ_ONLY = enum.auto()


Tag = SettingsTag


@dataclass
class SettingsSectionColumnInfo:
    """Contains information on how to place a settings section.

    Attributes
    ----------
    position : int, optional
        Decides in which column of the SettingsDialog the
        SettingsSection will be displayed. Default is 0, which
        means the section will be added to the first column.
    alignment : Qt.AlignTop or Qt.AlignBottom or Qt.AlignVCenter
                or None, optional
        Influences vertical position of the section. Default is None.
        Qt.AlignTop: Moves section to the top with a bottom stretch.
        Qt.AlignBottom: Moves section to the bottom with a top stretch.
        Qt.AlignVCenter: Moves section to the middle with top
                         and bottom stretches.
        None: The section will not have any additional stretches.

    Raises
    ------
    TypeError
        If position is not an int, or if alignment is given, but
        not among the allowed values listed above.
    """

    position: int = 0
    alignment: qtc.Qt.Alignment = None

    def __post_init__(self):
        """Check that we have the correct attribute types."""
        if not isinstance(self.position, int):
            raise TypeError('Column position must be an int.')
        allowed = (qtc.Qt.AlignTop, qtc.Qt.AlignBottom, qtc.Qt.AlignVCenter)
        if self.alignment and self.alignment not in allowed:
            raise TypeError('Alignment must be an allowed alignment or None.')


class SettingsTagHandler:
    """This class can return whether it has certain tags or not."""

    def __init__(self, **kwargs):
        self._tags = kwargs.pop('tags', Tag.NONE)

    @property
    def tags(self):
        """Return tags."""
        return self._tags

    @tags.setter
    def tags(self, new_tags):
        """Set tags."""
        self._tags = new_tags or Tag.NONE

    def has_tag(self, tag):
        """Return whether the instance has this tag."""
        return bool(tag & self.tags)


class _QContainerMeta(type(collections.abc.Container), type(qtc.QObject)):
    """Meta-class for a Container and QObject."""


class SettingsHandler(collections.abc.MutableMapping, qtc.QObject,
                      metaclass=_QContainerMeta):
    """Class that takes care of syncing handler widgets and settings."""

    settings_changed = qtc.pyqtSignal()
    redraw_needed = qtc.pyqtSignal()  # Should trigger redraw of dialog
    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, config, parent=None, show_path_to_config=False):
        """Initialize instance.

        Parameters
        ----------
        config : ViPErLEEDSettings
            The settings from which the handler will be generated.
        parent : QObject
            The parent QObject of this handler.
        show_path_to_config : bool, optional
            Whether the name of the settings file (and path to the file)
            should be displayed. Default is False.

        Returns
        -------
        None.
        """
        super(qtc.QObject, self).__init__(parent)
        self.__dict = collections.defaultdict(dict)
        self.__sub_handlers = []
        self.__sections = {}
        self.__complex_sections = []
        self.__config = config
        self.__widgets = []

        # Use a timer to delay a little the emission of .updated.
        # This way we should trigger a redraw only once if multiple
        # complex sections update within a few ms time span
        self.__updated_timer = qtc.QTimer(self)
        self.__updated_timer.setInterval(10)
        self.__updated_timer.setSingleShot(True)
        self.__updated_timer.timeout.connect(self.redraw_needed)

        if show_path_to_config:
            widget = qtw.QLabel()
            file = config.last_file
            widget.setText(file.stem if file else 'None')
            self.add_static_option(
                'File', 'config', widget,
                display_name='Settings file',
                tooltip=str(file) if file else '',
                tags=Tag.REGULAR,
                )

    def __delitem__(self, item):
        """Delete a section."""
        del self.__dict[item]
        self.__widgets.remove(self.__sections.pop(item))

    def __getitem__(self, item):
        """Return a section."""
        return self.__dict[item]

    def __iter__(self):
        """Return an iterable version of self."""
        return iter(self.__dict)

    def __len__(self):
        """Return the number of sections in self."""
        return len(self.__dict)

    def __setitem__(self, key, value):
        """Set a section of self."""
        self.__dict[key] = value

    @staticmethod
    def handler_from_type(type_):
        """Return a suitable handler widget given the type of entry."""
        if type_ == str:
            return qtw.QLineEdit
        if type_ == float:
            return qtw.QDoubleSpinBox
        if type_ == int:
            return qtw.QSpinBox
        if type_ == bool:
            return qtw.QCheckBox
        return None

    @property
    def settings(self):
        """Return the ViPErLEEDSettings managed."""
        return self.__config

    @property
    def widgets(self):
        """Return handler widgets."""
        return tuple(self.__widgets)

    def add_from_handler(self, handler):
        """Add section/options from another handler."""
        if not isinstance(handler, SettingsHandler):
            raise TypeError("Not a SettingsHandler")
        if self.settings is not handler.settings:
            raise ValueError("Cannot add a handler for different settings")

        # Keep a reference, otherwise we may risk loosing it, as well
        # as all the internal connections to the managed settings
        self.__sub_handlers.append(handler)

        # Then fill up self. The order is the one of .widgets,
        # but we have to fill in other parts as well
        for widget in handler.widgets:
            if isinstance(widget, SettingsDialogSection):
                self.__sections[widget.section_name] = widget
                self.__widgets.append(widget)
                for option in widget.options:
                    self[widget.section_name][option.option_name] = option
            else:
                self.__widgets.append(widget)
        handler.settings_changed.connect(self.settings_changed)

    def add_option(self, section_name, option_name, *args,
                   handler_widget=None, option_type=None, **kwargs):
        """Add a handler for a section/option pair."""
        if not self.__config.has_option(section_name, option_name):
            raise ValueError(
                "Configuration file does not contain a section"
                f"/option pair {section_name}/{option_name}"
                )
        if section_name not in self:
            warnings.warn(
                f"{self.__class__.__name__}: section {section_name} "
                f"was not added with add_section(). Option {option_name}"
                " will appear without a bounding frame."
                )

        if handler_widget is None and option_type is not None:
            handler_widget = self.handler_from_type(option_type)
        if handler_widget is None:
            value = self.__config[section_name][option_name]
            handler_widget = self.__guess_handler_from_value(value)
        if handler_widget is None:
            raise TypeError(
                "Could not infer a default handler widget for section/"
                f"option pair {section_name}/{option_name} with type "
                f"{option_type}. Pass an appropriate widget explicitly."
                )
        option = SettingsDialogOption(option_name, handler_widget,
                                      *args, **kwargs)
        self[section_name][option_name] = option

        if not option.has_tag(Tag.READ_ONLY):
            option.value_changed.connect(
                self.__option_setter(section_name, option_name)
                )

        if self.has_section(section_name):
            self.__sections[section_name].add_option(option)
        else:
            self.__widgets.append(option)

    def add_complex_section(self, section):
        """Add a section that internally handles writing to settings.

        Parameters
        ----------
        section : SettingsDialogSectionBase
            The section to be added. Typically a subclass of
            SettingsDialogSectionBase that implements handling
            the settings internally. This can be used to
            access and modify multiple settings in a complex
            manner.

        Raises
        ------
        TypeError
            If section is not a subclass of SettingsDialogSectionBase
        """
        if not isinstance(section, SettingsDialogSectionBase):
            raise TypeError(
                "add_complex_section requires a SettingsDialogSectionBase"
                )
        section.settings_changed.connect(self.settings_changed)
        self.__widgets.append(section)
        self.__complex_sections.append(section)
        try:
            section.error_occurred.connect(self.error_occurred)                 # TODO: should regular options and sections also have an error_occurred signal?
        except AttributeError:
            pass
        section.updated.connect(self.__updated_timer.start)

    def add_static_option(self, section_name, option_name, handler_widget,
                          *args, **kwargs):
        """Add a static option that may not come from settings."""
        if section_name not in self:
            warnings.warn(
                f'{self.__class__.__name__}: section {section_name} '
                f'was not added with add_section(). Option {option_name}'
                ' will appear without a bounding frame.'
                )
        # pylint: disable-next=unsupported-binary-operation
        tags = kwargs.pop('tags', Tag.NONE) | Tag.READ_ONLY
        option = StaticSettingsDialogOption(option_name, handler_widget,
                                            *args, tags=tags, **kwargs)
        self[section_name][option_name] = option
        if self.has_section(section_name):
            self.__sections[section_name].add_option(option)
        else:
            self.__widgets.append(option)

    def add_section(self, section_name, **kwargs):
        """Add a titled section to self."""
        if section_name not in self.__config:
            raise ValueError(f"No section {section_name} in config file")
        self._add_section(section_name, **kwargs)

    def add_static_section(self, section_name):
        """Add a static section that may not come from settings."""
        self._add_section(section_name)

    def get_widgets_with_tags(self, tags):
        """Return all widgets with a specific tag."""
        return tuple(wid for wid in self.widgets if tags in wid.tags)

    def has_advanced_options(self):
        """Return whether self contains any advanced option."""
        _adv = any(o.has_tag(Tag.A) for s in self.values() for o in s.values())
        return _adv or any(s.has_tag(Tag.A) for s in self.__complex_sections)

    def has_section(self, section):
        """Return whether a section is directly handled."""
        return section in self.__sections

    def make_from_config(self):
        """Create appropriate handlers for entries in self.config."""
        for section_name in self.__config.sections():
            if not self.__config[section_name]:
                # Empty section
                continue
            self.add_section(section_name, tags=Tag.REGULAR)
            for option_name, value in self.__config[section_name].items():
                handler = self.__guess_handler_from_value(value)
                if handler is None:
                    handler = qtw.QLineEdit()
                self.add_option(section_name, option_name,
                                handler_widget=handler)

    def update_widgets(self):
        """Update values in handlers from the config."""
        for sec_name, options in self.items():
            for opt_name, option in options.items():
                if isinstance(option, StaticSettingsDialogOption):
                    continue
                option.set_(self.__config[sec_name][opt_name])
        for section in self.__complex_sections:
            section.update_widgets()

    def _add_section(self, section_name, **kwargs):
        """Add section to self."""
        section = SettingsDialogSection(section_name, **kwargs)
        self.__sections[section_name] = section
        self.__widgets.append(section)

    def __guess_handler_from_value(self, value):
        """Return a handler widget guessed from a config value."""
        if not value:   # We cannot guess from empty values
            return None

        as_path = Path(value)
        if as_path.is_file(): # A file?
            return PathSelector()
        if as_path.exists():  # A directory?
            return PathSelector(select_file=False)

        # See if the value can be interpreted
        try:
            value = ast.literal_eval(value)
        except (TypeError, ValueError, SyntaxError,
                MemoryError, RecursionError):
            return None

        return self.handler_from_type(type(value))

    def __option_setter(self, section, option):
        """Return a pyqtSlot setter for a given section/option of config."""
        def _setter(new_value):
            old_value = self.__config[section][option]
            if old_value != new_value:
                self.__config[section][option] = new_value
                self.settings_changed.emit()
        return qtc.pyqtSlot(str)(_setter)


class SettingsDialogOption(qtc.QObject, SettingsTagHandler):
    """Class for handling a single settings option."""

    value_changed = qtc.pyqtSignal(str)
    handler_widget_changed = qtc.pyqtSignal()

    def __init__(self, option_name, handler_widget, *args, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        option_name : str
            The name of this option in the settings file.
        handler_widget : QWidget or type(QWidget)
            The widget instance or widget class to be used to
            display the option. If a class, an instance is created.
        *args : object
            Other positional arguments  passed on to handler_widget
            if only a QWidget class is given.
        display_name : str, optional
            The name of this section when displayed in a SettingsDialog.
            If not given or False-y, the name is taken from option_name:
            underscores are removed, and all words are capitalized.
        tooltip : str, optional
            A descriptive text that will be used as tooltip, displayed
            when hovering over, or clicking on the info icon. If it
            is an empty string, no tooltip is shown. Default is an
            empty string.
        tags : SettingsTag, optional
            Contains additional tags of this option. Possible tags are
            a bitwise-or of:
                SettingsTag.REGULAR
                    Present if the option contains regular settings.
                SettingsTag.ADVANCED
                    Present if the option contains advanced settings.
                SettingsTag.MEASUREMENT
                    Present if the option contains settings
                    related to the measurement.
                SettingsTag.READ_ONLY
                    Present if the option contains read-only settings.
            Default is None.
        label_alignment : {'top', 'centre', 'bottom'}
            Vertical alignment of label field relative to handler_widget.
            Only the first character matters. Any character other than 'c'
            or 'b' will be treated as 'top'. Default is 'top'.
        **kwargs : object
            Other keyword arguments passed on to handler_widget if
            only a QWidget class is given.

        Returns
        -------
        None.
        """
        display_name = kwargs.pop('display_name', None)
        tooltip = kwargs.pop('tooltip', '')
        v_align = kwargs.pop('label_alignment', 't')

        super().__init__(*args, **kwargs)
        self.option_name = option_name

        if isinstance(handler_widget, type(qtw.QWidget)):
            handler_widget = handler_widget(*args, **kwargs)

        self._handler_widget = handler_widget
        self._label = qtw.QLabel()
        self._info = None
        self._check_handler()
        self._connect_handler()
        self._update_handler_from_read_only()

        self.display_name = self._make_label_widget(display_name, tooltip,
                                                    v_align)
        if not self.has_tag(Tag.REGULAR):
            self.setVisible(False)

    def __iter__(self):
        """Return the label and the handler for this option."""
        return iter((self.display_name, self.handler_widget))

    def __repr__(self):
        """Return a string representation of self."""
        return f"SettingsDialogOption({self.option_name})"

    @property
    def handler_widget(self):
        """Return the handler widget for this option."""
        return self._handler_widget

    @handler_widget.setter
    def handler_widget(self, new_widget):
        """Set a new handler widget."""
        if new_widget is self.handler_widget:
            return
        self._handler_widget = new_widget
        self._check_handler()
        self._connect_handler()
        self._update_handler_from_read_only()
        self.handler_widget_changed.emit()

    @property
    def label(self):
        """Return the QLabel widget for this option's label."""
        return self._label

    def get_(self):
        """Return the value of this option as a string."""
        return self.handler_widget.get_()

    def set_(self, value):
        """Set value displayed by this handler."""
        self.handler_widget.set_(value)

    def set_enabled(self, enabled):
        """Enable or disable this option."""
        for child in self:
            child.setEnabled(enabled)

    def set_info_text(self, text):
        """Set informative text."""
        self._info.set_info_text(text)
        self._info.setVisible(bool(text))

    def setVisible(self, visible):   # pylint: disable=invalid-name
        """Set the visibility of all children."""
        for child in self:
            child.setVisible(visible)

    def _check_handler(self):
        """Check that the handler widget can be used."""
        to_have = ('get_', 'set_')
        if not self.has_tag(Tag.READ_ONLY):
            to_have += ('notify_',)
        handler = self.handler_widget
        missing = {m: True for m in to_have if not hasattr(handler, m)}
        if any(missing) and not _bind_default_methods(handler, **missing):
            raise ValueError(
                f"Invalid handler_widget {handler.__class__.__name__} "
                f"has no {'/'.join(missing)} method(s)"
                )

    def _connect_handler(self):
        """Connect self.handler_widget if not read only."""
        if self.has_tag(Tag.READ_ONLY):
            return

        signal = self.handler_widget.notify_
        if not isinstance(signal, qtc.pyqtBoundSignal):
            # Probably a method returning the signal itself
            signal = signal()
        signal.connect(self._notify_change)

    def _make_label_widget(self, label_text, info_text, v_align):
        """Return a QWidget to act as option label."""
        # Get a reasonable label_text
        if not label_text:
            label_text = self.option_name.replace('_', ' ').title()
        label_text = label_text.strip()
        if not label_text.endswith(':'):
            label_text += ':'
        self.label.setText(label_text)

        # Prepare a container widget and its layout
        container = qtw.QWidget()
        container.setLayout(qtw.QVBoxLayout())
        v_align_layout = container.layout()
        h_align_layout = qtw.QHBoxLayout()
        v_align_layout.setContentsMargins(0, 0, 0, 0)
        h_align_layout.setContentsMargins(0, 0, 0, 0)

        self._info = FieldInfo.for_widget(self.label, tooltip=info_text)

        # Fill layout
        h_align_layout.addWidget(self.label)
        h_align_layout.addWidget(self._info)
        v_align_layout.addLayout(h_align_layout)

        # Sort out vertical alignment, using stretches. "Top"
        # is the default for QFormLayout, i.e., nothing to do
        if v_align.startswith('c'):
            v_align_layout.insertStretch(0, 1)
            v_align_layout.insertStretch(-1, 1)
        elif v_align.startswith('b'):
            # Mimic the Qt implementation for top alignment,
            # i.e., we leave a little space at the bottom
            v_align_layout.insertStretch(0, 7)
            v_align_layout.insertStretch(-1, 1)

        # Now horizontal alignment: Decide where to
        # place a stretch to keep text & info together
        is_left_align = qtw.QFormLayout().labelAlignment() == qtc.Qt.AlignLeft
        h_align_layout.insertStretch(-1 if is_left_align else 0, 1)

        self.set_info_text(info_text)
        return container

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(object)
    def _notify_change(self, _=None):
        """Emit a value_changed signal with the new value of this option."""
        self.value_changed.emit(self.get_())

    def _update_handler_from_read_only(self):
        """Update the state of handler_widget according to READ_ONLY tag."""
        handler = self.handler_widget

        # See if there is a readOnly method
        try:
            handler.setReadOnly(self.has_tag(Tag.READ_ONLY))
        except AttributeError:
            pass
        else:
            return

        # See if there is a writeable read_only property
        try:
            handler.read_only = self.has_tag(Tag.READ_ONLY)
        except (AttributeError, TypeError):
            # Just disable it instead
            handler.setEnabled(not self.has_tag(Tag.READ_ONLY))


class SettingsDialogSectionBase(qtw.QGroupBox, SettingsTagHandler):
    """A base class for handling groups of settings.

    Use this class only as the parent class for "advanced"
    sections, i.e., those that handle themselves all the
    sync with the configuration file, and which may want
    to have a specialized way of showing the options they
    handle. You can do that by preparing a layout and set
    it as the layout of self.central_widget.
    """

    # The next signal can be used to notify when any of the
    # editable settings are changed. Emit this signal in
    # "advanced" sections subclasses, as it is used by the
    # SettingsHandler to notify of changes.
    settings_changed = qtc.pyqtSignal()

    # This signal can be used to trigger a check whether the
    # settings in a SettingsDialog are ok to be accepted. For
    # this to work, the are_settings_ok() method must be
    # reimplemented to check whether the settings of the
    # section are ok.
    settings_ok_changed = qtc.pyqtSignal()

    # The next signal is emitted when this section undergoes
    # an automatic update of its widgets that should trigger
    # a redraw of the dialog
    updated = qtc.pyqtSignal()

    def __init__(self, *__args, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        display_name : str, keyword only
            The name of this section when displayed in a SettingsDialog.
        tooltip : str, optional
            A descriptive text that will be used as tooltip, displayed
            when the mouse cursor hovers over the section title. If an
            empty string no tooltip is shown. Default is an empty string.
        column_info : SettingsSectionColumnInfo, optional
            Information on how to display the section. Check out the
            SettingsSectionColumnInfo class for more information. If
            not given, the section will be displayed in the main column.
        parent : QWidget, optional
            The parent widget of this SettingsDialogSection. Default
            is None.
        tags : SettingsTag, optional
            Contains additional tags of this section. Possible tags are
            a bitwise-or of:
                SettingsTag.REGULAR
                    Present if the section contains regular settings.
                SettingsTag.ADVANCED
                    Present if the section contains advanced settings.
                SettingsTag.MEASUREMENT
                    Present if the section contains settings
                    related to the measurement.
                SettingsTag.READ_ONLY
                    Present if the section contains read-only settings.
            Default is None.

        Raises
        -------
        TypeError
            If no display_name is given.
        TypeError
            If column_info is given, but not an instance
            of SettingsSectionColumnInfo.
        """
        display_name = kwargs.get('display_name', '')
        if not display_name:
            raise TypeError("Missing display_name")

        self.column_info = kwargs.get('column_info',
                                      SettingsSectionColumnInfo())
        if not isinstance(self.column_info, SettingsSectionColumnInfo):
            raise TypeError('Column info must be a SettingsSectionColumnInfo.')

        tooltip = kwargs.get('tooltip', '')
        self._info = qtw.QLabel()
        self.central_widget = qtw.QWidget()
        super().__init__(display_name, **kwargs)

        self.__compose()
        self.set_info(tooltip)

    def __repr__(self):
        """Return a string representation of self."""
        return f"{self.__class__.__name__}(display_name='{self.title()}')"

    def are_settings_ok(self):
        """Return whether the section settings are acceptable.

        Returns
        -------
        settings_ok : bool
            Whether the settings selected in the widget are
            acceptable or not.
        reason : str
            A descriptive string elaborating why the settings
            are not acceptable.
        """
        return True, ''

    def set_info(self, info_text):
        """Add informative text in a QLabel."""
        if not info_text:
            self._info.hide()
            return
        self._info.setText(info_text)
        self._info.show()

    def __compose(self):
        """Place children widgets."""
        _policy = self.sizePolicy()

        self._info.setWordWrap(True)
        # The next line suppresses annoying QWarnings about
        # QWindowsWindow being unable to ::setGeometry
        self._info.setSizePolicy(_policy.Preferred, _policy.Minimum)

        # Use a box layout to add some space left and right so the
        # text is a bit narrower than other option widgets below
        info_layout = qtw.QHBoxLayout()
        info_layout.setContentsMargins(0, 0, 0, 0)
        info_layout.addSpacing(10)
        info_layout.addWidget(self._info)
        info_layout.setSizeConstraint(info_layout.SetMinimumSize)
        info_layout.addSpacing(10)

        layout = qtw.QVBoxLayout()
        layout.addLayout(info_layout)

        # Other customization
        self.setFlat(True)     # Remove left, right, bottom frame edges
        self.central_widget.setSizePolicy(_policy.Minimum, _policy.Minimum)

        layout.addWidget(self.central_widget)

        self.setSizePolicy(_policy.Minimum, _policy.Minimum)
        self.setLayout(layout)


class SettingsDialogSection(SettingsDialogSectionBase):
    """A class for handling groups of settings."""

    def __init__(self, section_name, *__args, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        section_name : str
            The name of this section in the settings file
        display_name : str, optional
            The name of this section when displayed in a SettingsDialog.
            If not given or False-y, the name is take from section_name:
            underscores are removed, and all words are capitalized.
        tooltip : str, optional
            A descriptive text that will be used as tooltip, displayed
            when the mouse cursor hovers over the section title. If an
            empty string no tooltip is shown. Default is an empty string.
        tags : SettingsTag, optional
            Contains additional tags of this section. Possible tags are
            a bitwise-or of:
                SettingsTag.REGULAR
                    Present if the section contains regular settings.
                SettingsTag.ADVANCED
                    Present if the section contains advanced settings.
                SettingsTag.MEASUREMENT
                    Present if the section contains settings
                    related to the measurement.
                SettingsTag.READ_ONLY
                    Present if the section contains read-only settings.
            Default is None.
        options : Sequence, optional
            Elements are SettingsDialogOption instances. Options can
            also be added later with the .add_option/.add_options
            method. Default is no options.
        parent : QWidget, optional
            The parent widget of this SettingsDialogSection. Default
            is None.

        Returns
        -------
        None.
        """
        self.section_name = section_name

        if not kwargs.get('display_name', ''):
            kwargs['display_name'] = section_name.replace('_', ' ').title()

        super().__init__(**kwargs)
        self.__options = []

        self.__layout = qtw.QFormLayout()       # For inserting options
        self.__layout.setContentsMargins(0, 0, 0, 0)
        self.central_widget.setLayout(self.__layout)

        # Fill up widget
        self.add_options(kwargs.get('options', ()))

    def __repr__(self):
        """Return a string representation of self."""
        return f"{self.__class__.__name__}({self.section_name})"

    @property
    def options(self):
        """Return a tuple of options."""
        return tuple(self.__options)

    def add_option(self, option):
        """Add one SettingsDialogOption to self."""
        self.__options.append(option)
        self.__layout.addRow(*option)
        option.tags |= self.tags
        if self.has_tag(Tag.REGULAR):
            option.setVisible(True)
        option.handler_widget_changed.connect(self.__on_option_widget_changed)

    def add_options(self, options):
        """Add more SettingsDialogOptions to this section from a Sequence."""
        for option in options:
            self.add_option(option)

    def __on_option_widget_changed(self):
        """Replace the correct widget in the layout."""
        option = self.sender()
        row_idx = self.__options.index(option)

        # Remove the old widget and delete it
        old_row = self.__layout.takeRow(row_idx)
        old_row.fieldItem.widget().deleteLater()

        # Add back the new widget
        self.__layout.insertRow(row_idx, *option)


class SettingsDialog(qtw.QDialog):
    """A dialog to display settings."""

    # This signal is emitted whenever the OK, Cancel, or Apply
    # buttons are pressed, and only in case settings changed since
    # the last time this signal was emitted (or since the dialog
    # was shown). Users can .connect to this signal and use the
    # slot to update the object whose settings are being edited.
    settings_changed = qtc.pyqtSignal()

    # This signal is emitted every time the dialog finishes and if
    # any change occurred to the settings. It carries True if edited
    # settings were saved to file, False otherwise. Users can connect
    # to his signal, and may want to restore the original settings if
    # False.
    settings_saved = qtc.pyqtSignal(bool)

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, handled_obj=None, settings=None, title=None, **kwargs):
        """Initialize dialog instance.

        Parameters
        ----------
        handled_obj : object, optional
            The object whose settings are handled by this dialog.
            It should have a .settings attribute returning an instance
            of a ViPErLEEDSettings. If no handled_obj is given,
            settings should be given instead.
        settings : ViPErLEEDSettings, optional
            A ViPErLEEDSettings to be displayed in full. This
            argument is ignored if a handled_obj is passed.
        **kwargs : object
            Other keyword arguments passed on to QDialog

        Raises
        -------
        TypeError
            If neither a handled_obj nor a settings are given, or
            if a handled_obj is given, but it has no .settings
        """
        super().__init__(**kwargs)
        if handled_obj and not hasattr(handled_obj, 'settings'):
            raise TypeError(
                f"{self.__class__.__name__}: handled object has no .settings"
                )
        if not handled_obj and not settings:
            raise TypeError(
                f"{self.__class__.__name__}: need either a "
                "handled_obj or a settings to be displayed."
                )

        self._handled_obj = handled_obj
        if handled_obj:
            settings = handled_obj.settings
            handler = handled_obj.get_settings_handler()
        else:
            handler = SettingsHandler(settings)
            handler.make_from_config()
        self.handler = handler

        self._settings = {'current': settings,
                          'applied': copy.deepcopy(settings),
                          'original': copy.deepcopy(settings)}

        # Set up children widgets and self
        self._ctrls = {
            'accept': None, # Will be an 'Accept' button
            'apply': None,  # Will be an 'Apply' button
            'advanced': QNoDefaultPushButton("Show less"),
            }
        self._compose_and_connect()

        # And finally the window properties
        self.update_title(title, settings)
        self.setWindowFlags(self.windowFlags()        # Remove '?'
                            & ~qtc.Qt.WindowContextHelpButtonHint)

    @property
    def adv_button(self):
        """Return the advanced settings button."""
        return self._ctrls['advanced']

    @property
    def settings(self):
        """Return the settings currently displayed."""
        return self._settings['current']

    @property
    def handled_object(self):
        """Return the object whose settings are shown."""
        return self._handled_obj

    @qtc.pyqtSlot()
    def accept(self):
        """Notify if settings changed, decide whether to save, then close."""
        self.__on_apply_pressed()
        # Ask to save the settings to file.
        if self.settings != self._settings['original']:
            action = self._save_edited_settings(self._ask_to_save())
            if action == self.Rejected:
                super().reject()
                return
        super().accept()

    @qtc.pyqtSlot()
    def reject(self):
        """Load back the original settings, then close."""
        if not self.isVisible():
            super().reject()
            return

        # Ask confirmation if settings changed
        _changed = self._ctrls['apply'].isEnabled()
        _changed |= self._settings['original'] != self._settings['applied']
        if _changed:
            reply = _MSGBOX.question(
                self, "Discard changes?",
                f"{self.windowTitle()} were edited.\nAre you"
                " sure you want to discard changes?",
                _MSGBOX.Discard | _MSGBOX.Cancel
                )
            if reply == _MSGBOX.Cancel:
                return

        self.settings.read_dict(self._settings['original'])
        self.__on_apply_pressed()
        if self._settings['original'] != self._settings['applied']:
            self.settings_saved.emit(False)
        super().reject()

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Update widgets content from settings when shown."""
        if not event.spontaneous():
            # i.e., not a show after minimized
            self.settings.read_again()
            self.handler.update_widgets()
            self.adv_button.setChecked(False)
            if self.handled_object:
                self.update_title()
            # Update all settings with the current ones, and
            # fix the enabled state of the "Apply" button
            for key in ('applied', 'original'):
                self._settings[key] = copy.deepcopy(self.settings)
            self.__update_apply_enabled()
        super().showEvent(event)

    def update_title(self, title='', settings=None):
        """Update title from handled_obj, if given."""
        if not title:
            handled_obj = self._handled_obj
            if hasattr(handled_obj, "name"):
                title = handled_obj.name
            elif hasattr(handled_obj, "display_name"):
                title = handled_obj.display_name
            elif settings:
                title = settings.__class__.__name__
            else:
                title = handled_obj.__class__.__name__
            title += " settings"
        self.setWindowTitle(title)

    def _ask_to_save(self):
        """Ask the users whether to save edits or not.

        Returns
        -------
        reply : QMessageBox.Constant
            The action selected by the user.
        """
        message_box = self._get_ask_to_save_dialog()
        return message_box.exec()

    def _compose_and_connect(self):
        """Place and update children widgets."""
        self.handler.error_occurred.connect(self.error_occurred)
        self.handler.update_widgets()  # Fill widgets from settings

        buttons = self._compose_dialog_buttons()
        columns = self._compose_columns()

        columns_layout = qtw.QHBoxLayout()
        for column in columns:
            columns_layout.addLayout(column)

        outer_layout = qtw.QVBoxLayout()
        outer_layout.addLayout(columns_layout)
        outer_layout.addWidget(buttons)
        outer_layout.setSizeConstraint(outer_layout.SetMinimumSize)
        self.setLayout(outer_layout)
        outer_layout.setSpacing(round(outer_layout.spacing()*1.4))

        _policy = self.sizePolicy()
        self.setSizePolicy(_policy.Minimum, _policy.Minimum)

        # Connect signals
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self._ctrls['apply'].clicked.connect(self.__on_apply_pressed)
        self.adv_button.toggled.connect(self.__on_show_advanced_toggled)
        self.handler.settings_changed.connect(self.__update_apply_enabled)
        self.handler.redraw_needed.connect(self.__update_advanced_btn)

        self.__update_advanced_btn()

    def _compose_dialog_buttons(self):
        """Compose the buttons of the dialog.

        Returns
        -------
        buttons : QNoDefaultDialogButtonBox
            Contains the buttons necessary to handle the dialog.
        """
        _bbox = QNoDefaultDialogButtonBox
        buttons = _bbox(_bbox.Ok | _bbox.Cancel | _bbox.Apply)

        self._ctrls['apply'] = buttons.buttons()[-1]
        self._ctrls['accept'] = buttons.buttons()[0]
        self._ctrls['apply'].setEnabled(False)
        self.adv_button.setCheckable(True)

        # Use a ResetRole to have the button placed in a
        # different spot than all others on every platform
        buttons.addButton(self.adv_button, _bbox.ResetRole)
        return buttons

    def _compose_columns(self):
        """Compose the columns of the layout and add widgets.

        Returns
        -------
        columns : list
            Contains columns of the dialog layout.
        """
        columns = [qtw.QVBoxLayout(),]
        for widg in self.handler.widgets:
            if widg.has_tag(Tag.MEASUREMENT) and not widg.has_tag(Tag.REGULAR):
                continue
            if isinstance(widg, SettingsDialogSectionBase):
                position = widg.column_info.position
                alignment = widg.column_info.alignment
                while position >= len(columns):                                 # TODO: We could make a defaultlist class similar to defaultdict.
                    columns.append(qtw.QVBoxLayout())
                if alignment in (qtc.Qt.AlignVCenter, qtc.Qt.AlignBottom):
                    columns[position].addStretch(1)
                columns[position].addWidget(widg)
                if alignment in (qtc.Qt.AlignVCenter, qtc.Qt.AlignTop):
                    columns[position].addStretch(1)
                continue
            # An option. Need to add it in a QFormLayout
            _form = qtw.QFormLayout()
            _form.addRow(*widg)
            columns[0].addLayout(_form)
        columns[0].addStretch(1)
        return columns

    def _get_ask_to_save_dialog(self):
        """Get the QMessageBox that asks if settings can be saved."""
        message_box = _MSGBOX(
            _MSGBOX.Question, "Save settings to file?",
            (f"{self.windowTitle()} were edited.\n\n"
             "Would you like to save changes to file?"),
            parent=self,
            )
        message_box.addButton(_MSGBOX.Save)
        message_box.addButton(_MSGBOX.Discard)
        return message_box

    @qtc.pyqtSlot(bool)
    def __on_apply_pressed(self, _=False):
        """React to the user pressing 'Apply'."""
        self._ctrls['apply'].setEnabled(False)
        if self.settings != self._settings['applied']:
            self.settings_changed.emit()
        self._settings['applied'].read_dict(self.settings)

    @qtc.pyqtSlot(bool)
    def __on_show_advanced_toggled(self, visible):
        """Show or hide advanced options."""
        if not self.handler.has_advanced_options():
            return
        if visible:
            btn_text = "Show less"
        else:
            btn_text = "Show all"
        self.adv_button.setText(btn_text)
        for widg in self.handler.widgets:
            if not widg.has_tag(Tag.A | Tag.R):
                continue
            if isinstance(widg, SettingsDialogSection):
                _section_visible = False
                for option in widg.options:
                    option_visible = visible or (not option.has_tag(Tag.A)
                                                 and option.has_tag(Tag.R))
                    _section_visible |= option_visible
                    option.setVisible(option_visible)
                widg.setVisible(_section_visible)
            else:
                widg.setVisible(visible or (not widg.has_tag(Tag.A)
                                            and widg.has_tag(Tag.R)))
        self.adjustSize()   # TODO: does not always adjust when going smaller?

    def _save_edited_settings(self, reply):
        """Save changes to the current settings to file.

        Parameters
        ----------
        reply : qtw.QMessageBox.Constant
            The response seleced by the user. Decides
            whether the settings should be saved or not.

        Returns
        -------
        action : QDialog.Accepted
            The action that is to be performed on the dialog.
        """
        if reply == _MSGBOX.Save:
            try:
                self.settings.update_file()
            except FileNotFoundError:
                # The file must have been moved before
                # the settings could be saved.                                  # TODO: open QMessageBox to ask the user how to proceed from here. Ask whether the user wants to save the settings and if yes, ask where and under which name.
                pass
            finally:
                self._settings['original'].read_dict(self.settings)
        self.settings_saved.emit(reply == _MSGBOX.Save)
        return self.Accepted

    @qtc.pyqtSlot()
    def __update_advanced_btn(self):
        """Update visibility of button and options."""
        self.adv_button.setVisible(self.handler.has_advanced_options())
        self.__on_show_advanced_toggled(self.adv_button.isChecked())

    def __update_apply_enabled(self):
        """Enable/disable 'Apply' depending on whether anything changed."""
        settings_changed = self.settings != self._settings['applied']
        self._ctrls['apply'].setEnabled(settings_changed)


class MeasurementSettingsDialog(SettingsDialog):
    """A dialog to display measurement settings."""

    @qtc.pyqtSlot()
    def accept(self):
        """Store device settings, then notify if settings changed."""
        for widget in self.handler.widgets:
            try:
                widget.store_lower_level_settings()
            except AttributeError:
                pass
        super().accept()

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Check if the updated widgets have faulty settings."""
        super().showEvent(event)
        self._check_if_settings_ok()

    def _check_if_settings_ok(self):
        """Check if settings are ok and enable/disable accept button."""
        settings_ok = True
        self._ctrls['accept'].setToolTip('')
        for widget in self.handler.widgets:
            try:
                settings_ok, reason = widget.are_settings_ok()
            except AttributeError:
                pass
            if not settings_ok:
                self._ctrls['accept'].setToolTip(reason)
                break
        self._ctrls['accept'].setEnabled(settings_ok)

    def _compose_and_connect(self):
        """Place and update children widgets."""
        super()._compose_and_connect()
        self._ctrls['accept'].setText('Start measurement')
        self._ctrls['apply'].setVisible(False)
        for widget in self.handler.widgets:
            try:
                widget.settings_ok_changed.connect(self._check_if_settings_ok)
            except AttributeError:
                pass

    def _get_ask_to_save_dialog(self):
        """Ask whether to save or not.

        Returns
        -------
        reply : QMessageBox.Constant
            The action selected by the user.
        """
        message_box = super()._get_ask_to_save_dialog()
        for button in message_box.buttons():
            if message_box.buttonRole(button) == _MSGBOX.DestructiveRole:
                button.setText('Discard settings')
                break
        message_box.addButton(_MSGBOX.Abort)
        return message_box

    def _save_edited_settings(self, reply):
        """Save changes to the current settings to file.

        Parameters
        ----------
        reply : qtw.QMessageBox.Constant
            The response seleced by the user. Decides
            whether the settings should be saved or not.

        Returns
        -------
        action : QDialog.DialogCode
            The action that is to be performed on the dialog.
        """
        action = super()._save_edited_settings(reply)
        if reply == _MSGBOX.Abort:
            action = self.Rejected
        return action


class StaticSettingsDialogOption(SettingsDialogOption):
    """Read only settings dialog option."""

    def __init__(self, option_name, handler_widget, *args, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        option_name : str
            The name of this option in the settings file.
        handler_widget : QWidget or type(QWidget)
            The widget instance or widget class to be used to
            display the option. If a class, an instance is created.
        *args : object
            Other positional arguments passed on to `handler_widget`
            if only a QWidget class is given.
        display_name : str, optional
            The name of this section when displayed in a SettingsDialog.
            If not given or False-y, the name is taken from option_name:
            underscores are removed, and all words are capitalized.
        tooltip : str, optional
            A descriptive text that will be used as tooltip, displayed
            when hovering over, or clicking on the info icon. If it
            is an empty string, no tooltip is shown. Default is an
            empty string.
        tags : SettingsTag, optional
            Contains additional tags of this option. Possible tags are
            a bitwise-or of:
                SettingsTag.REGULAR
                    Present if the option contains regular settings.
                SettingsTag.ADVANCED
                    Present if the option contains advanced settings.
                SettingsTag.MEASUREMENT
                    Present if the option contains settings
                    related to the measurement.
                SettingsTag.READ_ONLY
                    Option contains read-only settings.
            SettingsTag.READ_ONLY must be present, otherwise a
            ValueError is raised.
        label_alignment : {'top', 'centre', 'bottom'}
            Vertical alignment of label field relative to handler_widget.
            Only the first character matters. Any character other than 'c'
            or 'b' will be treated as 'top'. Default is 'top'.
        **kwargs : object
            Other keyword arguments passed on to handler_widget if
            only a QWidget class is given.

        Returns
        -------
        None.

        Raises
        ------
        ValueError
            read_only is not True. Indicates wrong
            use of StaticSettingsDialogOption.
        """
        tags = kwargs.get('tags', Tag.NONE)
        if Tag.READ_ONLY not in tags:
            raise ValueError(
                'A StaticSettingsDialogOption may only be read only. You '
                'are seeing this message due to a faulty implementation.'
                )
        super().__init__(option_name, handler_widget, *args, **kwargs)

    def _check_handler(self):
        """Disable handler check."""
