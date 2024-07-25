"""Module settingsdialog of viperleed.guilib.measure.dialogs.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2022-09-02
Author: Michele Riva

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

import ast
import copy
import collections
from pathlib import Path
from types import MethodType as _bound_method
import warnings

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure.widgets.pathselector import PathSelector
from viperleed.guilib.measure.widgets.fieldinfo import FieldInfo
from viperleed.guilib.widgets.basewidgets import QNoDefaultDialogButtonBox

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


class _QContainerMeta(type(collections.abc.Container), type(qtc.QObject)):
    """Meta-class for a Container and QObject."""


class SettingsHandler(collections.abc.MutableMapping, qtc.QObject,
                      metaclass=_QContainerMeta):
    """Class that takes care of syncing handler widgets and settings."""

    settings_changed = qtc.pyqtSignal()
    redraw_needed = qtc.pyqtSignal()  # Should trigger redraw of dialog

    def __init__(self, config, parent=None, display_config=False):
        """Initialize instance.

        Parameters
        ----------
        display_config : bool, optional
            Whether the name of the settings file (and path to the file)
            should be displayed. Default is False.

        Returns
        -------
        None
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

        if display_config:
            widget = qtw.QLabel()
            file = config.last_file
            widget.setText(file.stem if file else 'None')
            self.add_static_option(
                'File', 'config', widget,
                display_name='Settings file',
                tooltip=str(file) if file else None,
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

        if not option.read_only:
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
        option = StaticSettingsDialogOption(option_name, handler_widget, *args,
                                            read_only=True, **kwargs)
        self[section_name][option_name] = option
        if self.has_section(section_name):
            self.__sections[section_name].add_option(option)
        else:
            self.__widgets.append(option)

    def add_static_section(self, section_name):
        """Add a static section that may not come from settings."""
        section = SettingsDialogSection(section_name)
        self.__sections[section_name] = section
        self.__widgets.append(section)

    def add_section(self, section_name, **kwargs):
        """Add a titled section to self."""
        if section_name not in self.__config:
            raise ValueError(f"No section {section_name} in config file")
        section = SettingsDialogSection(section_name, **kwargs)
        self.__sections[section_name] = section
        self.__widgets.append(section)

    def has_advanced_options(self):
        """Return whether self contains any advanced option."""
        _adv = any(o.advanced for s in self.values() for o in s.values())
        return _adv or any(s.advanced for s in self.__complex_sections)

    def has_section(self, section):
        """Return whether a section is directly handled."""
        return section in self.__sections

    def make_from_config(self):
        """Create appropriate handlers for entries in self.config."""
        for section_name in self.__config.sections():
            if not self.__config[section_name]:
                # Empty section
                continue
            self.add_section(section_name)
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


class SettingsDialogOption(qtc.QObject):
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
        read_only : bool, optional
            Whether the user is allowed to change the ViPErLEEDSettings
            by interacting with the handler of this SettingsDialogOption.
            Default is False.
        is_advanced : bool, optional
            Whether this option is to be displayed only in "Advanced"
            mode. Default is False.
        label_alignment : {'top', 'centre', 'bottom'}
            Vertical alignment of label field relative to handler_widget.
            Only the first character matters. Any character other than 'c'
            or 'b' will be treated as 'top'. Default is 'top'.
        relevant_for_meas: bool, optional
            Whether this option is to be displayed when
            editing measurement settings. Default is False.
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

        super().__init__(kwargs.pop('parent', None))
        self.option_name = option_name
        self._advanced = kwargs.pop('is_advanced', False)
        self._read_only = bool(kwargs.pop('read_only', False))
        self._relevant_for_meas = kwargs.pop('relevant_for_meas', False)

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

    def __iter__(self):
        """Return the label and the handler for this option."""
        return iter((self.display_name, self.handler_widget))

    def __repr__(self):
        """Return a string representation of self."""
        return f"SettingsDialogOption({self.option_name})"

    @property
    def advanced(self):
        """Return whether this option is considered advanced."""
        return self._advanced

    @advanced.setter
    def advanced(self, is_advanced):
        """Set whether this option is considered advanced."""
        self._advanced = bool(is_advanced)

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

    @property
    def read_only(self):
        """Return whether this option can be modified."""
        return self._read_only

    @property
    def relevant_for_meas(self):
        """Return whether this option is relevant for measurements."""
        return self._relevant_for_meas

    @relevant_for_meas.setter
    def relevant_for_meas(self, relevant):
        """Set whether this option is relevant for measurements."""
        self._relevant_for_meas = bool(relevant)

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
        if not self.read_only:
            to_have += ('notify_',)
        handler = self.handler_widget
        missing = {m: True for m in to_have if not hasattr(handler, m)}
        if any(missing) and not _bind_default_methods(handler, **missing):
            raise ValueError(
                f"Invalid handler_widget {handler.__class__.__name__} "
                f"has no {'/'.join(missing)} method(s)"
                )

    def _connect_handler(self):
        """Connect self.handler_widget if not self.read_only."""
        if self.read_only:
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

        # Get an appropriate size for the info object
        info_size = self.label.fontMetrics().boundingRect(label_text).height()
        self._info = FieldInfo(info_text, size=info_size)

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
        """Update the state of handler_widget according to self.read_only."""
        handler = self.handler_widget

        # See if there is a readOnly method
        try:
            handler.setReadOnly(self.read_only)
        except AttributeError:
            pass
        else:
            return

        # See if there is a writeable read_only property
        try:
            handler.read_only = self.read_only
        except (AttributeError, TypeError):
            # Just disable it instead
            handler.setEnabled(not self.read_only)


class SettingsDialogSectionBase(qtw.QGroupBox):
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
        is_advanced : bool, optional
            Whether this section contains only settings that are to be
            displayed if a user selects the 'Advanced' display.
            Default is False.
        relevant_for_meas: bool, optional
            Whether this this section contains only settings that are to
            be displayed when editing measurement settings.
            Default is False.
        parent : QWidget, optional
            The parent widget of this SettingsDialogSection. Default
            is None.

        Raises
        -------
        TypeError
            If no display_name is given
        """
        display_name = kwargs.get('display_name', '')
        if not display_name:
            raise TypeError("Missing display_name")

        tooltip = kwargs.get('tooltip', '')

        self._advanced = kwargs.get('is_advanced', False)
        self._relevant_for_meas = kwargs.pop('relevant_for_meas', False)

        self._info = qtw.QLabel()
        self.central_widget = qtw.QWidget()
        super().__init__(display_name, kwargs.get('parent', None))

        self.__compose()
        self.set_info(tooltip)

    def __repr__(self):
        """Return a string representation of self."""
        return f"{self.__class__.__name__}(display_name='{self.title()}')"

    @property
    def advanced(self):
        """Return whether this section contains only advanced settings."""
        return self._advanced

    @advanced.setter
    def advanced(self, advanced):
        """Mark this section as containing only advanced settings."""
        self._advanced = bool(advanced)

    @property
    def relevant_for_meas(self):
        """Return whether this section contains only measurement settings."""
        return self._relevant_for_meas

    @relevant_for_meas.setter
    def relevant_for_meas(self, relevant):
        """Mark this section as containing only measurment settings."""
        self._relevant_for_meas = bool(relevant)

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
        layout.addStretch(1)                # Keep all options together

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
        is_advanced : bool, optional
            Whether this section contains only settings that are to be
            displayed if a user selects the 'Advanced' display.
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
        if self.advanced:
            option.advanced = True
        if self.relevant_for_meas:
            option.relevant_for_meas = True
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

        self.__settings = {'current': settings,
                           'applied': copy.deepcopy(settings),
                           'original': copy.deepcopy(settings)}
        self._measurement_settings_only = False

        # Set up children widgets and self
        self.__ctrls = {
            'apply': None,  # Will be an 'Apply' button
            'advanced': qtw.QPushButton("Show less"),
            }
        self.__compose_and_connect()

        # And finally the window properties
        self.update_title(title, settings)
        self.setWindowFlags(self.windowFlags()        # Remove '?'
                            & ~qtc.Qt.WindowContextHelpButtonHint)

    @property
    def adv_button(self):
        """Return the advanced settings button."""
        return self.__ctrls['advanced']

    @property
    def settings(self):
        """Return the settings currently displayed."""
        return self.__settings['current']

    @property
    def handled_object(self):
        """Return the object whose settings are shown."""
        return self._handled_obj

    @property
    def measurement_settings_only(self):
        """Return whether only measurement settings should be shown."""
        return self._measurement_settings_only

    @measurement_settings_only.setter
    def measurement_settings_only(self, meas_only):
        """Set whether only measurement settings should be shown."""
        self._only_show_measurement_settings()
        self._measurement_settings_only = meas_only

    @qtc.pyqtSlot()
    def accept(self):
        """Notify if settings changed, then close."""
        self.__on_apply_pressed()

        # Ask to save the settings to file.
        if self.settings != self.__settings['original']:
            reply = _MSGBOX.question(
                self, "Save settings to file?",
                f"{self.windowTitle()} were edited.\n\n"
                "Would you like to save changes to file?",
                _MSGBOX.Save | _MSGBOX.Discard
                )
            _saved = reply == _MSGBOX.Save
            if _saved:
                try:
                    self.settings.update_file()
                except FileNotFoundError:
                    # The file must have been moved before
                    # the settings could be saved.                              # TODO: open QMessageBox to ask the user how to proceed from here. Ask whether the user wants to save the settings and if yes, ask where and under which name.
                    pass
                finally:
                    self.__settings['original'].read_dict(self.settings)
            self.settings_saved.emit(_saved)
        super().accept()

    @qtc.pyqtSlot()
    def reject(self):
        """Load back the original settings, then close."""
        if not self.isVisible():
            super().reject()
            return

        # Ask confirmation if settings changed
        _changed = self.__ctrls['apply'].isEnabled()
        _changed |= self.__settings['original'] != self.__settings['applied']
        if _changed:
            reply = _MSGBOX.question(
                self, "Discard changes?",
                f"{self.windowTitle()} were edited.\nAre you"
                " sure you want to discard changes?",
                _MSGBOX.Discard | _MSGBOX.Cancel
                )
            if reply == _MSGBOX.Cancel:
                return

        self.settings.read_dict(self.__settings['original'])
        self.__on_apply_pressed()
        if self.__settings['original'] != self.__settings['applied']:
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
                self.__settings[key] = copy.deepcopy(self.settings)
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

    def __compose_and_connect(self):
        """Place and update children widgets."""
        self.handler.update_widgets()  # Fill widgets from settings

        # Dialog buttons
        _bbox = QNoDefaultDialogButtonBox
        buttons = _bbox(_bbox.Ok | _bbox.Cancel | _bbox.Apply)

        self.__ctrls['apply'] = apply_btn = buttons.buttons()[-1]
        apply_btn.setEnabled(False)

        self.adv_button.setDefault(False)
        self.adv_button.setAutoDefault(False)
        self.adv_button.setCheckable(True)

        # Use a ResetRole to have the button placed in a
        # different spot than all others on every platform
        buttons.addButton(self.adv_button, _bbox.ResetRole)
        # self.adv_button.setVisible(self.handler.has_advanced_options())

        layout = qtw.QVBoxLayout()
        for widg in self.handler.widgets:
            if isinstance(widg, SettingsDialogSectionBase):
                layout.addWidget(widg)
                continue
            # An option. Need to add it in a QFormLayout
            _form = qtw.QFormLayout()
            _form.addRow(*widg)
            layout.addLayout(_form)
        layout.addWidget(buttons)
        layout.setSizeConstraint(layout.SetMinimumSize)
        self.setLayout(layout)
        layout.setSpacing(round(layout.spacing()*1.4))

        _policy = self.sizePolicy()
        self.setSizePolicy(_policy.Minimum, _policy.Minimum)

        # Connect signals
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        apply_btn.clicked.connect(self.__on_apply_pressed)
        self.adv_button.toggled.connect(self.__on_show_advanced_toggled)
        self.handler.settings_changed.connect(self.__update_apply_enabled)
        self.handler.redraw_needed.connect(self.__update_advanced_btn)

        self.__update_advanced_btn()

    @qtc.pyqtSlot(bool)
    def __on_apply_pressed(self, _=False):
        """React to the user pressing 'Apply'."""
        self.__ctrls['apply'].setEnabled(False)
        if self.settings != self.__settings['applied']:
            self.settings_changed.emit()
        self.__settings['applied'].read_dict(self.settings)

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
            if isinstance(widg, SettingsDialogSection):
                _section_visible = False
                for option in widg.options:
                    _section_visible |= visible or not option.advanced
                    option.setVisible(visible or not option.advanced)
                widg.setVisible(_section_visible)
            else:
                widg.setVisible(visible or not widg.advanced)
        self.adjustSize()   # TODO: does not always adjust when going smaller?

    def _only_show_measurement_settings(self):
        """Hide all settings that are not relevant for the measurement."""
        for widg in self.handler.widgets:
            if isinstance(widg, SettingsDialogSection):
                _section_visible = False
                for option in widg.options:
                    _section_visible |= option.relevant_for_meas
                    option.setVisible(option.relevant_for_meas)
                widg.setVisible(_section_visible)
            else:
                widg.setVisible(widg.relevant_for_meas)
        self.adjustSize()

    @qtc.pyqtSlot()
    def __update_advanced_btn(self):
        """Update visibility of button and options."""
        self.adv_button.setVisible(self.handler.has_advanced_options())
        self.__on_show_advanced_toggled(self.adv_button.isChecked())

    def __update_apply_enabled(self):
        """Enable/disable 'Apply' depending on whether anything changed."""
        settings_changed = self.settings != self.__settings['applied']
        self.__ctrls['apply'].setEnabled(settings_changed)


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
        read_only : bool
            Whether the user is allowed to change the ViPErLEEDSettings
            by interacting with the handler of this SettingsDialogOption.
            Must be True, otherwise a RuntimeError is raised.
        is_advanced : bool, optional
            Whether this option is to be displayed only in "Advanced"
            mode. Default is False.
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
        RuntimeError
            read_only is not True. Indicates wrong
            use of StaticSettingsDialogOption.
        """
        if kwargs.get('read_only') != True:
            raise RuntimeError(
                'A StaticSettingsDialogOption may only be read_only. You '
                'are seeing this message due to a faulty implementation.'
                )
        super().__init__(option_name, handler_widget, *args, **kwargs)

    def _check_handler(self):
        """Disable handler check."""
