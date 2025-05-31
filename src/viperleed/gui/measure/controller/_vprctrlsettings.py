"""Private module _vprctrlsettings of viperleed.gui.measure.controller.

This module defines widgets used in the settings dialog for
ViPErinoController instances. It is exclusively used internally
by controller.viperinocontroller.ViPErinoController.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-09-02'
__license__ = 'GPLv3+'

import functools
from pathlib import Path
import random
import string

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw

from viperleed.gui.helpers import resources_path
from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.datapoints import QuantityInfo
from viperleed.gui.measure.classes.settings import NotASequenceError
from viperleed.gui.measure.classes import thermocouple
from viperleed.gui.measure.dialogs.settingsdialog import (
    FieldInfo,
    SettingsDialogSectionBase,
    )
from viperleed.gui.measure.serial.viperleedserial import ViPErLEEDHardwareError
from viperleed.gui.measure.serial.viperleedserial import ExtraSerialErrors
from viperleed.gui.measure.widgets.spinboxes import TolerantCommaSpinBox
from viperleed.gui.widgets.lib import change_control_text_color
from viperleed.gui.widgets.lib import move_to_front

# pylint: disable=too-many-lines
# Makes sense to keep all the widgets in a single module


_INVOKE = qtc.QMetaObject.invokeMethod
_NOT_SET = "** NOT SET **"
_UNKNOWN = QuantityInfo.UNKNOWN
_UNIQUE = qtc.Qt.UniqueConnection
_QMSG = qtw.QMessageBox


class FWVersionViewer(qtw.QLabel):
    """Simple class to show a read-only firmware version.

    Normally only one value is shown (as it should be the same in
    hardware and config). However, should the two differ, a warning
    text is displayed.
    """

    def __init__(self, controller, **kwargs):
        """Initialise instance."""
        self.controller = controller
        super().__init__(**kwargs)
        self.__normal_color = self.palette().color(self.palette().WindowText)

    def get_(self):
        """Return the value of the firmware version."""
        ctrl = self.controller
        with ctrl.lock:
            hardware_v = ctrl.hardware.get('firmware', None)
        if not hardware_v:
            return ctrl.settings.get('controller', 'firmware_version')
        return str(hardware_v)

    def set_(self, config_v):
        """Set the value displayed."""
        hardware_v = self.controller.firmware_version
        txt = f"v{hardware_v}"
        color = self.__normal_color
        tip_txt = ''
        if hardware_v > config_v:
            color = 'red'
            tip_txt = ("ERROR: your settings are for an earlier version "
                       f"(v{config_v})! It may be impossible to perform "
                       "a measurement")
        elif hardware_v < config_v:
            color = 'orange'
            tip_txt = ("WARNING: your settings are for a later version "
                       f"(v{config_v})! You may want to upgrade your "
                       "firmware version in the Tools menu")
        change_control_text_color(self, color)
        self.setText(txt)
        self.setToolTip(tip_txt)


class SerialNumberEditor(qtw.QWidget):
    """A class that allows viewing and setting a serial number."""

    serial_number_changed = qtc.pyqtSignal()

    def __init__(self, controller, **kwargs):
        """Initialise instance."""
        super().__init__(**kwargs)
        self.__ctrl = controller

        self.__edit = qtw.QLineEdit()
        self.__rand_btn = qtw.QPushButton("Generate randomly")
        self.__set_btn = qtw.QPushButton("Set")
        self.__old_serial = ''
        self.notify_ = self.serial_number_changed

        self.__compose()
        self.__connect()

    @staticmethod
    def get_random_serial():
        """Return a valid, random serial number."""
        valid_chars = string.ascii_letters + string.digits
        return ''.join(random.choices(valid_chars.upper(), k=4))

    def get_(self):
        """Return the value to be stored in the config."""
        txt = self.__ctrl.name
        *_, serial = txt.split()
        if not self.valid_serial(serial):
            txt = "ViPErLEED ____"
        return txt

    def set_(self, device_name):
        """Set displayed value."""
        *_, serial = device_name.split()
        text = serial.upper()
        if not self.valid_serial(text):
            # Probably serial number is not set
            text = _NOT_SET
        self.__edit.setText(text)

    def valid_serial(self, serial):
        """Return whether serial is a valid serial number."""
        validator = self.__edit.validator()
        return validator.validate(serial, 0) != validator.Invalid

    def __compose(self):
        """Place children widgets."""
        self.__edit.setValidator(
            qtg.QRegExpValidator(qtc.QRegExp("[a-zA-Z0-9]"*4))
            )
        width = self.__edit.fontMetrics().boundingRect(_NOT_SET).width()
        width += 10
        self.__edit.setMinimumWidth(width)
        self.__edit.setMaximumWidth(width)

        policy = self.__set_btn.sizePolicy()
        self.__set_btn.setSizePolicy(policy.Fixed, policy.Fixed)
        self.__set_btn.setEnabled(False)

        layout = qtw.QHBoxLayout()
        layout.addWidget(self.__edit)
        layout.addWidget(self.__rand_btn)
        layout.addWidget(self.__set_btn)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

    def __connect(self):
        """Connect children signals."""
        self.__edit.editingFinished.connect(self.__on_edit_finished)
        self.__edit.textChanged.connect(self.__on_serial_text_changed)
        self.__rand_btn.clicked.connect(self.__on_make_random_clicked)
        self.__set_btn.clicked.connect(self.__on_set_serial_clicked)

    def __on_edit_finished(self):
        """Make serial number uppercase after editing is over."""
        self.__edit.setText(self.__edit.text().upper())

    def __on_hardware_info(self):
        """Update after hardware information arrived."""
        base.safe_disconnect(self.__ctrl.hardware_info_arrived,
                             self.__on_hardware_info)
        _INVOKE(self.__ctrl, 'disconnect_')
        with self.__ctrl.lock:
            new_serial = self.__ctrl.hardware.get('serial_nr', _NOT_SET)
        self.__edit.setText(new_serial)
        if new_serial == self.__old_serial:
            return
        self.serial_number_changed.emit()

        # The config name must change to keep consistency:
        settings = self.__ctrl.settings
        old_name = settings.last_file
        new_name = old_name.with_name(f"{self.__ctrl.name_clean}.ini")
        with new_name.open('w', encoding='utf-8') as fproxy:
            settings.write(fproxy)

        # What to do with the old one?
        msg = _QMSG(_QMSG.Question, "Delete old configuration file?",
                    f"Serial number changed from {self.__old_serial} to "
                    f"{new_serial}. The configuration file stored will be "
                    "renamed. Would you like to keep the old file too?",
                    parent=self)
        msg.addButton("Keep both files", _QMSG.YesRole)
        delete = msg.addButton("Delete old file", _QMSG.NoRole)
        msg.exec_()

        if msg.clickedButton() == delete:
            old_name.unlink(missing_ok=True)

    def __on_make_random_clicked(self):
        """Generate a random number."""
        self.__edit.setText(self.get_random_serial())

    def __on_serial_text_changed(self):
        """Activate "Set" if serial is changed."""
        txt = self.__edit.text().upper()
        with self.__ctrl.lock:
            hw_serial = self.__ctrl.hardware.get('serial_nr', None)
        if hw_serial is None or len(txt) != 4:
            self.__set_btn.setEnabled(False)
            return
        active = self.valid_serial(txt)
        active &= txt != hw_serial
        self.__set_btn.setEnabled(active)

    def __on_set_serial_clicked(self):
        """Set serial number, then update hardware info."""
        with self.__ctrl.lock:
            self.__old_serial = self.__ctrl.hardware.get('serial_nr', _NOT_SET)
        # NB: New and old are different, otherwise "Set" is disabled
        new_serial = self.__edit.text().upper()

        if self.__old_serial == _NOT_SET:
            reply = _QMSG.Yes
        else:
            reply = _QMSG.question(
                self, "Set new serial number?",
                "Are you sure you want to replace your serial number "
                f"({self.__old_serial}) with {new_serial}?",
                _QMSG.Yes | _QMSG.Cancel
                )
        if reply == _QMSG.Cancel:
            # Set back the old value
            self.__edit.setText(self.__old_serial)
            return
        _QMSG.information(self, "About to set new serial number",
                          "Setting new serial number. Do not "
                          "disconnect the USB cable to the unit.",
                          _QMSG.Ok)

        # ViPErinoController.set_serial_number also gets hardware
        # info after successfully setting the new serial number
        base.safe_connect(self.__ctrl.hardware_info_arrived,
                          self.__on_hardware_info, type=_UNIQUE)
        _INVOKE(self.__ctrl, 'set_serial_number', qtc.Q_ARG(str, new_serial))


class HardwareConfigurationEditor(SettingsDialogSectionBase):
    """Class for viewing and setting ADC inputs.

    Provides one line per each ADC present, one entry per each channel.
    When the user clicks on an editable quantity, a dialog pops up.
    Editable quantities are those for which one can potentially select
    a gain factor or a calibration curve. If a channel can only measure
    a single non-editable quantity, user interaction with it is disabled.

    This widget should be used as a special section of a settings
    dialog, as it internally handles all the relevant values stored
    in the configuration file.
    """

    def __init__(self, controller=None, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        controller : ViPErinoController
            The controller whose hardware configuration
            settings are managed.
        **kwargs : dict
            Optional arguments passed on to SettingsDialogSectionBase

        Returns
        -------
        None.
        """
        self.__ctrl = controller

        # Modify arguments for the following super() call
        kwargs['display_name'] = "Hardware configuration"
        kwargs['is_advanced'] = True
        kwargs['tooltip'] = (
            "This section lists what each hardware channel "
            "is capable of measuring. Quantities appearing "
            "on the same line cannot be measured simultaneously."
            )
        super().__init__(**kwargs)

        # Set up some static children widgets
        self.__no_power_warning = qtw.QWidget()
        self.__try_again = qtw.QPushButton("Retry")
        self.__hw_info = qtw.QWidget()  # just a container
        self.__adcs = []  # Will be a list of {channel: combo}

        self.__dialogs = {
            QuantityInfo.I0: _I0EditDialog(self.__ctrl),
            QuantityInfo.AUX: _AUXEditDialog(self.__ctrl),
            QuantityInfo.TEMPERATURE: _TemperatureEditDialog(self.__ctrl),
            }

        self.__compose_static_children()
        self.__connect()

        # QueuedConnection prevents serial port not open errors
        self.__ctrl.error_occurred.connect(self.__on_ctrl_error,
                                           type=qtc.Qt.QueuedConnection)

    @property
    def advanced(self):
        """Return whether this section contains only advanced settings."""
        # We want this section to be normally visible (i.e., not
        # advanced) only if there are some problems with the settings.
        # We want to see this if
        # (1) we don't have info yet
        if not self.__adcs:
            return False
        # (2) something is wrong with the current selection
        if any('?' in combo.currentText()
               for adc in self.__adcs for combo in adc.values()):
            return False
        return super().advanced

    @property
    def hardware_info(self):
        """Return a dictionary of information suitable for a settings file."""
        info = {}
        devices = []
        for adc in self.__adcs:
            devices.append(())
            for channel, combo in adc.items():
                quantity = combo.current_quantity
                if quantity is _UNKNOWN:
                    continue
                info[quantity.label] = str(channel)
                devices[-1] += (quantity.label,)

        with self.__ctrl.lock:
            if self.__ctrl.hardware['lm35']:
                cjc = QuantityInfo.COLD_JUNCTION.label
                devices.append((cjc,))
                info[cjc] = '0'  # unused anyway

        info['measurement_devices'] = str(tuple(devices))
        return info

    @qtc.pyqtSlot()
    def update_widgets(self):
        """Place children widgets."""
        with self.__ctrl.lock:
            detected_adcs = [k for k, present in self.__ctrl.hardware.items()
                             if 'adc' in k and present]
        has_power = any(detected_adcs)
        if has_power:
            self.__show_adcs(detected_adcs)

        self.__no_power_warning.setVisible(not has_power)
        self._info.setVisible(has_power)
        self.__hw_info.setVisible(has_power)
        self.__connect()
        self.adjustSize()
        self.updated.emit()

    def __clear_hw_info(self):
        """Remove all widgets from __hw_info."""
        layout = self.__hw_info.layout()
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()

    def __compose_static_children(self):
        """Set up some static children that are shown under some conditions."""
        lay = qtw.QHBoxLayout()
        _label = qtw.QLabel("No hardware information. Power\n"
                            "supply may be disconnected.")
        _policy = _label.sizePolicy()
        _label.setSizePolicy(_policy.Fixed, _policy.Preferred)
        lay.addWidget(_label)
        lay.addWidget(self.__try_again)
        lay.setContentsMargins(0, 0, 0, 0)
        self.__no_power_warning.setLayout(lay)

        self.__hw_info.setLayout(qtw.QGridLayout())
        self.__hw_info.layout().setContentsMargins(0, 0, 0, 0)
        self.__hw_info.layout().setSizeConstraint(lay.SetMinimumSize)
        self.__hw_info.setSizePolicy(_policy.Minimum, _policy.Minimum)

        # Now layout of the central widget: just a box. Only one
        # of its contents will ever be shown at a time
        central_layout = qtw.QHBoxLayout()
        central_layout.addWidget(self.__no_power_warning)
        central_layout.addWidget(self.__hw_info)
        central_layout.setContentsMargins(0, 0, 0, 0)
        central_layout.setSizeConstraint(lay.SetMinimumSize)
        self.central_widget.setLayout(central_layout)

        self.__try_again.clicked.connect(self.__ctrl.get_hardware)
        self.__ctrl.hardware_info_arrived.connect(self.update_widgets)

    def __connect(self):
        """Connect (only once) relevant signals and slots."""
        _connect = base.safe_connect
        for dialog in self.__dialogs.values():
            try:
                _connect(self.window().finished, dialog.reject, type=_UNIQUE)
            except AttributeError:
                pass
            _connect(dialog.settings_changed, self.settings_changed,
                     type=_UNIQUE)
            _connect(dialog.accepted, self.__on_settings_changed, type=_UNIQUE)

    def __get_quantity_combo(self, quantity):
        """Return an _ADCChannelCombo for a quantity."""
        combo = _ADCChannelCombo(ctrl=self.__ctrl)
        combo.add_quantity(quantity)
        combo.current_quantity = quantity
        combo.quantity_selected.connect(self.__on_quantity_selected)

        for i in range(combo.count()):
            dialog = self.__dialogs.get(combo.itemData(i), None)
            if dialog:
                dialog.settings_changed.connect(combo.update_items)
        return combo

    @qtc.pyqtSlot(tuple)
    def __on_ctrl_error(self, error_info):
        """React to to controller errors."""
        # pylint: disable=confusing-consecutive-elif
        # This is triggered by the inner "if reply ==" in the first
        # case. This would be an ideal use case for a switch, but
        # we want to maintain compatibility with Win7, i.e., py<=3.8
        *_, err_msg = error_info
        if error_info is ViPErLEEDHardwareError.ERROR_NO_HARDWARE_DETECTED:
            # No power. Suggest to try again.
            reply = _QMSG.warning(self, "No hardware info", err_msg,
                                  _QMSG.Retry | _QMSG.Ignore)
            if reply == _QMSG.Retry:
                self.__try_again.click()
        elif error_info is ViPErLEEDHardwareError.ADC_POWER_FAULT:
            # Info about available ADCs is inconsistent. Probably
            # something is wrong with the hardware itself.
            self.__on_hardware_fault(err_msg)
        elif error_info is ExtraSerialErrors.PORT_NOT_OPEN:
            _QMSG.critical(
                self, "Device lost",
                err_msg + "\n\nThe port may be in use or disconnected.",
                _QMSG.Ok
                )
            self.window().close()

    def __on_hardware_fault(self, err_msg):
        """React to a hardware problem."""
        msg = _QMSG(_QMSG.Critical, "Hardware Fault", err_msg,
                    _QMSG.Ok, self)
        if self.__ctrl.firmware_version >= 1.0:
            # Only for public releases of ViPErLEED!
            # TODO: Is there a decent way to pick the file name?
            _path = Path(resources_path(
                'hardware/schematics/viperLEED_HW_v8 - basic_configuration.pdf'
                )).resolve()
            show_schem = msg.addButton("Open schematics...", _QMSG.ActionRole)
            show_schem.disconnect()  # Box stays open when clicked
            show_schem.clicked.connect(functools.partial(
                qtg.QDesktopServices.openUrl,
                qtc.QUrl.fromLocalFile(str(_path))
                ))
        msg.exec_()
        self.window().close()

    def __on_quantity_selected(self, quantity):
        """Show a dialog for editing a measurable quantity."""
        this_combo = self.sender()
        previous_quantity = this_combo.previous_quantity
        dialog = self.__dialogs.get(quantity, None)
        if dialog:
            move_to_front(dialog)
        if previous_quantity is quantity:
            return

        if not dialog:
            self.__on_settings_changed()

        for combo in (c for adc in self.__adcs for c in adc.values()):
            combo.enable_quantity(previous_quantity)
            if combo is not this_combo and quantity is not _UNKNOWN:
                combo.enable_quantity(quantity, False)

    def __on_settings_changed(self):
        """React to a change of settings in one of the dialogs."""
        # Update settings with the info from the combos
        new_settings = {'controller': self.hardware_info}
        settings = self.__ctrl.settings

        # Remove all current channel/quantity pairs, then update
        for quantity in QuantityInfo:
            settings.remove_option('controller', quantity.label)
        settings.read_dict(new_settings)
        self.settings_changed.emit()

    def __show_adcs(self, detected_adcs):
        """Add one line per active ADC."""
        self.__adcs = []
        try:
            all_adcs = self.__ctrl.available_adcs()
        except (NotASequenceError, KeyError,
                TypeError, ValueError, RuntimeError):
            # Settings are wrong. User will have to fix this.
            all_adcs = ((_UNKNOWN, _UNKNOWN),)*len(detected_adcs)
            with self.__ctrl.lock:
                if self.__ctrl.hardware['lm35']:
                    all_adcs += ((QuantityInfo.COLD_JUNCTION,),)

        # Now populate self.__hw_info, adding one row per ADC,
        # one column per channel. Clicking on combo item allows
        # to edit.
        self.__clear_hw_info()
        layout = self.__hw_info.layout()
        for row, adc in enumerate(all_adcs, start=1):
            if QuantityInfo.COLD_JUNCTION in adc:
                continue
            self.__adcs.append({})
            _label = qtw.QLabel(f"ADC#{row - 1}:")
            _label.setSizePolicy(_label.sizePolicy().Fixed,
                                 _label.sizePolicy().Preferred)
            layout.addWidget(_label, row, 0)
            for col, quantity in enumerate(adc, start=1):
                if quantity is not _UNKNOWN:  # Use the known channel
                    col = adc[quantity] + 1
                combo = self.__get_quantity_combo(quantity)
                self.__adcs[-1][col - 1] = combo
                layout.addWidget(combo, row, col)

        # Now header line with channel numbers
        for col in range(layout.columnCount()):
            if not col:  # skip the one with "ADC#x" labels
                continue
            layout.addWidget(qtw.QLabel(f"CH{col - 1}"), 0, col)

        # Finally "+" and "-" buttons for adding/removing channels              # TODO. Also, have one button to start fresh from ?? everywhere (in case of fucked up settings)


class _ADCChannelCombo(qtw.QComboBox):
    """A combo for showing/editing quantities measured at an ADC channel.

    The combo is automatically disabled if the user is not supposed
    to change anything. Items have associated tooltips.
    """

    quantity_selected = qtc.pyqtSignal(QuantityInfo)

    def __init__(self, ctrl=None, parent=None):
        """Initialise instance."""
        self.__ctrl = ctrl
        self.__quantities = {}
        self.__last_idx = -1
        super().__init__(parent)
        self.activated.connect(self.__on_qty_clicked)

        red_palette = self.palette()
        red_palette.setColor(red_palette.Normal,  # only if not disabled
                             red_palette.Text, qtg.QColor('red'))
        red_palette.setColor(red_palette.Inactive,
                             red_palette.Text, qtg.QColor('red'))
        self.__palettes = {True: self.palette(),  # selection is OK
                           False: red_palette}    # selection is not OK

    @property
    def current_quantity(self):
        """Return the QuantityInfo currently selected."""
        return self.itemData(self.currentIndex())

    @current_quantity.setter
    def current_quantity(self, quantity):
        """Set quantity as the quantity currently selected."""
        _name, _tooltip, _ok = self.__qty_name_and_tooltip(quantity)
        self.setCurrentText(_name)
        self.setToolTip(_tooltip)

        # Decide if self should be enabled or not:
        # Enabled if: the current selection has problems, or
        # there are more options to pick from, or it is one
        # of those quantities that the user can edit
        _editable = (QuantityInfo.I0, QuantityInfo.TEMPERATURE,
                     QuantityInfo.AUX)
        enabled = not _ok
        if self.count() != 1:
            enabled = True
        elif quantity in _editable:
            enabled = True
        self.setEnabled(enabled)
        self.setPalette(self.__palettes[_ok])
        self.__last_idx = self.currentIndex()

    @property
    def previous_quantity(self):
        """Return the quantity that was selected before the last change."""
        return self.itemData(self.__last_idx)

    def add_quantity(self, quantity):
        """Add a quantity to this combo.

        This may also add other quantities that are known to be
        measurable by the same ADC channel.

        Parameters
        ----------
        quantity : QuantityInfo
            The quantity to be added. Text and tooltip is
            automatically generated.
        """
        self.__add_one_qty(quantity)

        # Now decide if we should add other quantities that may be
        # measured at the same channel. In the future, we can check
        # self.__ctrl.firmware_version, to decide if others should
        # be added.
        if quantity is QuantityInfo.TEMPERATURE:
            self.__add_one_qty(QuantityInfo.AUX)
        elif quantity is QuantityInfo.AUX:
            self.__add_one_qty(QuantityInfo.TEMPERATURE)
        elif quantity is _UNKNOWN:
            # Here add all measurable quantities
            measurable = (
                QuantityInfo.HV, QuantityInfo.I0, QuantityInfo.ISAMPLE,
                QuantityInfo.TEMPERATURE, QuantityInfo.AUX
                )
            for extra in measurable:
                self.__add_one_qty(extra)

    def enable_quantity(self, quantity, enabled=True):
        """Disable a quantity."""
        if quantity not in self.__quantities:
            return
        _name = self.__quantities[quantity]
        item = self.model().findItems(_name)[0]
        if enabled:
            flags = item.flags() | qtc.Qt.ItemIsEnabled
        else:
            flags = item.flags() & ~qtc.Qt.ItemIsEnabled
        item.setFlags(flags)

    def update_items(self):
        """Update text and tooltips of all items."""
        for i in range(self.count()):
            quantity = self.itemData(i)
            _name, _tooltip, _ = self.__qty_name_and_tooltip(quantity)
            self.setItemData(i, _name, qtc.Qt.DisplayRole)
            self.__quantities[quantity] = _name
            self.setItemData(i, _tooltip, qtc.Qt.ToolTipRole)
        _ok = '?' not in self.currentText()
        self.setPalette(self.__palettes[_ok])

    def __add_one_qty(self, quantity):
        """Add one quantity to self."""
        _name, _tooltip, _ = self.__qty_name_and_tooltip(quantity)
        self.addItem(_name, userData=quantity)
        self.setItemData(self.findText(_name), _tooltip, qtc.Qt.ToolTipRole)
        self.__quantities[quantity] = _name

    def __on_qty_clicked(self, index):
        """React to a user clicking on a quantity."""
        quantity = self.itemData(index)
        name = self.itemText(index)
        tooltip = self.itemData(index, qtc.Qt.ToolTipRole)
        self.setToolTip(tooltip)
        self.quantity_selected.emit(quantity)

        _ok = '?' not in name
        self.setPalette(self.__palettes[_ok])
        self.__last_idx = index

    def __qty_name_and_tooltip(self, quantity):
        """Return name and tooltip for quantity, and whether quantity is OK."""
        with self.__ctrl.lock:
            hardware = self.__ctrl.hardware.copy()
        _tooltip = quantity.description
        _name = quantity.label.replace("_", " ")
        if quantity.units:
            _name += f" ({quantity.units})"
        _qty_ok = quantity is not _UNKNOWN

        if quantity is QuantityInfo.I0:
            _name = f"I\u2080 ({quantity.units})"
            _tooltip += f". Range: {hardware['i0_range']}"
        elif quantity is QuantityInfo.ISAMPLE:
            # Unicode chars are for "sample" as subscripts
            _name = f"I\u209b\u2090\u2098\u209a\u2097\u2091 ({quantity.units})"
        elif quantity is QuantityInfo.TEMPERATURE:
            _tc_type = "??"
            _qty_ok = False
            if self.__ctrl.thermocouple:
                _tc_type = self.__ctrl.thermocouple.type_
                _qty_ok = True
            _name += f" \u2013 {_tc_type} thermocouple"
            _tooltip += f"Thermovoltage range: {hardware['aux_range']}"
        elif quantity is QuantityInfo.AUX:
            _tooltip += f"Range: {hardware['aux_range']}"

        return _name, _tooltip, _qty_ok


class _InputRangeSelector(qtw.QWidget):
    """A widget to show/edit input ranges."""

    range_changed = qtc.pyqtSignal(str)  # The new range

    def __init__(self, **kwargs):
        """Initialise instance.

        Parameters
        ----------
        parent : QWidget, optional
            The parent widget of self
        ranges : Sequence, optional
            The ranges among which this selector allows to choose.
            If not given, the keys in tooltips are used as ranges,
            if any. Otherwise fall back to ("0 \u2013 2.5 V",
            "0 \u2013 10 V").
        tooltips : dict, optional
            keys are ranges, values as strings to be used as help text
        help_file : pathlib.Path, optional
            Path to the file to be opened when range switching can
            only be performed manually. The file is supposed to help
            users perform the range switch. The user is presented with
            the option of showing this file only if not self.editable

        Returns
        -------
        None.
        """
        super().__init__(kwargs.get('parent', None))

        if not kwargs.get('tooltips', None):
            kwargs['tooltips'] = {}

        ranges = kwargs.get('ranges', tuple())
        if not ranges:
            if kwargs['tooltips']:
                ranges = kwargs['tooltips'].keys()
            else:
                ranges = ("0 \u2013 2.5 V", "0 \u2013 10 V")

        self.__range_options = {r: qtw.QRadioButton(r) for r in ranges}
        self.__tooltips = []
        self.__btn_group = qtw.QButtonGroup()
        self.__range_change_info = qtw.QPushButton(
            "Help me to\nswitch range..."
            )
        self.__range_change_info.setEnabled(False)
        self.__range_change_info.setAutoDefault(False)
        self.__range_change_info.setDefault(False)

        help_file = kwargs.get('help_file', None)
        if help_file is not None and help_file.exists():
            self.__range_change_info.clicked.connect(functools.partial(
                qtg.QDesktopServices.openUrl,
                qtc.QUrl.fromLocalFile(str(help_file))
                ))
            self.__range_change_info.setEnabled(True)

        self.editable = False
        self.__compose(kwargs)

        self.__btn_group.buttonClicked.connect(self.__on_range_changed)

    @property
    def range_(self):
        """Return the selected range."""
        for range_, btn in self.__range_options.items():
            if btn.isChecked():
                return range_
        return None

    @range_.setter
    def range_(self, new_range):
        """Set the selected range."""
        _btn = self.__range_options.get(new_range, None)
        if _btn is not None:
            _btn.setChecked(True)
            return

        # Uncheck all
        self.__btn_group.setExclusive(False)
        for _btn in self.__range_options.values():
            _btn.setChecked(False)
        self.__btn_group.setExclusive(True)

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Show widget."""
        self.update_widgets()
        super().showEvent(event)

    def update_widgets(self):
        """Update the state of children."""
        for widg in (*self.__range_options.values(), *self.__tooltips):
            widg.setEnabled(self.editable)
        _enabled = self.__range_change_info.isEnabled()
        self.__range_change_info.setVisible(_enabled and not self.editable)

    def __compose(self, kwargs):
        """Place children widgets."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.__make_ranges_group(kwargs))
        layout.addWidget(self.__range_change_info)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        # Make help button a bit larger
        _btn = self.__range_change_info
        width, height = _btn.sizeHint().width(), _btn.sizeHint().height()
        _btn.setMinimumSize(round(width*1.2), round(height*1.4))

    def __make_ranges_group(self, kwargs):
        """Prepare a group box with as many buttons as ranges."""
        quantity = kwargs.get('quantity', '')
        tooltips = kwargs['tooltips']
        _range_btns = qtw.QGroupBox(f"{quantity} input range")
        _btns_lay = qtw.QVBoxLayout()
        for range_, btn in self.__range_options.items():
            _lay = qtw.QHBoxLayout()
            _lay.addWidget(btn)
            self.__btn_group.addButton(btn)
            _btns_lay.addLayout(_lay)
            tooltip = tooltips.get(range_, "")
            if not tooltip:
                continue
            _size = btn.fontMetrics().boundingRect(btn.text()).height()
            info = FieldInfo(tooltip, _size)
            _lay.addWidget(info)
            _lay.addStretch(1)
            self.__tooltips.append(info)
        _range_btns.setLayout(_btns_lay)
        return _range_btns

    def __on_range_changed(self, _):
        """React to the user picking a different range."""
        self.range_changed.emit(self.range_)


class _EditDialogBase(qtw.QDialog):
    """Base class for dialogs for editable quantities."""

    settings_changed = qtc.pyqtSignal()

    def __init__(self, controller, *args, quantity='', input_quantity='',
                 raw_quantity='', tooltips=None, help_file=None, **kwargs):
        """Initialise dialog."""
        super().__init__(*args, **kwargs)
        self._ctrl = controller
        self.__raw_quantity = raw_quantity

        if not quantity:
            quantity = raw_quantity
        if not input_quantity:
            input_quantity = raw_quantity

        # Input range selector + info on how to change in
        # case no relay is present for electronic switching
        self._input_range = _InputRangeSelector(quantity=input_quantity,
                                                tooltips=tooltips,
                                                help_file=help_file)
        self.central_widget = qtw.QWidget()
        self.__compose_and_connect()

        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)
        self.setWindowTitle(f"Edit {quantity} measurement")

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Show dialog and its children."""
        self.update_widgets()
        super().showEvent(event)

    def update_widgets(self):
        """Update widget contents from the controller."""
        with self._ctrl.lock:
            hardware = self._ctrl.hardware
            range_ = hardware.get(f"{self.__raw_quantity}_range".lower(), None)
        self._input_range.range_ = range_
        self._input_range.update_widgets()

    def __compose_and_connect(self):
        """Place children and connect signals."""
        # Dialog buttons
        buttons = qtw.QDialogButtonBox()
        _done = buttons.addButton("Done", buttons.AcceptRole)

        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = qtw.QVBoxLayout()
        layout.addWidget(self._input_range)
        layout.addWidget(self.central_widget)
        layout.addStretch(1)
        layout.addWidget(buttons)
        self.setLayout(layout)

        _done.setAutoDefault(True)
        _done.setDefault(True)
        _done.setFocus()


class _I0EditDialog(_EditDialogBase):
    """Dialog for editing I0 input: gain and input range."""

    _suffixes = {"0 \u2013 2.5 V": " mV/µA",
                 "0 \u2013 10 V": " V/µA",
                 None: ""}
    _gains = {"0 \u2013 2.5 V": "kohm",
              "0 \u2013 10 V": "Mohm",
               None: "kohm/Mohm"}
    _gain_tip_base = (
        "<nobr>Change this value only if you are not using the default</nobr> "
        "<nobr>1 {0}</nobr> current-to-voltage conversion impedance. "
        "E.g., use gain 0.5 if you use a <nobr>0.5 {0}</nobr> impedance."
        )

    def __init__(self, controller, *args, **kwargs):
        """Initialise dialog."""
        kwargs['quantity'] = "I0"
        kwargs['input_quantity'] = "I\u2080"
        kwargs['raw_quantity'] = "i0"
        kwargs['tooltips'] = {
            "0 \u2013 2.5 V": (
                "<nobr>This range is most commonly used when</nobr> "
                "measuring I0 as the voltage drop across a 1 kohm "
                "resistor. <b>Use this for a SPECS ErLEED optics</b>"
                ),
            "0 \u2013 10 V": (
                "<nobr>This range is most commonly used when</nobr> "
                "measuring I0 as the voltage drop across a 1 Mohm "
                "resistor. <b>Use this for an Omicron SPECTALEED optics<b>"
                )
            }
        kwargs['help_file'] = Path(resources_path(
            'hardware/schematics/viperLEED_HW_v8 - jumpers.pdf'
            ))
        super().__init__(controller, *args, **kwargs)

        self.__gain = TolerantCommaSpinBox()
        self.__gain_info = None

        self.__compose()
        self.__connect()

    def update_widgets(self):
        """Update widget contents from the controller."""
        with self._ctrl.lock:
            self._input_range.editable = self._ctrl.hardware.get('relay',
                                                                 False)
        super().update_widgets()
        self.__gain.setSuffix(self._suffixes[self._input_range.range_])

        settings = self._ctrl.settings
        try:
            gain = settings.getfloat("conversions", "i0_gain", fallback=None)
        except (TypeError, ValueError):
            gain = None
        if gain is None:
            gain = 1.0
            settings.set("conversions", "i0_gain", "1.0")
        self.__gain.setValue(gain)
        self.__update_gain_info()

    def __compose(self):
        """Place children widgets."""
        _label = qtw.QLabel("I\u2080 gain")
        _policy = _label.sizePolicy()
        _label.setSizePolicy(_policy.Fixed, _policy.Preferred)
        height = _label.fontMetrics().boundingRect(_label.text()).height()
        self.__gain_info = FieldInfo(size=height)
        self.__gain_info.setSizePolicy(_policy.Fixed, _policy.Fixed)
        self.__update_gain_info()

        layout = qtw.QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(_label)
        layout.addWidget(self.__gain_info)
        layout.addWidget(self.__gain)

        self.central_widget.setLayout(layout)

        self.__gain.setRange(float('-inf'), float('inf'))
        self.__gain.setDecimals(8)

    def __connect(self):
        """Connect appropriate signals."""
        self.__gain.valueChanged.connect(self.__on_gain_changed)
        self._input_range.range_changed.connect(self.__on_range_changed)

    def __on_gain_changed(self, new_gain):
        """React to a change in the gain factor."""
        try:
            old_gain = self._ctrl.settings.getfloat(
                "conversions", "i0_gain", fallback=None
                )
        except (TypeError, ValueError):
            old_gain = None
        if new_gain != old_gain:
            self._ctrl.settings.set("conversions", "i0_gain",
                                     str(round(new_gain, 8)))
            self.settings_changed.emit()

    def __on_range_changed(self, new_range):  # pylint: disable=unused-argument
        """Set new ADC range in the box."""
        # This is currently not possible, as there is no code for
        # electronically switching the range by using an external
        # relay. Here one would place an invocation to .connect_
        # and to switch_relay (or whatever the method is called)
        self.__update_gain_info()

    def __update_gain_info(self):
        """Update tooltip text for the gain widgets."""
        range_ = self._input_range.range_
        _tip = self._gain_tip_base.format(self._gains[range_])
        self.__gain_info.set_info_text(_tip)


class _TemperatureEditDialog(_EditDialogBase):                                  # TODO: after check of TC, add offset correction
    """Dialog for editing the TEMPERATURE input.

    Allows to select which thermocouple is used. Displays
    the input range used for reading the thermovoltage.
    """

    def __init__(self, controller, *args, **kwargs):
        """Initialise dialog."""
        kwargs['quantity'] = "Temperature"
        kwargs['input_quantity'] = "AUX"
        kwargs['raw_quantity'] = "aux"
        kwargs['tooltips'] = {
            "0 \u2013 2.5 V": (
                "<nobr>This range is the only range that make sense </nobr> "
                "for measuring the thermoelectric voltage of a thermocouple"
                ),
            "0 \u2013 10 V": (
                "<nobr>It is <b>not a good idea</b> to use this range for "
                "measuring the</nobr> thermoelectric voltage of a "
                "thermocouple. Temperature measurements in this range are "
                "likely inaccurate (thermovoltages are only a few millivolts)"
                )
            }
        kwargs['help_file'] = Path(resources_path(
            'hardware/schematics/viperLEED_HW_v8 - jumpers.pdf'
            ))
        super().__init__(controller, *args, **kwargs)
        self.__thermocouples = qtw.QComboBox()
        self.__cjc = qtw.QCheckBox()
        self.__cjc_not_found = qtw.QWidget()

        self.__compose(kwargs)
        self.__connect()

    def update_widgets(self):
        """Update widget contents from the controller."""
        super().update_widgets()
        idx = self.__thermocouples.findText(str(self._ctrl.thermocouple))
        self.__thermocouples.setCurrentIndex(idx)

        with self._ctrl.lock:
            cjc_available = self._ctrl.hardware['lm35']
        self.__cjc.setChecked(cjc_available)
        self.__cjc_not_found.setVisible(not cjc_available)
        self.adjustSize()

    def __compose(self, kwargs):
        """Place children widgets."""
        self.__thermocouples.addItem(_NOT_SET, None)
        for t_couple in thermocouple.THERMOCOUPLES:
            self.__thermocouples.addItem(str(t_couple), t_couple)

        _policy = self.sizePolicy()
        _cjc_text = qtw.QLabel(
            "No LM35 device detected for measuring a cold-junction "
            "compensation temperature. Temperatures measured with "
            "your thermocouple will be inaccurate. Consider installing one."
            )
        _cjc_text.setWordWrap(True)
        _cjc_text.setSizePolicy(_policy.Preferred, _policy.Minimum)
        change_control_text_color(_cjc_text, 'red')
        _cjc_help = qtw.QPushButton("Show me\nhow...")
        _cjc_help.setSizePolicy(_policy.Preferred, _policy.Expanding)
        _cjc_help.clicked.connect(functools.partial(
            qtg.QDesktopServices.openUrl,
            qtc.QUrl.fromLocalFile(str(kwargs['help_file'].resolve()))
            ))
        _no_cjc_lay = qtw.QHBoxLayout()
        _no_cjc_lay.setContentsMargins(0, 0, 0, 0)
        _no_cjc_lay.addWidget(_cjc_text)
        _no_cjc_lay.addWidget(_cjc_help)
        self.__cjc_not_found.setLayout(_no_cjc_lay)
        self.__cjc.setEnabled(False)

        form_layout = qtw.QFormLayout()
        form_layout.setContentsMargins(0, 0, 0, 0)
        form_layout.addRow("Thermocouple:", self.__thermocouples)
        form_layout.addRow("Use cold-junction compensation:", self.__cjc)

        layout = qtw.QVBoxLayout()
        layout.addLayout(form_layout)
        layout.addWidget(self.__cjc_not_found)
        self.central_widget.setLayout(layout)

    def __connect(self):
        """Connect relevant signals."""
        self.__thermocouples.currentIndexChanged.connect(
            self.__on_thermocouple_changed
            )

    def __on_thermocouple_changed(self, new_idx):
        """React to a change in the thermocouple selection."""
        _thermocouple = self.__thermocouples.itemData(new_idx)
        self._ctrl.thermocouple = _thermocouple
        self.settings_changed.emit()


class _AUXEditDialog(_EditDialogBase):
    """Dialog for editing the AUX input.

    Display the input range, and instructions on how to change it.
    """

    def __init__(self, controller, *args, **kwargs):
        """Initialise dialog."""
        kwargs['raw_quantity'] = "AUX"
        kwargs['help_file'] = Path(resources_path(
            'hardware/schematics/viperLEED_HW_v8 - jumpers.pdf'
            ))
        super().__init__(controller, *args, **kwargs)
