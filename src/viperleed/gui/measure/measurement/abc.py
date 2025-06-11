"""Module abc of viperleed.gui.measure.measurement.

This module contains the definition of the MeasurementABC class
which gives commands to the controller classes and its associated
ViPErLEEDErrorEnum class MeasurementErrors.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-19'
__license__ = 'GPLv3+'

from abc import abstractmethod
from collections.abc import Sequence
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import time

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.classes.abc import QObjectWithSettingsABC
from viperleed.gui.measure.classes.datapoints import DataPoints
from viperleed.gui.measure.classes.datapoints import QuantityInfo
from viperleed.gui.measure.classes.settings import NoSettingsError
from viperleed.gui.measure.classes.settings import NotASequenceError
from viperleed.gui.measure.classes.settings import SystemSettings
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.controller.abc import MeasureControllerABC
from viperleed.gui.measure.dialogs.settingsdialog import (
    SettingsHandler,
    SettingsSectionColumnInfo,
    SettingsTag,
    )
from viperleed.gui.measure.measurement import _meassettings as _settings
from viperleed.gui.measure.widgets.spinboxes import CoercingDoubleSpinBox


_INVOKE = qtc.QMetaObject.invokeMethod
_QUEUED = qtc.Qt.QueuedConnection
_UNIQUE = qtc.Qt.UniqueConnection


class MeasurementErrors(base.ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""

    RUNTIME_ERROR = (300, 'Runtime error. Info: {}')
    WRONG_CONTROLLER_CLASS = (
        301,
        'The secondary controller at address {!r} is not a subclass '
        'of MeasureControllerABC. All secondary controllers need '
        'to be a subclass of MeasureControllerABC.'
        )
    MISSING_PRIMARY = (
        302,
        'No primary controller available to set the beam energy. Check '
        'both the measurement and the primary controller configuration '
        'files.'
        )
    MISSING_CAMERA = (
        303,
        'No camera available for the measurement. Check both the '
        'measurement and the camera configuration files.'
        )
    TOO_MUCH_DATA = (
        304,
        'The devices {} returned more data than expected.'
        )


class MeasurementException(Exception):
    """Base exception of instances of MeasurementABC."""


class MeasurementDevicesConnectedError(MeasurementException):
    """Raised if measurement devices are connected while setting settings."""


class MeasurementIsRunningError(MeasurementException):
    """The same MeasurementABC was started again before it finished."""


class MeasurementReusedError(MeasurementException):
    """Raised if a measurement instance is used more than once."""


# Progression:
# Entry: .start()
# * ._continue_preparation()
# * ._check_preparation_finished()
# * BEGIN MEASURING by auto-call to ._begin_next_energy_step()
# * END OF CURRENT STEP: _on_controller_data_ready/_on_camera_busy_changed,
#   which call ._ready_for_next_measurement(). This decides to call
# * ._begin_next_energy_step() or ._prepare_finalization()

# too-many-instance-attributes
class MeasurementABC(QObjectWithSettingsABC):                                   # TODO: doc about inner workings
    """Generic measurement class."""

    # All cameras and controllers have been disconnected
    devices_disconnected = qtc.pyqtSignal()

    # Whole measurement is over
    finished = qtc.pyqtSignal()

    # New data can be processed
    new_data_available = qtc.pyqtSignal(dict)

    # Preparation is finished
    prepared = qtc.pyqtSignal()

    # Abort current tasks on all devices
    _request_stop_devices = qtc.pyqtSignal()                                    # TODO: Could use QMetaObject.invokeMethod

    # _preparation_started: emitted in ._begin_preparation right
    # before the first energy is set. Carries pairs of energies and
    # settle times, passed on to .primary_controller.set_energy()
    _preparation_started = qtc.pyqtSignal(tuple)                                # TODO: Could use QMetaObject.invokeMethod

    # _preparation_continued: emitted after all controllers have
    # completed the first segment of their preparation
    _preparation_continued = qtc.pyqtSignal()                                   # TODO: Could use QMetaObject.invokeMethod

    _mandatory_settings = (
        ('devices', 'primary_controller'),
        ('energies', 'start_energy'),
        ('energies', 'delta_energy'),
        ('energies', 'end_energy'),
        ('energies', 'step_profile'),
        )

    # Backwards compatibility fix                                               # TODO: #242
    _settings_synonyms = (
        (('energies', 'start_energy'),
         ('measurement_settings', 'start_energy'),),
        (('energies', 'delta_energy'),
         ('measurement_settings', 'delta_energy'),),
        (('energies', 'end_energy'),
         ('measurement_settings', 'end_energy'),),
        (('energies', 'step_profile'),
         ('measurement_settings', 'step_profile'),),
        )

    def __init__(self, settings):
        """Initialise measurement instance."""
        super().__init__(settings=settings)
        self._other_mandatory_settings = [('measurement_settings',
                                           'measurement_class',
                                           (self.__class__.__name__,))]
        self._current_energy = 0
        self._previous_energy = 0
        self._primary_controller = None
        self._secondary_controllers = []
        self._cameras = []
        self._aborted = False
        self._has_been_used_before = False # Used to stop reuse of object.      # TODO: We may want to modify the measurement to allow this behaviour.
        self._temp_dir = None   # Directory for saving files

        # Keep track of which data of which controller was
        # stored in self.data_points at this energy step
        self._missing_data = {}

        self.threads = []
        self.running = False     # Used for aborting from outside
        self.data_points = DataPoints(parent=self)
        self.data_points.error_occurred.connect(self.error_occurred)

        self.error_occurred.connect(self._store_delayed_error)
        self.error_occurred.connect(self._on_hardware_error, type=_QUEUED)

        self._camera_timer = qtc.QTimer(parent=self)
        self._camera_timer.setSingleShot(True)

        # _force_end_timer ensures that no-matter-what we will
        # wait at most 4.5 sec to clean up the measurement at its
        # end. Using a timer we can handle a loss (or timeout) of
        # the primary controller.
        self._force_end_timer = qtc.QTimer(parent=self)
        self._force_end_timer.setSingleShot(True)
        self._force_end_timer.setInterval(4500)
        self._force_end_timer.timeout.connect(self._cleanup_and_end)

        # We do not use the delayed_errors context here because the             # TODO: perhaps we can also store the sender together with the errors? Then we could sender.emit_error again.
        # device.error_occurred signals are not connected yet and
        # therefore the MeasurementABC collects those errors and
        # reports them with its own errors.
        self.set_settings(self._settings_to_load)

        if self._delayed_errors:
            self._delay_errors_timer.start(20)
        self.error_occurred.disconnect(self._store_delayed_error)

    @property
    def aborted(self):
        """Return whether the measurement was aborted."""
        return self._aborted

    @property
    def cameras(self):
        """Return the cameras used by this class."""
        return self._cameras

    @cameras.setter
    def cameras(self, new_cameras):
        """Set the cameras which should be used and handle signals.

        Parameters
        ----------
        new_cameras : list
            List of camera class objects.
        """
        self._disconnect_cameras()
        self._cameras = new_cameras
        self._connect_cameras()

    @property
    def controllers(self):
        """Return a tuple of all controllers."""
        if not self.primary_controller:
            return tuple()
        return (self.primary_controller, *self.secondary_controllers)

    @property
    def current_energy(self):
        """Return the current energy in electronvolts."""
        return self._current_energy

    @current_energy.setter
    def current_energy(self, new_energy):
        """Set a new value of the current_energy.

        Notice that this DOES NOT actually set the energy
        in the LEED electronics. Call set_leed_energy for
        that purpose.

        Parameters
        ----------
        new_energy : float
            The new current energy
        """
        self._previous_energy = self.current_energy
        self._current_energy = new_energy

    @property
    def current_step_nr(self):
        """Return an incremental number for the current energy step."""
        return len(self.data_points)

    @property
    def devices(self):
        """Return all controllers and cameras."""
        return (*self.controllers, *self.cameras)

    @property
    def hv_settle_time(self):
        """Return the time interval for the settling of energies."""
        if not self.primary_controller:
            raise RuntimeError("Cannot return a voltage-settling time "
                               "when no primary controller was set.")
        return self.primary_controller.hv_settle_time

    @property
    def primary_controller(self):
        """Return the primary controllers used by this class."""
        return self._primary_controller

    @primary_controller.setter
    def primary_controller(self, new_controller):
        """Set the primary controller.

        Parameters
        ----------
        new_controller : ControllerABC
            Controller that sets the energy.
        """
        self._disconnect_primary_controller()
        self._primary_controller = new_controller
        self._connect_primary_controller()

    @property
    def secondary_controllers(self):
        """Return the controllers used by this class."""
        return self._secondary_controllers

    @secondary_controllers.setter
    def secondary_controllers(self, new_controllers):
        """Set secondary controllers and move them to threads.

        Use the set_secondary_controllers method to set secondary
        controllers from another thread via a signal or via
        invokeMethod.

        Parameters
        ----------
        new_controllers : list
            List of MeasureControllerABC objects.
        """
        self.set_secondary_controllers(new_controllers)

    @qtc.pyqtSlot(object)
    def set_secondary_controllers(self, new_controllers):
        """Set the controllers which should be used and handle signals.

        Disconnect old secondary controllers and quit and clear their
        threads. After that, move new secondary controllers to newly
        created threads, start them and connect their signals.

        Parameters
        ----------
        new_controllers : list
            List of MeasureControllerABC objects.
        """
        self._disconnect_secondary_controllers()
        self.stop_threads()
        qtw.qApp.processEvents()
        self.threads.clear()
        for controller in new_controllers:
            thread = qtc.QThread()
            controller.moveToThread(thread)
            self.threads.append(thread)
        for thread in self.threads:
            thread.start(priority=thread.TimeCriticalPriority)
        self._secondary_controllers = new_controllers
        self._connect_secondary_controllers()

    @property
    def start_energy(self):
        """Return the first energy for the energy ramp."""
        if not self.settings:
            return 0.0
        try:
            return self.settings.getfloat('energies', 'start_energy',
                                          fallback=0)
        except (TypeError, ValueError):
            # Not a float
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/start_energy', '')
            return 0.0

    @property
    def step_profile(self):                                                     # TODO: probably move to generator class?
        """Return a list of energies and times for setting the next energy.

        The returned value excludes the very last step, i.e.,
        self.current_energy and the settling time for it.
        A typical call to set_leed_energy would be
            self.set_leed_energy(*self.step_profile,
                                 self.current_energy,
                                 last_settle_time,
                                 ...)

        Returns
        -------
        step_profile : tuple
            Sequence of energies and waiting intervals.
        """
        try:
            profile = self.settings.getsequence('energies', 'step_profile',
                                                fallback=("abrupt",))
        except NotASequenceError:
            profile = self.settings['energies']['step_profile']

        if isinstance(profile, str):
            profile = (profile,)

        # Now we have two acceptable cases:
        # (1) the sequence can be cast to (float, int, float, int...)
        # (2) the first entry is a known profile shape
        try:
            return self._step_profile_from_strings(profile)
        except (ValueError, TypeError):
            pass

        shape, *params = profile
        if shape.lower() == 'abrupt':
            values = tuple()
        elif shape.lower() == 'linear':
            values = self._get_linear_step(*params)
        else:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            f'Unknown profile shape {shape}')
            values = tuple()
        return values

    @abstractmethod
    @qtc.pyqtSlot()
    def abort(self):
        """Interrupt measurement, save partial data, set energy to zero.

        This function needs to be extended in subclasses.
        Implementation needs to reset all variables used in loop
        operation. Then call super().abort().

        Returns
        -------
        None.
        """
        if self._aborted:
            return
        self._aborted = True
        self._camera_timer.stop()
        if self.settings:
            self.settings.set('measurement_info', 'was_aborted', 'True')
        self._force_end_timer.start()
        # _force_end_timer performs a delayed call to _cleanup_and_end
        # which will stop all threads. The quitting of threads must be
        # delayed in case abort is called while secondary controllers
        # are still busy.
        if not self.running:
            return
        self._prepare_finalization()

    def are_runtime_settings_ok(self):
        """Return whether runtime settings are ok.

        This method is used to check if the runtime settings that result
        from the given settings are enough to run a measurement. The
        base implementation returns False if no primary controller
        exists.

        Returns
        -------
        settings_ok : bool
            True if the runtime settings are
            sufficient to start a measurement.

        Emits
        -----
        error_occurred
            If the primary controller is missing.
        """
        if not self.primary_controller:
            self.emit_error(MeasurementErrors.MISSING_PRIMARY)
            return False
        return True

    def are_settings_invalid(self, settings):
        """Check if there are any invalid settings.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            The new settings.

        Returns
        -------
        invalid_settings : list of tuples
            Invalid _mandatory_settings of self as a list of tuples.
            The first entry in each tuple can be either '<section>',
            '<section>/<option>', or
            '<section>/<option> not one of <value1>, <value2>, ...'.
            Further optional entries may be added by subclasses. They
            specify additional information on what is wrong with each
            invalid setting.
        """
        invalid_settings = settings.has_settings(
            *self._mandatory_settings,
            *self._other_mandatory_settings
            )

        # Backwards compatibility fix                                           # TODO: #242
        if not settings.has_section('energies'):                                # TODO: Auto-generate new sections if necessary!
            settings.add_section('energies')
        for new_setting, old_setting in self._settings_synonyms:
            if '/'.join(new_setting) in invalid_settings:
                old_missing = settings.has_settings(old_setting)
                if not old_missing:
                    settings.set(*new_setting, settings.get(*old_setting))
                    settings.remove_option(*old_setting)
                    settings.update_file()
                    invalid_settings.remove('/'.join(new_setting))

        return [(invalid,) for invalid in invalid_settings]

    def disconnect_devices_and_notify(self):
        """Disconnect devices and emit devices_disconnected."""
        for device in self.devices:
            base.safe_connect(device.connection_changed,
                              self._check_if_all_devices_disconnected,
                              type=_UNIQUE | qtc.Qt.QueuedConnection)
        self._disconnect_devices()
        if not any(device.connected for device in self.devices):
            self.devices_disconnected.emit()

    @abstractmethod
    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings.

        This method should be extended in subclasses, i.e., do
        handler = super().get_settings_handler(), and then add
        appropriate sections and/or options to it using the
        handler.add_section, and handler.add_option methods.

        Use the QNoDefaultPushButton from the widgets.buttons
        module in order to prevent any button from being set
        as the default button of the dialog.

        The base-class implementation returns a handler that
        already contains the following settings:
            'energies'/'start_energy'
            'energies'/'delta_energy'
            'energies'/'end_energy'
            'energies'/'step_profile'

        Returns
        -------
        handler : SettingsHandler
            The handler used in a SettingsDialog to display the
            settings of this measurement to users.
        """
        self.check_creating_settings_handler_is_possible()
        handler = SettingsHandler(self.settings, show_path_to_config=True)
        sys_config = SystemSettings()
        settings_path = sys_config.paths['configuration']

        handler.add_section('measurement_info', tags=SettingsTag.REGULAR)
        tip = ('<nobr>This string will be appended </nobr>'
               'to the file name of the data.')
        # Backwards compatibility fix                                           # TODO: #242
        has_option = self.settings.has_option('measurement_info', 'suffix')
        option = 'suffix' if has_option else 'tag'
        handler.add_option('measurement_info', option,
                           handler_widget=qtw.QLineEdit,
                           display_name='File suffix',
                           tooltip=tip,)
        # Backwards compatibility fix                                           # TODO: #242
        has_option = self.settings.has_option('measurement_info', 'comments')
        option = 'comments' if has_option else 'info'
        handler.add_option('measurement_info', option,
                           handler_widget=qtw.QTextEdit,
                           display_name='Comments',)

        handler.add_section('measurement_settings', tags=SettingsTag.REGULAR)
        type_display = qtw.QLabel(type(self).__name__)
        handler.add_static_option(
            'measurement_settings', 'measurement_class',
            type_display, display_name='Measurement type',
            )

        handler.add_section('energies', tags=SettingsTag.REGULAR)
        info = (
            ('start_energy', _settings.START_E_NAME,
             '<nobr>The energy at which the measurement starts.</nobr>'),
            ('delta_energy', _settings.DELTA_E_NAME,
             '<nobr>The energy difference between two measurement '
             'steps.</nobr>'),
            ('end_energy', _settings.END_E_NAME,
             '<nobr>The energy value at which </nobr>'
             'the measurement will finish.'),
            )
        for option_name, display_name, tip in info:
            widget = CoercingDoubleSpinBox(decimals=1, soft_range=(0, 1000),
                                           suffix=' eV')
            handler.add_option(
                'energies', option_name, handler_widget=widget,
                display_name=display_name, tooltip=tip
                )
        delta_energy = handler['energies']['delta_energy']
        delta_energy.handler_widget.soft_minimum = -1000
        delta_energy.handler_widget.setSingleStep(0.5)

        widget = _settings.StepProfileViewer()
        tip = '<nobr>How to move from </nobr>one energy to the next one.'
        handler.add_option('energies', 'step_profile',
                           handler_widget=widget, display_name='Step profile',
                           tooltip=tip)

        second_column = SettingsSectionColumnInfo(position=1)
        device_section = _settings.DeviceEditor(
            self.settings, default_folder=settings_path,
            may_have_cameras=True, column_info=second_column,
            )
        handler.add_complex_section(device_section)

        return handler

    @classmethod
    def is_matching_default_settings(cls, obj_info, config, match_exactly):
        """Determine if the default `config` file is for a measurement.

        Parameters
        ----------
        obj_info : SettingsInfo or None
            The information that should be used to check `config`.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether obj_info should be matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used to sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the items in the tuple is the
            order of their significance. This return value is used
            to determine the best-matching settings files when
            multiple files are found. An empty tuple signifies that
            `config` does not match the requirements.
        """
        return (1,)

    @classmethod
    def is_matching_user_settings(cls, obj_info, config, match_exactly):
        """Determine if the `config` file is for a measurement.

        Parameters
        ----------
        obj_info : SettingsInfo
            The information that should be used to check `config`.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether obj_info should be matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used to sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the items in the tuple is the
            order of their significance. This return value is used
            to determine the best-matching settings files when
            multiple files are found. An empty tuple signifies that
            `config` does not match the requirements.
        """
        return (1,)

    @classmethod
    def is_settings_for_this_class(cls, config):
        """Determine if the `config` file is for this measurement.

        Parameters
        ----------
        config : ConfigParser
            The settings to check.

        Returns
        -------
        is_suitable : bool
            True if the settings file is for this measurement.
        """
        meas_class = config.get('measurement_settings', 'measurement_class',
                                fallback=None)
        return cls.__name__ == meas_class

    def moveToThread(self, thread):  # pylint: disable=invalid-name
        """Move self and primary controller to a new thread."""
        # Notice that this explicit extension is necessary as
        # relying on parent-child relationships is not enough.
        # In fact, child.moveToThread(thread) is not called,
        # as is any overridden version of the method.
        super().moveToThread(thread)
        if self.primary_controller:
            self.primary_controller.moveToThread(thread)

    def set_leed_energy(self, energy, settle_time, *more_steps,
                        trigger_meas=True):
        """Set the electron energy used for LEED.

        In order to achieve quicker settling times for the LEED
        electronics one can do quick steps after one another
        whose effects cancel each other out. E.g., a small energy
        overshoot quickly followed by a correction.

        Parameters
        ----------
        energy : float
            Nominal electron energy in electronvolts.
        settle_time : integer
            Interval in milliseconds that the controller will wait
            before deeming the LEED optics stable at set energy.
        *more_steps : Number
            If given, it should be an even number of elements.
            Odd elements are energies, even ones settle-time intervals.
            Multiple steps can be executed quickly after each other.
            The last step will be the final energy that is set and
            should ensure stabilization of the electronics.
        trigger_meas : bool, optional
            True if the controllers are supposed to take
            measurements after the energy has been set.
            Default is True.

        Returns
        -------
        None.
        """
        primary = self.primary_controller

        if trigger_meas and not primary.measures():
            # If the primary controller cannot measure we anyway need a
            # timestamp as a reference for all other measurements.
            primary.time_stamp = time.perf_counter()
            dummy_dict = {QuantityInfo.TIMESTAMPS: []}
            dummy_dict[QuantityInfo.TIMESTAMPS].append(primary.time_stamp)
            self.data_points.add_data(dummy_dict, primary)

        if self.current_energy != self._previous_energy:
            primary.set_energy(energy, settle_time, *more_steps,
                               trigger_meas=trigger_meas)
            return

        # When the energy did not change, we can avoid setting
        # it again and go straight to measuring, if desired
        if not trigger_meas:
            return

        try:
            primary.measure_now()
        except AttributeError:
            pass  # Controller cannot measure

        primary.about_to_trigger.emit()

    @qtc.pyqtSlot(object)
    def set_settings(self, new_settings):
        """Change settings of the measurement.

        Settings are loaded only if they are valid. Otherwise
        the previous settings stay in effect. If the settings
        have been accepted, controller and camera objects as
        specified in the settings will be instantiated, told
        what they will be measuring, moved to their respective
        properties and connected to all necessary signals.

        In order to prevent attempting to connect to already
        connected devices, it is only possible to set settings
        after all devices have been disconnected. This means
        disconnect_devices_and_notify has to be executed and the
        devices_disconnected signal has to be received before
        calling set_settings on a measurement object that was
        used before.

        Parameters
        ----------
        new_settings : dict or ConfigParser or str or Path or ViPErLEEDSettings
            Configuration of the measurement.

        Returns
        -------
        settings_valid : bool
            True if the new settings given were accepted.

        Raises
        ------
        MeasurementDevicesConnectedError
            If any devices are still connected while attempting to
            set new settings. Call disconnect_devices_and_notify
            and wait for emission of devices_disconnected to
            prevent this error.
        TypeError
            If new_settings is neither a dict, ConfigParser, string
            or path and if an element of the mandatory_settings is
            None or has a length greater than 3.

        Emits
        -----
        QObjectSettingsErrors.MISSING_SETTINGS
            If new_settings is missing.
        QObjectSettingsErrors.INVALID_SETTINGS
            If any element of the new_settings does not fit the
            mandatory_settings.
        """
        if any(device.connected for device in self.devices):
            raise MeasurementDevicesConnectedError(
                'Setting settings is only allowed after all devices '
                'have been disconnected. Make sure to disconnect them '
                'before attempting to set new settings. See '
                'help(measurement.set_settings).'
                )

        self._aborted = False  # Set False in case of abort through settings

        # Notice that we clear data even if the settings are not accepted.
        # When the primary controller is set in the primary_controller
        # setter it is handed to self.data_points. This happens in
        # _make_primary_ctrl below.
        self.data_points = DataPoints(parent=self)
        self.data_points.error_occurred.connect(self.error_occurred)

        if not super().set_settings(new_settings):
            return False

        if not self._make_primary_ctrl():
            # Something went wrong (already reported in _make_primary).
            self.primary_controller = None                                      # TODO: probably good to clean up secondaries and cameras!
            return False
        self.data_points.primary_controller = self.primary_controller
        self._make_secondary_ctrls()
        self._make_cameras()
        self._make_tmp_directory_tree()

        self._missing_data = {c: 1
                              for c in self.controllers
                              if c.measures()}
        for camera in self.cameras:
            self._missing_data[camera] = 1
        return True

    @qtc.pyqtSlot()
    def start(self):
        """Check runtime settings and start measurement."""
        self._aborted = False
        if not self.are_runtime_settings_ok():
            return
        if self._has_been_used_before:
            raise MeasurementReusedError(
                'Cannot reuse measurements to run multiple cycles. '
                'Create a new instance instead.'
                )
        if self.running:
            raise MeasurementIsRunningError(
                'Attempted to start a new measurement on an instance '
                'that was already performing a measurement.'
                )

        self._has_been_used_before = True
        self.settings.set('measurement_info', 'was_aborted', 'False')
        self.settings.set('measurement_info', 'started',
                          time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime()))
        self.running = True
        self._begin_preparation()

    def stop_threads(self):
        """Quit all threads immediately."""
        for thread in self.threads:
            thread.quit()

    @abstractmethod
    def _begin_next_energy_step(self):
        """Set next energy (if required) and measure.

        This method does not increment self.current_energy, as
        this is done as part of the end of the previous "step".

        This method should be extended in subclasses. The base
        implementation only prepares a new 'data point'.

        Reimplementations should set the new energy by calling
        self.set_leed_energy(..., self.current_energy, ...).

        Returns
        -------
        None.
        """
        self._missing_data = dict.fromkeys(self._missing_data.keys(), 1)
        self.data_points.new_data_point(self.current_energy, self.controllers,
                                        self.cameras)

        image_name = (f"{{__count__:0>{self._n_digits}}}_"
                      f"{self.current_energy:.1f}eV.tiff")
        for camera in self.cameras:
            camera.process_info.filename = image_name
            camera.process_info.energy = self.current_energy

    def _begin_preparation(self):
        """Start preparation for measurements.

        Prepare the controllers and cameras for a measurement.

        This method triggers the first part (i.e., everything that is
        done before the starting energy is set), then automatically
        moves on to the second part. Finally, the actual measurement
        loop is triggered. Users can detect when the whole preparation
        is over by connecting to the .prepared() signal.

        Emits
        -----
        _preparation_started
            Starts the measurement preparation and carries
            a tuple of energies and times with it.
        """
        self.current_energy = self.start_energy

        primary = self.primary_controller
        about_to_trigger = primary.about_to_trigger
        for ctrl in self.controllers:
            # Make sure no controller can be triggered for
            # measurement during the whole preparation
            try:
                base.safe_disconnect(about_to_trigger, ctrl.measure_now)
            except AttributeError:
                pass
            # Measurements that may still be delivered during
            # preparation should be discarded:
            base.safe_disconnect(ctrl.data_ready,
                                 self._on_controller_data_ready)

            # When controllers will turn "not busy" at the end of
            # this first preparation segment, go to second segment
            ctrl.busy_changed.connect(self._continue_preparation)

        # Disconnect the camera.busy_changed signal here, and reconnect
        # it later in ._check_preparation_finished(). This prevents
        # early calls to _on_camera_busy_changed.
        for camera in self.cameras:
            base.safe_disconnect(camera.busy_changed,
                                 self._on_camera_busy_changed)

        # Notice that we have to handle only controllers:
        # ._preparation_started is already connected to camera.start
        self._preparation_started.emit(
            (self.start_energy, primary.long_settle_time)
            )

    @qtc.pyqtSlot(bool)
    def _check_if_all_devices_disconnected(self, device_connected):
        """Check and notify when all devices are disconnected.

        Emits
        -----
        devices_disconnected
            After all devices have been disconnected.
        """
        device = self.sender()
        if device_connected or device.connected:
            return
        base.safe_disconnect(self.sender().connection_changed,
                             self._check_if_all_devices_disconnected)
        if any(device.connected for device in self.devices):
            return
        self.devices_disconnected.emit()

    @qtc.pyqtSlot(bool)
    def _check_preparation_finished(self, _):
        """Check if measurement preparation is done.

        Whenever a device is done with its preparation
        this function will be called. The busy attribute is
        used to check if all devices are done with their
        preparation. If this is True, the busy signal going to
        this function will be disconnected and the first
        measurement will immediately be started.

        Emits
        -------
        prepared
            After all devices are done with the preparation,
            and right before entering the measurement loop.
        """
        if any(device.busy for device in self.devices):
            return
        for device in self.devices:
            base.safe_disconnect(device.busy_changed,
                                 self._check_preparation_finished)

        # Finally, reconnect all devices to
        # be ready to actually take measurements
        self._connect_primary_controller()
        self._connect_secondary_controllers()
        self._connect_cameras()

        self.prepared.emit()           # Signal that we're done.
        self.current_energy = self.start_energy
        self._begin_next_energy_step()  # And start the measurement loop

    @qtc.pyqtSlot(bool)
    @qtc.pyqtSlot()
    def _cleanup_and_end(self, *__args):
        """Conclude measurement and clean up.

        Quit threads and emit finished().

        Emits
        -----
        finished
            Emitted right before this method returns. All
            operations are done when this happens.
        """
        self._disconnect_devices()
        self._force_end_timer.stop()
        self.stop_threads()
        self.running = False
        self.finished.emit()

    @qtc.pyqtSlot(bool)
    def _continue_preparation(self, _):
        """Continue preparation for measurements.

        Do nothing till all controllers are done with te first part
        of the preparation, then move on to the second segment (i.e.,
        everything that is done after the starting energy is set).
        Finally, move on to trigger the beginning of the measurement
        loop.

        Emits
        -----
        _preparation_continued
            Starts the second part of the measurement
            preparation.
        """
        if any(controller.busy for controller in self.controllers):
            return

        # Use the controller.busy_changed to move from this segment
        # of the preparation to the exit point of the preparation
        # that will later start the measurement loop.
        for ctrl in self.controllers:
            base.safe_disconnect(ctrl.busy_changed,
                                 self._continue_preparation)
        # The camera.busy_changed signals are connected only now, rather
        # than during _begin_preparation. This prevents early calls
        # to the _check_preparation_finished method, should cameras
        # be ready early.
        for device in self.devices:
            device.busy_changed.connect(self._check_preparation_finished,
                                        type=_UNIQUE)
        self._preparation_continued.emit()

    def _connect_cameras(self):
        """Connect necessary camera signals."""
        # It is not necessary to call .connect_(), as it is called
        # already in the settings setter of the camera. Not calling
        # it again saves some time, as loading camera settings can be
        # slow. We are thus supposing that the camera objects already
        # have the right settings.
        for camera in self.cameras:
            base.safe_connect(camera.busy_changed,
                              self._on_camera_busy_changed,
                              type=_UNIQUE | qtc.Qt.QueuedConnection)
            base.safe_connect(self._request_stop_devices, camera.stop,
                              type=_UNIQUE)
            base.safe_connect(self._camera_timer.timeout, camera.trigger_now,
                              type=_UNIQUE)
            base.safe_connect(self._preparation_started, camera.start,
                              type=_UNIQUE)
            base.safe_connect(camera.image_saved, self._on_image_saved,
                              type=_UNIQUE)

    def _connect_controller(self, ctrl):
        """Connect serial and signals of a controller."""
        with ctrl.errors_delayed():
            ctrl.connect_()
        base.safe_connect(ctrl.data_ready, self._on_controller_data_ready,
                          type=_UNIQUE)
        base.safe_connect(self._request_stop_devices, ctrl.stop, type=_UNIQUE)
        base.safe_connect(self._preparation_started, ctrl.begin_preparation,
                          type=_UNIQUE)
        base.safe_connect(self._preparation_continued,
                          ctrl.continue_preparation, type=_UNIQUE)

    def _connect_primary_controller(self):
        """Connect signals of the primary controller."""
        if not self.primary_controller:
            return
        self._connect_controller(self.primary_controller)

    def _connect_secondary_controllers(self):
        """Connect necessary controller signals."""
        about_to_trigger = self.primary_controller.about_to_trigger
        for ctrl in self.secondary_controllers:
            base.safe_connect(about_to_trigger, ctrl.measure_now, type=_UNIQUE)
            self._connect_controller(ctrl)

    def _disconnect_cameras(self):
        """Disconnect necessary camera signals."""
        disconnect = base.safe_disconnect
        for camera in self.cameras:
            disconnect(camera.busy_changed, self._on_camera_busy_changed)
            disconnect(self._camera_timer.timeout, camera.trigger_now)
            disconnect(self._preparation_started, camera.start)
            disconnect(camera.stopped, self._finalize)
            disconnect(camera.image_saved, self._on_image_saved)
            _INVOKE(camera, 'disconnect_')

    def _disconnect_controller(self, ctrl):
        """Disconnect a generic controller."""
        disconnect = base.safe_disconnect
        _INVOKE(ctrl, 'disconnect_')
        disconnect(ctrl.data_ready, self._on_controller_data_ready)
        disconnect(self._request_stop_devices, ctrl.stop)
        disconnect(self._preparation_started, ctrl.begin_preparation)
        disconnect(self._preparation_continued, ctrl.continue_preparation)
        busy_slots = (self._continue_preparation,
                      self._check_preparation_finished,
                      self._cleanup_and_end, self._finalize)
        for slot in busy_slots:
            disconnect(ctrl.busy_changed, slot)

    def _disconnect_devices(self):
        """Disconnect all connected devices but keep them alive."""
        self._disconnect_secondary_controllers()
        self._disconnect_cameras()
        self._disconnect_primary_controller()

    def _disconnect_primary_controller(self):
        """Disconnect serial and signals of the primary controller."""
        primary = self.primary_controller
        if primary is None:
            return
        about_to_trigger = primary.about_to_trigger
        for ctrl in self.secondary_controllers:
            base.safe_disconnect(about_to_trigger, ctrl.measure_now)
        self._disconnect_controller(primary)

    def _disconnect_secondary_controllers(self):
        """Disconnect serials and signals of the secondary controllers."""
        about_to_trigger = getattr(self.primary_controller,
                                   'about_to_trigger',
                                   None)  # No primary yet
        for ctrl in self.secondary_controllers:
            if about_to_trigger:
                base.safe_disconnect(about_to_trigger, ctrl.measure_now)
            self._disconnect_controller(ctrl)

    @qtc.pyqtSlot(bool)
    @qtc.pyqtSlot()
    def _finalize(self, *_):
        """Finish the measurement: save data, then set energy to zero.

        This is the slot connected to the busy_changed signal of
        all controllers, and to the .stopped signal of all cameras,
        after ._is_finished() returns True, or after .abort()ing the
        measurement. Signals are disconnected again once all devices
        turn "not busy".

        Returns
        -------
        None.
        """
        if any(device.busy for device in self.devices):
            return
        if any(camera.is_running for camera in self.cameras):
            # Not busy, but some images are still being processed
            return
        self._missing_data.clear()
        try:
            self._save_data()
        finally:
            # Disconnect all devices and their signals
            self._disconnect_devices()

            # Keep only the primary controller connected, so we can
            # set the LEED energy to zero (and detect it has been set)
            primary = self.primary_controller
            primary.connect_()
            primary.busy_changed.connect(self._cleanup_and_end, type=_UNIQUE)
            self.current_energy = 0
            self.set_leed_energy(self.current_energy, 50, trigger_meas=False)

    def _get_device_settings(self, configname):
        """Return a ViPErLEEDSettings for a device given its path.

        This method reads the correct settings, distinguishing
        the two cases in which self.settings was read from an
        actual file or from within an archive.

        Returns
        -------
        device_settings : ViPErLEEDSettings
            The loaded settings

        Raises
        ------
        RuntimeError
            In case any error occurs while loading the settings.
        """
        configname = configname.strip()
        device_cfg = ViPErLEEDSettings()

        # Treat now different cases:
        # (1) configname does not begin with './' --> read file
        if not configname.startswith('./'):
            try:
                device_cfg = device_cfg.from_settings(configname)
            except (ValueError, NoSettingsError):
                raise RuntimeError() from None
            return device_cfg

        # (2) configname begins with './' --> has to be replaced
        configname = configname[2:]

        # (2.1) self.settings was read from an archive
        # --> attempt reading config from the same archive
        if self.settings.base_dir.endswith('.zip'):
            with ZipFile(self.settings.base_dir, 'r') as arch:
                try:
                    cfg_lines = arch.read(configname).decode()
                except KeyError:
                    raise RuntimeError(
                        f"No config file '{configname}' in "
                        f"archive {self.settings.base_dir}"
                        ) from None
            device_cfg.read_string(cfg_lines)
            device_cfg.base_dir = self.settings.base_dir
            return device_cfg

        # (2.2) self.settings was read from a folder
        # --> attempt reading config from the same folder
        configname = Path(self.settings.base_dir) / configname
        try:
            device_cfg = device_cfg.from_settings(configname)
        except (ValueError, NoSettingsError):
            raise RuntimeError(f"No config file '{configname}' in "
                               f"folder {self.settings.base_dir}") from None
        return device_cfg

    def _get_linear_step(self, *params):
        """Return energies and times for a simple linear step."""
        if len(params) != 2:
            # Too many/too few
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            'Too many/few parameters for linear profile. '
                            f'Expected 2, found {len(params)}')
            return tuple()

        try:
            n_steps, tot_time = (int(p) for p in params)
        except (TypeError, ValueError):
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            'Could not convert to integer the '
                            'parameters for linear profile')
            return tuple()

        if n_steps <= 0 or tot_time < 0:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'energies/step_profile',
                            'Linear-step parameters should be '
                            'positive integers')
            return tuple()

        delta_t = tot_time // n_steps
        delta_e = self.current_energy - self._previous_energy
        if not delta_t or abs(delta_e) < 1e-4:
            return tuple()

        # Make a line of the form f(t) = t/tot_time, with t == 0
        # the time at which self.current_energy is set, and choose
        # an (almost) equally-spaced time grid, with the first
        # interval compensating for non integer-divisibility
        times = [-tot_time,
                 *(-(i-1)*delta_t for i in range(n_steps, 0, -1))]
        intervals = (tj - ti for ti, tj in zip(times, times[1:]))

        # The best way to approximate a function with a piecewise
        # constant signal is to have values fk equal to the average
        # of f over the k-th interval. For our line, this means
        # f[k] = (t[k] + t[k+1]) / (2*tot_time)
        slope = delta_e / (2 * tot_time)
        energies = (slope*(ti + tj) + self.current_energy
                    for ti, tj in zip(times, times[1:]))

        # Finally interleave energies and times
        return tuple(v for tup in zip(energies, intervals) for v in tup)

    @abstractmethod
    def _is_finished(self):
        """Check if the full measurement cycle is done.

        This method must be overridden in subclasses. It should
        check if the measurement cycle is done via the settings.

        Returns
        -------
        finished : bool
            True when the measurement is finished.
        """
        return True

    def _make_camera(self, camera_settings):
        """Instantiate camera class object.

        Take camera settings and generate a camera object from it.

        Parameters
        ----------
        camera_settings : dict, ConfigParser, str, path
            Settings used for the instantiated camera.
            Has to contain the name of the camera class
            to be instantiated.

        Returns
        -------
        instance : camera class object
            The camera class used for the connected camera.

        Raises
        ------
        RuntimeError
            If anything goes wrong with instantiation.

        Emits
        -----
        QObjectSettingsErrors.INVALID_SETTINGS
            If the camera_settings given could not be found,
            or if failed to make a camera instance.
        """
        try:
            config = self._get_device_settings(camera_settings)
        except RuntimeError as err:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'devices/path to camera configuration', err)
            raise

        invalid = config.has_settings(('camera_settings', 'class_name'))
        if invalid:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'camera_settings/class_name',
                            f'No class_name in {config.last_file}')
            raise RuntimeError

        camera_cls_name = config.get('camera_settings', 'class_name')
        try:
            camera_class = base.class_from_name('camera', camera_cls_name)
        except ValueError:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'camera_settings/class_name', '')
            raise RuntimeError from None

        camera = camera_class(settings=config)
        if camera.mode != 'triggered':
            # Force mode to be triggered
            camera.settings.set("camera_settings", "mode", "triggered")
            # self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,           # TODO: Should we make this a non-critical warning?
                            # 'camera_settings/mode', 'Camera {camera.name}: '
                            # 'Cannot measure in live mode.')
        return camera

    def _make_cameras(self):
        """Make cameras from self.settings."""
        try:
            cam_settings = self.settings.getsequence('devices', 'cameras',
                                                     fallback=())
        except NotASequenceError:                                               # TODO: probably report the error
            cam_settings = tuple()

        if isinstance(cam_settings, str):
            # Only one camera, without brackets, or missing comma
            cam_settings = (cam_settings,)

        cameras = []
        for settings in cam_settings:
            try:
                cam = self._make_camera(settings)
            except RuntimeError:
                continue
            cam.error_occurred.connect(self._store_delayed_error)
            cam.error_occurred.connect(self._on_hardware_error)
            cam.process_info.count = 1  # Image counter, start at 1
            cameras.append(cam)
        self.cameras = cameras

    def _make_controller(self, configname, measurements, is_primary=False):
        """Instantiate controller class object.

        Take controller settings and generate a controller object
        from it. This method is used exclusively in .set_settings
        to generate controllers.

        Parameters
        ----------
        configname : str or path
            Path to the file of the settings to be given to the
            returned controller. In addition to the mandatory
            settings required by the controller, the file pointed
            by this path must contain the name of the controller
            class to be instantiated.
        measurements : Sequence of str
            The quantities measured by this controller. Expected to
            be compatible with the argument of .set_measurements in
            MeasureControllerABC.
        is_primary : boolean
            True if the controller is the primary controller.

        Returns
        -------
        controller : ControllerABC
            The controller that is going to handle the hardware.
            May be a subclass of ControllerABC.

        Raises
        ------
        RuntimeError
            If anything goes wrong with the creation of the controller

        Emits
        -----
        MeasurementErrors.MISSING_CLASS_NAME
            If the controller class name is missing.
        QObjectSettingsErrors.INVALID_SETTINGS
            If the controller could not be instantiated
            from the given name.
        """
        try:
            config = self._get_device_settings(configname)
        except RuntimeError as err:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'devices/path to controller configuration', err)
            raise

        invalid = config.has_settings(('controller', 'controller_class'))
        if invalid:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'controller/controller_class',
                            f'No controller_class in {config.last_file}')
            raise RuntimeError

        # For now, check that the address is in the settings.                   # TODO: add getting address from device list
        # Later on, this check will only happen if the unique name
        # of the controller in the settings file that was passed
        # is not found in the device list.
        # Backwards compatibility fix for port_name:                            # TODO: #242
        address = 'address'
        invalid = config.has_settings(('controller', 'address'))
        if invalid:
            address = 'port_name'
            invalid = config.has_settings(('controller', 'port_name'))
        if invalid:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'controller/address',
                            f'No address in {config.last_file}')
            raise RuntimeError
        address = config.get('controller', address)
        if not address:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'controller/address',
                            f'No address in {config.last_file}')
            raise RuntimeError

        cls_name = config['controller']['controller_class']
        try:
            cls = base.class_from_name('controller', cls_name)
        except ValueError:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'controller/controller_class',
                            f'Unknown class {cls_name} in {config.last_file}')
            raise RuntimeError from None

        if isinstance(measurements, str):
            # Only one quantity, with no brackets or no comma
            measurements = (measurements,)

        if not isinstance(measurements, Sequence):
            section = ('primary_controller' if is_primary
                       else 'secondary_controllers')
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            f"devices/{section}",
                            "Measured quantities is not a sequence "
                            f"in {self.settings.last_file}")

        controller = cls(settings=config, address=address,
                         sets_energy=is_primary)
        # Connect error signal. The connection with _store_delayed_error
        # has any impact only during initialization of self, but it
        # does not hurt to leave it connected. _on_hardware_error
        # causes abortion of the measurement.
        controller.error_occurred.connect(self._store_delayed_error)
        controller.error_occurred.connect(self._on_hardware_error)
        controller.set_measurements(measurements)
        return controller

    def _make_primary_ctrl(self):
        """Make primary controller from self.settings."""
        try:
            info = self.settings.getsequence('devices', 'primary_controller')
        except NotASequenceError:
            info = ()

        if len(info) != 2:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'devices/primary_controller', '')
            return False

        try:
            ctrl = self._make_controller(*info, is_primary=True)
        except RuntimeError:
            # Something went wrong with instantiating, and is
            # already reported by emitting in _make_controller.
            return False

        self.primary_controller = ctrl
        return True

    def _make_secondary_ctrls(self):
        """Make secondary controllers from self.settings."""
        # infos should be a sequence with elements of
        # the form (path, (stuff, to, measure, ...))
        try:
            infos = self.settings.getsequence(
                'devices', 'secondary_controllers', fallback=()
                )
        except NotASequenceError:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'devices/secondary_controllers',
                            'Could not be converted to a sequence')
            infos = tuple()

        secondary_controllers = []
        for info in infos:
            if len(info) != 2:
                self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                                'devices/secondary_controllers', '')
                continue
            try:
                ctrl = self._make_controller(*info, is_primary=False)
            except RuntimeError:
                continue
            if not isinstance(ctrl, MeasureControllerABC):
                self.emit_error(MeasurementErrors.WRONG_CONTROLLER_CLASS,
                                ctrl.address)
                continue
            secondary_controllers.append(ctrl)
        # The next one also quits the old threads
        self.secondary_controllers = secondary_controllers

    def _make_tmp_directory_tree(self):
        """Prepare temporary folder tree where data will be saved."""
        sys_config = SystemSettings()
        base_path = sys_config.paths["measurements"]
        if not base_path:                                                       # TODO: emit a SystemSettingsError of some kind
            raise RuntimeError("Invalid path in system config file")

        base_path = Path(base_path).resolve() / "__tmp__"
        base_path.mkdir(parents=True, exist_ok=True)

        self._temp_dir = base_path

        for camera in self.cameras:
            cam_dir = base_path / camera.name_clean
            camera.process_info.base_path = str(cam_dir)
            cam_dir.mkdir(exist_ok=True)

        # Create a blank zip archive (or delete
        # the contents if one is already present)
        with ZipFile(self._temp_dir.parent / '__tmp__.zip', 'w'):
            pass

    @property
    def _n_digits(self):
        """Return the appropriate number of digits for padding image names."""
        # With 4 digits we can store 9999 images with the same prefix           # TODO: nicer implementation: ask n_steps to the generator and use a large value for _n_digits if < 0 (i.e., unknown no. steps)
        # size. More images will have longer file names.
        return 4

    @qtc.pyqtSlot(bool)
    def _on_camera_busy_changed(self, busy):
        """Receive not busy signal from camera.

        Receive the camera.busy_changed signal and check if all cameras
        and controllers are ready if the camera that emitted the
        signal was not busy.

        Returns
        -------
        None.
        """
        if busy:
            # Just triggered
            return

        # Collected all frames, and will save a processed image later
        camera = self.sender()
        self.data_points.add_image(camera)
        camera.process_info.count += 1

        self._missing_data[camera] -= 1
        self._ready_for_next_measurement()

    @qtc.pyqtSlot(dict)
    def _on_controller_data_ready(self, data):
        """Receive measurement data from the controller.

        Append data to the internal dictionary, then check if all the
        controllers are ready for the next measurement. This is the
        slot connected to the .data_ready signal of each controller.

        Parameters
        ----------
        data : dict
            A dictionary containing the measurements. Keys are
            QuantityInfo, values are Sequences (most commonly
            1-element long).

        Returns
        -------
        None.
        """
        controller = self.sender()
        self.data_points.add_data(data, controller)
        if controller.measures():
            self._missing_data[controller] -= 1
        self._ready_for_next_measurement()

    @qtc.pyqtSlot(tuple)
    def _on_hardware_error(self, *_):
        """Abort if a hardware error occurs."""
        self.abort()

    @qtc.pyqtSlot(str)
    def _on_image_saved(self, img_name):
        """Archive the latest image saved."""
        if not self._temp_dir or not self._temp_dir.exists():
            return
        img_name = Path(img_name).resolve()
        arch_name = self._temp_dir.with_suffix('.zip')
        with ZipFile(arch_name, 'a', compression=ZIP_DEFLATED,
                     compresslevel=2) as archive:
            archive.write(img_name, img_name.relative_to(self._temp_dir))

        return  # TODO: remove after testing

        # Remove the image just appended to the archive
        img_name.unlink()

    def _prepare_finalization(self):
        """Prepare for finalization.

        This method is called both when ._is_finished() returns
        True, and while .abort()ing a measurement.

        This method may need to be extended in subclasses.
        Ensure that finalization is prepared properly and that
        self._finalize() gets called. Can be only a call to
        self._finalize() if no preparation is needed.

        Returns
        -------
        None.
        """
        if self.primary_controller:
            # This check is to prevent raising AttributeError
            # when this method is called because of errors in
            # setting up the primary controller.
            about_to_trigger = self.primary_controller.about_to_trigger
            for ctrl in self.secondary_controllers:
                base.safe_disconnect(about_to_trigger, ctrl.measure_now)

        for ctrl in self.controllers:
            base.safe_disconnect(ctrl.busy_changed,
                                 self._continue_preparation)
            base.safe_disconnect(ctrl.busy_changed,
                                 self._check_preparation_finished)
            # Force all controllers to busy, such that ._finalize()
            # is called for all when they turn "not busy" anymore
            ctrl.busy = True
            base.safe_connect(ctrl.busy_changed, self._finalize, type=_UNIQUE)

        for camera in self.cameras:
            base.safe_connect(camera.stopped, self._finalize, type=_UNIQUE)

        self._request_stop_devices.emit()

    def _ready_for_next_measurement(self):                                      # TODO: > _continue_measuring_when_ready?
        """Proceed after all measurements have been received.

        After all measurements have been received, a check if
        the loop is done will be called. This decides whether
        to go to the next energy step, or if we should wrap up.

        Returns
        -------
        None.
        """
        if self.aborted:
            return
        if any(device.busy for device in self.devices):
            return
        if any(miss < 0 for miss in self._missing_data.values()):
            # Check if any device returned more data than expected.
            _too_many = {d.name: v for d, v in self._missing_data.items()
                         if v < 0}
            self.emit_error(MeasurementErrors.TOO_MUCH_DATA, _too_many)
        if any(self._missing_data.values()):
            # Don't go to the next energy if we haven't processed the
            # data from all devices. Notice that this solution prevents
            # a race condition in which one device may have returned
            # data while another one has turned not busy (i.e., all are
            # not busy), but the signal of the other devices was not
            # processed yet.
            return

        self.data_points.calculate_times()
        self.data_points.nr_steps_done += 1
        self.new_data_available.emit(self.data_points[-1].copy())

        if self._is_finished():
            self._prepare_finalization()
        else:
            self._begin_next_energy_step()

    # too-many-locals, too-complex
    def _save_data(self):                                                       # TODO: clean up after testing zip
        """Save data."""
        if not self._temp_dir:
            return

        tmp_path = self._temp_dir
        base_path = tmp_path.parent
        tag = self.settings.get('measurement_info', 'tag', fallback='')
        name_tag = (' ' + tag) if tag else ''
        current_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
        self.settings.set('measurement_info', 'ended', current_time)
        final_path = base_path / (current_time + name_tag)
        move_to_archive = []

        if self.data_points:
            fname = tmp_path / 'measurement.csv'
            self.data_points.save_data(fname)
            move_to_archive.append(fname)

        ctrl_locations = []
        for ctrl in self.controllers:
            fname = "controller_" + ctrl.name_clean + ".ini"
            with open(tmp_path / fname, 'w', encoding='utf-8') as fproxy:
                ctrl.settings.write(fproxy)
            move_to_archive.append(tmp_path / fname)
            ctrl_locations.append("./" + fname)

        cam_locations = []
        for camera in self.cameras:
            fname = "camera_" + camera.name_clean + ".ini"
            with open(tmp_path / fname, 'w', encoding='utf-8') as fproxy:
                camera.settings.write(fproxy)
            move_to_archive.append(tmp_path / fname)
            cam_locations.append("./" + fname)

        if cam_locations:                                                       # TODO: do the same for all controllers (without editing the measured stuff!)
            self.settings.set("devices", "cameras", str(tuple(cam_locations)))

        file_name = tmp_path / "measurement.ini"
        move_to_archive.append(file_name)
        with open(file_name, 'w', encoding='utf-8') as configfile:
            self.settings.write(configfile)

        arch_name = tmp_path.with_suffix('.zip')
        with ZipFile(arch_name, 'a', compression=ZIP_DEFLATED,
                     compresslevel=2) as archive:
            for fname in move_to_archive:
                archive.write(fname, fname.relative_to(tmp_path))

        # Move all files by renaming the archive and temp folder
        try:
            arch_name.rename(arch_name.with_name(final_path.name + '.zip'))
        except OSError as err:
            self.emit_error(MeasurementErrors.RUNTIME_ERROR, err)
        try:
            tmp_path.rename(final_path)
        except OSError as err:
            self.emit_error(MeasurementErrors.RUNTIME_ERROR, err)

        self._temp_dir = None  # Prevents re-saving

        return  # TODO: remove

        for fname in move_to_archive:
            fname.unlink()
        for fname in tuple(final_path.glob('*')):                               # TODO: not sure I want to do this cleanup. Maybe just remove camera folders higher up?
            if fname.is_file():
                fname.unlink()
                continue
            try:
                fname.rmdir()
            except OSError:
                pass

    def _step_profile_from_strings(self, profile):                              # TODO: Warn for .ini files created before 23/05/2025.
        """Return a tuple of energies and times from strings."""
        delta = self.current_energy - self._previous_energy
        if abs(delta) < 1e-4:
            return tuple()

        energies_times = [0]*len(profile)
        for i, fraction in enumerate(profile[::2]):
            # We shift by -1 here in order to display to the user that
            # the 'current_energy' (the energy before the energy step)
            # is equivalent to a step fraction of 0 and the next energy
            # is equivalent to 1. We have to do this as the
            # current_energy is already incremented to the next energy.
            this_delta = (float(fraction) - 1) * delta
            energies_times[2*i] = this_delta + self.current_energy
        for i, time_ in enumerate(profile[1::2]):
            time_ = int(time_)
            if time_ < 0:
                self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                                'energies/step_profile',
                                '\nInfo: Time intervals must be non-negative')
                return tuple()
            energies_times[2*i+1] = time_
        return tuple(energies_times)
