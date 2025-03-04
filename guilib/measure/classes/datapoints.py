"""Module datapoints of viperleed.guilib.measure.classes.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-11-12
Author: Michele Riva
Author: Florian Doerr

Defines the QuantityInfo and DataPoints classes.
"""
import csv
import re
from copy import deepcopy
from collections.abc import MutableSequence, Sequence
from collections import defaultdict
import enum

# ViPErLEED modules
from viperleed.guilib.measure.classes.abc import QMetaABC
from viperleed.guilib.measure.classes.abc import QObjectWithError
from viperleed.guilib.measure.hardwarebase import ViPErLEEDErrorEnum


_ALIASES = {
    'Energy': ('nominal_energy', 'energy',),
    'Measured_Energy': ('hv', 'measured_energy',),
    'Times': ('times', 'measurement_t',),
    'I0': ('i0',),
    'I_Sample': ('isample', 'i_sample',),
    'Temperature': ('temperature',),
    'Aux': ('aux',),
    'Cold_Junction': ('cold_junction',),
    'Timestamps': ('timestamps',),
    }


NAN = float('nan')


class DataErrors(ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""

    INVALID_MEASUREMENT = (400,
                           'The returned data dictionary contained a key '
                           'that was not specified in the DataPoints class.')
    UNKNOWN_QUANTITIES = (401,
                         'Unknown quantity/quantities {} will be ignored')
    NO_DATA_FOR_CONTROLLER = (402,
                              'Controller at {} did not return any data. '
                              'Consider increasing energy_step_duration.')


class QuantityInfo(enum.Enum):
    """Measurement quantities with unit, scaling, type, name, etc.

    New measurement quantities have to be added here.
    """

    # Info:   units, scale, dtype, label, axis, common_label, tooltip
    IMAGES = ('Number', None, str, 'Images', None, None, "")
    ENERGY = ('eV', 'lin', float, 'Energy', 'x', None,
              "The nominal value of the primary electron energy")
    HV = ('eV', 'lin', float, 'Measured_Energy', 'y', 'Voltage',
          "<nobr>The actual value of the primary electron energy"
          "</nobr> measured on the LEED optics at high voltage")
    TIMES = ('s', 'lin', float, 'Times', 'x', None, "")
    I0 = ('µA', 'lin', float, 'I0', 'y', 'Current',
          "<nobr>The total electron current emitted by the "
          "electron gun,</nobr> measured on the LEED optics")
    ISAMPLE = ('µA', 'lin', float, 'I_Sample', 'y', 'Current',
               "<nobr>The total electron current emitted by the electron "
               "gun,</nobr> measured by biasing the sample to +33 V via the "
               "'I_target' BNC connector. This is an alternative to "
               "I<sub>0</sub> in case your LEED optics does not provide an "
               "I<sub>0</sub> output. LEED-IV videos should not be acquired "
               "at the same time to avoid electric-field-induced distortions")
    TEMPERATURE = ('°C', 'lin', float, 'Temperature', 'y', 'Temperature',
                   "")
    AUX = ('mV', 'lin', float, 'Aux', 'y', 'Aux', "")
    COLD_JUNCTION = ('°C', 'lin', float, 'Cold_Junction', 'y', 'Temperature',
                     "Reference temperature measured internally in the "
                     "ViPErLEED unit to convert the measured thermocouple "
                     "voltage to a temperature")
    TIMESTAMPS = ('s', None, str, 'Timestamp', None, None, "")
    UNKNOWN = ('', None, str, '??', None, None, "")

    @classmethod
    def from_label(cls, label):
        """Return the QuantityInfo member associated with label.

        Parameters
        ----------
        label : str
            The label to search for. It may be the .name of the
            attribute or one of the known aliases.

        Returns
        -------
        member : QuantityInfo
            The member of QuantityInfo associated to label.

        Raises
        ------
        TypeError
            If label is not a str
        ValueError
            If label does not match any of the known QuantityInfo(s).
        """
        if not isinstance(label, str):
            raise TypeError(f'Unexpected type {type(label).__name__} for '
                            'QuantityInfo.from_label. Expected str.')
        try:
            return getattr(cls, label.upper())
        except AttributeError:
            pass
        for quantity, aliases in _ALIASES.items():
            if label.lower() in aliases:
                label = quantity
                break
        try:
            return cls.get_labels()[label]
        except KeyError as err:
            raise ValueError(f'{cls.__name__}: {label!r} unknown') from err

    @classmethod
    def get_axis_labels(cls, axis):
        """Return all members of QuantityInfo with .axis == axis.

        Parameters
        ----------
        axis : str
            The axis to be searched for. Can only be one of
            'x' or 'y'.

        Returns
        -------
        members : list
            Each element is the label of the QuantityInfo(s)
            whose .axis is equal to axis.
        """
        return [q.label for q in cls if q.axis == axis]

    @classmethod
    def get_labels(cls):
        """Return a dict {label: enum} of all members of QuantityInfo."""
        return {l.label: l for l in cls}

    @property
    def axis(self):
        """Return the default axis on which self is plotted.

        Typically, only TIMES and ENERGY are used as
        abscissae ('x'), while all the other quantities
        are represented as ordinates ('y').

        Returns
        -------
        axis : {'x', 'y'}
            The default plot axis for this quantity.
        """
        return self.value[4]

    @property
    def common_label(self):                                                     # TODO: rename generic_label? denomination?
        """Return the generic name of self (e.g., "Current", "Voltage")."""
        return self.value[5]

    @property
    def dtype(self):
        """Return the data type of self as a callable (e.g., float)."""
        return self.value[2]

    @property
    def description(self):
        """Return a descriptive text for this quantity."""
        return self.value[6]

    @property
    def label(self):
        """Return the unique label of self as a str (e.g., "Energy")."""
        return self.value[3]

    @property
    def plot_scale(self):
        """Return the default plotting scale ('lin', 'log') for self."""
        return self.value[1]

    @property
    def units(self):
        """Return the units of measure for self as a string."""
        return self.value[0]


# _EXCEPTIONAL contains quantities that are treated differently
# while saving data. Either processed separately, or not saved
_EXCEPTIONAL = (QuantityInfo.IMAGES, QuantityInfo.ENERGY,
                QuantityInfo.TIMESTAMPS, QuantityInfo.UNKNOWN)


# too-many-instance-attributes
class DataPoints(QObjectWithError, MutableSequence, metaclass=QMetaABC):
    """Data storage class."""

    def __init__(self, *args, primary_controller=None, time_resolved=None,
                 continuous=None, parent=None):
        """Initialise data class.

        Parameters
        ----------
        *args : Sequence of dict, optional
            The contents of this DataPoints
        primary_controller : ControllerABC, or str, optional
            The primary controller that is used for the measurement
            in this DataPoints object. It must be set (also later)
            before any operation can be performed.
        time_resolved : bool, optional
            True if this DataPoints object represents a time-resolved
            series. Can also be set ONLY ONCE later by using the
            .time_resolved attribute. It must be set before any
            operation can be done on this object. The value passed
            may be ignored depending on the <continue> argument.
        continuous : bool, optional
            True if this is a time-resolved data series originating
            from a "continuous" measurement, where data arrives at
            maximum possible speed. If True, it overrides the value
            of time_resolved (to True). It can be set ONCE later
            using the .continuous attribute. It must be set before
            any operation can be done on this object if it also
            .is_time_resolved.
        parent : QtCore.QObject
            The parent of this DataPoints object

        Returns
        -------
        None.
        """
        super().__init__()
        self.setParent(parent)
        for element in args:
            self.__check_data(element)
        self.__list = list(args)
        self.__time_resolved = True if continuous else time_resolved
        self.__continuous = continuous
        self.primary_controller = primary_controller

        # Here some counters (set from outside):
        # .nr_steps_done
        #       number of steps that have been completed, i.e.,
        #       excluding steps that already have data but new
        #       measurements can still arrive.
        # .nr_steps_total
        #       Total number of steps to be done. This value is
        #       immutable for one DataPoints instance. Once all
        #       data has been acquired, len(self) == .nr_steps_total
        #       Currently this is only used for the colors of
        #       plotted curves in separate_steps mode.
        #       TODO: is there a better way to do this?
        self.nr_steps_done = 0
        self.nr_steps_total = 1

    @property
    def nominal_energies(self):
        """Return energies as a list."""
        return [d[QuantityInfo.ENERGY] for d in self]

    @property
    def cameras(self):
        """Return cameras as a tuple."""
        if QuantityInfo.IMAGES not in self[0]:
            return tuple()
        return tuple(self[0][QuantityInfo.IMAGES].keys())

    @property
    def controllers(self):
        """Return controllers as a tuple."""
        if not self:
            return tuple()
        return tuple(self[0][QuantityInfo.TIMES].keys())

    @property
    def has_data(self):
        """Check if there is already data in the class."""
        if not self:
            return False
        if QuantityInfo.TIMESTAMPS not in self[0]:  #  Loaded from file
            return any(self[0][QuantityInfo.TIMES].values())
        return any(self[0][QuantityInfo.TIMESTAMPS].values())

    @property
    def continuous(self):
        """Return whether the data comes from a continuous measurement.

        This information should only be used internally. It will
        not be correct for data that has been read from file.

        Returns
        -------
        continuous : bool or None
            True if self contains data originating from
            a continuous time-resolved measurement. None
            if this information is not available.
        """
        return False if not self.__time_resolved else self.__continuous

    @continuous.setter
    def continuous(self, continuous):
        """Set whether data came from a continuous, t-resolved measurement."""
        if self.__continuous is not None:
            raise RuntimeError(
                f"Cannot set {self.__class__.__name__}"
                ".continuous more than once"
                )
        if self.has_data:
            raise RuntimeError(
                f"Cannot set {self.__class__.__name__}"
                ".continuous after there is already data"
                )
        # pylint: disable=compare-to-zero
        # Complains about the "is False", but we need to check for
        # "False" and not False-y. None means "undefined"
        if self.__time_resolved is False and continuous:
            raise ValueError(
                f"{self.__class__.__name__} cannot be "
                ".continuous and not .time_resolved"
                )
        self.__continuous = bool(continuous)
        if continuous:
            self.__time_resolved = True

    @property
    def is_time_resolved(self):
        """Check if the contained data is time-resolved."""
        return self.time_resolved

    @property
    def time_resolved(self):
        """Return whether data is time-resolved."""
        return True if self.__continuous else self.__time_resolved

    @time_resolved.setter
    def time_resolved(self, resolved):
        """Set whether data is time-resolved if not present."""
        if self.__time_resolved is not None:
            raise RuntimeError(
                f"Cannot set {self.__class__.__name__}"
                ".time_resolved more than once"
                )
        if self.has_data:
            raise RuntimeError(
                f"Cannot set {self.__class__.__name__}"
                ".time_resolved after there is already data"
                )
        if self.__continuous and not resolved:
            raise ValueError(
                f"{self.__class__.__name__} cannot be "
                ".continuous and not .time_resolved"
                )
        self.__time_resolved = bool(resolved)
        if not resolved:
            self.__continuous = False

    def __deepcopy__(self, memo):
        """Return a deep copy of self."""
        cls = self.__class__
        # The primary controller object is passed on because it is used
        # for indexing and its times are used as a reference for other
        # controllers in time-resolved measurements.
        kwargs = {'primary_controller' : self.primary_controller,
                  'time_resolved' : self.__time_resolved,
                  'continuous' : self.__continuous,
                  'parent' : self.parent(),}
        result = cls(*deepcopy(self.__list, memo), **kwargs)
        result.nr_steps_done = self.nr_steps_done
        result.nr_steps_total = self.nr_steps_total
        return result

    def __str__(self):
        """Return a string representation of self."""
        return str(self.__list)

    def __repr__(self):
        """Return a string representation of self."""
        return f"DataPoints({self.__list})"

    def __getitem__(self, index):
        """Return element(s) at index or slice."""
        return self.__list[index]

    def __setitem__(self, index, value):
        """Set element(s) at index or slice."""
        if isinstance(value, Sequence):
            for element in value:
                self.__check_data(element)
        else:
            self.__check_data(value)
        self.__list[index] = value

    def __delitem__(self, index):
        """Delete element(s) at index or slice."""
        del self.__list[index]

    def __len__(self):
        """Return the length of self."""
        return len(self.__list)

    def insert(self, index, value):
        """Insert one element at index."""
        self.__check_data(value)
        self.__list.insert(index, value)

    def add_data(self, new_data, controller):
        """Add new data to the currently active data point.

        When not self.continuous, new_data must contain also
        a QuantityInfo.TIMESTAMPS key, with as many values as
        there are in the other quantities.

        Parameters
        ----------
        new_data : dict
            keys : QuantityInfo
                The quantity for which data should be added.
            values : Sequence of float
                The data to be added for this quantity.
        controller : ControllerABC
            The controller that performed the measurement.

        Emits
        -----
        DataErrors.INVALID_MEASUREMENT
            If new_data contains unexpected quantities.

        Raises
        ------
        ValueError
            If controller is not one of self.controllers. This means
            that add_data cannot be called before new_data_point().
        """
        if controller not in self.controllers:
            raise ValueError(
                f"Controller {controller} is not in {self.__class__.__name__}"
                )
        for quantity, values in new_data.items():
            if quantity not in self[-1]:
                self.emit_error(DataErrors.INVALID_MEASUREMENT)
                continue
            # In continuous-time-resolved mode, store only the first
            # timestamp (it is strictly needed only for the primary)
            if (quantity is QuantityInfo.TIMESTAMPS
                and self.continuous
                    and self[-1][quantity][controller]):
                continue
            self[-1][quantity][controller].extend(values)

    def add_image(self, camera):
        """Add the name of the last image of camera."""
        # TODO: camera name added to file name in a 'path'
        # format. Unclear if it isn't just simpler to use
        # the information in the header to decide from
        # which folder images should be opened (in ImageJ)
        self[-1][QuantityInfo.IMAGES][camera].append(
            f"{camera.name_clean}/{camera.process_info.filename}"
            )

    def calculate_times(self, complain=True):  # too-complex
        """Calculate times for the last data point."""
        err = ""
        if not self.has_data:
            err += "Cannot calculate_times without data. "
        if not self.primary_controller:
            err += "Cannot calculate_times without a primary controller. "
        if self and not QuantityInfo.TIMESTAMPS in self[0]:
            err += "Should not calculate_times for data read from file. "

        if err and complain:
            raise RuntimeError(err)
        if err:
            return

        start = self[0][QuantityInfo.TIMESTAMPS][self.primary_controller][0]
        for ctrl, timestamps in self[-1][QuantityInfo.TIMESTAMPS].items():
            if not ctrl.measures():
                # No times needed for a non-measuring controller
                continue
            if not timestamps:
                # May happen in time-resolved when the energy step
                # duration is so small that a controller did not yet
                # return any measurement.
                if complain:
                    self.emit_error(DataErrors.NO_DATA_FOR_CONTROLLER,
                                    ctrl.address)
                continue

            if self.continuous:
                self.__calculate_times_continuous(ctrl, timestamps, start)
                continue
            self.__calculate_times_non_continuous(ctrl, timestamps, start)

    def get_energy_resolved_data(self, *quantities):
        """Return energy resolved data for each controller.

        This method can only be called on non-time-resolved
        data.

        Parameters
        ----------
        *quantities : QuantityInfo
            The measured quantities to be returned. Only those
            quantities that have been measured are returned.

        Returns
        -------
        extracted : dict
            Contains all requested data.
            keys : ControllerABC or str
                The controller that measured data. keys are
                strings if data was read from a saved file.
            values : dict
                keys are QuantityInfo objects (from quantities),
                values are lists of the measurements for the
                quantities (one value per each energy). Only those
                requested data that were measured by this controller
                are returned.
        energies : list
            Electron energies (in electronvolts).

        Raises
        ------
        RuntimeError
            If this method is called on a time-resolved data series.
        """
        # Keep only the quantities that were measured
        quantities = [q for q in quantities if q in self[0]]
        if not quantities:
            return {}, []

        self.__check_attributes()
        if self.is_time_resolved:
            raise RuntimeError(
                f"{self.__class__.__name__}: cannot return "
                "energy-resolved data from a time-resolved series."
                )

        # Prepare the structure of the dictionary to be returned:
        #   {controller:
        #        {quantity: list of measurements (one per energy)}
        #   }
        # including only those controllers that measured (at least
        # one of) the requested quantities.
        extracted = {controller: defaultdict(list)
                     for quantity in quantities
                     for controller in self[0][quantity]}

        # Now fill in the dictionary, looping through the energies.
        # There should only be one measurement for each energy.
        for data_point in self:
            for quantity in quantities:
                for controller, measurements in data_point[quantity].items():
                    if not measurements:
                        # No data for this controller (yet).
                        continue
                    extracted[controller][quantity].append(measurements[0])

        return extracted, self.nominal_energies

    def get_time_resolved_data(self, *quantities, separate_steps=True):
        """Get time resolved data for each controller.

        Separate_steps cannot be False if absolute_times is False.
        If separate_steps is True and absolute_times is False,
        then the steps are justified.

        Parameters
        ----------
        *quantities : QuantityInfo
            The measured quantities to be returned. Only those
            quantities that have been measured are returned.
        separate_steps : bool, optional
            If True, return one list per nominal energy for each
            controller. Additionally, each energy step has its time
            axis starting at zero seconds when the energy was changed.
            If False, return only a single list for each controller,
            containing data from all energies. In this case, times
            are all relative to the moment the first energy was set.
            Default is True.

        Returns
        -------
        extracted : dict
            Contains all requested data. Times are always included.
            keys : ControllerABC or str
                The controller that measured data. keys are
                strings if data was read from a saved file.
            values : dict
                keys are QuantityInfo objects (from quantities),
                values are lists of the measurements for the
                quantities. Only those requested data that were
                measured by this controller are returned. Values
                are 1D lists of floats if separate_steps is False,
                otherwise lists of lists, one per energy step.
        energies : list
            One element per each energy step, i.e.,
            len(energies) == len(self). Notice that energies
            may be longer than the data extracted if the last
            step underway is not finished yet.
        """
        # Keep only the quantities that were measured
        quantities = [q for q in quantities if q in self[0]]
        if not quantities:
            return {}, []

        self.__check_attributes()

        # Prepare the structure of the dictionary to be returned:
        #       {controller: {quantity: measurements}, }
        # including only those controllers that measured (at least
        # one of) the requested quantities. measurements may be
        # a single, "flat" list (if not separate_steps) or a list
        # of lists (one per energy step).
        extracted = {controller: defaultdict(list)
                     for quantity in quantities
                     for controller in self[0][quantity]}

        # Loop through each energy step and fill in the data
        q_time = QuantityInfo.TIMES
        for data_point in self[:self.nr_steps_done]:
            # First, fill in all the quantities requested
            for quantity in quantities:
                for ctrl, measurements in data_point[quantity].items():
                    if not separate_steps:
                        extracted[ctrl][quantity].extend(measurements)
                    else:
                        extracted[ctrl][quantity].append(measurements)

            # Then, take care of the times
            try:
                start_time = data_point[q_time][self.primary_controller][0]
            except IndexError:
                start_time = 0
            for ctrl in extracted:
                ctrl_times = data_point[q_time][ctrl]
                if not separate_steps:
                    extracted[ctrl][q_time].extend(ctrl_times)
                else:
                    ctrl_times = [t - start_time for t in ctrl_times]
                    extracted[ctrl][q_time].append(ctrl_times)

        return extracted, self.nominal_energies

    def new_data_point(self, energy, controllers, cameras):
        """Append a blank data point at a new energy.

        Parameters
        ----------
        energy : float
            Energy that is currently set.
        controllers : list of controllers
            All controllers used.
        cameras : list of cameras
            All cameras used.

        Raises
        -------
        RuntimeError
            If this method is called before .primary_controller,
            .time_resolved, and (for time-resolved only) .continuous
            were set
        ValueError
            If the new_data_point created does not fit with
            already present data (i.e., if it contains different
            controllers and different cameras than already there).
        """
        if not self.primary_controller:
            raise RuntimeError(
                f"{self.__class__.__name__}: cannot add "
                "a new data point before .primary_controller "
                "has been set."
                )
        self.__check_attributes()
        if self:
            if set(controllers) != set(self.controllers):
                raise ValueError(
                    f"{self.__class__.__name__}: inconsistent "
                    "controllers passed to new_data_point. Must "
                    "always contain the same controllers."
                    )
            if set(cameras) != set(self.cameras):
                raise ValueError(
                    f"{self.__class__.__name__}: inconsistent "
                    "cameras passed to new_data_point. Must "
                    "always contain the same cameras."
                    )
            controllers = self.controllers
            cameras = self.cameras

        quantities = (q for q in QuantityInfo
                      if any(c.measures(q) for c in controllers))
        # Create a new data point, keeping only those controllers
        # that measured a certain quantity. We cannot simply use
        # a defaultdict(list) for the inner dictionary because
        # we need to maintain the same order of controllers.
        data_point = {q: {c: [] for c in controllers if c.measures(q)}
                      for q in quantities}
        data_point[QuantityInfo.TIMES] = {c: [] for c in controllers}
        data_point[QuantityInfo.TIMESTAMPS] = {c: [] for c in controllers}
        data_point[QuantityInfo.ENERGY] = energy
        if cameras:
            data_point[QuantityInfo.IMAGES] = {c: [] for c in cameras}
        self.append(data_point)

    def read_csv(self, csv_name):
        """Read data from csv file.

        Parameters
        ----------
        csv_name : str or Path
            Full path (incl. file name) to the file
            from which data is to be read.

        Raises
        -------
        RuntimeError
            If the file to read does not contain a valid energy
            column, or it contains no valid data.
        """
        with open(csv_name, 'r', encoding='UTF8', newline='') as csv_file:
            self.read_lines(csv_file, source=csv_name)

    def read_lines(self, lines, source=''):   # too-many-locals
        """Read data from an iterable returning lines of data.

        Parameters
        ----------
        lines : iterable
            Each element is a string of comma-separated values.
            The first line is expected to contain the "header".
        source : string, optional
            The source from which  data is read. Used exclusively
            for error-reporting purposes.

        Raises
        ------
        RuntimeError
            If the lines to read do not contain a valid energy
            column, or they contain no valid data.
        """
        reader = csv.DictReader(lines, delimiter=',')
        try:
            first_row = next(reader)
        except StopIteration:
            raise RuntimeError(f'{source} contains no data.') from None

        # Prepare a dictionary mapping of the form
        #    {column_label: (QuantityInfo, device)}
        # that we will later use to populate the data points.
        (column_map,
         egy_col,
         invalid) = self.__make_column_header_map(first_row)

        if invalid:
            self.emit_error(DataErrors.UNKNOWN_QUANTITIES, ", ".join(invalid))

        _egy = QuantityInfo.ENERGY
        # Prepare a dummy data point dictionary
        # {energy: value, other_quantity: {device: [values]}, ...}
        # that will be appended to self each time the energy
        # changes, and then populated.
        empty_data_point = {q: defaultdict(list)
                            for q, _ in column_map.values()}
        empty_data_point[_egy] = None

        # Start filling up self
        self.__list = []
        self.append(deepcopy(empty_data_point))
        for row in (first_row, *reader):
            # Handle energy before other quantities
            curr_energy = _egy.dtype(row[egy_col])
            last_energy = self[-1][_egy]
            if last_energy is None:
                # Energy not set yet for this data point.
                # Effectively this happens only for the
                # first data point. All others are handled
                # by the next condition.
                self[-1][QuantityInfo.ENERGY] = curr_energy
            elif abs(curr_energy - last_energy) > 1e-5:
                # Energy changed --> add a new entry
                self.append(deepcopy(empty_data_point))
                self[-1][_egy] = curr_energy

            for column, value in row.items():
                quantity, device = column_map.get(column, (None, None))
                if not quantity or quantity is _egy:
                    # Unknown or energy --> skip
                    continue
                value = quantity.dtype(value)
                self[-1][quantity][device].append(value)

        # Finally define a primary controller
        self.primary_controller = self.controllers[0]
        self.nr_steps_done = len(self)

    def save_data(self, csv_name):
        """Save data to file.

        Parameters
        ----------
        csv_name : str or Path
            Full path and file name to
            which the data will be saved.

        Returns
        -------
        None.
        """
        # TODO: left-pad columns to make it look nice in a text editor.
        #       CANNOT BE DONE WITH csv, have to write manually.
        with open(csv_name, 'w', encoding='UTF8', newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            writer.writerow(self.__make_header_line())

            # Now write each energy step consecutively.
            # To write to file, we un-rag possibly ragged lists
            # of data (each controller/camera may have returned a
            # different number of values/images) by 'padding' with
            # not-a-number/'--' toward the end of the energy step.
            # Each block will be as long as the longest list of data.
            for data_point in self[:self.nr_steps_done]:
                max_length = self.__block_size(data_point)
                data = []
                if self.cameras:
                    for images in data_point[QuantityInfo.IMAGES].values():
                        extra_length = max_length - len(images)
                        data.append(images + ['--']*extra_length)
                data.append([data_point[QuantityInfo.ENERGY]]*max_length)
                for quantity, ctrl_dict in data_point.items():
                    if quantity in _EXCEPTIONAL:
                        # images and energies, already processed
                        continue
                    for ctrl, measurement in ctrl_dict.items():
                        if not ctrl.measures():
                            continue
                        extra_length = max_length - len(measurement)
                        data.append(measurement + [NAN]*extra_length)
                for line in zip(*data):
                    writer.writerow(line)

    def __block_size(self, block):
        """Return the longest of the data entries in a block."""
        max_ctrl = max(
            len(times)
            for times in block[QuantityInfo.TIMES].values()
            )
        max_ctrl = max(1, max_ctrl)  # at least one entry, even if NaN
        if not self.cameras:
            return max_ctrl

        max_cam = max(
            len(images)
            for images in block[QuantityInfo.IMAGES].values()
            )
        return max(max_ctrl, max_cam)

    def __calculate_times_continuous(self, ctrl, timestamps, start):
        """Calculate time axis for the last step."""
        # In a continuous, time-resolved measurement we treat the
        # primary and secondary controllers differently:
        # - primary: is stopped at each energy step, and has a
        #            different timestamp that has to be used
        # - secondary: returns measurements at measurement_interval
        #            intervals. Its timestamp makes any sense only
        #            for the first step, and is the same value since.
        #            Otherwise, we just take the last time of the
        #            previous step and go from there.
        is_first_step = len(self) == 1
        if ctrl == self.primary_controller or is_first_step:
            first_time = timestamps[0] - start
            first_time += ctrl.initial_delay / 1000
        else:
            first_time = (self[-2][QuantityInfo.TIMES][ctrl][-1] +
                          ctrl.measurement_interval / 1000)

        quantity = ctrl.measured_quantities[0]
        n_measurements = len(self[-1][quantity][ctrl])
        self[-1][QuantityInfo.TIMES][ctrl] = [
            first_time + ctrl.measurement_interval * i / 1000
            for i in range(n_measurements)
            ]

    def __calculate_times_non_continuous(self, ctrl, timestamps, start):
        """Calculate time axis for the last step."""
        # This may be an energy-resolved measurement (which has
        # only one timestamp and one measurement) or a triggered
        # time-resolved one (which has as many timestamps as there
        # are measurements)
        self[-1][QuantityInfo.TIMES][ctrl] = times = []
        for timestamp in timestamps:
            timestamp -= start
            timestamp += ctrl.initial_delay / 1000
            times.append(timestamp)

    def __check_attributes(self):
        """Check that .time_resolved (and .continuous) were set."""
        if self.__time_resolved is None:
            raise RuntimeError(
                f"{self.__class__.__name__}.time_resolved "
                "should be set before data can be added"
                )
        if self.__time_resolved and self.__continuous is None:
            raise RuntimeError(
                f"{self.__class__.__name__}.continuous "
                "should be set before data can be added"
                )

    def __check_data(self, value):
        """Check if an element is of acceptable type."""
        if not isinstance(value, dict):
            raise TypeError(
                f"{self.__class__.__name__}: invalid data "
                f"type {type(value).__name__!r}. Expected 'dict'"
                )

    @staticmethod
    def __make_column_header_map(headers):
        """Return a dict of QuantityInfo from text headers.

        Used for reading in data from file.

        Parameters
        ----------
        headers : iterable
            Sequence of strings that will be used for
            constructing the map.

        Returns
        -------
        mapping : dict
            keys: items in headers that could be converted
            to QuantityInfo objects.
            values: (QuantityInfo, str), where str is whatever
            else was found in header that did not match and
            acceptable QuantityInfo.
        egy_col : str
            The column label for the only column containing
            energies
        skipped : list
            Entries in headers that could not be parsed/converted.

        Raises
        -------
        RuntimeError
            If the headers do not contain a valid energy column,
            or they contain no valid data.
        """
        # All column headers, except for the "energy" one are
        # expected to be of the form "quantity(measuring_device_id)".
        extractor = re.compile(
            r'\s*(?P<quantity>.*)\s*\(\s*(?P<ctrl>.*)\s*\)\s*'
            )

        # Keep track of non-matching or unknown
        # quantities for reporting to the user.
        skipped = []

        column_map = {}
        for label in headers:
            if any(label.lower() == e for e in _ALIASES['Energy']):
                column_map[label] = (QuantityInfo.ENERGY, None)
                continue
            extracted = extractor.match(label)
            if not extracted:     # Invalid format
                skipped.append(label)
                continue
            try:
                info = QuantityInfo.from_label(extracted.group("quantity"))
            except ValueError:    # Unknown quantity
                skipped.append(label)
                continue
            column_map[label] = (info, extracted.group("ctrl"))

        # Make sure there is only one energy column, at
        # least one times column, and report any problems
        # encountered with other column headers
        _egy = QuantityInfo.ENERGY
        egy_col = [col
                   for col, (q, _) in column_map.items()
                   if q is _egy]
        if not egy_col:
            raise RuntimeError(
                f'Could not find an "{_egy.label}" column. '
                'File may not be not a ViPErLEED measurement.'
                )
        if len(egy_col) > 1:
            raise RuntimeError(
                f'Multiple columns labeled "{_egy.label}". '
                'File may not be not a ViPErLEED measurement.'
                )
        egy_col = egy_col[0]

        if not any(q is QuantityInfo.TIMES
                   for q, _ in column_map.values()):
            raise RuntimeError('File is a ViPErLEED measurement, '
                               'but it contains no data.')

        return column_map, egy_col, skipped

    def __make_header_line(self):
        """Return a header line that can be written to file.

        Used for saving data to file.

        Returns
        -------
        header : list
            One entry per 'column'. Default order is:
            <all camera names> <energy>
            <quantities for each controller>
            <times for each controller>
        """
        _img = QuantityInfo.IMAGES.label
        header = [f"{_img}({camera.name})" for camera in self.cameras]
        header.append(QuantityInfo.ENERGY.label)
        for quantity, controllers in self[0].items():
            if quantity in _EXCEPTIONAL:
                # cameras and energy, already processed
                continue
            for ctrl in controllers:
                if not ctrl.measures():
                    continue
                header.append(f"{quantity.label}({ctrl.name})")
        return header
