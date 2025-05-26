"""Module datapoints of viperleed.gui.measure.

Defines the QuantityInfo and DataPoints classes.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-11-12'
__license__ = 'GPLv3+'

import csv
from collections import defaultdict
from collections.abc import MutableSequence
from collections.abc import Sequence
from copy import deepcopy
import enum
import re

from PyQt5 import QtCore as qtc

from viperleed.gui.measure.hardwarebase import QMetaABC
from viperleed.gui.measure.hardwarebase import ViPErLEEDErrorEnum
from viperleed.gui.measure.hardwarebase import emit_error


_ALIASES = {
    'Energy': ('nominal_energy', 'energy',),
    'Measured_Energy': ('hv', 'measured_energy',),
    'Times': ('times', 'measurement_t',),
    'I0': ('i0',),
    'I_Sample': ('isample', 'i_sample',),
    'Temperature': ('temperature',),
    'Cold_Junction': ('cold_junction',),
    }


NAN = float('nan')


class DataErrors(ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""
    INVALID_MEASUREMENT = (400,
                           "The returned data dictionary contained a section "
                           "that was not specified in the DataPoints class.")
    UNKOWN_QUANTITIES = (401,
                         "Unkown quantite(s) {} will be ignored")
    NO_DATA_FOR_CONTROLLER = (402,
                              "Controller at {} did not return any data. "
                              "Consider increasing energy_step_duration.")


class QuantityInfo(enum.Enum):
    """Measurement types with unit, scaling, type and name.

    New measurement types have to be added here.
    """
    IMAGES = ('Number', None, str, 'Images', None, None)
    ENERGY = ('eV', 'lin', float, 'Energy', 'x', None)
    HV = ('eV', 'lin', float, 'Measured_Energy', 'y', 'Voltage')
    TIMES = ('s', 'lin', float, 'Times', 'x', None)
    I0 = ('uA', 'lin', float, 'I0', 'y', 'Current')
    ISAMPLE = ('uA', 'lin', float, 'I_Sample', 'y', 'Current')
    TEMPERATURE = ('°C', 'lin', float, 'Temperature', 'y', 'Temperature')
    COLD_JUNCTION = ('°C', 'lin', float, 'Cold_Junction', 'y', 'Temperature')

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


# too-many-instance-attributes
class DataPoints(qtc.QObject, MutableSequence, metaclass=QMetaABC):
    """Data storage class."""
    # Is emitted when an error occurs
    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, *args, time_resolved=None):
        """Initialise data class."""
        super().__init__()
        self.__list = list(args)
        self.__delimiter = ','
        self.__primary_first_time = None
        self.__exceptional_keys = (QuantityInfo.IMAGES, QuantityInfo.ENERGY)
        self.__time_resolved = time_resolved

        # Keep track of whether times were calculated already for
        # the most-recent data point, to prevent multiple calls
        # to functions that would screw up stuff.
        self.__times_calculated = True

        # primary_controller needs to be given to this class
        # before starting measurements.
        self.primary_controller = None

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
        return bool(self and any(self[0][QuantityInfo.TIMES].values()))

    @property
    def is_time_resolved(self):
        """Check if the contained data is time-resolved."""
        if self.time_resolved is not None:
            return self.time_resolved
        if self.has_data:
            # Only if there are times saved already it is possible to
            # decide if the measurement is a time- or energy-resolved
            # measurement.
            if any(len(t) > 1 for t in self[0][QuantityInfo.TIMES].values()):
                return True
        return False

    @property
    def time_resolved(self):
        """Return whether data is time-resolved.

        Returns
        -------
        time_resolved : bool or None
            None if it has not been set yet.
        """
        return self.__time_resolved

    @time_resolved.setter
    def time_resolved(self, resolved):
        """Set whether data is time-resolved if not present."""
        if self.__time_resolved is not None:
            raise RuntimeError(
                "time_resolved can only be set once, "
                "but it was attempted to set it again."
                )
        self.__time_resolved = bool(resolved)

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
            If new_data contains unexpected quantitites.
        """
        for quantity, values in new_data.items():
            if quantity not in self[-1]:
                emit_error(self, DataErrors.INVALID_MEASUREMENT)
            else:
                self[-1][quantity][controller].extend(values)

        if not self[-1][QuantityInfo.TIMES][controller]:
            time = (controller.serial.time_stamp
                    + controller.time_to_trigger/1000)
            self[-1][QuantityInfo.TIMES][controller].append(time)

    def add_image_names(self, name):                                            # TODO: this will not work correctly for a time-resolved with multiple "acquisitions" per energy (i.e., not continuous) -- right now we take only one image, which is inadequate.
        """Add image names to currently active data point."""
        for camera in self[-1][QuantityInfo.IMAGES]:
            # TODO: camera name added to file name in a 'path'
            # format. Unclear if it isn't just simpler to use
            # the information in the header to decide from
            # which folder images should be opened (in ImageJ)
            self[-1][QuantityInfo.IMAGES][camera] = f"{camera.name}/{name}"

    def calculate_times(self, continuous=False):
        """Calculate times for the last data point."""
        if self.__times_calculated:
            raise RuntimeError(
                "Times were already calculated for this data "
                "point. Use .recalculate_last_step_times()."
                )

        self.__times_calculated = True
        time = QuantityInfo.TIMES
        if self.__primary_first_time is None:
            self.__primary_first_time = (
                    self[0][time][self.primary_controller][0]
                    )
        for ctrl, ctrl_times in self[-1][time].items():
            if not ctrl_times:
                # May happen in time-resolved when the energy step              # TODO: this happens especially while replotting. May be solved by moving measurement (and primary) to its own thread.
                # duration is so small that (typically the primary)
                # controller did not yet return any measurement.
                emit_error(self, DataErrors.NO_DATA_FOR_CONTROLLER,
                           ctrl.port_name)
                continue
            first_time = ctrl_times[0]
            if (not continuous
                or len(self) == 1  # First energy step
                    or ctrl == self.primary_controller):
                first_time -= self.__primary_first_time
                first_time += ctrl.initial_delay / 1000
            else:
                first_time = (self[-2][time][ctrl][-1] +
                              ctrl.measurement_interval / 1000)
            if ctrl.measured_quantities:
                quantity = ctrl.measured_quantities[0]
                n_measurements = len(self[-1][quantity][ctrl])
                self[-1][time][ctrl] = [
                    first_time + ctrl.measurement_interval * i / 1000           # TODO: incorrect for non-continuous. Should probably use the timestamps from triggering.
                    for i in range(n_measurements)
                    ]
            else:
                self[-1][time][ctrl][0] = first_time

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
                    if len(measurements) > 1:
                        raise RuntimeError(
                            f"{self.__class__.__name__}: cannot return "
                            "energy-resolved data from a time-resolved "
                            "series."
                            )
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
            may be longer then the data extracted if the last
            step underway is not finished yet.
        """
        # Keep only the quantities that were measured
        quantities = [q for q in quantities if q in self[0]]
        if not quantities:
            return {}, []

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
            If this method is called before .primary_controller
            was set.
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
        data_point[QuantityInfo.ENERGY] = energy
        if cameras:
            data_point[QuantityInfo.IMAGES] = {c: "" for c in cameras}
        self.append(data_point)
        self.__times_calculated = False

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

    def read_lines(self, lines, source=''):
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
        reader = csv.DictReader(lines, delimiter=self.__delimiter)
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
            emit_error(self, DataErrors.UNKOWN_QUANTITIES,
                       ", ".join(invalid))

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

    def recalculate_last_step_times(self):
        """Recalculate the whole time axis for the last energy step.

        This is necessary because the secondaries keep returning
        measurements before they are turned off.
        """
        if not self.has_data:
            return
        for ctrl, ctrl_times in self[-1][QuantityInfo.TIMES].items():
            if not ctrl_times:
                # Controller has not yet finished measuring this step
                continue
            first_time = ctrl_times[0]
            if ctrl.measured_quantities:
                quantity = ctrl.measured_quantities[0]
                n_measurements = len(self[-1][quantity][ctrl])
                self[-1][QuantityInfo.TIMES][ctrl] = [
                    first_time + ctrl.measurement_interval * i / 1000
                    for i in range(n_measurements)
                    ]
            else:
                self[-1][QuantityInfo.TIMES][ctrl][0] = first_time

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
            writer = csv.writer(csv_file, delimiter=self.__delimiter)
            writer.writerow(self.__make_header_line())

            # Now write each energy step consecutively.
            # To write to file, we un-rag possibily ragged lists
            # of data (each controller may have retured a different
            # number of values) by 'padding' with not-a-number toward
            # the end of the energy step. Each block will be as long
            # as the longest list of data.
            for data_point in self[:self.nr_steps_done]:
                max_length = max(
                    len(times)
                    for times in data_point[QuantityInfo.TIMES].values()
                    )
                data = [[data_point[QuantityInfo.IMAGES][camera]]*max_length
                        for camera in self.cameras]
                data.append([data_point[QuantityInfo.ENERGY]]*max_length)
                for quantity, ctrl_dict in data_point.items():
                    if quantity in self.__exceptional_keys:
                        # images and energies, already processed
                        continue
                    for measurement in ctrl_dict.values():
                        extra_length = max_length - len(measurement)
                        data.append(measurement + [NAN]*extra_length)
                for line in zip(*data):
                    writer.writerow(line)

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
            or they contains no valid data.
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
            if quantity in self.__exceptional_keys:
                # cameras and energy, already processed
                continue
            for controller in controllers:
                header.append(
                    f"{quantity.label}({controller.name})"
                    )
        return header
