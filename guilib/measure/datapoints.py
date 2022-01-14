"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-11-12
Author: Michele Riva
Author: Florian Doerr

Defines the DataPoints class.
"""
import csv
import re
import copy
from collections.abc import MutableSequence, Sequence
import enum

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.hardwarebase import (
    ViPErLEEDErrorEnum, emit_error, QMetaABC
    )


class DataErrors(ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""
    INVALID_MEASUREMENT = (400,
                           "The returned data dictionary contained a section "
                           "that was not specified in the DataPoints class.")


_ALIASES = {'Energy': ('nominal_energy',),
            'Measured_Energy': ('hv',),
            'I_Sample': ('isample',),
           }


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
        if not isinstance(label, str):
            raise TypeError(f'Unexpected type {type(label).__name__} for '
                            'QuantityInfo.from_label. Expected str.')
        try:
            return getattr(cls, label)
        except AttributeError:
            pass
        for quantity, aliases in _ALIASES.items():
            if label.lower() in aliases:
                label = quantity
                break
        try:
            return cls.get_labels()[label]
        except AttributeError as err:
            raise AttributeError(f'Unknown {label=} for quantity in '
                                 'QuantityInfo')

    @property
    def common_label(self):
         return self.value[5]

    @property
    def axis(self):
        return self.value[4]

    @property
    def label(self):
        return self.value[3]

    @property
    def dtype(self):
        return self.value[2]

    @property
    def plot_scale(self):
        return self.value[1]

    @property
    def units(self):
        return self.value[0]

    @classmethod
    def get_labels(cls):
        return {l.label: l for l in cls}

    @classmethod
    def get_axis_labels(cls, axis):
        return [q.label for q in cls if q.axis == axis]


class DataPoints(qtc.QObject, MutableSequence, metaclass=QMetaABC):
    """Data storage class."""
    # Is emitted when an error occurs
    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, *args):
        """Initialise data class."""
        super().__init__()
        self.__list = list(args)
        # primary_controller needs to be given to this class
        # before starting measurements.
        self.primary_controller = None
        self.controllers = None
        self.delimiter = ','
        self.primary_first_time = None
        self.__exceptional_keys = (QuantityInfo.IMAGES, QuantityInfo.ENERGY )

    @property
    def nominal_energies(self):
        """Return energies."""
        return [d[QuantityInfo.ENERGY] for d in self]

    @property
    def cameras(self):
        if not QuantityInfo.IMAGES in self[0]:
            return tuple()
        return tuple(self[0][QuantityInfo.IMAGES].keys())

    def __str__(self):
        return str(self.__list)

    def __repr__(self):
        return f"DataPoints({self.__list})"

    def __getitem__(self, index):
        return self.__list[index]

    def __setitem__(self, index, value):
        if isinstance(value, Sequence):
            for element in value:
                self.__check_data(element)
        else:
            self.__check_data(value)
        self.__list[index] = value

    def __delitem__(self, index):
        del self.__list[index]

    def __len__(self):
        return len(self.__list)

    def insert(self, index, value):
        self.__check_data(value)
        self.__list.insert(index, value)

    def __check_data(self, value):
        """Check if received data is acceptable."""
        if not isinstance(value, dict):
            raise ValueError(f'{self.__class__.__name__}: received '
                             'unexpected data type')

    def add_data(self, new_data, controller, primary_delay=None):
        """Add received data to currently active data point."""
        for quantity, value in new_data.items():
            if quantity not in self[-1]:
                emit_error(self, DataErrors.INVALID_MEASUREMENT)
            else:
                self[-1][quantity][controller].append(value)

        if not self[-1][QuantityInfo.TIMES][controller]:
            time = controller.serial.time_stamp
            if primary_delay:
                time += primary_delay/1000
            self[-1][QuantityInfo.TIMES][controller].append(time)

    def add_image_names(self, name):
        """Add image names to currently active data point."""
        for camera in self[-1][QuantityInfo.IMAGES]:
            # TODO: may want to add camera name to image name
            self[-1][QuantityInfo.IMAGES][camera] = name

    def calculate_times(self, continuous=False):
        """Calculate times for the last data point."""
        if not self.primary_first_time:
            self.primary_first_time = (
                    self[0][QuantityInfo.TIMES][self.primary_controller][0]
                    )
        for ctrl in self[-1][QuantityInfo.TIMES]:
            if not len(self[-1][QuantityInfo.TIMES][ctrl]):
                continue
            if not continuous:
                self[-1][QuantityInfo.TIMES][ctrl][0] -= self.primary_first_time
                self[-1][QuantityInfo.TIMES][ctrl][0] += ctrl.initial_delay
            else:
                if len(self) > 1 and not ctrl == self.primary_controller:
                    self[-1][QuantityInfo.TIMES][ctrl][0] = self[-2][QuantityInfo.TIMES][ctrl][-1] + ctrl.measurement_interval
            quantity = ctrl.measured_quantities[0]
            length = len(self[-1][quantity][ctrl])
            for i in range(length - 1):
                self[-1][QuantityInfo.TIMES][ctrl].append(
                    self[-1][QuantityInfo.TIMES][ctrl][0] +
                    ctrl.measurement_interval * (i + 1)
                    )

    def get_energy_resolved_data(self, key, include_energies=False):
        """Get energy resolved data for each controller.

        Parameters
        ----------
        key : str
            The key associated with the requested
            measured quantity.
        include_energies : bool, optional
            Whether energy should be returned. Default is False.

        Returns
        -------
        data : list
            Contains measurements.
        self.nominal_energies : list
            Contains set energies.
        """
        data = []
        if key not in self[0]:
            raise ValueError(f'{self.__class__.__name__}: invalid key {key}')
        for measurement in self[0][key]:
            if measurement:
                data.append([])

        for data_point in self:
            i = 0
            for ctrl, measurement in data_point[key].items():
                if len(measurement) == 1:
                    data[i].append(*measurement)
                else:
                    raise ValueError('DataPoints contains time-resolved data '
                                     'but it was attempted to return energy '
                                     'resolved data.')
                i += 1
        if include_energies:
            return data, self.nominal_energies
        return data

    def get_time_resolved_data(self, key, include_energies=False,
                               separate_steps=True, absolute_times=False):
        """Get time resolved data for each controller.

        Separate_steps cannot be False if absolute_times is False.
        If separate_steps is True and absolute_times is False,
        then the steps are justified.

        Parameters
        ----------
        key : str
            The key associated with the requested
            measured quantity.
        include_energies : bool, optional
            Whether energy should be returned. Default is False.
        separate_steps : bool, optional
            If True, return one list per nominal energy for each
            controller. If False, return only a single list for
            each controller, containing data from all energies.
            Default is True.
        absolute_times: bool, optional
            If True, the times returned are 'absolute' (but relative
            to the beginning of the first step), otherwise, each energy
            step starts at time == 0. This argument is ignored if
            include_times is False. Default is False.

        Returns
        -------
        data : list
            Contains measurements and associated
            starting times.
        times : list
            Contains times at which measurements
            were taken.
        nominal_energies : list
            Contains set energies. This is only returned if
            include_energies == True.
        """
        if not separate_steps and not absolute_times:
            raise ValueError('Cannot return flat data '
                             'with relative times.')
        data = []
        times = []
        if key not in self[0]:
            raise ValueError(f'{self.__class__.__name__}: invalid key {key}')

        # Prepare data and times to contain as many empty lists as
        # there are controllers that measured key during the first
        # (and all the other) energy step.
        for controller in self[0][key]:
            if controller:
                data.append([])
                times.append([])

        # Run over each energy step
        # For each energy step we add information from the
        # ctrl_idx-th controller into data[ctrl_idx]. At the
        # end, data[ctrl_idx] will contain all the measurements
        # done by the ctrl_idx-th controller at all energy steps,
        # with the energy being the 'innermost' index.
        for data_point in self:
            # ctrl_idx is an index running over controllers that
            # measured key. We cannot simply enumerate, because
            # there may be controllers that did not measure the
            # stuff.
            ctrl_idx = 0
            for ctrl, measurements in data_point[key].items():
                if len(measurements) == 1:
                    raise ValueError('DataPoints contains energy-resolved '
                                     'data but it was attempted to return '
                                     'time resolved data.')
                if not separate_steps:
                    data[ctrl_idx].extend(measurements)
                    times[ctrl_idx].extend(data_point[QuantityInfo.TIMES][ctrl])
                else:
                    data[ctrl_idx].append(measurements)
                    times[ctrl_idx].append(data_point[QuantityInfo.TIMES][ctrl])
                    for i, time_set in enumerate(times[ctrl_idx]):
                        start_time = times[ctrl_idx][i][0]
                        times[ctrl_idx][i] = [t - start_time
                                              for t in times[ctrl_idx][i]]
                ctrl_idx += 1

        if include_energies:
            return data, times, self.nominal_energies
        return data, times

    def new_data_point(self, energy, controllers, cameras):
        """Create a new data_point.

        Parameters
        ----------
        energy : float
            Energy that is currently set.
        controllers : list of controllers
            All controllers used.
        cameras : list of cameras
            All cameras used.

        Returns
        -------
        None.
        """
        data_points_dict = {}
        for quantity in QuantityInfo:
            if any(ctrl.measures(quantity) for ctrl in controllers):
                data_points_dict[quantity] = {
                    controller: []
                    for controller in controllers
                    if controller.measures(quantity)
                    }
        data_points_dict[QuantityInfo.TIMES] = {}
        data_points_dict[QuantityInfo.IMAGES] = {}
        for controller in controllers:
            data_points_dict[QuantityInfo.TIMES][controller] = []
        for camera in cameras:
            data_points_dict[QuantityInfo.IMAGES][camera] = None
        data_points_dict[QuantityInfo.ENERGY] = energy
        self.append(data_points_dict)

    def read_data(self, csv_name):
        """Read data from csv file.

        Parameters
        ----------
        csv_name : Path
            Full path and file name to
            which the data will be saved.

        Returns
        -------
        None.
        """
        # TODO: fix reading!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Does not check for aliases yet.
        with open(csv_name, 'r', encoding='UTF8', newline='') as csv_file:
            reader = csv.DictReader(csv_file, delimiter=self.delimiter)
            first_row = next(reader)
            if QuantityInfo.ENERGY.label not in first_row:
                # Read file is not a ViPErLEED measurement.
                raise RuntimeError(
                    f'{csv_name} is not a ViPErLEED measurement.'
                    )
            key_dict = {}
            extractor = re.compile(r'(.*)\((.*)\)')
            for key in first_row:
                if key == QuantityInfo.ENERGY.label:
                    key_dict[key] = (key, None)
                    continue
                extracted = extractor.match(key)
                if not extracted:
                    continue
                key_dict[key] = list(extracted.groups())

            data_points_dict = {QuantityInfo.ENERGY: -1}
            for key in key_dict:
                quantity, device = key_dict[key]
                quantity_obj = QuantityInfo.from_label(quantity)
                if quantity_obj == QuantityInfo.ENERGY:
                    continue
                if quantity_obj not in data_points_dict:
                    data_points_dict[quantity_obj] = {}
                data_points_dict[quantity_obj][device] = []
            self.append(copy.deepcopy(data_points_dict))

            for row in (first_row, *reader):
                for column, value in row.items():
                    quantity, device = key_dict[column]
                    quantity_obj = QuantityInfo.from_label(quantity)
                    for allowed_quantity in QuantityInfo:
                        if quantity_obj == allowed_quantity:
                            converter = allowed_quantity.dtype
                            break
                    value = converter(value)
                    # if QuantityInfo.from_label(column) == QuantityInfo.ENERGY:
                    if quantity_obj == QuantityInfo.ENERGY:
                        if self[-1][QuantityInfo.ENERGY] == -1:
                            self[-1][QuantityInfo.ENERGY] = value
                        if value != self[-1][QuantityInfo.ENERGY]:
                            self.append(copy.deepcopy(data_points_dict))
                            self[-1][QuantityInfo.ENERGY] = value
                    else:
                        self[-1][quantity_obj][device].append(value)
            for point in self:
                print(point)
                print('####################################')

    def save_data(self, csv_name):
        """Save data.

        Parameters
        ----------
        csv_name : Path
            Full path and file name to
            which the data will be saved.

        Returns
        -------
        None.
        """
        first_line = []
        if self.cameras:
            for camera in self.cameras:
                first_line.append(camera.name)
        first_line.append(QuantityInfo.ENERGY.label)
        for quantity, ctrl_dict in self[0].items():
            if quantity not in self.__exceptional_keys:
                for controller in ctrl_dict:
                    first_line.append(
                        f"{quantity.label}({controller.serial.port_name})"
                        )

        with open(csv_name, 'w', encoding='UTF8', newline='') as file_name:
            writer = csv.writer(file_name, delimiter=self.delimiter)
            writer.writerow(first_line)
            for data_point in self:
                data = []
                max_length = max(
                    len(mt)
                    for ctrl, mt in data_point[QuantityInfo.TIMES].items()
                    )
                if self.cameras:
                    for camera in self.cameras:
                        images = [
                            data_point[QuantityInfo.IMAGES][camera]
                            ]*max_length
                        data.append(images)
                energy_list = [data_point[QuantityInfo.ENERGY]]*max_length
                data.append(energy_list)
                for quantity, ctrl_dict in data_point.items():
                    if quantity not in self.__exceptional_keys:
                        for ctrl, measurement in ctrl_dict.items():
                            extra_length = max_length - len(measurement)
                            measurement += ['Nan']*extra_length
                            data.append(measurement)
                for line in zip(*data):
                    writer.writerow(line)

    def is_time_resolved(self):
        """Check if the contained data is time-resolved."""
        if self.has_data():
            # Only if there are times saved already it is possible to
            # decide if the measurement is a time or energy resolved
            # measurement.
            if any(len(t) > 1 for t in self[0][QuantityInfo.TIMES].values()):
                return True
        return False

    def has_data(self):
        """Check if there is already data in the class."""
        if any(t for t in self[0][QuantityInfo.TIMES].values()):
            return True
        return False
