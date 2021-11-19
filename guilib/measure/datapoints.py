"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-11-12
Author: Michele Riva
Author: Florian Doerr

Defines the DaDataPointsta class.
"""
import csv
from collections.abc import MutableSequence, Sequence

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.hardwarebase import (
    ViPErLEEDErrorEnum, QMetaABC, emit_error)


class DataErrors(ViPErLEEDErrorEnum):
    """Errors that might occur during a measurement cycle."""
    INVALID_MEASUREMENT = (400,
                           "The returned data dictionary contained a section "
                           "that was not specified in the DataPoints class.")


class DataPoints(qtc.QObject, MutableSequence, metaclass=QMetaABC):
    """Data storage class
    .
    The plot_info dictionary in this class may be edited.
    Each section (label) needs to contain the unit and
    the scaling ('lin' for linear and 'log' for logarithmic
    scaling.)
    """
    # Is emitted when an error occurs
    error_occurred = qtc.pyqtSignal(tuple)

    plot_info = {}
    plot_info['images'] = ['Number']
    plot_info['nominal_energy'] = ['eV', 'lin']
    plot_info['measurement_t'] = ['s', 'lin']
    plot_info['I0'] = ['uA', 'lin']
    plot_info['HV'] = ['eV', 'lin']
    plot_info['Isample'] = ['V', 'log']
    plot_info['temperature'] = ['°C', 'lin']
    plot_info['cold_junction'] = ['°C', 'lin']

    def __init__(self, *args):
        """Initialise data class."""
        super().__init__()
        self.__list = list(args)
        # self.__exceptional_keys = ('images', 'nominal_energy', 'measurement_t')
        self.__exceptional_keys = ('images', 'nominal_energy')
        # primary_controller needs to be given to this class
        # before starting measurements.
        self.primary_controller = None
        self.controllers = None

    @property
    def nominal_energies(self):
        """Return energies."""
        energies = []
        for data_point in self:
            energies.append(data_point['nominal_energy'])
        return energies

    @property
    def cameras(self):
        return tuple(camera for camera, _ in self[0]['images'])

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

    def add_data(self, receive, controller, primary_delay=None):
        """Add received data to currently active data point."""
        for quantity, value in receive.items():
            if quantity not in self[-1]:
                emit_error(self, DataErrors.INVALID_MEASUREMENT)
            else:
                self[-1][quantity][controller].append(value)

        if not self[-1]['measurement_t'][controller]:
            time = controller.serial.time_stamp
            if primary_delay:
                time += primary_delay/1000
            self[-1]['measurement_t'][controller].append(time)

    def add_image_names(self, name):
        """Add image names to currently active data point."""
        for camera in self[-1]['images']:
            # TODO: may want to add camera name to image name
            self[-1]['images'][camera] = name

    def calculate_times(self):
        """Calculate times for the last data point."""
        primary_time = (
                self[-1]['measurement_t'][self.primary_controller][0]
                )
        for ctrl in self[-1]['measurement_t']:
            self[-1]['measurement_t'][ctrl][0] -= primary_time
            self[-1]['measurement_t'][ctrl][0] += ctrl.initial_delay
            quantity = ctrl.measured_quantities[0]
            length = len(self[-1][quantity][ctrl])
            for i in range(length - 1):
                self[-1]['measurement_t'][ctrl].append(
                    self[-1]['measurement_t'][ctrl][0] +
                    ctrl.measurement_interval * (i + 1)
                    )

    def get_energy_resolved_data(self, key, include_energies=False):
        """Get energy resolved data for each controller.

        Parameters
        ----------
        key : str
            The key associated with the requested
            measured quantity.
        include_energies : bool
            Whether energy should be returned.

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
                    raise ValueError('DataPoints contains timeresolved data '
                                     'but it was attempted to return energy '
                                     'resolved data.')
                i += 1
        if include_energies:
            return data, self.nominal_energies
        return data

    def get_time_resolved_data(self, key, include_energies=False,
                               include_times=False):
        """Get time resolved data for each controller.

        Parameters
        ----------
        key : str
            The key associated with the requested
            measured quantity.
        include_energies : bool
            Whether energy should be returned.
        include_times : bool
            Whether measurement times should be returned.

        Returns
        -------
        data : list
            Contains measurements and associated
            starting times.
        self.nominal_energies : list
            Contains set energies.
        times : list
            Contains times at which measurements
            were taken.
        """
        data = []
        times = []
        if key not in self[0]:
            raise ValueError(f'{self.__class__.__name__}: invalid key {key}')
        for measurement in self[0][key]:
            if measurement:
                data.append([])
                times.append([])

        for data_point in self:
            i = 0
            for ctrl, measurement in data_point[key].items():
                if len(measurement) == 1:
                    raise ValueError('DataPoints contains energy resolved '
                                     'data but it was attempted to return '
                                     'time resolved data.')
                else:
                    data[i].append(measurement)
                    times[i].append(data_point['measurement_t'][ctrl])
                i += 1
        if include_energies and include_times:
            return data, times, self.nominal_energies
        elif include_energies:
            return data, self.nominal_energies
        elif include_times:
            return data, times
        return data

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
        for quantity in self.plot_info:
            if any(ctrl.measures(quantity) for ctrl in controllers):
                data_points_dict[quantity] = {
                    controller: []
                    for controller in controllers
                    if controller.measures(quantity)
                    }
        data_points_dict['measurement_t'] = {}
        data_points_dict['images'] = {}
        for controller in controllers:
            data_points_dict['measurement_t'][controller] = []
        for camera in cameras:
            data_points_dict['images'][camera] = None
        data_points_dict['nominal_energy'] = energy
        self.append(data_points_dict)

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
        for camera in self.cameras:
            first_line.append(camera.name)
        first_line.append('nominal_energy')
        for quantity, ctrl_dict in self[0].items():
            if quantity not in self.__exceptional_keys:
                for controller in ctrl_dict:
                    first_line.append(
                        f"{quantity}_{controller.serial.port_name}"
                        )

        with open(csv_name, 'w', encoding='UTF8', newline='') as file_name:
            writer = csv.writer(file_name, delimiter = ';')
            writer.writerow(first_line)
            for data_point in self:
                data = []
                max_length = max(
                    len(mt)
                    for ctrl, mt in data_point['measurement_t'].items()
                    )
                for camera in self.cameras:
                    images = [data_points_dict['images'][camera]]*max_length
                    data.append(images)
                energy_list = [data_point['nominal_energy']]*max_length
                data.append(energy_list)
                for quantity, ctrl_dict in data_point.items():
                    if quantity not in self.__exceptional_keys:
                        for ctrl, measurement in ctrl_dict.items():
                            extra_length = max_length - len(measurement)
                            measurement += ['Nan']*extra_length
                            data.append(measurement)
                for line in zip(*data):
                    writer.writerow(line)

