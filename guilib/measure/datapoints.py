"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-11-12
Author: Michele Riva
Author: Florian Doerr

Defines the DaDataPointsta class.
"""

from collections.abc import MutableSequence, Sequence


class DataPoints(MutableSequence):
    """Data storage class."""

    def __init__(self, *args):
        """Initialise data class"""
        self.__list = list(args)
        self.__controller_delays = [] #TODO: append delays to this list after instantiating all controllers

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

    def get_energy_resolved_data(self, key):
        """Get energy resolved data for each controller."""
        # TODO: only for energy resolved data!
        data = []
        if key not in self[0]:
            raise ValueError(f'{self.__class__.__name__}: invalid key {key}')
        for measurement in self[0][key]:
            if measurement:
                data.append([])

        for data_point in self:
            i = 0
            for measurement in data_point[key]:
                if len(measurement) == 1:
                    data[i].append(*measurement)
                else:
                    data[i].append(measurement)
                i += 1
        return data
        
    def get_time_resolved_data(self, key):
        """Get time resolved data for each controller."""
        # TODO: only for time resolved data!
        data = []


if __name__ == '__main__':
    tmp = DataPoints()
    a = dict()
    a['energy'] = 10
    a['times'] = [0, 1]
    a['measured_quantity'] = [[5,4,], [6,7,]]
    # a['measured_quantity'] = [[5,], [6,]]
    tmp[0:2] = a, a, a                 #uses __setitem__
    tmp.append(a)                      #uses insert
    tmp[0] = a                         #uses __setitem__
    print(tmp.get_energy_resolved_data('measured_quantity'))
    # print(tmp.__list)