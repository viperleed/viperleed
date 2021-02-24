"""
========================================
   ViPErLEED Graphical User Interface
========================================
*** module guilib.leedsim.LEEDParameters ***

Created: 2021-01-23
Author: Michele Riva

This module contains the definition of the LEEDParameters class and
LEEDParametersList, that fully represent a viperleed.LEEDPattern and
viperleed.LEEDPatternList.
"""

from collections import MutableMapping, MutableSequence
from configparser import ConfigParser
import copy
import ast
from warnings import warn as warning   # eventually will replace with logging

import numpy as np

from viperleed import guilib as gl
# import guilib as gl


class LEEDParameters(MutableMapping):
    """
    Dictionary-like container of parameters needed to represent a
    viperleed.LEEDPattern. Keys cannot be deleted, and only acceptable keys
    can be set.
    """
    mandatory_keys = {'eMax': None,
                      'surfBasis': None,
                      'SUPERLATTICE': None,
                      'surfGroup': None,
                      'bulkGroup': None}
    # optional keys with their defaults
    optional_keys = {'screenAperture' : 110,
                     'bulk3Dsym': None,
                     'beamIncidence': (0, 0)  # (theta, phi)
                    }
    
    def __init__(self, data):
        if isinstance(data, LEEDParameters):
            # probably can do it better, just copying data.__data into
            # self.__data and skipping all the rest, with the assumption that
            # a LEEDParameters instance is already fine as is
            # data = {k: v for k, v in data.items()}
            self.__data = data.__data
            return None
        if not isinstance(data, (dict, ConfigParser)):
            raise TypeError("LEEDParameters: 'dict' or 'ConfigParser' needed. "
                            f"Found {type(data).__name__!r} instead")

        self.__data = {**self.mandatory_keys, **self.optional_keys}
        if isinstance(data, dict):
            self.__from_dict(data)
        else:
            self.__from_file(data)

    def __len__(self):
        return len(self.__data)

    def __iter__(self):
        return iter(self.__data)

    def __setitem__(self, k, v):
        if not isinstance(k, str):
            raise TypeError("LEEDParameters: keys need to be strings")
        try:
            idx = [key.lower() for key in self.__data].index(k.lower())
        except ValueError:
            warning(f"Unknown LEED parameter {k!r} will be ignored.")
        else:
            # key is there, but may be in the wrong case. Get the right one.
            k = tuple(self.__data.keys())[idx]
            data = {k: v}
            self._check_and_process_all(data)
            self.__data[k] = data[k]

    def __delitem__(self, k):
        raise NotImplementedError

    def __getitem__(self, k):
        return self.__data[k]

    def __contains__(self, k):
        return k in self.__data
    
    def __repr__(self):
        return repr(self.__data)

    def __from_dict(self, data):
        """
        Create LEED parameters from a dictionary. Checks that the dictionary
        contains the mandatory keys, and sets the optional ones to their default
        values, issuing some warnings if the warn keyword is True
        """
        # make all keys lowercase so that duplicates are removed,
        # and such that we don't modify the original data passed
        data = {k.lower(): v for k, v in data.items()}
        self._check_mandatory(data)
        self._check_and_process_all(data)
        for k, v in data.items():
            # find the key in self.__data that matches the one in data
            # when compared case-insensitive. There should be at most one.
            # All those that do not appear are neglected.
            keys = [kk for kk in self.__data if kk.lower() == k.lower()]
            if len(keys):
                self.__data[keys[0]] = v
            else:
                warning(f"Unknown LEED parameter {k!r} will be ignored.")
    
    def __from_file(self, parser):
        if len(parser.sections()) > 1:
            raise TypeError("LEEDParameters: multiple structures found. Use "
                            "LEEDParametersList instead")
        data = dict(parser[parser.sections()[0]])
        self.__from_dict(data)

    def _check_mandatory(self, data):
        """
        Check that all necessary keys are present in the dictionary
        """
        if any(key not in data for key in self.mandatory_keys):
            # perhaps they are there, just with the wrong case
            low = [key.lower() for key in data]
            missing = ', '.join(key for key in self.mandatory_keys
                                if key.lower() not in low)
            if missing:
                raise NameError("Not enough LEED parameters. "
                                f"Missing: {missing}")
    
    def _check_and_process_all(self, data):
        self._process_input_data(data)
        self._check_acceptable(data)
    
    def _process_input_data(self, data):
        """
        This is called after _check_mandatory has been already executed to get
        the data dictionary with the correct types
        """
        for k, v in data.items():
            # skip those that can be strings:
            if k.lower() in ('bulk3Dsym'.lower(),):
                continue
            # and process all the others
            if k.lower() == 'eMax'.lower():
                try:
                    data[k] = float(v)
                except ValueError:
                    raise ValueError("LEEDParameters: invalid maximum energy "
                                     f"{v!r} found")
            elif k.lower() == 'screenAperture'.lower():
                try:
                    data[k] = float(v)
                except ValueError:
                    raise ValueError("LEEDParameters: invalid screen aperture "
                                     f"{v!r} found")
            elif k.lower() == 'surfBasis'.lower():
                m = None
                if isinstance(v, str):
                    m = gl.string_matrix_to_numpy(v, dtype=float,
                                                  needs_shape=(2, 2))
                else:
                    try:
                        m = np.asarray(v)
                    except TypeError:
                        raise ValueError("LEEDParameters: invalid matrix "
                                         "defining the surface periodicity")
                if m is None:
                    raise ValueError("LEEDParameters: invalid shape of the "
                                     "matrix defining the surface periodicity")
                data[k] = m
            elif k.lower() == 'SUPERLATTICE'.lower():
                m = None
                if isinstance(v, str):
                    m = gl.string_matrix_to_numpy(v, dtype=int,
                                                  needs_shape=(2, 2))
                else:
                    try:
                        m = np.asarray(v)
                    except TypeError:
                        raise ValueError("LEEDParameters: invalid SUPERLATTICE"
                                         " matrix.")
                if m is None:
                    raise ValueError("LEEDParameters: invalid shape of the "
                                     "SUPERLATTICE matrix")
                data[k] = m
            elif k.lower() == 'surfGroup'.lower():
                data[k] = gl.PlaneGroup(v)
            elif k.lower() == 'bulkGroup'.lower():
                data[k] = gl.PlaneGroup(v)
            elif k.lower() == 'beamIncidence'.lower():                          # check theta < 0 and update phi accordingly
                if isinstance(v, str):
                    v = ast.literal_eval(v)
                if not isinstance(v, (list, tuple, np.ndarray)) or len(v) != 2:
                    raise ValueError("LEEDParameters: invalid beamIncidence "
                                     f"entry found: {v}")
                data[k] = v
    
    def _check_acceptable(self, data_dict):
        """
        Check that the parameters dictionary contains acceptable input to be a
        LEED parameter. Raises ValueError if not fine.
        """
        for key, v in data_dict.items():
            if key.lower() == 'eMax'.lower() and v <= 0:
                raise ValueError("Maximum LEED energy should be positive.")
            if (key.lower() == 'screenAperture'.lower()
                and (v <= 0 or v > 180)):
                raise ValueError("screenAperture should be between "
                                 f"0 and 180. Found {v} instead.")
            if (key.lower() in ('surfBasis'.lower(), 'SUPERLATTICE'.lower())
                and np.shape(v) != (2, 2)):
                raise ValueError("Lattice basis and SUPERLATTICE need to have a"
                                 f" (2, 2) shape. Found {np.shape(v)} instead")
            if (key.lower() == 'SUPERLATTICE'.lower()
                and any(np.abs(np.asarray(v)%1).ravel() > 1e-4)):
                raise ValueError("SUPERLATTICE needs to be an integer-valued "
                                 "matrix")


class LEEDParametersList(MutableSequence):
    """
    1-D only list of LEEDParameters that are suitable for producing LEED
    patterns, i.e., only structures with the same bulk lattice and symmetry can
    be added after instantiation. During instantiation checks that all are
    consistent
    """
    def __init__(self, data=None):
        super().__init__()
        self.__list = []
        data = self._to_leed_parameters(data)
        self._check_acceptable(*data)
        self.__list = data  # the list that underlies this container

    def __delitem__(self, el):
        del self.__list[el]

    def __getitem__(self, el):
        # do this type-preserving. If a slice, should return same type,
        # if index, just an element
        if isinstance(el, slice):
            return LEEDParametersList(self.__list[el])
        return self.__list[el]
    
    def __setitem__(self, idx, data):
        data = self._to_leed_parameters(data)
        self._check_acceptable(*data)
        if not data:
            return None
        self.__list[idx] = data

    def __len__(self):
        return len(self.__list)

    def __repr__(self):
        return repr(self.__list)

    def insert(self, idx, data):
        data = self._to_leed_parameters(data)
        self._check_acceptable(*data)
        # when inserting, we should make sure that the data is actually a single
        # LEEDParameters to avoid creating a list with entries of different type
        # at different indices
        if not data:
            # do nothing if trying to insert nothing
            return None
        if len(data) > 1:
            # more than one parameter passed. Not acceptable for insertion.
            raise ValueError("Cannot insert multiple LEEDParameters at one "
                             "index. Use LEEDParametersList.append(params) "
                             "or a list of LEEDParametersList")
        self.__list.insert(idx, data[0])

    def append(self, data):
        self.insert(len(self), data)
    
    # def extend(self, data):  # not sure if I want this
        # # data needs to be an unpackable iterable
        # self.insert(len(self), *data)
    
    def _to_leed_parameters(self, data):
        """
        Processes the data given as input, and returns an appropriate list
        of LEEDParameters
        """
        if data is None:
            ret = []
        elif isinstance(data, dict):
            # delegate to LEEDParameters the check on the dictionary
            ret = [LEEDParameters(data)]
        elif isinstance(data, ConfigParser):
            ret = [LEEDParameters(dict(data[section])) for section in data]
        elif isinstance(data, LEEDParameters):
            ret = [data]    # was [LEEDParameters], but it must have been wrong!
        elif isinstance(data, LEEDParametersList):
            ret = data
        elif isinstance(data, (list, tuple, np.ndarray)):
            ret = [LEEDParameters(d) for d in data]
        else:
            raise TypeError("LEEDParametersList: invalid data type. Expected "
                            "dict, ConfigParser, LEEDParameters, "
                            "LEEDParametersList, list, tuple, or np.ndarray. "
                            f"Found {type(data).__name__!r} instead")
        return ret
    
    def _check_acceptable(self, *data):
        # First store the largest eMax and screenAperture. Will be set
        # as the eMax and screenAperture values if the data is acceptable
        all_params = (*self.__list, *data)
        if len(all_params) == 0:
            return None
        emax = max(param['eMax'] for param in all_params)
        aperture = max(param['screenAperture'] for param in all_params)
        
        leeds = []
        for params in copy.deepcopy(all_params):
            # set up LEEDPatterns with eMax = 1 to make initialization faster,
            # as this does not calculate all the beams.
            params['eMax'] = 1  # this is the reason for the deepcopy
            leeds.append(gl.LEEDPattern(params))
        
        # check consistency of bulk
        bulk = leeds[0].reciprocal_lattices['bulk']
        b_ops = set(bulk.group.group_ops(include_3d=True))
        for leed in leeds[1:]:
            this_bulk = leed.reciprocal_lattices['bulk']
            # bulk lattices should be the same
            if not np.allclose(bulk.basis, this_bulk.basis):
                raise ValueError("LEEDParametersList: Inconsistent bulk bases "
                                 "found in the input parameters")
            # bulk groups also should be the same. Not sure how to handle the
            # case in which some of the parameters does not contain the
            # bulk3Dsym but others do. Right now it raises errors.
            same_bulk_group = (
                bulk.group == this_bulk.group
                and b_ops == set(this_bulk.group.group_ops(include_3d=True)))
            if not same_bulk_group:
                raise ValueError("LEEDParametersList: Inconsistent symmetry "
                                 "operations of bulk lattices in the "
                                 "input parameters")
            # TODO: check also that all elements give unique LEED patterns, 
            #       which includes rotational domains. Remove those that do not.
            #       May be easier to do by adding first __eq__ and __neq__
            #       to LEEDPattern
        
        # after passing the check, eMax and screenAperture are set to the values
        # found above
        for param in all_params:
            param['eMax'] = emax
            param['screenAperture'] = aperture



