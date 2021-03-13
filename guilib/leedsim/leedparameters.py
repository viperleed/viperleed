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
    # keys that are not given as input, but are calculated at instantiation
    calculated_keys = {'bulkReciprocalBasis': None,}
    
    def __init__(self, data):
        # self.__data is the underlying dictionary that is accessed
        # when issuing self[key]
        if isinstance(data, LEEDParameters):
            # Just copy data into self.__data and skip all the rest
            self.__data = data
            return None
        if not isinstance(data, (dict, ConfigParser)):
            raise TypeError("LEEDParameters: 'dict' or 'ConfigParser' needed. "
                            f"Found {type(data).__name__!r} instead")

        self.__data = {**self.mandatory_keys,
                       **self.optional_keys,
                       **self.calculated_keys}
        if isinstance(data, dict):
            self.__from_dict(data)
        else:
            self.__from_file(data)

        # Calculate the parameters in self.calculated_keys
        self._calculate_missing()
        
        # And fix the 3D symmetry operations of the bulk group by embedding
        # them into the PlaneGroup instance directly
        self._3d_ops_to_bulk_group()

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
        return 'LEEDParameters(' + repr(self.__data).replace('\n','') + ')'
    
    def __eq__(self, other):
        """
        Equality method for LEEDParameters. Instances are considered equal if
        they produce the very same LEED pattern, i.e., they have:
        - the same beamIncidence
        - the same surface group
        - the same bulk basis, expressed in the same coordinate system, but
          accounting for possible sign changes of one of the unit vectors
        - the same surface lattices, i.e., the same domains. Also in this case,
          sign changes are acceptable
          TODO: probably this thing with the sign change needs revisiting.

        The algorithm used assumes that both self and other would produce
        lattices with highest symmetry

        Returns
        -------
        True if LEEDParameters are equal, NotImplemented otherwise
        """
        if not isinstance(other, LEEDParameters):
            # Try to see if other can be used as an argument to the
            # LEEDParameters constructor
            try:
                other = LEEDParameters(other)
            except (TypeError, NameError, ValueError):
                # If instantiation raises one of the known exceptions
                return NotImplemented

        # (1) first a basic check
        if self is other:
            return True

        # (2) beamIncidence
        if self['beamIncidence'] != other['beamIncidence']:
            return NotImplemented

        # (3) surface plane group
        if self['surfGroup'] != other['surfGroup']:
            return NotImplemented

        # (4) reciprocal bulk bases, accounting for sign changes
        bulk_transform = np.dot(np.linalg.inv(other['bulkReciprocalBasis']),
                                self['bulkReciprocalBasis']).T
        if not np.allclose(np.abs(bulk_transform), gl.PlaneGroup.E, atol=1e-4):
            return NotImplemented

        # (5) check whether the real-space surface bases are the same,
        #     accounting for sign changes
        surf_transform = np.dot(other['surfBasis'],
                                 np.linalg.inv(self['surfBasis']))
        if np.allclose(np.abs(surf_transform), gl.PlaneGroup.E):
            # (5.1) in this case, for the two patterns to be the same, we need
            #       also the bulk groups to be equal, including the operations
            #       that generate domains
            if ((self['bulkGroup'] == other['bulkGroup'])
                and (self['bulkGroup'].same_operations(other['bulkGroup'],
                                                       include_3d=True))):
                return True
            return NotImplemented

        # At this point we have:
        # * same beamIncidence
        # * same bulk basis (except sign changes)
        # * same surface group
        # * somewhat different surface bases

        # (6) The surface lattices could still be two symmetry-related domains.
        #     If this is the case, the two superlattice matrices should be
        #     related to one another by one of the symmetry operations of the
        #     bulk group
        # First we have to take into account that the BULK basis of the other
        # instance may have some sign change. This requires transforming its
        # superlattice to a common basis:
        super_other = np.dot(other['SUPERLATTICE'], bulk_transform.round())
        inv_super_self = np.linalg.inv(self['SUPERLATTICE'])

        # get the matrices that transform the superlattice of self into the
        # superlattice of other. Here we also have to account
        # for sign changes of the SURFACE basis
        sign_changes = (gl.PlaneGroup.E,   # no sign change
                        gl.PlaneGroup.C2,  # both vectors change sign
                        gl.PlaneGroup.Mx,  # b changes sign
                        gl.PlaneGroup.My)  # a changes sign
        # super_transforms is a "list" where the i-th element equals
        #    np.dot(inv_super_self, np.dot(sign_changes[i], super_other))
        super_transforms = np.einsum("ilm,jl,mk->ijk",
                                     sign_changes, inv_super_self, super_other)

        # (6.1) if the transform is a bulk operation, it certainly needs to
        #       have only integer entries. Keep only the ones that do:
        int_transforms = [t.round().astype(int)
                          for t in super_transforms
                          if np.allclose(t, t.round())]
        if len(int_transforms) == 0:
            return NotImplemented

        # (6.2) check if indeed one of the transforms equals one of the bulk
        #       symmetry operations
        bulk_ops = self['bulkGroup'].operations(include_3d=True)
        for transform in int_transforms:
            transform = gl.two_by_two_array_to_tuple(transform)
            if transform in bulk_ops:
                return True
        return NotImplemented

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
            # find the key in self that matches the one in data
            # when compared case-insensitive. There should be at most one.
            # All those that do not appear are neglected.
            keys = [kk for kk in self if kk.lower() == k.lower()]
            if len(keys):
                self[keys[0]] = v
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
            elif k.lower() == 'beamIncidence'.lower():
                if isinstance(v, str):
                    v = ast.literal_eval(v)
                if not isinstance(v, (list, tuple, np.ndarray)) or len(v) != 2:
                    raise ValueError("LEEDParameters: invalid beamIncidence "
                                     f"entry found: {v}")
                data[k] = gl.conventional_angles(*v)
    
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

    def _calculate_missing(self):
        """
        Calculates the extra parameters that are not in the original input
        but are useful to keep
        """
        # (1) bulkReciprocalBasis
        # build a dummy real-space surface lattice just to get its reciprocal
        # basis, from which one can find the bulk reciprocal basis. Also,
        # round a bit the result (8 decimal digits).        
        dummy_surf = gl.Lattice(self['surfBasis'])
        bulk_basis = np.dot(self['SUPERLATTICE'].T, dummy_surf.reciprocal_basis)
        self['bulkReciprocalBasis'] = bulk_basis.round(decimals=8)

    def _3d_ops_to_bulk_group(self):
        """
        Edit the 'bulkGroup' entry, which should be a PlaneGroup instance,
        such that it contains the correct glide and screw symmetry operations
        """
        # do something only if it makes sense
        if self['bulk3Dsym'] is None:
            return None

        # get dummy bulk to get the lattice shape
        dummy_bulk = gl.Lattice(self['bulkReciprocalBasis'], space='reciprocal',
                                group=self['bulkGroup'])
        bulk_shape = dummy_bulk.cell_shape

        # And apply the bulk operations to the bulk group
        self['bulkGroup'].screws_glides = (self['bulk3Dsym'], bulk_shape)


class LEEDParametersList(MutableSequence):
    """
    1-D only list of LEEDParameters that are suitable for producing LEED
    patterns, i.e., only structures with the same bulk lattice and symmetry can
    be added after instantiation. During instantiation checks that all are
    consistent
    """
    def __init__(self, data=None, duplicates=False):
        super().__init__()
        # self.__list is the list that underlies this container, i.e., the
        # one that is accessed when calling self[index]. It needs initializing
        # to empty because it is accessed in _process_input_data
        self.__list = []
        self.__list = self._process_input_data(data, duplicates)

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
        txt = 'LEEDParametersList([\n'
        txt += '\n'.join(repr(d) for d in self)
        txt += '\n])'
        return txt 

    def insert(self, idx, data, duplicates=False):
        data = self._process_input_data(data, duplicates)
        # when inserting, we should make sure that the data is actually a single
        # LEEDParameters to avoid creating a list with entries of different type
        # at different indices
        if not data:
            # do nothing if trying to insert nothing
            return None
        if len(data) > 1:
            # more than one parameter passed. Not acceptable for insertion.
            raise ValueError("Cannot insert multiple LEEDParameters at one "
                             "index. Use LEEDParametersList.extend(params) "
                             "or a list of LEEDParametersList")
        self.__list.insert(idx, data[0])

    def append(self, data, duplicates=False):
        self.insert(len(self), data, duplicates)

    def extend(self, data, duplicates=False):
        for d in data:
            self.append(d, duplicates)

    def _process_input_data(self, data, duplicates=False):
        """
        Do all the processing needed to accept the input data, i.e.,
        converting each entry to a LEEDParameters instance, verify that the
        the data passed is consistent, and remove duplicates if needed
        """
        data = self._to_leed_parameters(data)
        self._check_acceptable(*data)
        if not duplicates:
            # first remove duplicates from data only
            data = gl.remove_duplicates(data)
            # then keep only those that are not already in the list
            # this way, if, at earlier calls we wanted to keep duplicates,
            # they remain in self
            data = [d for d in data if d not in self]
        return data
    
    def _to_leed_parameters(self, data):
        """
        Processes the data given as input, and returns an appropriate list
        of LEEDParameters
        """
        if not data:
            ret = []
        elif isinstance(data, dict):
            # delegate to LEEDParameters the check on the dictionary
            ret = [LEEDParameters(data)]
        elif isinstance(data, ConfigParser):
            ret = [LEEDParameters(dict(data[section])) for section in data]
        elif isinstance(data, LEEDParameters):
            ret = [data]
        elif isinstance(data, LEEDParametersList):
            ret = data
        elif isinstance(data, (list, tuple, np.ndarray)):
            # here we check if any of the entries in the list is a
            # LEEDParametersList first, in which case we replace the entry
            # with as many LEEDParameters instances as there are in the list
            ret = []
            for d in data:
                if isinstance(d, LEEDParametersList):
                    ret.extend(d)
                else:
                    ret.append(LEEDParameters(d))
        else:
            raise TypeError("LEEDParametersList: invalid data type. Expected "
                            "dict, ConfigParser, LEEDParameters, "
                            "LEEDParametersList, or list, tuple, or "
                            "np.ndarray of the previous. "
                            f"Found {type(data).__name__!r} instead")
        return ret
    
    def _check_acceptable(self, *data):
        # First store the largest eMax and screenAperture. Will be set
        # as the eMax and screenAperture values if the data is acceptable
        all_params = (*self, *data)
        if len(all_params) == 0:
            return None
        emax = max(param['eMax'] for param in all_params)
        aperture = max(param['screenAperture'] for param in all_params)
        
        # check consistency of bulk
        inv_bulk_basis = np.linalg.inv(all_params[0]['bulkReciprocalBasis'])
        bulk_group = all_params[0]['bulkGroup']
        bulk_ops = set(bulk_group.operations(include_3d=False))
        for param in all_params[1:]:
            # bulk lattices should be the same, except, possibly, for a sign
            # change. Compute the matrix that transforms one into the other.
            # Then check that the transform is, except for a sign change, the
            # identity matrix. Since the bulk bases are the reciprocal ones,
            # one would need to calculate
            #        T = ((B'*)^-1 @ B*)^T,
            # but it is more convenient to calculate
            #        T* = B'* @ B*^(-1)
            # since abs(T) == E <--> abs(T*) == E
            transform = np.dot(param['bulkReciprocalBasis'], inv_bulk_basis)
            if not np.allclose(np.abs(transform), gl.PlaneGroup.E):
                raise ValueError("LEEDParametersList: Inconsistent bulk bases "
                                 "found in the input parameters")
            # bulk groups also should be the same, but excluding the 3D symmetry
            # as this only affects how many symmetry domains will be generated
            this_group = param['bulkGroup']
            same_bulk_group = (
                bulk_group == this_group
                and bulk_ops == set(this_group.operations(include_3d=False)))
            if not same_bulk_group:
                raise ValueError("LEEDParametersList: Inconsistent symmetry "
                                 "operations of bulk lattices in the "
                                 "input parameters")

        # after passing the check, eMax and screenAperture are set to the values
        # found above
        for param in all_params:
            param['eMax'] = emax
            param['screenAperture'] = aperture


