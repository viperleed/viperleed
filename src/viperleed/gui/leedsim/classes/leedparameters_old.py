"""Module leedparameters of viperleed.guilib.leedsim.

========================================
   ViPErLEED Graphical User Interface
========================================
*** module guilib.leedsim.leedparameters ***

Created: 2021-01-23
Author: Michele Riva

This module contains the definition of the LEEDParameters class and
LEEDParametersList, that fully represent a viperleed.LEEDPattern and
viperleed.LEEDPatternList.
"""

from collections.abc import MutableMapping, MutableSequence
from configparser import ConfigParser
# import copy
import ast
from warnings import warn as warning   # eventually will replace with logging

import numpy as np

from viperleed.gui.classes import planegroup
from viperleed.gui.classes.lattice2d import Lattice2D as Lattice
from viperleed.gui.helpers import conventional_angles
from viperleed.gui.helpers import remove_duplicates
from viperleed.gui.helpers import single_spaces_only
from viperleed.gui.helpers import string_matrix_to_numpy
from viperleed.gui.helpers import two_by_two_array_to_tuple
from viperleed.gui.leedsim.classes.leedparser import LEEDParser


class LEEDParameters(MutableMapping):
    """Dictionary-like container of LEED pattern parameters.

    Dictionary-like container of parameters needed to represent a
    viperleed.LEEDPattern. Keys cannot be deleted, and only acceptable
    keys can be set.
    """

    mandatory_keys = {'eMax': None,
                      'surfBasis': None,
                      'SUPERLATTICE': None,
                      'surfGroup': None,
                      'bulkGroup': None}

    # optional keys with their defaults:
    #     'name' will be the empty string unless it is defined in the
    #     dictionary used for instantiation, or if the LEEDParameters
    #     is loaded from a file (in which case it can be defined either
    #     as a section header or as a key of the configuration file;
    #     the latter takes precedence)
    optional_keys = {'screenAperture': 110,
                     'bulk3Dsym': None,
                     'beamIncidence': (0, 0),  # (theta, phi)
                     'name': '', }

    # keys that are not given as input, but are calculated at instantiation
    calculated_keys = {'bulkReciprocalBasis': None, }

    def __init__(self, data):
        # self.__data is the underlying dictionary that is accessed
        # when issuing self[key]
        if isinstance(data, LEEDParameters):
            # Just copy data into self.__data and skip all the rest
            self.__data = dict(data)
            return
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
        """len(self)."""
        return len(self.__data)

    def __iter__(self):
        """Return iterator of self."""
        return iter(self.__data)

    def __setitem__(self, k, value):
        """Item setter.

        Parameters
        ----------
        k : str
            Name of parameter to set. Case-insensitive. Only
            acceptable parameters will be set.
        value : depends on k
            Value of the parameter to set.

        Raises
        ------
        Warning
            If k is not an acceptable parameter name.
        """
        if not isinstance(k, str):
            raise TypeError("LEEDParameters: keys need to be strings")
        try:
            idx = [key.lower() for key in self.__data].index(k.lower())
        except ValueError:
            warning(f"Unknown LEED parameter {k!r} will be ignored.")
        else:
            # key is there, but may be in the wrong case. Get the right one.
            k = tuple(self.__data.keys())[idx]
            data = {k: value}
            self._check_and_process_all(data)
            self.__data[k] = data[k]

    def __delitem__(self, k):
        """Not possible to delete items."""
        raise NotImplementedError

    def __getitem__(self, k):
        """Get item with key k."""
        return self.__data[k]

    def __contains__(self, k):
        """Return whether k is in self."""
        return k in self.__data

    def __repr__(self):
        """Return string representation of self."""
        return 'LEEDParameters(' + repr(self.__data).replace('\n', '') + ')'

    def __eq__(self, other):
        """Equality method for LEEDParameters.

        Instances are considered equal if they produce
        the very same LEED pattern, i.e., they have:
        - the same beamIncidence
        - the same surface group
        - the same bulk basis, expressed in the same coordinate
          system, but accounting for possible sign changes of
          one of the unit vectors
        - the same surface lattices, i.e., the same domains.
          Also in this case, sign changes are acceptable
          TODO: probably this thing with the sign change needs revisiting.

        The algorithm used assumes that both self and other
        would produce lattices with highest symmetry

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
        if not np.allclose(np.abs(bulk_transform), planegroup.E, atol=1e-4):
            return NotImplemented

        # (5) check whether the real-space surface bases are the same,
        #     accounting for sign changes
        surf_transform = np.dot(other['surfBasis'],
                                np.linalg.inv(self['surfBasis']))
        if np.allclose(np.abs(surf_transform), planegroup.E):
            # (5.1) in this case, for the two patterns to be the same,
            #       we need also the bulk groups to be equal, including
            #       the operations that generate domains
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

        # (6) The surface lattices could still be two symmetry-related
        #     domains. If this is the case, the two superlattice
        #     matrices should be related to one another by one of the
        #     symmetry operations of the bulk group
        # First we have to take into account that the BULK basis of
        # the other instance may have some sign change. This requires
        # transforming its superlattice to a common basis:
        super_other = np.dot(other['SUPERLATTICE'], bulk_transform.round())
        inv_super_self = np.linalg.inv(self['SUPERLATTICE'])

        # get the matrices that transform the superlattice of self
        # into the superlattice of other. Here we also have to
        # account for sign changes of the SURFACE basis
        sign_changes = (planegroup.E,   # no sign change
                        planegroup.C2,  # both vectors change sign
                        planegroup.Mx,  # b changes sign
                        planegroup.My)  # a changes sign
        # super_transforms is a "list" where the i-th element equals
        # np.dot(inv_super_self, np.dot(sign_changes[i], super_other))
        super_transforms = np.einsum("ilm,jl,mk->ijk",
                                     sign_changes,
                                     inv_super_self,
                                     super_other)

        # (6.1) if the transform is a bulk operation, it certainly
        #       needs to have only integer entries. Keep only the
        #       ones that do:
        int_transforms = [t.round().astype(int)
                          for t in super_transforms
                          if np.allclose(t, t.round())]
        if len(int_transforms) == 0:
            return NotImplemented

        # (6.2) check if indeed one of the transforms equals one
        #       of the bulk symmetry operations
        bulk_ops = self['bulkGroup'].operations(include_3d=True)
        for transform in int_transforms:
            transform = two_by_two_array_to_tuple(transform)
            if transform in bulk_ops:
                return True
        return NotImplemented

    def __from_dict(self, data):
        """Create LEEDParameters from a dictionary.

        Checks that the dictionary contains the mandatory keys, and
        sets the optional ones to their default values, issuing some
        warnings if the warn keyword is True

        Parameters
        ----------
        data : dict or LEEDParameters
            Dictionary from which to create self.

        Returns
        -------
        None.
        """
        # Make all keys lowercase so that duplicates are removed,
        # and such that we don't modify the original data passed
        data = {k.lower(): v for k, v in data.items()}
        self._check_mandatory(data)
        self._check_and_process_all(data)
        for k, value in data.items():
            # Find the key in self that matches the one in data when
            # compared case-insensitive.  There should be at most
            # one.  All those that do not appear are neglected.
            keys = [kk for kk in self if kk.lower() == k.lower()]
            if len(keys):
                self[keys[0]] = value
            else:
                warning(f"Unknown LEED parameter {k!r} will be ignored.")

    def __from_file(self, parser):
        """Create LEEDParameters from a configuration file.

        Parameters
        ----------
        parser : ConfigParser
            The parser for the configuration file from which to
            initialize. It should contain only one section.

        Returns
        -------
        None.
        """
        if len(parser.sections()) > 1:
            raise TypeError("LEEDParameters: multiple structures found. Use "
                            "LEEDParametersList instead")
        # Take the section name as a 'name' attribute
        section = parser.sections()[0]
        self.__from_dict({**{'name': section}, **parser[section]})

    def _check_mandatory(self, data):
        """Check that all necessary keys are present in the dictionary."""
        if any(key not in data for key in self.mandatory_keys):
            # Perhaps they are there, just with the wrong case
            low = [key.lower() for key in data]
            missing = ', '.join(key for key in self.mandatory_keys
                                if key.lower() not in low)
            if missing:
                raise NameError("Not enough LEED parameters. "
                                f"Missing: {missing}")

    def _check_and_process_all(self, data):
        """Process data and check that the input is OK.

        Parameters
        ----------
        data : dict or LEEDParameters
            dictionary of parameters. Modified in place.

        Returns
        -------
        None.
        """
        self._process_input_data(data)
        self._check_acceptable(data)

    @staticmethod
    def _process_input_data(data):
        """Process data to correct format.

        This is called after _check_mandatory has been already
        executed to get the data dictionary with the correct types.

        Parameters
        ----------
        data : dict or LEEDParameters
            The data to be processed to get the correct types.
            Modified in place.

        Returns
        -------
        None.
        """
        for k, value in data.items():
            # Process all the values that should not be strings
            if k.lower() in ('emax', 'screenaperture'):
                try:
                    data[k] = float(value)
                except ValueError as err:
                    raise ValueError(f"LEEDParameters: invalid {k} "
                                     f"{value!r} found") from err

            elif k.lower() in ('surfbasis', 'superlattice'):
                dtype = float if k.lower() == 'surfbasis' else int
                try:
                    data[k] = LEEDParameters.__array_from_input(value, dtype)
                except ValueError as err:
                    raise ValueError(f"{err}: {k}") from err

            elif k.lower() in ('surfgroup', 'bulkgroup'):
                data[k] = planegroup.PlaneGroup(value)

            elif k.lower() == 'beamincidence':
                data[k] = LEEDParameters.__beam_incidence_from_input(value)

    @staticmethod
    def __array_from_input(input_matrix, dtype):
        """Convert input to 2x2 numpy.array.

        Parameters
        ----------
        input_matrix : str or iterable
        dtype : dtype
            Data type to use for the conversion from string

        Returns
        -------
        numpy.ndarray
            shape (2, 2)

        Raises
        ------
        ValueError
            if not possible to convert to (2, 2) array
        """
        matrix = None
        if isinstance(input_matrix, str):
            matrix = string_matrix_to_numpy(input_matrix, dtype=dtype,
                                            needs_shape=(2, 2))
        else:
            try:
                matrix = np.asarray(input_matrix)
            except TypeError as err:
                raise ValueError("LEEDParameters: invalid "
                                 "matrix input") from err
        if matrix is None or matrix.shape != (2, 2):
            raise ValueError("LEEDParameters: invalid matrix shape")
        return matrix

    @staticmethod
    def __beam_incidence_from_input(angles):
        """Get beamIncidence tuple from input.

        Parameters
        ----------
        angles : iterable or str
            When an iterable, it should have only two elements.
            When a str, it should be possible to eval-it to a
            2-elements iterable.

        Raises
        ------
        ValueError
            If not possible to convert to beamIncidence angles

        Returns
        -------
        theta, phi
            Theta in [0, 90] range, phi in [0,360).
        """
        if isinstance(angles, str):
            angles = ast.literal_eval(angles)
        if (not isinstance(angles, (list, tuple, np.ndarray))
                or np.shape(angles) != (2,)
                or not -90 <= angles[0] <= 90):
            raise ValueError("LEEDParameters: invalid beamIncidence "
                             f"entry found: {angles}")
        return conventional_angles(*angles)

    @staticmethod
    def _check_acceptable(data_dict):
        """Check that the dict input is acceptable.

        Parameters
        ----------
        data_dict : dict or LEEDParameters

        Raises
        ------
        ValueError if not data_dict contains unacceptable data.
        """
        for key, value in data_dict.items():
            if key.lower() == 'emax' and value <= 0:
                raise ValueError("Maximum LEED energy should be positive.")
            if (key.lower() == 'screenaperture'
                    and (value <= 0 or value > 180)):
                raise ValueError("screenAperture should be between "
                                 f"0 and 180. Found {value} instead.")
            if (key.lower() in ('surfbasis', 'superlattice')
                    and np.shape(value) != (2, 2)):
                raise ValueError("Lattice basis and SUPERLATTICE need to have "
                                 f"a (2, 2) shape. Found {np.shape(value)} "
                                 "instead")
            if (key.lower() == 'superlattice'
                    and any(np.abs(np.asarray(value) % 1).ravel() > 1e-4)):
                raise ValueError("SUPERLATTICE needs to be an integer-valued "
                                 "matrix")

    def _calculate_missing(self):
        """Calculate extra parameters that used only internally.

        Calculates the extra parameters that are not in the original
        input but are useful to keep
        """
        # (1) bulkReciprocalBasis
        # build a dummy real-space surface lattice just to get its reciprocal
        # basis, from which one can find the bulk reciprocal basis. Also,
        # round a bit the result (8 decimal digits).
        bulk_basis = np.dot(self['SUPERLATTICE'].T,
                            Lattice(self['surfBasis']).reciprocal_basis)
        self['bulkReciprocalBasis'] = bulk_basis.round(decimals=8)

    def _3d_ops_to_bulk_group(self):
        """Edit the 'bulkGroup' entry to contain screws/glides.

        Edit the 'bulkGroup' entry, which should be a PlaneGroup
        instance, such that it contains the correct glide and screw
        symmetry operations.
        """
        # do something only if it makes sense
        if self['bulk3Dsym'] is None:
            return

        # get dummy bulk to get the lattice shape
        dummy_bulk = Lattice(self['bulkReciprocalBasis'],
                             space='reciprocal',
                             group=self['bulkGroup'])
        bulk_shape = dummy_bulk.cell_shape

        # And apply the bulk operations to the bulk group
        self['bulkGroup'].set_screws_glides(self['bulk3Dsym'], bulk_shape)


class LEEDParametersList(MutableSequence):
    """1-D only list of LEEDParameters.

    1-D only list of LEEDParameters that are suitable for producing
    LEED patterns, i.e., only structures with the same bulk lattice
    and symmetry can be added after instantiation.

    During instantiation checks that all are consistent.
    """

    def __init__(self, data=None, keep_duplicates=False):
        """Initialize LEEDParametersList.

        Parameters
        ----------
        data : iterable or LEEDParametersList
            Each element is either dict, ConfigParser, or
            LEEDParameters. Each element will be checked for
            compatibility, i.e., all should have the same bulk.
        keep_duplicates : bool (default=False)
            If False, only one of the entries in data that
            would produce identical LEED patterns (i.e., same
            spots, same symmetry domains, same symmetry relations
            between spots) is kept.

        Returns
        -------
        None.
        """
        super().__init__()
        # self.__list is the list that underlies this container,
        # i.e., the one that is accessed when calling self[index].
        # It needs initializing to empty because it is accessed
        # in _process_input_data
        self.__list = []
        self.__list = self._process_input_data(data, keep_duplicates)

    def __delitem__(self, index):
        """Remove item."""
        del self.__list[index]

    def __getitem__(self, index):
        """Get item at index."""
        # Do this type-preserving. If a slice, should
        # return same type, if index, just an element
        if isinstance(index, slice):
            return LEEDParametersList(self.__list[index])
        return self.__list[index]

    def __setitem__(self, idx, data):
        """Set item at idx to data.

        The item is set only if it is an acceptable LEEDParameter,
        given the data already present in self.

        Parameters
        ----------
        idx : int or slice
        data : dict, ConfigParser or LEEDParameters
        """
        data = self._to_leed_parameters(data)
        self._check_acceptable(*data)
        if not data:
            return
        self.__list[idx] = data

    def __len__(self):
        """Return len(self)."""
        return len(self.__list)

    def __repr__(self):
        """Return string representiation of self."""
        txt = 'LEEDParametersList([\n'
        txt += '\n'.join(repr(d) for d in self)
        txt += '\n])'
        return txt

    def insert(self, idx, data, keep_duplicates=False):
        """Insert data at index.

        The data is inserted only if it is compatible with the
        data already present in self.

        Parameters
        ----------
        idx : int
            Index at which to insert
        data : dict, ConfigParser, or LEEDParameters
            The data to be inserted. Will be checked for compatibility.
        keep_duplicates : bool (default=False)
            If False, data is inserted only if it is not a duplicate
            of elements already present in self.
        """
        data = self._process_input_data(data, keep_duplicates)
        # when inserting, we should make sure that the data is actually
        # a single LEEDParameters to avoid creating a list with entries
        # of different type at different indices
        if not data:
            # Do nothing if trying to insert nothing
            return
        if len(data) > 1:
            # More than one parameter passed.
            # Not acceptable for insertion.
            raise ValueError("Cannot insert multiple LEEDParameters at one "
                             "index. Use LEEDParametersList.extend(params) "
                             "or a list of LEEDParametersList")
        self.__list.insert(idx, data[0])

    def append(self, data, keep_duplicates=False):
        """Insert a single element at the end.

        The data is appended only if it is compatible with the
        data already present in self.

        Parameters
        ----------
        data : dict, ConfigParser, or LEEDParameters
            The data to be inserted. Will be checked for compatibility.
        keep_duplicates : bool (default=False)
            If False, data is appended only if it is not a duplicate
            of elements already present in self.
        """
        self.insert(len(self), data, keep_duplicates)

    def extend(self, data, keep_duplicates=False):
        """Insert a elements from an iterable at the end.

        Only those entries that are compatible with the
        data already present in self are appended.

        Parameters
        ----------
        data : iterable or LEEDParametersList
            Each element is either dict, ConfigParser, or LEEDParameters
            The data to be inserted. Each element will be checked for
            compatibility.
        keep_duplicates : bool (default=False)
            If False, append only those elements in data that are not
            a duplicate of elements already present in self.
        """
        for datum in data:
            self.append(datum, keep_duplicates)

    def set_beam_incidence(self, theta=None, phi=None):
        """Set the 'beamIncidence' key to the given angles.

        The same angles will be set for all the entries in self.

        Parameters
        ----------
        theta : float, default=None
            Polar incidence angle. Should be in the [-90, 90] range.
            Normal icidence is theta=0. No modification if not given
            or if None.
        phi : float
            Azimuthal incidence angle. Positive counterclockwise,
            measured with respect to the positive x axis of the
            Cartesian coordinates in which the real-space bases
            are expressed. Will be taken modulo 360.
        """
        if theta is None and phi is None:
            return
        old_angles = self[0]['beamIncidence']
        if theta is None:
            theta = old_angles[0]
        elif phi is None:
            phi = old_angles[1]
        try:
            theta = float(theta)
            phi = float(phi)
        except TypeError as err:
            raise TypeError("Invalid angles. "
                            "Expected numbers.") from err
        if (theta, phi) == old_angles:
            return
        new_angles = conventional_angles(theta, phi)
        if new_angles == old_angles:
            return
        for param in self:
            param['beamIncidence'] = new_angles

    def _process_input_data(self, data, keep_duplicates=False):
        """Prcess the given data.

        Do all the processing needed to accept the input data, i.e.,
        converting each entry to a LEEDParameters instance, verify
        that the data passed is consistent, and remove duplicates
        if needed.

        Parameters
        ----------
        data : iterable or LEEDParametersList
            Each element is either dict, ConfigParser, or
            LEEDParameters. Each element will be checked for
            compatibility.
        keep_duplicates : bool, default=False
            If False, discard those elements in data that are
            a duplicate of elements already present in self,
            as well as duplicates within data.
        """
        data = self._to_leed_parameters(data)
        self._check_acceptable(*data)
        if not keep_duplicates:
            # first remove duplicates from data only
            data = remove_duplicates(data)
            # then keep only those that are not already in the list.
            # This way, if at earlier calls we wanted to keep
            # duplicates, they remain in self
            data = [d for d in data if d not in self]
        return data

    @staticmethod
    def _to_leed_parameters(data):
        """Convert data to a list of LEEDParameters.

        Parameters
        ----------
        data : dict, ConfigParser, LEEDParameters, or iterable of the same

        Returns
        -------
        list of LEEDParameters
        """
        if not data:
            ret = []
        elif isinstance(data, dict):
            # delegate to LEEDParameters the check on the dictionary
            ret = [LEEDParameters(data)]
        elif isinstance(data, ConfigParser):
            ret = [LEEDParameters({**{'name': section}, **data[section]})
                   for section in data]
        elif isinstance(data, LEEDParameters):
            ret = [data]
        elif isinstance(data, LEEDParametersList):
            ret = data
        elif isinstance(data, (list, tuple, np.ndarray)):
            # Here we check if any of the entries in the list
            # is a LEEDParametersList first, in which case we
            # replace the entry with as many LEEDParameters
            # instances as there are in the LEEDParametersList
            ret = []
            for elem in data:
                if isinstance(elem, LEEDParametersList):
                    ret.extend(elem)
                else:
                    ret.append(LEEDParameters(elem))
        else:
            raise TypeError("LEEDParametersList: invalid data type. Expected "
                            "dict, ConfigParser, LEEDParameters, "
                            "LEEDParametersList, or list, tuple, or "
                            "np.ndarray of the previous. "
                            f"Found {type(data).__name__!r} instead")
        return ret

    def _check_acceptable(self, *data):
        """Check that the data is acceptable.

        Data is considered acceptable if the bulk lattices have the
        same bases (includng possibly sign changes), and if the
        bulk groups contain the same symmetry operations.  At the
        same time, this method replaces all 'eMax' and 'screenAperture'
        keys with the maxima among all the data (and of self).

        Parameters
        ----------
        *data : dict, ConfigParser, or LEEDParameters

        Raises
        ------
        ValueError
            If any of the data is inconsistent

        Returns
        -------
        None.

        """
        # First store the largest eMax and screenAperture.
        # Will be set as the eMax and screenAperture values
        # if the data is acceptable
        all_params = (*self, *data)
        if len(all_params) == 0:
            return
        emax = max(param['eMax'] for param in all_params)
        aperture = max(param['screenAperture'] for param in all_params)

        # check consistency of bulk
        inv_bulk_basis = np.linalg.inv(all_params[0]['bulkReciprocalBasis'])
        bulk_group = all_params[0]['bulkGroup']
        for param in all_params[1:]:
            # Bulk lattices should be the same, except, possibly, for a
            # sign change. Compute the matrix that transforms one into
            # the other. Then check that the transform is, except for a
            # sign change, the identity matrix. Since the bulk bases
            # are the reciprocal ones, one would need to calculate
            #        T = ((B'*)^-1 @ B*)^t,
            # but it is more convenient to calculate
            #        T* = B'* @ B*^(-1)
            # since abs(T) == E <--> abs(T*) == E
            transform = np.dot(param['bulkReciprocalBasis'], inv_bulk_basis)
            if not np.allclose(np.abs(transform), planegroup.E):
                raise ValueError("LEEDParametersList: Inconsistent bulk bases "
                                 "found in the input parameters")

            # Bulk groups also should be the same, but excluding the
            # 3D symmetry as this only affects how manysymmetry
            # domains will be generated
            this_group = param['bulkGroup']
            same_bulk_group = (
                bulk_group == this_group
                and bulk_group.same_operations(this_group, include_3d=False)
                )
            if not same_bulk_group:
                raise ValueError("LEEDParametersList: Inconsistent symmetry "
                                 "operations of bulk lattices in the "
                                 "input parameters")

        # Check passed. Set eMax and screenAperture to the max
        for param in all_params:
            param['eMax'] = emax
            param['screenAperture'] = aperture
