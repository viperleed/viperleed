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

from collections.abc import MutableMapping, MutableSequence, Sequence
from configparser import ConfigParser
from warnings import warn as warning   # eventually will replace with logging

import numpy as np

from viperleed.guilib.classes import planegroup
from viperleed.guilib.classes.lattice2d import Lattice2D
from viperleed.guilib.helpers import conventional_angles
from viperleed.guilib.helpers import remove_duplicates
from viperleed.guilib.helpers import single_spaces_only
from viperleed.guilib.helpers import two_by_two_array_to_tuple
from viperleed.guilib.leedsim.classes.leedparser import LEEDParser

non_string_keys = ('emax', 'surfbasis', 'superlattice', 'bulkgroup',
                   'surfgroup', 'beamincidence', 'screenaperture')


# TODO: may need to pull the old LEEDParameters from master

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
            self.__from_parser(data)

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
        TypeError
            If k is not a string.
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
        txt = single_spaces_only(repr(self.__data)).replace('\n', '')
        txt = txt.replace('[ ', '[').replace(' ]', ']').replace(' ,', ',')
        return f"LEEDParameters({txt})"

    # TODO: pylint too many return (11/6)
    # TODO: pylint too complex mccabe=14
    def __eq__(self, other):
        """Return whether self is equal to other.

        Instances are considered equal if they produce
        the very same LEED pattern, i.e., they have:
        - the same 'name'
        - the same 'beamIncidence'
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
        equal : bool or NotImplementedType
            Returns True if LEEDParameters are equal,
            NotImplemented otherwise
        """
        if not isinstance(other, LEEDParameters):
            # Try to see if other can be used as an argument to the
            # LEEDParameters constructor
            try:
                other = LEEDParameters(other)
            except (TypeError, NameError, ValueError):
                # If instantiation raises one of the known exceptions
                return NotImplemented

        # (1) first a couple of basic checks:
        #     identity and difference of names (if given)
        if self is other:
            return True
        if self['name'] != other['name']:
            return NotImplemented

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
        if not int_transforms:
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

    def __from_parser(self, parser):
        """Create LEEDParameters from a configuration file.

        Parameters
        ----------
        parser : ConfigParser
            The parser for the configuration file from which to
            initialize. It should contain only one section.

        Returns
        -------
        None.

        Raises
        ------
        RuntimeError
            When parser contains multiple sections.
        """
        if len(parser.sections()) > 1:
            raise RuntimeError("LEEDParameters: multiple structures "
                               "found. Use LEEDParametersList instead")
        leed_parser = LEEDParser()
        leed_parser.read_structures(parser)
        self.__from_dict(leed_parser.as_dict())

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
        # Check whether any of the values that are not supposed to
        # be strings, are in fact strings.  If they are, use a
        # LEEDParser to convert them, and write them back into data
        if any(isinstance(v, str)
               for k, v in data.items()
               if k in non_string_keys):
            name = data.get('name', 'S1')
            parser = LEEDParser()
            parser.read_dict({name: data})
            for k, value in parser.as_dict().items():
                data[k.lower()] = value

    # TODO: pylint too complex mccabe=11
    @staticmethod
    def _check_acceptable(data_dict):
        """Check that the dict input is acceptable.

        Parameters
        ----------
        data_dict : dict or LEEDParameters
            The data to be checked

        Raises
        ------
        ValueError
            if data_dict contains unacceptable data.
        """
        for key, value in data_dict.items():
            if key.lower() == 'emax' and value <= 0:
                raise ValueError("Maximum LEED energy should be positive.")

            if (key.lower() == 'screenaperture'
                    and (value <= 0 or value > 180)):
                raise ValueError("screenAperture should be between "
                                 f"0 and 180. Found {value} instead.")

            if key.lower() in ('surfbasis', 'superlattice'):
                try:
                    value = np.asarray(value)
                except TypeError as err:
                    raise ValueError(f"Invalid {key}"
                                     "matrix input") from err
                if value.shape != (2, 2):
                    raise ValueError("Invalid shape "
                                     f"{value.shape} for {key}. "
                                     "Expected (2, 2)")
                data_dict[key] = value

            if (key.lower() == 'superlattice'
                    and any(np.abs(np.asarray(value) % 1).ravel() > 1e-4)):
                raise ValueError("SUPERLATTICE must be an "
                                 "integer-valued matrix")

            if key.lower() == 'beamincidence':
                if (not isinstance(value, (list, tuple, np.ndarray))
                    or np.shape(value) != (2,)
                        or not -90 <= value[0] <= 90):
                    raise ValueError(f"Invalid beamIncidence: {value}")
                data_dict[key] = conventional_angles(*value)

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
                            Lattice2D(self['surfBasis']).reciprocal_basis)
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
        dummy_bulk = Lattice2D(self['bulkReciprocalBasis'],
                               space='reciprocal',
                               group=self['bulkGroup'])
        bulk_shape = dummy_bulk.cell_shape

        # And apply the bulk operations to the bulk group
        self['bulkGroup'].set_screws_glides(self['bulk3Dsym'], bulk_shape)


# Probably bug in pylint. Should not complain about too-many-ancestors
# when subclassing an collections.abc
# Should be solved in pylint 2.9.0, but it isn't as of pylint 2.9.3
# pylint: disable=too-many-ancestors
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
        data : Sequence
            Each element is either dict, ConfigParser,
            LEEDParameters, or LEEDParametersList. Each element
            will be checked for compatibility, i.e., all should
            have the same bulk. A single LEEDParametersList is
            also an acceptable Sequence.
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

    def __eq__(self, other):
        """Return whether self is equal to other."""
        if not isinstance(other, (LEEDParametersList, Sequence)):
            return NotImplemented
        # if isinstance(other, LEEDParametersList):
            # return self.__list == other.__list
        return self.__list == other

    def __getitem__(self, index):
        """Get item at index."""
        # Do this type-preserving. If a slice, should
        # return same type, if index, just an element
        if isinstance(index, slice):
            return LEEDParametersList(self.__list[index])
        return self.__list[index]

    def __setitem__(self, idx, data):
        """Set item at idx to data.

        The items are set only if they are acceptable LEEDParameter(s),
        given the data already present in self. This behaves somewhat
        differently than normal list.__setitem__, as it prevents length
        changes as well as increases in the number of duplicates.

        Parameters
        ----------
        idx : int or slice
        data : dict, ConfigParser, LEEDParameters (or
               iterables of the same if idx is a slice)

        Raises
        ------
        ValueError
            * When a slice is passed, and the data passed contains
              duplicates
            * When the data passed would increase the number of
              duplicates in self
            * When replacement of items would change len(self)
        """
        # data = self._process_input_data(data)
        # At first, just check whether there are duplicates within
        # the data themselves
        new_data = self._process_input_data(data, keep_duplicates=True)
        new_data = remove_duplicates(new_data)
        if isinstance(idx, slice) and len(new_data) < len(data):
            raise ValueError("Data contains duplicates")

        # Make sure we're replacing the right number of entries
        if isinstance(idx, slice):
            to_replace = self[idx]
        else:
            to_replace = [self[idx]]
        if len(new_data) != len(to_replace):
            raise ValueError("Inconsistent number of replaced LEEDParameters: "
                             f"Cannot replace {len(to_replace)} "
                             f"elements with {len(new_data)}")

        # Now see if the data contains duplicates of elements in self:
        # (1) We can't replace something with a parameter that
        #     is already in self and won't be replaced, as
        #     this would increase the number of duplicates
        if any(d in self and d not in to_replace for d in new_data):
            raise ValueError("Not possible to add duplicates "
                             "of already stored data")

        # (2) We can still accept the data if each
        #     duplicate will replace its own copy
        duplicated = [d in to_replace for d in new_data]
        if (any(duplicated)
            and any(d_new != d_old if new_is_duplicate else False
                    for d_new, d_old, new_is_duplicate in zip(new_data,
                                                              to_replace,
                                                              duplicated))):
            raise ValueError("Not possible to add duplicates "
                             "of already stored data")
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

    # Disable pylint check. It's not great to make the signature
    # different than the one of list.insert, but it seems the best
    # option for this very specific case of 'list without duplicates'
    # pylint: disable=arguments-differ
    def insert(self, idx, data, keep_duplicates=False):
        """Insert data at index.

        The data is inserted only if it is compatible with the
        data already present in self.

        Parameters
        ----------
        idx : int
            Index at which to insert
        data : dict, ConfigParser, LEEDParameters
            The data to be inserted. Will be checked for compatibility.
        keep_duplicates : bool, optional
            If False, data is inserted only if it is not a duplicate
            of elements already present in self. Default is False.

        Raises
        ------
        ValueError
            When trying to insert multiple LEEDParameters at one index,
            as this has the potential of making the LEEDParametersList
            multi-dimensional
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

    # See comment at .insert()
    # pylint: disable=arguments-differ
    def append(self, data, keep_duplicates=False):
        """Insert a single element at the end.

        The data is appended only if it is compatible with the
        data already present in self.

        Parameters
        ----------
        data : dict, ConfigParser, LEEDParameters
            The dictionary-like LEED parameters to be inserted.
            Will be checked for compatibility.
        keep_duplicates : bool, optional
            If False, data is appended only if it is not a duplicate
            of elements already present in self. Default is False.
        """
        self.insert(len(self), data, keep_duplicates)
    # pylint: enable=arguments-differ

    # See comment at .insert()
    # pylint: disable=arguments-differ
    def extend(self, iterable, keep_duplicates=False):
        """Insert a elements from an iterable at the end.

        Only those entries that are compatible with the
        data already present in self are appended.

        Parameters
        ----------
        iterable : iterable or LEEDParametersList
            Each element is either dict, ConfigParser, or LEEDParameters
            The data to be inserted. Each element will be checked for
            compatibility.
        keep_duplicates : bool, optional
            If False, append only those elements in data that are not
            a duplicate of elements already present in self. Duplicates
            are tested with the equality method for LEEDParameters.
            Default is False.
        """
        for elem in iterable:
            self.append(elem, keep_duplicates)
    # pylint: enable=arguments-differ

    def set_beam_incidence(self, theta=None, phi=None):
        """Set the 'beamIncidence' key to the given angles.

        The same angles will be set for all the entries in self.

        Parameters
        ----------
        theta : float, optional
            Polar incidence angle. Should be in the [-90, 90] range.
            Normal incidence is theta=0. No modification if not given
            or if None. Default is None
        phi : float, optional
            Azimuthal incidence angle. Positive counterclockwise,
            measured with respect to the positive x axis of the
            Cartesian coordinates in which the real-space bases
            are expressed. Will be taken modulo 360. Not modified
            if not given or if None. Default is None.

        Raises
        ------
        TypeError
            If either theta or phi are given, but cannot
            be converted to floating point numbers.
        """
        if theta is None and phi is None:
            return
        old_angles = self[0]['beamIncidence']
        if theta is None:
            theta = old_angles[0]
        elif phi is None:
            phi = old_angles[1]
        # Disable pylint check as this is done on purpose:
        # Both conversions to float must be possible.
        try:  # pylint: disable=too-many-try-statements
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
        """Process the given data.

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

        Returns
        -------
        list
            The checked and processed input. Each element is a
            LEEDParameters instance.
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
        data : dict, ConfigParser, LEEDParameters, or Sequence
            The data to be converted to a list of LEEDParameters.
            When a Sequence, each element should be dict, ConfigParser,
            or LEEDParameters. A LEEDParametersList is also acceptable.

        Returns
        -------
        list
            Each element is a LEEDParameters whatever the input.

        Raises
        ------
        TypeError
            If any of the data is not one of the acceptable types.
        """
        if not data:
            ret = []
        elif isinstance(data, ConfigParser):
            leed_parser = LEEDParser()
            leed_parser.read_dict(data)
            ret = [LEEDParameters(leed_parser.as_dict(section))
                   for section in data.sections()]
        elif isinstance(data, dict):
            # delegate to LEEDParameters the check on the dictionary
            ret = [LEEDParameters(data)]
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
        if not all_params:
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
                raise ValueError("Inconsistent bulk bases "
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
                raise ValueError("Inconsistent symmetry "
                                 "operations of bulk lattices "
                                 "in the input parameters")

        # Check passed. Set eMax and screenAperture to the max
        for param in all_params:
            param['eMax'] = emax
            param['screenAperture'] = aperture
