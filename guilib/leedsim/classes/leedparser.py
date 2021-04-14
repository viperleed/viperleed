"""Module leedparser of viperleed.guilib.leesdim.classes.

======================================
  ViPErLEED Graphical User Interface
======================================

Created: 2021-04-11
Author: Michele Riva

This module defines the LEEDParser class that is used to read, write
and process the information needed to build a LEEDPattern to/from a
text file. The file is structured like a .ini file, where each
'section' (except 'DEFAULT') contains information about a single
structure. Section 'DEFAULT' has the default values for parameters
that may be missing from the other structures. NB: the contents of
the DEFAULT section may be overwritten with up-to-date defaults,
as defined in this module in the <defaults> dictionary.
"""

import configparser
import ast
import warnings
from collections.abc import Sequence, Mapping
from io import IOBase

import numpy as np

from viperleed import guilib as gl

warnings.simplefilter('always', DeprecationWarning)

defaults = {'screenAperture': '110',
            'bulk3Dsym': 'None',
            'beamIncidence': '(0, 0)'}

param_names = ('eMax', 'surfBasis', 'SUPERLATTICE', 'bulkGroup', 'surfGroup',
               'bulk3Dsym', 'beamIncidence', 'screenAperture', 'name')
param_names_map = {k.lower(): k for k in param_names}

comment_chars = ('#', '!', '%')


class LEEDParser(  # pylint: disable=too-many-ancestors
        configparser.ConfigParser):
    """Parse an input file containing one or more structures."""

    def __init__(self):
        """Initialize LEEDParser instance."""
        super().__init__(comment_prefixes=comment_chars,
                         inline_comment_prefixes=comment_chars,
                         strict=False)
        self['DEFAULT'] = defaults

    @staticmethod
    def write_structures(structures, fname):
        """Write structures to file.

        A DEFAULT section will also be created.

        Parameters
        ----------
        structures : Sequence, dict or ConfigParser
            Structures to be written to file. If a single
            dictionary, or a sequence of dictionaries
            the 'name' keys will be used as section headers,
            falling back on 'S{i+1}' if not present. If a
            dictionary of dictionaries, the outer keys are
            used as section headers.
        fname : str or path
            File to write the dictionary to

        Returns
        -------
        None.
        """
        parser = LEEDParser()
        parser.read_structures(structures)
        parser.write(fname)

    def as_dict(self, structure=None):
        """Return a section as a dictionary.

        All the values that should not be strings are converted.

        Parameters
        ----------
        structure : str
            This parameter is optional only in case there is a
            single structure, in which case that structure will
            be returned.

        Returns
        -------
        dict

        Raises
        ------
        ValueError
            If structure is not passed when multiple are present
        RuntimeError
            If no structure is passed and self is empty
        RuntimeError
            If any of the values in self[structure] does not conform
            to the expected types/shapes/limiting values.
        """
        if structure is None:
            if len(self.sections()) > 1:
                raise ValueError("LEEDParser: need to specify a structure "
                                 "when multiple structures are present")
            if not self.sections():
                raise RuntimeError("LEEDParser: no structures in parser")
            structure = self.sections()[0]

        self.__check_structure(structure)

        out_dict = dict(self[structure])
        for k in out_dict:
            # Process all the values that should not be strings
            if k in ('eMax', 'screenAperture'):
                try:
                    out_dict[k] = self.getfloat(structure, k)
                except ValueError as err:
                    raise RuntimeError(f"LEEDParser: invalid {k} "
                                       f"{out_dict[k]!r} found") from err

            elif k in ('surfBasis', 'SUPERLATTICE'):
                out_dict[k] = self.getarray(structure, k)

            elif k in ('surfGroup', 'bulkGroup'):
                out_dict[k] = self.getgroup(structure, k)

            elif k == 'beamIncidence':
                out_dict[k] = self.getbeamincidence(structure)

            elif k == 'bulk3Dsym' and out_dict[k].lower() == 'none':
                out_dict[k] = None

        if 'name' not in out_dict:
            out_dict['name'] = structure

        return out_dict

    def clear(self):
        """Clear the parser, but keep the defaults."""
        super().clear()
        self['DEFAULT'] = defaults

    def getarray(self, structure, key):
        """Convert a (structure, key) to 2x2 numpy.ndarray.

        Parameters
        ----------
        structure : str
            One of the structures in self
        key : {'surfBasis', 'SUPERLATTICE'}
            Which matrix is to be returned

        Returns
        -------
        numpy.ndarray
            shape (2, 2)

        Raises
        ------
        KeyError
            if key is not one of 'surfBasis', 'SUPERLATTICE'
        RuntimeError
            if not possible to convert to (2, 2) array
        """
        self.__check_structure(structure)
        key = self.optionxform(key)

        if key not in ('surfBasis', 'SUPERLATTICE'):
            raise KeyError(f"LEEDParser: {key} cannot be converted to array")
        self.__check_key(structure, key)

        dtype = float if key == 'surfBasis' else int

        try:
            matrix = gl.string_matrix_to_numpy(self[structure][key],
                                               dtype=dtype,
                                               needs_shape=(2, 2))
        except RuntimeError as err:
            if matrix.shape == (2, 2):
                raise RuntimeError("LEEDParser: could not convert "
                                   f"{key} of {structure} into a "
                                   "numpy.ndarray.") from err
            raise RuntimeError(f"LEEDParser: invalid shape "
                               f"of {key} in {structure}") from err
        return matrix

    def getbeamincidence(self, structure):
        """Get beamIncidence tuple from (structure, key).

        Parameters
        ----------
        structure : str
            One of the structures in self

        Raises
        ------
        ValueError
            If not possible to convert to beamIncidence angles

        Returns
        -------
        theta, phi
            Theta in [0, 90] range, phi in [0,360).

        Raises
        ------
        RuntimeError
            If self[structure]['beamIncidence'] is not an acceptable
            beamIncidence
        """
        key = 'beamIncidence'
        self.__check_key(structure, key)

        angles = ast.literal_eval(self[structure][key])
        if (not isinstance(angles, (list, tuple, np.ndarray))
                or np.shape(angles) != (2,)
                or not -90 <= angles[0] <= 90):
            raise RuntimeError("LEEDParser: invalid beamIncidence in "
                               f"found in structure {structure}: {angles}")
        return gl.conventional_angles(*angles)

    def getgroup(self, structure, key):
        """Convert a (structure, key) to PlaneGroup.

        Parameters
        ----------
        structure : str
            One of the structures in self
        key : {'surfGroup', 'bulkGroup'}
            Which group is to be returned

        Returns
        -------
        viperleed.guilib.PlaneGroup
        """
        self.__check_structure(structure)
        key = self.optionxform(key)

        if key not in ('surfGroup', 'bulkGroup'):
            raise KeyError(f"LEEDParser: {key} cannot be "
                           "converted to PlaneGroup")
        self.__check_key(structure, key)

        return gl.PlaneGroup(self[structure][key])

    def optionxform(self, optionstr):
        """Overridden ConfigParser method.

        Given an 'option' (i.e., key), return the form in which this
        option should be stored in self. Here we convert one of the
        acceptable names to its standard format.

        Parameters
        ----------
        option : str

        Returns
        -------
        str
        """
        return param_names_map[optionstr.lower()]

    def read(self, filenames, encoding=None):
        """Read the contents of a single file into self.

        It also checks whether the file is in the 'old' format, i.e.,
        without 'section headers'. If this is the case, the standard
        section header '[S1]' is added, and the file is overwritten
        to comply with the new format.

        Parameters
        ----------
        filenames : str or list
            File names of the files to be read. If a single string,
            it is treated already as a single filename to read.
        encoding : str, optional
            See ConfigParser.read() documentation

        Returns
        -------
        None.

        Raises
        ------
        DeprecationWarning if the file is in the old format.
        """
        if isinstance(filenames, str):
            filenames = [filenames]
        n_open = 0
        for fname in filenames:
            try:
                open_file = open(fname, 'r')
            except FileNotFoundError:
                pass
            else:
                open_file.close()
                n_open +=1
        if not n_open:
            raise FileNotFoundError("None of the files passed exists!")
        for fname in filenames:
            try:
                super().read(fname, encoding)
            except configparser.MissingSectionHeaderError:
                warnings.warn(f"{filenames} is an old-style input file. "
                              "Will be reformatted and overwritten.",
                              DeprecationWarning)
                with open(fname, 'r') as open_file:
                    lines = ['[S1]\n', *open_file.readlines()]
                with open(fname, 'w') as open_file:
                    open_file.writelines(lines)
        super().read(filenames, encoding)

        self['DEFAULT'] = defaults

    def read_dict(self, dictionary, source='<dict>'):
        """Read a dictionary of structures."""
        # Cast to a proper dict of dict, as e.g.,
        # LEEDParameters does not support key deletion
        dictionary = {k: dict(v) for k, v in dictionary.items()}

        # keep only acceptable keys
        extras = []
        for structure in dictionary.values():
            extras.extend((structure, k)
                          for k in structure
                          if k.lower() not in param_names_map)
        for structure, k in extras:
            structure.pop(k)

        # Now process any numpy array, and can-be-None values.
        # All the others will be correctly converted with str().
        for structure in dictionary.values():
            for key, value in structure.items():
                if (self.optionxform(key) in ('SUPERLATTICE', 'surfBasis')
                        and isinstance(value, np.ndarray)):
                    structure[key] = gl.array2string(value)
                elif self.optionxform(key) == 'bulk3Dsym':
                    structure[key] = str(value)

        super().read_dict(dictionary, source)

    def read_structures(self, structures):
        """Read a Sequence or dict of structures in self.

        Parameters
        ----------
        structures : str, Sequence, dict or ConfigParser
            Structures to be read. If a single dictionary,
            or a sequence of dictionaries the 'name' keys
            will be used as section headers, falling back
            on 'S{i+1}' if not present. If a dictionary
            of dictionaries, the outer keys are used as
            section headers. If a string, it is interpreted
            as a file name to read from.

        Returns
        -------
        None.

        Raises
        ------
        TypeError : when structures is not an acceptabe type
        """
        if isinstance(structures, str):
            # Assume it's a file name
            self.read(structures)
            return

        # When a sequence of dict, or a single dict see if
        # any have a 'name' otherwise use a standard one
        if isinstance(structures, Sequence):
            input_dict = {}
            for i, structure in enumerate(structures):
                key = structure.get('name', f"S{i+1}")
                input_dict[key] = structure
        elif isinstance(structures, (Mapping, configparser.ConfigParser)):
            try:
                structures['bulkGroup']
            except KeyError:
                # It's a dict of dicts
                input_dict = structures
            else:
                # A single structure.
                if 'name' in structures:
                    key = structures['name']
                else:
                    key = 'S1'
                input_dict = {key: structures}
        else:
            raise TypeError("LEEDParser: invalid structures type "
                            f"{type(structures).__name__!r}")
        self.read_dict(input_dict)

    def write(self, fp, space_around_delimiters=True):
        """Save self to file with name fname."""
        file_proxy = fp
        if isinstance(file_proxy, IOBase):
            super().write(file_proxy, space_around_delimiters)
            return
        with open(file_proxy, 'w') as open_file:
            super().write(open_file, space_around_delimiters)

    def __check_key(self, structure, key):
        """Raise KeyError if key is not in self[structure]."""
        if not self.has_option(structure, key):
            raise KeyError(f"LEEDParser: {key} missing in {structure}.")

    def __check_structure(self, structure):
        """Raise ValueError if structure is not in self."""
        if not self.has_section(structure):
            raise ValueError(f"LEEDParser: {structure} not found.")
