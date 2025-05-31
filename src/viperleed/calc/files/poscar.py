"""Functions for reading and writing POSCAR files.

Also defines the POSCARError specific exception, as well as some subclasses.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

from collections import defaultdict
from contextlib import AbstractContextManager
from io import TextIOBase
import logging
from pathlib import Path

import numpy as np

from viperleed.calc.classes import atom as calc_atom
from viperleed.calc.classes import slab as calc_slab

_LOGGER = logging.getLogger(__name__)


class POSCARError(Exception):
    """A generic exception related to a POSCAR file."""


class POSCARSyntaxError(POSCARError):
    """An exceptions for syntax errors in POSCAR files."""


def read(filename='POSCAR'):
    """Return a Slab with the contents of a POSCAR file.

    Parameters
    ----------
    filename : str or Path or TextIOBase, optional
        The file to read from. May be an (open) file object.
        The default is 'POSCAR'.

    Returns
    -------
    slab : Slab
        The slab object. Notice that atoms are normally sorted the              # TODO: discuss if OK
        same way as they were in filename. The only exception is when
        the file contains atoms that are non-element-contiguous, i.e.,
        there are multiple "blocks" with the same chemical element. In
        the latter case, atoms are re-sorted to be element-contiguous,
        in order of appearance of the chemical elements.

    Raises
    ------
    FileNotFoundError
        If filename is not the path to a valid file.
    POSCARSyntaxError
        If unit cell vectors have invalid shape (not enough
        vectors or not enough components).
    POSCARSyntaxError
        If the file does not contain the (normally optional)
        line of chemical-element labels, or if this line is
        inconsistent with the one specifying the atom counts.
    POSCARSyntaxError
        If no (or too few) atomic coordinate is read.
    """
    if isinstance(filename, TextIOBase):
        poscar = POSCARStreamReader(filename)
        return poscar.read()

    with POSCARFileReader(filename) as poscar:
        return poscar.read()


def write(slab, filename='CONTCAR', comments='none', silent=False):
    """Write a POSCAR-style file from a slab.

    If a file named 'POSCAR' exists in the current folder, its first,
    i.e., comment line is copied over.

    Parameters
    ----------
    slab : Slab
        The Slab object to write.
    filename : str or Path or TextIOBase, optional
        The file(name) to write to. The default is 'CONTCAR'.
    comments : str, optional
        Defines whether additional information for each atom
        should be printed:
        - 'none' writes a normal POSCAR (Default)
        - 'all' all annotations
        - 'nodir' all annotations except directions relating
          to a or b (meant to be used with POSCAR_oricell)
        - 'bulk' is like 'none' but writes the space group
          (meant to be used with POSCAR_bulk).
    silent : bool, optional
        If True, will print less to log.

    Raises
    -------
    OSError
        If writing to filename fails.
    """
    slab.sort_by_element()   # This is the minimum that has to happen

    if isinstance(filename, TextIOBase):
        poscar = POSCARStreamWriter(filename, comments=comments)
        poscar.write(slab)
        return

    try:  # pylint: disable=too-many-try-statements
        with POSCARFileWriter(filename, comments=comments) as poscar:
            poscar.write(slab)
    except OSError:
        _LOGGER.error(f'Failed to write {filename}')
        raise
    if not silent:
        _LOGGER.debug(f'Wrote to {filename} successfully')


class POSCARReader:
    """Base class for reading POSCAR structures into a Slab."""

    def __init__(self, source, *, ucell_eps=1e-4):
        """Initialize instance.

        Parameters
        ----------
        source : object
            Positional-only. The source from which the POSCAR
            structure should be read.
        ucell_eps : float, optional
            Cartesian tolerance (Angstrom). All the unit cell
            Cartesian components smaller than ucell_eps will
            be zeroed exactly. Default is 1e-4.
        """
        self._source = source
        self.ucell_eps = ucell_eps

    @property
    def stream(self):
        """Return the stream from which to read structure information."""
        raise NotImplementedError

    def read(self):
        """Return a slab with info read from source."""
        slab = calc_slab.Slab()

        next(self.stream, None)  # Skip first line: comment
        slab.poscar_scaling = float(next(self.stream, '1').split()[0])          # TODO: POSCAR wiki says it's more complex!

        self._read_unit_cell(slab)    # Three lines of unit cell
        element_info = self._read_elements(slab)
        cartesian, _ = self._read_cartesian_and_group(slab)
        positions = self._read_atom_coordinates(slab, cartesian)
        self._make_atoms(slab, positions, *element_info)
        return slab

    @staticmethod
    def _make_atoms(slab, positions, elements, element_counts):
        """Populate slab with atoms from their fractional coordinates.

        Parameters
        ----------
        slab : Slab
            The slab to which atoms are to be added.
        positions : Sequence
            Each element is the (a, b, c) fractional coordinate of
            an atom. There should be exactly sum(element_counts)
            coordinates. See also element_counts.
        elements : Sequence
            Each element is a string, corresponding to the distinct
            chemical species to be added. Chemical species need not
            be unique.
        element_counts : Sequence
            Same number of items as elements. element_counts[i] is
            the number of atoms with element elements[i].

        Raises
        -------
        ValueError
            If the number of atomic coordinates given is not equal
            to the sum of element_counts, if any position is not
            a three-element Sequence, and if element_counts and
            elements have different number of items.
        """
        n_positions, n_coords = np.shape(positions)
        if n_positions != sum(element_counts):
            raise ValueError('Inconsistent number of atomic fractional '
                             f'coordinates ({n_positions}) and chemical '
                             f'element counts ({sum(element_counts)}).')
        if n_coords != 3:
            raise ValueError('Invalid fractional coordinates. Expected '
                             f'3 components, found {n_coords}.')
        if len(elements) != len(element_counts):
            raise ValueError('Inconsistent number of items in elements '
                             f'({len(elements)} and in element_counts '
                             f'({len(element_counts)}')

        positions = iter(positions)
        atoms = defaultdict(list)  # Merge same-element blocks
        for element, n_atoms in zip(elements, element_counts):
            for n_added, fractional_pos in enumerate(positions, start=1):
                atoms[element].append(
                    calc_atom.Atom(element, fractional_pos, 0, slab)
                    )
                if n_added >= n_atoms:
                    break
        slab.atlist.clear()
        slab.atlist.strict = False  # All atoms have the same .num
        slab.atlist.extend(at for el_ats in atoms.values() for at in el_ats)
        slab.update_atom_numbers()  # Also updates atlist map
        slab.atlist.strict = True
        slab.update_element_count()
        slab.update_cartesian_from_fractional()

    def _read_atom_coordinates(self, slab, cartesian):
        """Return an array of fractional atomic coordinates read from file_."""
        n_atoms = sum(slab.n_per_elem.values())
        positions = []
        for line in self.stream:
            coordinates = line.split()
            if not coordinates:
                _LOGGER.debug('POSCAR: Empty line found; '
                              'stopping position readout')
                break
            try:
                positions.append([float(c) for c in coordinates[:3]])
            except ValueError:  # Reached the optional "Lattice velocities"
                _LOGGER.debug(f'POSCAR: no coordinates in {line!r}; '
                              'stopping position readout')
                break
            if len(positions) >= n_atoms:
                break

        if not positions:
            raise POSCARSyntaxError('POSCAR: No atomic coordinates found.')
        if len(positions) < n_atoms:
            raise POSCARSyntaxError(
                'POSCAR: Too few atomic coordinates. Found '
                f'{len(positions)}, expected {n_atoms}'
                )

        # pylint: disable=redefined-variable-type
        # Disabled because used only locally, and clear enough.
        positions = np.array(positions)
        if cartesian:
            uc_inv = np.linalg.inv(slab.ucell.T)
            positions = positions.dot(uc_inv)
        positions[:, :2] %= 1.0  # Collapse coordinates in a and b
        return positions

    def _read_elements(self, slab):
        """Return chemical elements and counts, and update slab."""
        # Line of element labels, then their counts
        elements = [el.capitalize() for el in next(self.stream, '').split()]

        # Element symbols may be terminated by a '/' + some hash; remove it
        # See https://www.vasp.at/forum/viewtopic.php?t=19113
        elements = [el.split('/')[0] for el in elements]

        # Element labels line is optional, but we need it
        try:
            _ = [int(e) for e in elements]
        except ValueError:  # Some non-numeric text. OK
            pass
        else:  # Only numbers. Probably just the element counts
            elements = []
        if not elements:
            raise POSCARSyntaxError(
                'POSCAR: Element labels line not found. This is an optional '
                'line in the POSCAR specification, but ViPErLEED needs it.'
                )
        if any(not el for el in elements):
            raise POSCARSyntaxError(
                'POSCAR: Empty element label found. This is not allowed.'
                )

        element_counts = [int(c) for c in next(self.stream, '').split()]
        if len(element_counts) != len(elements):
            raise POSCARSyntaxError('Length of element list does not match '
                                    'length of atoms-per-element list')
        n_per_elem = defaultdict(int)
        for element, n_atoms in zip(elements, element_counts):
            n_per_elem[element] += n_atoms
        slab.n_per_elem = dict(n_per_elem)
        if len(slab.n_per_elem) != len(elements):
            logging.warning(
                'POSCAR: atoms are not element-contiguous, i.e., '
                'there are multiple blocks with the same chemical-'
                f'species label ({elements}). This structure will '
                'not be preserved, i.e., atoms will be re-sorted '
                f'(to {list(slab.n_per_elem.keys())}) so that only '
                'one block for each element is present.'
                )
        return elements, element_counts

    def _read_cartesian_and_group(self, slab):
        """Read the line containing the "Cartesian"/"Direct" information.

        slab : Slab
            The slab object in which information is to be stored.
            The only attribute modified is .preprocessed, depending
            on whether there is suitable plane-group information.

        Returns
        -------
        cartesian : bool
            Whether the atomic coordinates in the POSCAR file are
            expressed as fractional (False) or Cartesian (True).
        plane_group : str
            The plane group found in the POSCAR file. It will be
            'unknown' if no plane group information exists, or in
            case the information read is not a valid plane group.
        """
        # Skip the optional "S[elective Dynamics]" line
        cartesian_line = next(self.stream, '').lower()
        if cartesian_line.startswith('s'):
            _LOGGER.debug("POSCAR: skipping 'Selective dynamics' line (no. 9)")
            cartesian_line = next(self.stream, '').lower()
        cartesian = cartesian_line.startswith(('c', 'k'))

        # Check whether POSCAR was pre-processed, i.e., whether
        # the 'Plane group = ...' comment is already there
        plane_group = 'unknown'
        if all(s in cartesian_line for s in ('plane group = ', '  n')):
            _, plane_group_str = cartesian_line.split('plane group = ')
            plane_group_str = plane_group_str.split('  n')[0]
            plane_group_str = plane_group_str.split('(')[0].strip()
            if not plane_group_str.startswith('*'):
                slab.preprocessed = True
            plane_group = plane_group_str
        return cartesian, plane_group

    def _read_unit_cell(self, slab):
        """Read unit-cell vectors into slab from the next three lines."""
        ucell = []
        for line in self.stream:
            try:
                ucell.append([float(coord) for coord in line.split()])
            except ValueError:
                break
            if len(ucell) == 3:
                break
        slab.ucell = np.array(ucell).T * slab.poscar_scaling
        if slab.ucell.shape != (3, 3):
            try:
                n_vec, n_comp = slab.ucell.shape
            except ValueError:
                _err = ('Invalid unit-cell vectors. Expected shape (3, 3), '
                        f'found {slab.ucell.shape} instead')
            else:
                _err = ('Invalid unit-cell vectors: not enough vectors '
                        f'({n_vec}/3) or not enough Cartesian components '
                        f'({n_comp}/3).')
            _LOGGER.error(_err)
            raise POSCARSyntaxError(_err)

        slab.ucell_ori = slab.ucell.copy()
        slab.ucell[abs(slab.ucell) < self.ucell_eps] = 0.                       # TODO: was done only in x,y. OK to do all?


class POSCARFileReader(AbstractContextManager, POSCARReader):
    """A context manager for reading POSCAR files into a slab."""

    def __init__(self, filename, *, ucell_eps=1e-4):
        """Initialize instance.

        Parameters
        ----------
        filename : str or Path
            The path to the file to be read.
        ucell_eps : float, optional
            Cartesian tolerance (Angstrom). All the unit cell
            Cartesian components smaller than ucell_eps will
            be zeroed exactly. Default is 1e-4.
        """
        super().__init__(Path(filename), ucell_eps=ucell_eps)
        self._file_object = None

    @property
    def filename(self):
        """Return the name of the file from which to read."""
        return self._source  # From POSCARReader

    @property
    def stream(self):
        """Return the stream from which to read structure information."""
        return self._file_object

    def __enter__(self):
        """Enter context."""
        if not self.filename.is_file():
            _LOGGER.error('POSCAR not found.')
            raise FileNotFoundError(
                f'POSCAR file at {self.filename} not found.'
                )
        self._file_object = self.filename.open('r', encoding='utf-8')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close file object when exiting."""
        try:
            self._file_object.close()
        except AttributeError:
            pass
        return super().__exit__(exc_type, exc_val, exc_tb)

    def read(self):
        """Return a slab with info read from file."""
        if self._file_object is None or self._file_object.closed:
            raise ValueError(
                'I/O operation on closed file. Use\n'
                f'with {type(self).__name__}(filename) as poscar:\n'
                '    slab = poscar.read()'
                )
        return super().read()


class POSCARStreamReader(POSCARReader):
    """Class for reading a POSCAR structure from an open stream."""

    @property
    def stream(self):
        """Return the stream from which to read structure information."""
        return self._source

    def read(self):
        """Return a slab with info read from file."""
        if self.stream.closed:
            raise ValueError('I/O operation on closed stream.')
        return super().read()


class POSCARWriter:
    """Base class for writing a Slab to POSCAR format."""

    def __init__(self, target, *, comments='none', **__kwargs):
        """Initialize instance.

        Parameters
        ----------
        target : object
            Positional-only. An object to be used for writing the
            structure in POSCAR format. Subclasses are more specific
            as to which objects they require.
        comments : str, optional
            Defines whether additional information for each atom
            should be printed:
            - 'none' writes a normal POSCAR (Default)
            - 'all' all annotations
            - 'nodir' all annotations except directions relating
              to a or b (meant to be used with POSCAR_oricell)
            - 'bulk' is like 'none' but writes the space group
              (meant to be used with POSCAR_bulk).
        **__kwargs : object, optional
            Other unused keyword-only arguments.

        Returns
        -------
        None.
        """
        self._target = target
        self.comments = comments
        self.slab = None
        self.header = '\n'  # Default is no header

    @property
    def stream(self):
        """Return the stream to which the POSCAR should be written."""
        raise NotImplementedError

    def write(self, slab):
        """Write the info in slab to the self.stream as a POSCAR."""
        self.stream.write(self.header)
        self.stream.writelines(self._get_scaling(slab))
        self.stream.writelines(self._get_unit_cell(slab))
        self.stream.writelines(self._get_elements_and_numbers(slab))
        self._write_selective_dynamics_line()
        self.stream.writelines(self._get_cartesian_and_group_line(slab))
        self.stream.writelines(self._get_atom_coordinates(slab))

    @staticmethod
    def _get_scaling(slab):
        """Yield the 'poscar-scaling' line."""
        yield f'{slab.poscar_scaling:10.6f}\n'

    @staticmethod
    def _get_unit_cell(slab):
        """Yield one line per unit cell vector."""
        for vector in slab.ucell.T:
            vec = ''.join(f'{v / slab.poscar_scaling:22.16f}' for v in vector)
            yield vec + '\n'

    @staticmethod
    def _get_elements_and_numbers(slab):
        """Yield the two lines of chemical species and their counts."""
        yield ''.join(f'{el:>5}' for el in slab.elements) + '\n'
        yield ''.join(f'{n:>5}' for n in slab.n_per_elem.values()) + '\n'

    def _get_atom_coordinates(self, slab):
        """Yield a line of coordinates and other info for each atom in slab."""
        more_than_one = (link for link in slab.linklists if len(link) > 1)
        linklists = {id(link): i+1 for i, link in enumerate(more_than_one)}
        for i, atom in enumerate(slab):
            line = ''.join(f'{coord:20.16f}' for coord in atom.pos)
            line += self._get_comments_for_atom(atom, i, linklists)
            yield line + '\n'

    def _get_cartesian_and_group_line(self, slab):
        """Yield the line for "Direct", group and other column headers."""
        if self.comments == 'none':
            yield 'Direct\n'
            return

        line = f'{"Direct":<16}Plane group = '
        if self.comments == 'nodir':
            line += '*'

        group = slab.planegroup
        if (slab.planegroup in ["pm", "pg", "cm", "rcm", "pmg"]
                and self.comments != 'nodir'):
            line += f"{group}{slab.orisymplane.par}"
        else:
            line += group

        if slab.planegroup != slab.foundplanegroup.split('[')[0]:
            line += ' (found '
            if '[' in slab.foundplanegroup and self.comments == 'nodir':
                line += slab.foundplanegroup.split('[')[0]
            else:
                line += slab.foundplanegroup
            line += ')'
        line = line.ljust(60)
        if self.comments == 'all':
            line += (f'{"N":>5}{"SiteLabel":>12}{"Layer":>7}'
                     f'{"Linking":>9}{"FreeDir":>12}')
        elif self.comments == 'nodir':
            line += (f'{"N":>5}{"SiteLabel":>12}{"Layer":>7}'
                     f'{"Linking":>9}{"FreeDir":>12}: see POSCAR')
        yield line + '\n'

    def _get_comments_for_atom(self, atom, atom_index, linklists):
        """Return a line of comments for atom at atom_index."""
        if self.comments in ('none', 'bulk'):
            # Only coordinates
            return ''

        line = (f'{atom_index+1:>5}'                  # N
                # f'{atom.el}{NperElCount}'.rjust(9)  # NperEl
                f'{atom.site.label:>12}'              # SiteLabel
                f'{atom.layer.num + 1:>7}')           # Layer
        if len(atom.linklist) <= 1:                   # Linking
            linked = 'none'
        else:
            linked = linklists[id(atom.linklist)]
        line += f'{linked:>9}'
        if self.comments == 'nodir':
            return line
        #                                              # FreeDir
        _free_dir = ""
        if atom.is_bulk:
            _free_dir = 'bulk'
        elif isinstance(atom.freedir, np.ndarray):
            _free_dir = str(atom.freedir)  # has to be made into a string
        elif atom.freedir == 0:
            _free_dir = 'locked'
        else:
            _free_dir = 'free'
        line += f"{_free_dir:>12}"
        return line

    def _write_selective_dynamics_line(self):
        """Write the 'Selective dynamics' line needed for DFT relaxation."""


class POSCARFileWriter(AbstractContextManager, POSCARWriter):
    """A context manager for writing a slab to POSCAR."""

    def __init__(self, filename, *, comments='none'):
        """Initialize instance.

        Parameters
        ----------
        filename : str or Path
            The path to the file to be read.
        comments : str, optional
            Defines whether additional information for each atom
            should be printed:
            - 'none' writes a normal POSCAR (Default)
            - 'all' all annotations
            - 'nodir' all annotations except directions relating
              to a or b (meant to be used with POSCAR_oricell)
            - 'bulk' is like 'none' but writes the space group
              (meant to be used with POSCAR_bulk).

        Returns
        -------
        None.
        """
        super().__init__(Path(filename), comments=comments)
        self._file_object = None

        # Take header line from existing POSCAR, or use a dummy header
        poscar = Path('POSCAR')
        if poscar.is_file():
            with poscar.open('r', encoding='utf-8') as _file:
                self.header = _file.readline()
        else:
            self.header = 'unknown'
        if not self.header.endswith('\n'):
            self.header += '\n'

    @property
    def filename(self):
        """Return the path to the file to be written."""
        return self._target

    @property
    def stream(self):
        """Return the stream to which the POSCAR should be written."""
        return self._file_object

    def __enter__(self):
        """Enter context."""
        self._file_object = self.filename.open('w', encoding='utf-8')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close file object when exiting."""
        try:
            self._file_object.close()
        except AttributeError:
            pass
        super().__exit__(exc_type, exc_val, exc_tb)

    def write(self, slab):
        """Write the info contained in slab to file."""
        if self._file_object is None or self._file_object.closed:
            raise ValueError(
                'I/O operation on closed file. Use\n'
                f'with {type(self).__name__}(filename) as poscar:\n'
                '    poscar.write(slab)'
                )
        super().write(slab)


class POSCARStreamWriter(POSCARWriter):
    """Class for writing a Slab to an open stream in POSCAR format."""

    @property
    def stream(self):
        """Return the stream to which the POSCAR should be written."""
        return self._target

    def write(self, slab):
        """Write the info contained in slab to the current stream."""
        if self.stream.closed:
            raise ValueError('I/O operation on closed stream')
        super().write(slab)


class VASPPOSCARWriter(POSCARStreamWriter):
    """Class for writing a Slab to a POSCAR, ready for VASP relaxation."""

    def __init__(self, stream, *, comments='none', relax_info=None):
        """Initialize instance.

        Parameters
        ----------
        stream : TextIOBase
            The stream to which the POSCAR will be written.
        comments : str, optional
            Defines whether additional information for each atom
            should be printed:
            - 'none' writes a normal POSCAR (Default)
            - 'all' all annotations
            - 'nodir' all annotations except directions relating
              to a or b (meant to be used with POSCAR_oricell)
            - 'bulk' is like 'none' but writes the space group
              (meant to be used with POSCAR_bulk).
            This argument is ignored if relax_info is given.
        relax_info : dict or None, optional
            Contains optional information on how the slab should
            be relaxed. Keys in use:
            - 'above_c' : float, optional
                Only atoms above this c fraction are allowed to move.
                Should be between zero and one (both included). If not
                given, all atoms can move.
            - 'c_only' : bool, optional
                Are movements allowed only along c (True) or also in
                plane (False)? Default is True.

        Raises
        ------
        ValueError
            If the 'above_c' key of relax_info is not a
            floating-point number between 0 and 1 (included)
        """
        super().__init__(stream, comments=comments)
        if relax_info is None:
            relax_info = {}
        self.relax_info = relax_info
        self._check_relax_info()

    def _check_relax_info(self):
        """Raise if self.relax_info contains inappropriate values."""
        above_c = self.relax_info.get('above_c', 0)
        if not 0 <= above_c <= 1:
            raise ValueError('Invalid relax_info: "above_c" '
                             'must be in range [0, 1].')

    def _write_selective_dynamics_line(self):
        """Write the 'Selective dynamics' line needed for DFT relaxation."""
        if self.relax_info:
            self.stream.write('Selective dynamics\n')

    def _get_comments_for_atom(self, atom, atom_index, linklists):
        """Return a line of comments for atom at atom_index."""
        if not self.relax_info:
            return super()._get_comments_for_atom(atom, atom_index, linklists)

        c_cutoff = self.relax_info.get('above_c', -np.inf)
        relax_only_along_c = self.relax_info.get('c_only', True)

        if atom.pos[2] <= c_cutoff:
            return '   F   F   F'
        if relax_only_along_c:
            return '   F   F   T'
        return '   T   T   T'
