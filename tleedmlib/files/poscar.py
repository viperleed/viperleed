# -*- coding: utf-8 -*-
"""Module poscar of viperleed.tleedmlib.files.

Created on Tue Aug 18 17:27:46 2020

@author: Florian Kraushofer
@author: Michele Riva (refactor on 2023-01-16)
@author: Alexander M. Imre

Functions for reading and writing POSCAR files. Also defines the
POSCARError specific exception, as well as some subclasses.
"""

from io import TextIOBase
from collections import defaultdict
from contextlib import AbstractContextManager
import logging
from pathlib import Path

import numpy as np

from viperleed.tleedmlib.classes import atom as tl_atom, slab as tl_slab


_LOGGER = logging.getLogger("tleedm.files.poscar")


class POSCARError(Exception):
    """A generic exception related to a POSCAR file."""


class POSCARSyntaxError(POSCARError):
    """An exceptions for syntax errors in POSCAR files."""


def readPOSCAR(filename='POSCAR'):
    """Return a Slab with the contents of a POSCAR file.

    Parameters
    ----------
    filename : str or Path, optional
        The file to read from. The default is 'POSCAR'.

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
        filepath = filename
    else:
        filepath = Path(filename)
        if not filepath.is_file():
            _LOGGER.error("POSCAR not found.")
            raise FileNotFoundError(f"POSCAR file at {filepath} not found.")

    with POSCARReader(filepath) as poscar:
        return poscar.read()


def writePOSCAR(slab, filename='CONTCAR', reorder=False, comments='none',
                silent=False, relax_info=None):
    """Write a POSCAR-style file from a slab.

    If a file named 'POSCAR' exists in the current folder, its first,
    i.e., comment line is copied over.

    Parameters
    ----------
    slab : Slab
        The Slab object to write.
    filename : str or Path, optional
        The file name to write to. The default is 'CONTCAR'.
    reorder : bool, optional
        If True, atoms will be sorted by z coordinate; in all
        cases atoms are sorted by original-element order. The
        default is False.
    comments : str, optional
        Defines whether additional information for each atom
        should be printed:
        - 'none' writes a normal POSCAR (Default)
        - 'all' all annotations
        - 'nodir' all annotations except directions relating
          to a or b (meant to be used with POSCAR_oricell)
        - 'bulk' is like 'none' but writes the space group
          (meant to be used with POSCAR_bulk).
        - 'relax'
    silent : bool, optional
        If True, will print less to log.

    Raises
    -------
    OSError
        If writing to filename fails.
    """
    if reorder:
        slab.sort_by_z()
    slab.sort_by_element()   # this is the minimum that has to happen

    try:  # pylint: disable=too-many-try-statements
        with POSCARWriter(filename, comments, relax_info) as poscar:
            poscar.write(slab)
    except OSError:
        _LOGGER.error(f"Failed to write {filename}")
        raise
    if not silent:
        _LOGGER.log(1, f"writePOSCAR() comments: {comments}")
        _LOGGER.debug(f"Wrote to {filename} successfully")


def ensure_away_from_c_edges(positions, eps):
    """Ensure positions are not closer than eps to unit-cell edges along c.

    Parameters
    ----------
    positions : Sequence
        Shape (n, 3), n > 0. Fractional coordinates to be checked.
    eps : float
        Fractional distance from 0 and 1 edges. Should be in range
        (0, 1).

    Returns
    -------
    positions : numpy.ndarray
        Shape (n, 3). Modified in place if an array was passed.
        The new fractional coordinates. All elements are such
        that their z fractional coordinate is between eps and
        1 - eps.

    Raises
    ------
    ValueError
        If eps is out of range, or if positions cannot be shifted
        rigidly to fit the (eps, 1 - eps) range.
    """
    if not 0 < eps < 1:
        raise ValueError(f"ensure_away_from_c_edges: Invalid {eps=}. "
                         "Should be between zero and one.")
    positions = np.asarray(positions)
    min_c, max_c = positions[:, 2].min(), positions[:, 2].max()
    if max_c - min_c > 1 - 2*eps:                                               # TODO: couldn't we expand the unit cell?
        raise ValueError(
            "ensure_away_from_c_edges: Cannot shift positions "
            f"to fit between {eps} and {1-eps} along the c axis."
            )
    offset = np.zeros(3)
    if max_c < eps or max_c > 1 - eps:
        offset[2] = 1 - eps - max_c
    elif min_c < eps:
        offset[2] = eps - min_c
    positions += offset
    return positions


class POSCARReader(AbstractContextManager):
    """A context manager for reading POSCAR files into a slab."""

    def __init__(self, filename, ucell_eps=1e-4, min_frac_dist=1e-4):
        """Initialize instance.

        Parameters
        ----------
        filename : str or Path
            The path to the file to be read.
        ucell_eps : float, optional
            Cartesian tolerance (Angstrom). All the unit cell
            Cartesian components smaller than ucell_eps will
            be zeroed exactly. Default is 1e-4.
        min_frac_dist : float, optional
            Smallest fractional distance from the c=0 and c=1
            edges of the unit cell. If atoms are closer, they
            are shifted, if possible.
        """
        super().__init__()
        if isinstance(filename, TextIOBase):
            self.filename = filename
        else:
            self.filename = Path(filename)
        self.file_object = None
        self.ucell_eps = ucell_eps
        self.min_frac_dist = min_frac_dist

    def __enter__(self):
        """Enter context."""
        if not isinstance(self.filename, TextIOBase):
            self.file_object = self.filename.open('r', encoding='utf-8')
        else:
            self.file_object = self.filename
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close file object when exiting."""
        if not isinstance(self.filename, TextIOBase):
            self.file_object.close()
            super().__exit__(exc_type, exc_val, exc_tb)

    def read(self):
        """Return a slab with info read from file."""
        if self.file_object is None or self.file_object.closed:
            raise ValueError(
                "I/O operation on closed file. Use\n"
                "with POSCARReader(filename) as poscar:\n"
                "    slab = poscar.read()"
                )

        slab = tl_slab.Slab()

        next(self.file_object, None)  # Skip first line: comment
        slab.poscar_scaling = float(next(self.file_object, '1').split()[0])     # TODO: POSCAR wiki says it's more complex!

        self._read_unit_cell(slab)    # Three lines of unit cell
        element_info = self._read_elements(slab)
        cartesian, _ = self._read_cartesian_and_group(slab)
        positions = self._read_atom_coordinates(slab, cartesian)

        try:
            # Move atoms away from edges in z to avoid conflicts
            positions = ensure_away_from_c_edges(positions, self.min_frac_dist)
        except ValueError:
            _LOGGER.warning(
                "POSCAR contains atoms close to c=0 and atoms close "
                "to c=1. This cannot be corrected automatically and "
                "will likely cause problems with layer assignment!"
                )

        # And add them to the slab
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
            raise ValueError("Inconsistent number of atomic fractional "
                             f"coordinates ({n_positions}) and chemical "
                             f"element counts ({sum(element_counts)}).")
        if n_coords != 3:
            raise ValueError("Invalid fractional coordinates. Expected "
                             f"3 components, found {n_coords}.")
        if len(elements) != len(element_counts):
            raise ValueError("Inconsistent number of items in elements "
                             f"({len(elements)} and in element_counts "
                             f"({len(element_counts)}")

        positions = iter(positions)
        atoms = defaultdict(list)  # Merge same-element blocks
        for element, n_atoms in zip(elements, element_counts):
            for n_added, fractional_pos in enumerate(positions, start=1):
                atoms[element].append(
                    tl_atom.Atom(element, fractional_pos, 0, slab)
                    )
                if n_added >= n_atoms:
                    break
        slab.atlist = [at for el_ats in atoms.values() for at in el_ats]
        slab.updateAtomNumbers()
        slab.updateElementCount()
        slab.getCartesianCoordinates()

    def _read_atom_coordinates(self, slab, cartesian):
        """Return an array of fractional atomic coordinates read from file_."""
        n_atoms = sum(slab.n_per_elem.values())
        positions = []
        for line in self.file_object:
            coordinates = line.split()
            if not coordinates:
                _LOGGER.debug("POSCAR: Empty line found; "
                             "stopping position readout")
                break
            try:
                positions.append([float(c) for c in coordinates[:3]])
            except ValueError:  # Reached the optional "Lattice velocities"
                _LOGGER.debug(f"POSCAR: no coordinates in {line!r}; "
                             "stopping position readout")
                break
            if len(positions) >= n_atoms:
                break

        if not positions:
            raise POSCARSyntaxError("POSCAR: No atomic coordinates found.")
        if len(positions) < n_atoms:
            raise POSCARSyntaxError(
                "POSCAR: Too few atomic coordinates. Found "
                f"{len(positions)}, expected {n_atoms}"
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
        file_ = self.file_object

        # Line of element labels, then their counts
        elements = [el.capitalize() for el in next(file_, '').split()]

        # Element labels line is optional, but we need it
        try:
            _ = [int(e) for e in elements]
        except ValueError:  # Some non-numeric text. OK
            pass
        else:  # Only numbers. Probably just the element counts
            elements = []
        if not elements:
            raise POSCARSyntaxError(
                "POSCAR: Element labels line not found. This is an optional "
                "line in the POSCAR specification, but ViPErLEED needs it."
                )

        element_counts = [int(c) for c in next(file_, '').split()]
        if len(element_counts) != len(elements):
            raise POSCARSyntaxError('Length of element list does not match '
                                    'length of atoms-per-element list')
        n_per_elem = defaultdict(int)
        for element, n_atoms in zip(elements, element_counts):
            n_per_elem[element] += n_atoms
        slab.n_per_elem = dict(n_per_elem)
        if len(slab.n_per_elem) != len(elements):
            logging.warning(
                "POSCAR: atoms are not element-contiguous, i.e., "
                "there are multiple blocks with the same chemical-"
                f"species label ({elements}). This structure will "
                "not be preserved, i.e., atoms will be re-sorted "
                f"(to {list(slab.n_per_elem.keys())}) so that only "
                "one block for each element is present."
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
        cartesian_line = next(self.file_object, '').lower()
        if cartesian_line.startswith('s'):
            _LOGGER.debug("POSCAR: skipping 'Selective dynamics' line (no. 9)")
            cartesian_line = next(self.file_object, '').lower()
        cartesian = cartesian_line.startswith(('c', 'k'))

        # Check whether POSCAR was pre-processed, i.e., whether
        # the 'Plane group = ...' comment is already there
        plane_group = 'unknown'
        if all(s in cartesian_line for s in ("plane group = ", "  n")):
            _, plane_group_str = cartesian_line.split("plane group = ")
            plane_group_str = plane_group_str.split("  n")[0]
            plane_group_str = plane_group_str.split('(')[0].strip()
            if not plane_group_str.startswith("*"):
                slab.preprocessed = True
            plane_group = plane_group_str
        return cartesian, plane_group

    def _read_unit_cell(self, slab):
        """Read unit-cell vectors into slab from the next three lines."""
        ucell = []
        for line in self.file_object:
            try:
                ucell.append([float(coord) for coord in line.split()])
            except ValueError:
                break
            if len(ucell) == 3:
                break
        slab.ucell = np.array(ucell).T * slab.poscar_scaling
        if slab.ucell.shape != (3, 3):
            n_vec, n_comp = slab.ucell.shape
            _err = ("Invalid unit-cell vectors: not enough vectors "
                    f"({n_vec}/3) or not enough Cartesian components "
                    f"({n_comp}/3).")
            _LOGGER.error(_err)
            raise POSCARSyntaxError(_err)

        slab.ucell_ori = slab.ucell.copy()
        slab.ucell[abs(slab.ucell) < self.ucell_eps] = 0.                       # TODO: was done only in x,y. OK to do all?


class POSCARWriter(AbstractContextManager):
    """A context manager for writing a slab to POSCAR."""

    def __init__(self, filename, comments='none', relax_info=None):
        """Initialize instance."""
        super().__init__()
        if isinstance(filename, TextIOBase):
            self.filename = filename
        else:
            self.filename = Path(filename)
        self.comments = comments
        self.file_object = None
        self.slab = None
        self.relax_info = relax_info
        if self.comments == 'relax':
            self._check_relax_info()

        if not isinstance(self.filename, TextIOBase):
            # Take header line from existing POSCAR, or use a dummy header
            poscar = Path('POSCAR')
            if poscar.is_file():
                with poscar.open('r', encoding='utf-8') as _file:
                    self.header = _file.readline()
            else:
                self.header = 'unknown'
            if not self.header.endswith('\n'):
                self.header += '\n'
        else:
            self.header = '\n'

    def __enter__(self):
        """Enter context."""
        if isinstance(self.filename, TextIOBase):
            self.file_object = self.filename
        else:
            self.file_object = self.filename.open('w', encoding='utf-8')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close file object when exiting."""
        if isinstance(self.filename, TextIOBase):
            return
        self.file_object.close()
        super().__exit__(exc_type, exc_val, exc_tb)

    def write(self, slab):
        """Write the info contained in slab to file."""
        if self.file_object is None or self.file_object.closed:
            raise ValueError(
                "I/O operation on closed file. Use\n"
                "with POSCARWriter(filename) as poscar:\n"
                "    poscar.write(slab)"
                )

        self.file_object.writelines(self.header)
        self.file_object.writelines(self._get_scaling(slab))
        self.file_object.writelines(self._get_unit_cell(slab))
        self.file_object.writelines(self._get_elements_and_numbers(slab))
        if self.comments == 'relax':
            self.file_object.writelines("Selective dynamics\n")
        self.file_object.writelines(self._get_cartesian_and_group_line(slab))
        self.file_object.writelines(self._get_atom_coordinates(slab))

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

    def _check_relax_info(self):
        if self.relax_info is None and self.comments == 'relax':
            raise ValueError("Specified comments='relax' but missing "
                             "relaxation_info.")
        if self.relax_info is not None and self.comments != 'relax':
            _LOGGER.warning("Specified relaxation_info but comments is not"
                           "'relax'.")
        if (self.relax_info is not None
            and (self.relax_info['above_c'] <= 0.0
            or self.relax_info['above_c'] >= 1.0)):
            raise ValueError("Invalid relaxation_info: above_c must be in range"
                             "[0,1].")
        self.above_c = self.relax_info['above_c']
        self.relax_c_only = self.relax_info['c_only']

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

    def _get_atom_coordinates(self, slab):
        """Yield a line of coordinates and other info for each atom in slab."""
        more_than_one = (link for link in slab.linklists if len(link) > 1)
        linklists = {id(link): i+1 for i, link in enumerate(more_than_one)}
        for i, atom in enumerate(slab.atlist):
            line = ''.join(f"{coord:20.16f}" for coord in atom.pos)
            if self.comments in ('none', 'bulk'):
                # Only coordinates
                yield line + '\n'
                continue
            if self.comments == 'relax':
                if atom.pos[2] <= self.above_c:
                    relax_pos_flags = "   F   F   F"
                elif atom.pos[2] > self.above_c and self.relax_c_only:
                    relax_pos_flags = "   F   F   T"
                else:
                    relax_pos_flags = "   T   T   T"
                yield line + relax_pos_flags + '\n'
                continue

            line += (f"{i+1:>5}"                           # N
                     # f"{atom.el}{NperElCount}".rjust(9)  # NperEl
                     f"{atom.site.label:>12}"              # SiteLabel
                     f"{atom.layer.num + 1:>7}")           # Layer
            if len(atom.linklist) <= 1:                    # Linking
                linked = 'none'
            else:
                linked = linklists[id(atom.linklist)]
            line += f"{linked:>9}"
            if self.comments == 'nodir':
                yield line + '\n'
                continue
            #                                              # FreeDir
            _free_dir = ""
            if atom.layer.isBulk:
                _free_dir = 'bulk'
            elif isinstance(atom.freedir, np.ndarray):
                _free_dir = str(atom.freedir)  # has to be made into a string
            elif atom.freedir == 0:
                _free_dir = 'locked'
            else:
                _free_dir = 'free'
            line += f"{_free_dir:>12}"
            yield line + '\n'