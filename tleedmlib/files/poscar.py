# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 17:27:46 2020

@author: Florian Kraushofer

Functions for reading and writing POSCAR files
"""

import logging
import numpy as np

from viperleed.tleedmlib.classes.slab import Slab
from viperleed.tleedmlib.classes.atom import Atom

logger = logging.getLogger("tleedm.files.poscar")


def readPOSCAR(filename='POSCAR'):
    """
    Reads a POSCAR and returns a Slab object with the information.

    Parameters
    ----------
    filename : str, optional
        The file to read from. The default is 'POSCAR'.

    Raises
    ------
    ValueError
        Raised if unit cell vectors are invalid.

    Returns
    -------
    None.

    """

    def initAtomList(slab, poslist):
        """Creates a list of Atom objects based on the data read previously"""
        n = sum(slab.n_per_elem.values())
        slab.atlist = []
        elnum = 0
        atcount = 0
        for i in range(0, n):
            slab.atlist.append(Atom(slab.elements[elnum], poslist[i], i+1,
                                    slab))
            atcount += 1
            if atcount == slab.n_per_elem[slab.elements[elnum]]:
                elnum += 1
                atcount = 0
        slab.getCartesianCoordinates()

    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logger.error("POSCAR not found.")
        raise
    # start reading
    linenum = 1		# iterates the current line being read
    read = True		# set false to stop reading
    sl = Slab()
    c0 = False  # atom close to c=0
    c1 = False  # atom close to c=1
    eps = 1e-4
    poslist = []
    for line in rf:
        if linenum == 1:
            pass
        elif linenum == 2:
            scaling = float(line.split()[0])
        elif linenum <= 5:
            if linenum == 3:
                ucellList = []
            llist = [float(i) for i in line.split()]
            ucellList.append(llist)
            if linenum == 5:
                sl.ucell = scaling * np.transpose(np.array(ucellList))
                sl.ucell_ori = scaling * np.transpose(np.array(ucellList))
                if sl.ucell[2, 0] != 0.0 or sl.ucell[2, 1] != 0.0:
                    if sl.ucell[2, 0] < 1e-4 and sl.ucell[2, 1] < 1e-4:
                        sl.ucell[2, 0] = 0.
                        sl.ucell[2, 1] = 0.
                    else:
                        logger.error(
                            'Unit cell a and b vectors must not have an '
                            'out-of-surface (Z) component!')
                        raise ValueError('Invalid unit cell vectors in POSCAR')
        elif linenum == 6:  # element labels
            sl.elements = [v.capitalize() for v in line.split()]
        elif linenum == 7:
            il = [int(i) for i in line.split()]
            if len(il) != len(sl.elements):
                raise ValueError('Length of element list does not match '
                                 'length of atoms-per-element list')
            for (ind, val) in enumerate(il):
                sl.n_per_elem[sl.elements[ind]] = il[ind]
        elif linenum == 8:
            # may be 'Direct'/'Cartesian', or selective dynamics line
            # check whether POSCAR was pre-processed, ie whether the
            #    'Plane group = ...' comment is already there
            if "Plane group = " in line and "  N" in line:
                if line.split("Plane group = ")[1][0] != "*":
                    sl.preprocessed = True
        elif linenum == 9:
            # this line might already contain coordinates, or not, depending
            # on whether the "Selective dynamics" line was there
            llist = line.split()
            try:
                pos = np.array([float(llist[0]), float(llist[1]),
                                float(llist[2])])
                if abs(pos[2]) < eps:
                    c0 = True
                if abs(pos[2]-1) < eps:
                    c1 = True
                for i in range(0, 2):
                    # in a and b, make sure the values are between 0 and 1
                    pos[i] = pos[i] % 1.0
                poslist.append(pos)
            except Exception:
                logger.debug("POSCAR contains 'Selective dynamics' line, "
                             "skipping line 9")
                # exception was raised because of the 'Selective dynamics'
                # line; this is fine, move on.
        elif read:
            llist = line.split()
            if len(llist) == 0:
                read = False
                logger.debug("POSCAR: Empty line found; stopping position "
                             "readout")
            else:
                pos = np.array([float(llist[0]), float(llist[1]),
                                float(llist[2])])
                if abs(pos[2]) < eps:
                    c0 = True
                if abs(pos[2]-1) < eps:
                    c1 = True
                for i in range(0, 2):
                    # in a and b, make sure the values are between 0 and 1
                    pos[i] = pos[i] % 1.0
                poslist.append(pos)
        linenum += 1
    rf.close()
    initAtomList(sl, poslist)
    # if atoms are very close to c=0 or c=1, move all to avoid edge conflicts
    if c0 and not c1:
        # move up by epsilon
        for at in sl.atlist:
            at.pos[2] = (at.pos[2] + eps) % 1.0
        m = min([at.pos[2] for at in sl.atlist])
        if m < eps:  # move more
            for at in sl.atlist:
                at.pos[2] = (at.pos[2] - m + eps) % 1.0
    elif c1 and not c0:
        # move down by epsilon
        for at in sl.atlist:
            at.pos[2] = (at.pos[2] - eps) % 1.0
        m = max([at.pos[2] for at in sl.atlist])
        if m > 1-eps:  # move more
            for at in sl.atlist:
                at.pos[2] = (at.pos[2] - m + 1 - eps) % 1.0
    elif c0 and c1:
        logger.warning(
            "POSCAR contains atoms close to c=0 and atoms close to c=1. This "
            "cannot be corrected automatically and will likely cause problems "
            "with layer assignment!")
    sl.getCartesianCoordinates()
    return(sl)


def writeCONTCAR(sl, filename='CONTCAR', reorder=False, comments='none',
                 silent=False):
    """
    Takes a Slab object and writes a POSCAR-type file based on the atlist. If a
    POSCAR is in the folder, it will take header data from it, otherwise is
    will fill in a generic header.

    Parameters
    ----------
    sl : Slab
        The Slab object to write.
    filename : str, optional
        The file name to write to. The default is 'CONTCAR'.
    reorder : bool, optional
        If True, atoms will be sorted by element and by z coordinate;
        otherwise, the original order will be used. The default is False.
    comments : str, optional
        Defines whether additional information for each atom should be printed:
        - 'none' writes a normal POSCAR
        - 'all' all annotations
        - 'nodir' all annotations except directions relating to a or b
            (meant to be used with POSCAR_oricell)
        - 'bulk' is like 'none' but writes the space group
            (meant to be used with POSCAR_bulk).
    silent : bool, optional
        If True, will print less to log.

    Returns
    -------
    None.

    """

    if reorder:
        sl.sortByZ()
        sl.sortByEl()
    else:
        # this is the minimum that has to happen
        sl.sortByEl()
    output = '' 	# store output here temporarily before writing it
    # see if a POSCAR is there
    epos = True
    try:
        rf = open('POSCAR', 'r')
    except Exception:
        epos = False
    # header
    if epos:
        # copy the first line:
        output += rf.readline()
        rf.close()
    else:
        # fill dummy header:
        output += 'unknown\n'
    output += '   1.0\n'    # no extra scaling
    # unit cell
    m = np.transpose(sl.ucell)
    for vec in m:
        for val in vec:
            output += '{:22.16f}'.format(val)
        output += '\n'
    # atom types and numbers
    for el in sl.elements:
        output += '{:>5}'.format(el)
    output += '\n'
    for n in sl.n_per_elem.values():
        output += '{:>5}'.format(n)
    output += '\n'
    # header line 'Direct'
    ol = 'Direct'
    if comments != 'none':
        ol = ol.ljust(16)
        ol += 'Plane group = '
        if comments == 'nodir':
            ol += '*'
        if (sl.planegroup in ["pm", "pg", "cm", "rcm", "pmg"]
                and comments != 'nodir'):
            ol += sl.planegroup + str(sl.orisymplane.par)
        else:
            ol += sl.planegroup
        if sl.planegroup != sl.foundplanegroup.split('[')[0]:
            ol += ' (found '
            if '[' in sl.foundplanegroup and comments == 'nodir':
                ol += sl.foundplanegroup.split('[')[0]
            else:
                ol += sl.foundplanegroup
            ol += ')'
        ol = ol.ljust(60)
        if comments == 'all':
            ol += ('N'.rjust(5)+'SiteLabel'.rjust(12)
                   + 'Layer'.rjust(7)+'Linking'.rjust(9)+'FreeDir'.rjust(12))
        elif comments == 'nodir':
            ol += ('N'.rjust(5) + 'SiteLabel'.rjust(12)
                   + 'Layer'.rjust(7) + 'Linking'.rjust(9)
                   + 'FreeDir'.rjust(12) + ': see POSCAR')
    output += ol+'\n'
    # NperElCount = 0
    # lastEl = sl.atlist[0].el
    for i, at in enumerate(sl.atlist):
        # if lastEl != at.el:
        #     lastEl = at.el
        #     NperElCount = 1
        # else:
        #     NperElCount += 1
        ol = ''
        for coord in at.pos:
            ol += '{:20.16f}'.format(coord)
        if comments != 'none' and comments != 'bulk':
            ol += str(i+1).rjust(5)                         # N
            # ol += (at.el+str(NperElCount)).rjust(9)       # NperEl
            ol += at.site.label.rjust(12)                   # SiteLabel
            ol += str(at.layer.num+1).rjust(7)              # Layer
            if len(at.linklist) <= 1:                       # Linking
                ol += 'none'.rjust(9)
            else:
                ol += str((sl.linklists.index(at.linklist))+1).rjust(9)
            if comments != 'nodir':                         # FreeDir
                if at.layer.isBulk:
                    ol += 'bulk'.rjust(12)
                elif type(at.freedir) == int:
                    if at.freedir == 0:
                        ol += 'locked'.rjust(12)
                    else:
                        ol += 'free'.rjust(12)
                else:
                    ol += str(at.freedir).rjust(12)
        output += ol+'\n'
    # write output
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write "+filename)
        raise
    if not silent:
        logger.info("Wrote to "+filename+" successfully")
    return
