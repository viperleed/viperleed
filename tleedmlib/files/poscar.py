# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 17:27:46 2020

@author: Florian Kraushofer

Functions for reading and writing POSCAR files
"""

import logging
import numpy as np

from tleedmlib.classes.slab import Slab

logger = logging.getLogger("tleedm.files.poscar")

def readPOSCAR(filename='POSCAR'):
    """Reads a POSCAR and returns a Slab object with the information."""
    # open input file
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logger.error("POSCAR not found.")
        raise

    #read POSCAR:
    linenum = 1		# iterates the current line being read
    read = True		# set false to stop reading
    sl = Slab()
    c0 = False  # atom close to c=0
    c1 = False  # atom close to c=1
    eps = 1e-4
    for line in rf:
        if linenum == 1:
            pass
        elif linenum == 2:
            scaling = float(line.split()[0])
        elif linenum <= 5:
            if linenum == 3: ucellList = []
            llist = [float(i) for i in line.split()]
            ucellList.append(llist)
            if linenum == 5:
                sl.ucell = scaling * np.transpose(np.array(ucellList))
                sl.uCellOri = scaling * np.transpose(np.array(ucellList))
                if sl.ucell[2,0] != 0.0 or sl.ucell[2,1] != 0.0:
                    if sl.ucell[2,0] < 1e-4 and sl.ucell[2,1] < 1e-4:
                        sl.ucell[2,0] = 0.
                        sl.ucell[2,1] = 0.
                    else:
                        logger.error('ERROR: Unit cell a and b vectors must '
                                'not have an out-of-surface (Z) component!')
                        raise
        elif linenum == 6:
            sl.elements = line.split()       		# element labels
            sl.oriels = sl.elements[:]              # copy
            sl.nelem = len(sl.elements)			# number of different elements
        elif linenum == 7:
            sl.nperelem = [int(i) for i in line.split()]
            if len(sl.nperelem) != sl.nelem:
                logger.warning('\nPOSSIBLE PROBLEM: lenght of element list '
                  'does not match length of atoms-per-element list\n')
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
                for i in range(0,2):
                    # in a and b, make sure the values are between 0 and 1
                    pos[i] = pos[i] % 1.0
                sl.atpos.append(pos)
            except:
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
                for i in range(0,2):
                    # in a and b, make sure the values are between 0 and 1
                    pos[i] = pos[i] % 1.0
                sl.atpos.append(pos)
        linenum += 1
    rf.close()
    sl.initAtomList()
    # if atoms are very close to c=0 or c=1, move all to avoid edge conflicts
    if c0 == True and c1 == False:
        # move up by epsilon
        for at in sl.atlist:
            at.pos[2] = (at.pos[2] + eps) % 1.0
        m = min([at.pos[2] for at in sl.atlist])
        if m < eps: # move more
            for at in sl.atlist:
                at.pos[2] = (at.pos[2] - m + eps) % 1.0
    elif c1 == True and c0 == False:
        # move down by epsilon
        for at in sl.atlist:
            at.pos[2] = (at.pos[2] - eps) % 1.0
        m = max([at.pos[2] for at in sl.atlist])
        if m > 1-eps: # move more
            for at in sl.atlist:
                at.pos[2] = (at.pos[2] - m + 1 - eps) % 1.0
    elif c0 == True and c1 == True:
        logger.warning("POSCAR contains atoms close to c=0 and atoms close "
                "to c=1. This cannot be corrected automatically and will "
                "likely cause problems with layer assignment!")
    sl.getCartesianCoordinates()
    return(sl)

def writeCONTCAR(sl, filename='CONTCAR', reorder=False, comments='none',
                 silent=False):
    """Takes a Slab object and writes a CONTCAR based on the atlist. If a 
    POSCAR is in the folder, it will take header data from it, otherwise is 
    will fill in a generic header. Defining a filename other than CONTCAR is 
    optional. If the reorder tag is set as True, atoms will be sorted by 
    element and by z coordinate; otherwise, the original order will be used. 
    The 'comments' tag defines whether additional information for each atom 
    should be printed; 'none' writes a normal POSCAR, 'all' all annotations, 
    'nodir' all annotations except directions relating to a or b 
    (meant to be used with POSCAR_oricell), 'bulk' is like 'none' but 
    writes the space group (meant to be used with POSCAR_bulk)."""
    if reorder:
        sl.sortByZ()
        sl.sortByEl()
    else:
        sl.sortByEl()   #this is the minimum that has to happen, otherwise the 
                        #  output will be buggy
    output = '' 	# store output here temporarily before writing it
    # see if a POSCAR is there
    epos = True
    try:
        rf = open('POSCAR', 'r')
    except:
        epos = False
    #header
    if epos:
        #copy the first line:
        output += rf.readline()
        rf.close()
    else:
        #fill dummy header:
        output += 'unknown\n'
    output += '   1.0\n'    #scaling will be 1.0, since we apply the scaling 
                            #  to the unit cell
    #unit cell
    m = np.transpose(sl.ucell)
    for vec in m:
        for val in vec:
            s = '%.16f'%val
            output += '{:>22}'.format(s)
        output += '\n'
    #atom types and numbers
    for el in sl.oriels:
        output += '{:>5}'.format(el)
    output += '\n'
    for n in sl.nperelem:
        output += '{:>5}'.format(n)
    output += '\n'
    #header line 'Direct'
    ol = 'Direct'
    if comments != 'none':
        ol = ol.ljust(16)
        ol += 'Plane group = '
        if comments == 'nodir': ol += '*'
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
                   +'Layer'.rjust(7)+'Linking'.rjust(9)+'FreeDir'.rjust(12))
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
            s = '%.16f'%coord
            ol += '{:>20}'.format(s)
        if comments != 'none' and comments != 'bulk':
            ol += str(i+1).rjust(5)                         #N
            # ol += (at.el+str(NperElCount)).rjust(9)         #NperEl
            ol += at.site.label.rjust(12)                   #SiteLabel
            ol += str(at.layer.num+1).rjust(7)              #Layer
            if len(at.linklist) <= 1:                       #Linking
                ol += 'none'.rjust(9)
            else:
                ol += str((sl.linklists.index(at.linklist))+1).rjust(9)
            if comments != 'nodir':                         #FreeDir
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
    #write output
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    if not silent:
        logger.info("Wrote to "+filename+" successfully")
    return 0