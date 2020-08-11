# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Contains functions for reading the input files of the TensErLEED interface
"""

import logging
import numpy as np
import re
import os
import copy

import fortranformat as ff
import tleedmlib as tl
from tleedmlib import DEFAULT

logger = logging.getLogger("tleedm.read")

def collectFIN():
    """Combines AUXLATGEO, _BEAMLIST, AUXNONSTRUCT, _PHASESHIFTS, AUXBEAMS
    and AUXGEO into one string (input for refcalc), which it returns."""
    filenames = ["AUXLATGEO", "_BEAMLIST", "AUXNONSTRUCT", "_PHASESHIFTS",
                 "AUXBEAMS", "AUXGEO"]
    fin = ""
    for fn in filenames:
        with open(fn, "r") as rf:
            fin += rf.read()
        if fin[-1] != "\n":
            fin += "\n"
    return fin

def readBEAMLIST(filename="_BEAMLIST"):
    """Reads the _BEAMLIST file and returns the contents as a list of the
    lines as strings."""
    beamlist = []
    try:
        with open(filename, "r") as rf:
            beamlist = rf.readlines()
        logger.debug("_BEAMLIST file was read successfully")
    except:
        logger.error("Error opening _BEAMLIST file.")
        raise
    return beamlist

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
    sl = tl.Slab()
    c0 = False  # atom close to c=0
    c1 = False  # atom close to c=1
    eps = 1e-4
    for line in rf:
        if linenum == 1:
            pass
        elif linenum == 2:
            scaling = float(tl.linelist(line)[0])
        elif linenum <= 5:
            if linenum == 3: ucellList = []
            llist = [float(i) for i in tl.linelist(line)]
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
            sl.elements = tl.linelist(line)			# element labels
            sl.oriels = sl.elements[:]              # copy
            sl.nelem = len(sl.elements)			# number of different elements
        elif linenum == 7:
            sl.nperelem = [int(i) for i in tl.linelist(line)]
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
            llist = tl.linelist(line)
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
            llist = tl.linelist(line)
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

def updatePARAMETERS_searchOnly(rp, filename='PARAMETERS'):
    """
    Reads PARAMETERS file again, but ignores everything not concerning the 
    search. Updates the given Rparams object accordingly.
    
    Parameters
    ----------
    rp : Rparams
        Parameters for current run, as defined previously. Will be updated if 
        parameters have changed.
    filename : string, optional
        Read from which file. The default is 'PARAMETERS'.

    Returns
    -------
    None.

    """
    try:
        with open(filename, 'r') as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.error("PARAMETERS not found.")
        raise
    for line in lines:
        if "!" in line:
            line = line.split("!")[0].rstrip()
        if line.lstrip().startswith("SEARCH_KILL"):
            if not re.match(r"\s*SEARCH_KILL\s*=\s*[Ff](alse)?", line):
                rp.SEARCH_KILL = True
        if not "=" in line:
            continue #ignore all lines that don't have an "=" sign at all
        param = line.split('=')[0]        #parameter is defined left of "="
        if param:
            #get rid of spaces and check the leftmost entry.
            plist = param.split()
            if plist:
                param = plist[0]
        try:
            value = line.split('=', maxsplit=1)[1].rstrip()
            llist = tl.linelist(value)  #read the stuff to the right of "="
        except IndexError:
            llist = []
        if not llist:
            continue   
        if param == 'SEARCH_CONVERGENCE':
            flags = plist[1:]
            if flags[0].lower() not in ['dgen', 'gaussian']:
                continue
            fl = [None,None]
            for (i, s) in enumerate(llist[:2]):
                try:
                    fl[i] = float(s)
                except:
                    pass
            if flags[0].lower() == 'gaussian':
                if (fl[0] is not None and fl[0] > 0 
                    and fl[0] != rp.searchConvInit["gaussian"]):
                    rp.GAUSSIAN_WIDTH = fl[0]
                    rp.searchConvInit["gaussian"] = fl[0]
                if fl[1] is not None:
                    if 0 < fl[1] <= 1:
                        rp.GAUSSIAN_WIDTH_SCALING = fl[1]
            elif flags[0].lower() == 'dgen':
                if len(flags) == 1:
                    target = 'dec'
                elif flags[1].lower() in ['dec', 'best', 'all']:
                    target = flags[1].lower()
                else:
                    continue
                if fl[0] is not None and fl[0] > 0:
                    if fl[0] != rp.searchConvInit["dgen"][target]:
                        rp.SEARCH_MAX_DGEN[target] = fl[0]
                        rp.searchConvInit["dgen"][target] = fl[0]
                if fl[1] is not None:
                    if 1 <= fl[1]:
                        rp.SEARCH_MAX_DGEN_SCALING[target] = fl[1]


def readPARAMETERS(filename='PARAMETERS', slab=DEFAULT, silent=False):
    """Reads a PARAMETERS file and returns an Rparams object with the
    information"""
    # open input file
    loglevel = logger.level
    if silent:
        logger.setLevel(logging.ERROR)
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logger.error("PARAMETERS not found.")
        raise
    #read PARAMETERS:
    rpars = tl.Rparams()
    for line in rf:
        if "!" in line:
            line = line.split("!")[0].rstrip()
        if line.lstrip().startswith("SEARCH_KILL"):
            if not re.match(r"\s*SEARCH_KILL\s*=\s*[Ff](alse)?", line):
                logger.warning('PARAMETERS file: SEARCH_KILL is set at start '
                        'of program. This means the search will be stopped as '
                        'soon as it starts. Delete SEARCH_KILL from '
                        'PARAMETERS to avoid this.')
        if not "=" in line:
            continue     #ignore all lines that don't have an "=" sign at all
        param = line.split('=')[0]        #parameter is defined left of "="
        if not param:
            continue
        #get rid of spaces and check the leftmost entry.
        plist = param.split()
        if plist:
            param = plist[0]
        try:
            value = line.split('=', maxsplit=1)[1].rstrip()
            llist = tl.linelist(value)  #read the stuff to the right of "="
        except IndexError:
            llist = []
        if not llist:
            logger.warning('PARAMETERS file: ' + param + ' appears to '
                            'have no value')
            rpars.setHaltingLevel(1)
            continue
        if param == 'ATTENUATION_EPS':
            try:
                f = float(llist[0])
                if f >= 0.0001 and f < 1:
                    rpars.ATTENUATION_EPS = f
                else:
                    logger.warning('PARAMETERS file: Unexpected input for '
                                    'ATTENUATION_EPS. Input will be ignored.')
                    rpars.setHaltingLevel(1)
            except ValueError:
                logger.warning('PARAMETERS file: ATTENUATION_EPS: Could not '
                            'convert value to float. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'BEAM_INCIDENCE':
            if ',' in value:
                try:
                    sublists = tl.splitSublists(llist, ',')
                    for sl in sublists:
                        if sl[0].lower() == 'theta':
                            theta = float(sl[1])
                            if 0 <= theta and theta <= 90:
                                rpars.THETA = float(theta)
                            else:
                                logger.warning('PARAMETERS file: '
                                    'BEAM_INCIDENCE: Unexpected value for '
                                    'theta (should be between 0 and 90). '
                                    'Input will be ignored.')
                                rpars.setHaltingLevel(1)
                        elif sl[0].lower() == 'phi':
                            rpars.PHI = float(sl[1])%360
                        else:
                            logger.warning('PARAMETERS file: '
                                    'BEAM_INCIDENCE: Unknown flag found. '
                                    'Input will be ignored.')
                            rpars.setHaltingLevel(1)
                except ValueError:
                    logger.warning('PARAMETERS file: BEAM_INCIDENCE: '
                                    'Could not convert value to float. '
                                    'Input will be ignored.')
                    rpars.setHaltingLevel(1)
            else:
                try:
                    theta = float(llist[0])
                    if 0 <= theta and theta <= 90:
                        rpars.THETA = float(theta)
                    else:
                        logger.warning('PARAMETERS file: BEAM_INCIDENCE: '
                                'Unexpected value for theta (should be '
                                'between 0 and 90). Input will be ignored.')
                        rpars.setHaltingLevel(1)
                    rpars.PHI = float(llist[1])
                except ValueError:
                    logger.warning('PARAMETERS file: BEAM_INCIDENCE: '
                                    'Could not convert value to float. '
                                    'Input will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'BULKDOUBLING_EPS':
            try:
                rpars.BULKDOUBLING_EPS = float(llist[0])
            except ValueError:
                logger.warning('PARAMETERS file: BULKDOUBLING_EPS: '
                                'Could not convert value to float. Input '
                                'will be ignored.')
                rpars.setHaltingLevel(1)
            if rpars.BULKDOUBLING_EPS < 0.0001:
                rpars.BULKDOUBLING_EPS = 0.0001
                logger.warning('PARAMETERS file: BULKDOUBLING_EPS cannot be '
                    'smaller than 0.001 due to fortran reading it as an F7.4; '
                    'value was changed to 0.0001')
                rpars.setHaltingLevel(1)
        elif param == 'BULKDOUBLING_MAX':
            try:
                i = int(llist[0])
            except ValueError:
                logger.warning('PARAMETERS file: BULKDOUBLING_MAX: Could not '
                        'convert value to integer. Input will be ignored.')
                rpars.setHaltingLevel(1)
            if i > 0:
                rpars.BULKDOUBLING_MAX = i
            else:
                logger.warning('PARAMETERS file: BULKDOUBLING_MAX: '
                                'Unexpected input (0 or negative). Input '
                                'will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'BULK_REPEAT':
            s = value.lower()
            if not "[" in s:
                if not "(" in s:
                    try:
                        rpars.BULK_REPEAT = abs(float(llist[0]))
                    except ValueError:
                        logger.warning('PARAMETERS file: BULK_REPEAT: '
                                'Could not convert value to float. Input '
                                'will be ignored.')
                        rpars.setHaltingLevel(1)
                else:
                    m = re.match(r'\s*(c|z)\(\s*(?P<val>[0-9.]+)\s*\)', s)
                    if not m:
                        logger.warning('PARAMETERS file: BULK_REPEAT: '
                            'Could not parse input expression. Input will '
                            'be ignored.')
                    else:
                        try:
                            v = abs(float(m.group("val")))
                        except:
                            logger.warning('PARAMETERS file: '
                                'BULK_REPEAT: Could not conver value to '
                                'float. Input will be ignored.')
                            rpars.setHaltingLevel(1)
                        else:
                            if "z" in s:
                                rpars.BULK_REPEAT = v
                            else:  #c
                                rpars.BULK_REPEAT = slab.ucell[2,2] * v
            else:  # vector
                vec = tl.readVector(s, slab.ucell)
                if vec is None:
                    logger.warning('PARAMETERS file: BULK_REPEAT: '
                        'Could not parse input expression. Input will '
                        'be ignored.')
                else:
                    rpars.BULK_REPEAT = vec
        elif param == 'ELEMENT_MIX':
            ptl = [el.lower() for el in tl.periodic_table]
            found = False
            for el in llist:
                if el.lower() not in ptl:
                    logger.warning('PARAMETERS file: ELEMENT_MIX for '
                            +plist[1]+': '+el+' not found in periodic table. '
                            'ELEMENT_MIX for '+plist[1]+' will be ignored.')
                    rpars.setHaltingLevel(1)
                    found = True
            if not found:
                rpars.ELEMENT_MIX[plist[1]] = [el.capitalize()
                                            for el in llist]
        elif param == 'ELEMENT_RENAME':
            ptl = [el.lower() for el in tl.periodic_table]
            if llist[0].lower() not in ptl:
                logger.warning('PARAMETERS file: ELEMENT_RENAME for '
                        +plist[1]+': '+llist[0]+' not found in periodic table.'
                        ' ELEMENT_RENAME for '+plist[1]+' will be ignored.')
                rpars.setHaltingLevel(1)
            else:
                rpars.ELEMENT_RENAME[plist[1]] = llist[0].capitalize()
        elif param == 'FILAMENT_WF':
            if llist[0].lower() == 'w':
                rpars.FILAMENT_WF = 4.5
            elif llist[0].lower() == 'lab6':
                rpars.FILAMENT_WF = 2.65
            else:
                try:
                    rpars.FILAMENT_WF = float(llist[0])
                except:
                    logger.warning('PARAMETERS file: FILAMENT_WF parameter: '
                            'Error parsing values. Input will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'FORTRAN_COMP':
            if len(plist) <= 1 and llist[0].lower() in ["ifort","gfortran"]:
                rpars.getFortranComp(comp=llist[0].lower())
            elif (len(plist) > 1 and plist[1].lower() == "mpi" 
                      and llist[0].lower() in ["mpifort", "mpiifort"]):
                rpars.getFortranMpiComp(comp=llist[0].lower())
            else:
                delim = llist[0][0]     # should be quotation marks
                if delim not in ["'", '"']:
                    logger.warning('PARAMETERS file: FORTRAN_COMP '
                        'parameter: No valid shorthand and not delimited '
                        'by quotation marks. Value will be ignored.')
                    rpars.setHaltingLevel(1)
                else:
                    setTo = value.split(delim)[1]
                if len(plist) <= 1:
                    rpars.FORTRAN_COMP[0] = setTo
                elif plist[1].lower() == "post":
                    rpars.FORTRAN_COMP[1] = setTo
                elif plist[1].lower() == "mpi":
                    rpars.FORTRAN_COMP_MPI[0] = setTo
                elif plist[1].lower() == "mpipost":
                    rpars.FORTRAN_COMP_MPI[1] = setTo
                else:
                    logger.warning('PARAMETERS file: FORTRAN_COMP parameter: '
                        'Could not interpret flags. Value will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'HALTING':
            try:
                i = int(llist[0])
                if 1 <= i <= 3:
                    rpars.HALTING = i
                else:
                    logger.warning('PARAMETERS file: HALTING: Invalid '
                                'value given. Input will be ignored.')
                    rpars.setHaltingLevel(1)
            except ValueError:
                logger.warning('PARAMETERS file: HALTING: Could not convert '
                                'value to integer. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'IV_SHIFT_RANGE':
            if len(llist) in [2, 3]:
                fl = []
                try:
                    f = [float(s) for s in llist]
                except ValueError:
                    logger.warning('PARAMETERS file: Failed to convert '
                       'IV_SHIFT_RANGE input to floats Input will be '
                       'ignored')
                    rpars.setHaltingLevel(1)
                else:
                    if fl[1] >= fl[0]:
                        for i in range(0,2):
                            rpars.IV_SHIFT_RANGE[i] = fl[i]
                    else:
                        logger.warning('PARAMETERS file: IV_SHIFT_RANGE '
                            'end energy has to be greater than or equal '
                            'to start energy. Input will be ignored.')
                        rpars.setHaltingLevel(1)
                    if len(fl) == 3:
                        if fl[2] > 0:
                            rpars.IV_SHIFT_RANGE[2] = fl[2]
                        else:
                            logger.warning('PARAMETERS file: '
                                'IV_SHIFT_RANGE step has to be positive. '
                                'Input will be ignored.')
                            rpars.setHaltingLevel(1)
            else:
                logger.warning('PARAMETERS file: Unexpected number of '
                    'values for IV_SHIFT_RANGE. Input will be ignored.')
        elif param == 'LAYER_CUTS':
            # some simple filtering here, but leave as list of strings
            if "<" in value and ">" in value:
                logger.warning('PARAMETERS file: LAYER_CUTS parameter: '
                        'Cannot parse list containing both "<" and ">" '
                        'signs. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            elif "<" in value or ">" in value:
                newlist = []
                for s in llist:
                    s = s.replace("<", " < ")
                    s = s.replace(">", " > ")
                    newlist.extend(s.split())
                llist = newlist
            rgx = re.compile(r'\s*(dz|dc)\s*\(\s*(?P<cutoff>[0-9.]+)\s*\)')
            for (i,s) in enumerate(llist):
                if "dz" in s.lower() or "dc" in s.lower():
                    m = rgx.match(value.lower())
                    if m:
                        try:
                            float(m.group('cutoff'))
                            llist[i] = m.group(0)
                        except:
                            logger.warning('PARAMETERS file: LAYER_CUTS '
                                    'parameter: Could not parse function '
                                    + s + '. Input will be ignored.')
                            rpars.setHaltingLevel(1)
                            continue
                elif not (s == "<" or s == ">"):
                    try:
                        float(s)
                    except:
                        logger.warning('PARAMETERS file: LAYER_CUTS '
                                'parameter: Error parsing values. Input will '
                                'be ignored.')
                    rpars.setHaltingLevel(1)
                    continue
            rpars.LAYER_CUTS = llist
        elif param == 'LAYER_STACK_VERTICAL':
            s = llist[0].lower()
            if s in ['false','f','c']:
                rpars.LAYER_STACK_VERTICAL = False
            elif s in ['true','t','z']:
                rpars.LAYER_STACK_VERTICAL = True
            else:
                logger.warning('PARAMETERS file: LAYER_STACK_VERTICAL: Could '
                        'not parse given value. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'LMAX':
            try:
                i = int(llist[0])
            except ValueError:
                logger.warning('PARAMETERS file: LMAX: Could not convert '
                                'value to integer. Input will be ignored.')
                rpars.setHaltingLevel(1)
            if 1 <= i <= 15:
                if rpars.PHASESHIFT_EPS != 0:
                    logger.warning('PARAMETERS file: Both LMAX and '
                        'PHASESHIFT_EPS are being defined. PHASESHIFT_EPS '
                        'will be ignored.')
                rpars.LMAX = i
            else:
                logger.warning('PARAMETERS file: LMAX: Value has to be '
                                'between 1 and 15. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'LOG_DEBUG':
            s = llist[0].lower()
            if s in ['false','f']:
                rpars.LOG_DEBUG = False
            elif s in ['true','t']:
                rpars.LOG_DEBUG = True
            else:
                logger.warning('PARAMETERS file: LOG_DEBUG: Could not parse '
                                'given value. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'LOG_SEARCH':
            s = llist[0].lower()
            if s in ['false','f']:
                rpars.LOG_SEARCH = False
            elif s in ['true','t']:
                rpars.LOG_SEARCH = True
            else:
                logger.warning('PARAMETERS file: LOG_SEARCH: Could not parse '
                                'given value. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'N_BULK_LAYERS':
            try:
                rpars.N_BULK_LAYERS = int(llist[0])
            except:
                logger.warning('PARAMETERS file: N_BULK_LAYERS: could '
                                'not convert value to integer.')
            if rpars.N_BULK_LAYERS != 1 and rpars.N_BULK_LAYERS != 2:
                logger.warning('PARAMETERS file: N_BULK_LAYERS was set to '
                                +rpars.N_BULK_LAYERS+', 1 or 2 expected. '
                                'Value was set to 1 (default).')
                rpars.setHaltingLevel(2)
                rpars.N_BULK_LAYERS = 1
        elif param == 'N_CORES':
            i = 0
            try:
                i = int(llist[0])
            except:
                logger.warning('PARAMETERS file: N_CORES: could not '
                                'convert value to integer.')
            if i > 0:
                rpars.N_CORES = i
            else:
                logger.warning('PARAMETERS file: N_CORES should be a '
                                'positive integer.')
                rpars.setHaltingLevel(1)
        elif param == 'PHASESHIFT_EPS':
            if rpars.LMAX != 0:
                logger.warning('PARAMETERS file: Both LMAX and '
                    'PHASESHIFT_EPS are being defined. PHASESHIFT_EPS '
                    'will be ignored.')
            else:
                try:
                    f = float(llist[0])
                except ValueError:
                    s = llist[0].lower()[0]
                    if s == 'r':    #rough
                        f = 0.1
                    elif s == 'n':  #normal
                        f = 0.05
                    elif s == 'f':  #fine
                        f = 0.01
                    elif s == 'e':  #extrafine
                        f = 0.001
                    else:
                        logger.warning('PARAMETERS file: PHASESHIFT_EPS: '
                                'Could not convert value to float. '
                                'Input will be ignored.')
                        rpars.setHaltingLevel(1)
                if f > 0 and f < 1:
                    rpars.PHASESHIFT_EPS = f
                else:
                    logger.warning('PARAMETERS file: PHASESHIFT_EPS: '
                                'Unexpected value (should be between 0 '
                                'and 1). Input will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'PLOT_COLORS_RFACTOR':
            if len(llist) >= 2:
                if len(llist) > 2:
                    logger.warning('PARAMETERS file: PLOT_COLORS_RFACTOR '
                        'parameter: Expected two values, found {}. First two '
                        'values will be used.'.format(len(llist)))
                    rpars.setHaltingLevel(1)
                rpars.PLOT_COLORS_RFACTOR = (llist[0], llist[1])
            else:
                logger.warning('PARAMETERS file: PLOT_COLORS_RFACTOR '
                    'parameter: Expected two values, found {}. Input will be '
                    'ignored.'.format(len(llist)))
                rpars.setHaltingLevel(1)
        elif param == 'RUN':
            rl = []
            for s in llist:
                l = tl.readIntRange(s)
                if len(l) > 0:
                    rl.extend(l)
                else:
                    logger.warning('PARAMETERS file: RUN: Could not '
                            'interpret value '+s+', skipping value...')
                    rpars.setHaltingLevel(2)
            if len(rl) > 0:
                i = 0
                while i < len(rl):
                    if rl[i] not in [0,1,2,3,11,31]:
                        logger.warning('PARAMETERS file: RUN: Value '
                            +str(rl[i])+' does not correspond to a segment '
                            'and will be skipped.')
                        rl.pop(i)
                        rpars.setHaltingLevel(1)
                    else:
                        i += 1
            if len(rl) > 0:
                if rl[0] != 0:
                    rl.insert(0, 0)
                rpars.RUN = rl
            else:
                # could in principle try carrying on with some default, but
                #   this is probably a case in which a hard stop is
                #   appropriate
                logger.warning('PARAMETERS file: RUN was defined, but '
                        'no values were read. Execution will stop.')
                rpars.setHaltingLevel(3)
        elif param == 'R_FACTOR_SMOOTH':
            try:
                i = int(llist[0])
            except:
                logger.warning('PARAMETERS file: R_FACTOR_SMOOTH: '
                                'could not convert value to integer.')
                rpars.setHaltingLevel(1)
            else:
                if 1000 > i >= 0:
                    rpars.R_FACTOR_SMOOTH = i
                else:
                    if i < 0:
                        logger.warning('PARAMETERS file: R_FACTOR_SMOOTH '
                                        'should be >= 0')
                        rpars.setHaltingLevel(1)
                    else:
                        logger.warning('PARAMETERS file: R_FACTOR_SMOOTH '
                                        'should be < 1000')
                        rpars.setHaltingLevel(1)
                    rpars.setHaltingLevel(1)
        elif param == 'R_FACTOR_TYPE':
            try:
                rpars.R_FACTOR_TYPE = int(llist[0])
            except:
                logger.warning('PARAMETERS file: R_FACTOR_TYPE: could '
                                'not convert value to integer.')
            if rpars.R_FACTOR_TYPE != 1 and rpars.R_FACTOR_TYPE != 2:
                logger.warning('PARAMETERS file: R_FACTOR_TYPE was set to '
                                +rpars.R_FACTOR_TYPE+', 1 or 2 expected. '
                                'Value was set to 1 (default).')
                rpars.setHaltingLevel(1)
                rpars.R_FACTOR_TYPE = 1
        elif param == 'SEARCH_BEAMS':
            if llist[0][0].lower() in ["0","a"]:
                rpars.SEARCH_BEAMS = 0
            elif llist[0][0].lower() in ["1","i"]:
                rpars.SEARCH_BEAMS = 1
            elif llist[0][0].lower() in ["2","f"]:
                rpars.SEARCH_BEAMS = 2
            else:
                logger.warning('PARAMETERS file: SEARCH_BEAMS: value not '
                                'recognized. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'SEARCH_CONVERGENCE':
            if len(plist) == 1:
                if value.lower().strip() == 'off':
                    rpars.GAUSSIAN_WIDTH_SCALING = 1.
                else:
                    logger.warning('PARAMETERS file: SEARCH_CONVERGENCE: '
                        'no flag given, value '+value+' not recognized.')
                    rpars.setHaltingLevel(1)
                    continue
            flags = plist[1:]
            if flags[0].lower() not in ['dgen', 'gaussian']:
                logger.warning('PARAMETERS file: SEARCH_CONVERGENCE: flag "'
                    + flags[0] + '" not recognized.')
                rpars.setHaltingLevel(1)
                continue
            fl = [None,None]
            for (i, s) in enumerate(llist[:2]):
                try:
                    fl[i] = float(s)
                except:
                    logger.warning('PARAMETERS file: SEARCH_CONVERGENCE '
                            'gaussian: could not convert value to float.')
                    rpars.setHaltingLevel(1)
            if flags[0].lower() == 'gaussian':
                if fl[0] is not None and fl[0] > 0:
                    rpars.GAUSSIAN_WIDTH = fl[0]
                else:
                    logger.warning('PARAMETERS file: SEARCH_CONVERGENCE '
                        'gaussian should be a positive number.')
                    rpars.setHaltingLevel(1)
                if fl[1] is not None:
                    if 0 < fl[1] <= 1:
                        rpars.GAUSSIAN_WIDTH_SCALING = fl[1]
                    else:
                        logger.warning('PARAMETERS file: SEARCH_CONVERGENCE '
                            'gaussian: scaling value should be in range ]0, '
                            '1[')
                        rpars.setHaltingLevel(1)
            elif flags[0].lower() == 'dgen':
                if len(flags) == 1:
                    target = 'dec'
                elif flags[1].lower() in ['dec', 'best', 'all']:
                    target = flags[1].lower()
                else:
                    logger.warning('PARAMETERS file: SEARCH CONVERGENCE '
                        'dgen: flag "'+flags[1]+'" not recognized.')
                    rpars.setHaltingLevel(1)
                    continue
                if fl[0] is not None and fl[0] > 0:
                    rpars.SEARCH_MAX_DGEN[target] = fl[0]
                else:
                    logger.warning('PARAMETERS file: SEARCH_CONVERGENCE '
                        'dgen should be a positive number.')
                    rpars.setHaltingLevel(1)
                if fl[1] is not None:
                    if 1 <= fl[1]:
                        rpars.SEARCH_MAX_DGEN_SCALING[target] = fl[1]
                    else:
                        logger.warning('PARAMETERS file: SEARCH_CONVERGENCE '
                            'dgen '+target+': scaling value cannot be smaller '
                            'than 1.')
                        rpars.setHaltingLevel(1)
        elif param == 'SEARCH_CULL':
            try:
                f = float(llist[0])
            except:
                logger.warning('PARAMETERS file: SEARCH_CULL: could not '
                                'interpret value.')
                rpars.setHaltingLevel(1)
                continue
            if f >= 1:
                if f - int(f) < 1e-6:
                    rpars.SEARCH_CULL = int(f)
                else:
                    logger.warning('PARAMETERS file: SEARCH_CULL: Values '
                        'greater than 1 have to be integers.')
                    rpars.setHaltingLevel(1)
                    continue
            elif f >= 0:
                rpars.SEARCH_CULL = f
            else:
                logger.warning('PARAMETERS file: SEARCH_CULL cannot be '
                                'negative.')
                rpars.setHaltingLevel(1)
                continue
            if len(llist) > 1:
                if llist[1].lower() in ["clone", "genetic", "random"]:
                    rpars.SEARCH_CULL_TYPE = llist[1].lower()
                else:
                    logger.warning('PARAMETERS file: SEARCH_CULL type '
                                    'not recognized.')
                    rpars.setHaltingLevel(1)
        elif param == 'SEARCH_MAX_GEN':
            i = 0
            try:
                i = int(llist[0])
            except:
                logger.warning('PARAMETERS file: SEARCH_MAX_GEN: '
                                'could not convert value to integer.')
                rpars.setHaltingLevel(1)
            if i > 0:
                rpars.SEARCH_MAX_GEN = i
            else:
                logger.warning('PARAMETERS file: SEARCH_MAX_GEN should '
                                'be a positive integer.')
                rpars.setHaltingLevel(1)
        elif param == 'SEARCH_POPULATION':
            i = 0
            try:
                i = int(llist[0])
            except:
                logger.warning('PARAMETERS file: SEARCH_POPULATION: '
                                'could not convert value to integer.')
            if i > 0:
                rpars.SEARCH_POPULATION = i
            else:
                logger.warning('PARAMETERS file: SEARCH_POPULATION should '
                                'be a positive integer.')
                rpars.setHaltingLevel(1)
            if i < 16:
                logger.warning('SEARCH_POPULATION is very small. A '
                                'minimum value of 16 is recommended.')
        elif param == 'SEARCH_START':
            s = llist[0].lower()
            if s in ["random", "rand", "centered", "center", "control"]:
                if s in ["random, rand"]:
                    rpars.SEARCH_START = "random"
                elif s in ["centered", "center"]:
                    rpars.SEARCH_START = "centered"
                elif s == "control":
                    rpars.SEARCH_START = "control"
                elif s in ["cr", "crand", "crandom"]:
                    rpars.SEARCH_START = "crandom"
            else:
                logger.warning('PARAMETERS file: SEARCH_START: flag '
                                'not recognized.')
                rpars.setHaltingLevel(1)
        elif param == 'SITE_DEF':
            newdict = {}
            sublists = tl.splitSublists(llist, ',')
            for sl in sublists:
                atnums = []
                for i in range(1, len(sl)):
                    l = tl.readIntRange(sl[i])
                    if len(l) > 0:
                        atnums.extend(l)
                    elif "top(" in sl[i]:
                        if slab == DEFAULT:
                            logger.warning('PARAMETERS file: SITE_DEF '
                                'parameter contains a top() function, '
                                'but no slab was passed. The atoms '
                                'will be assigned the default site '
                                'type instead.')
                            rpars.setHaltingLevel(1)
                        else:
                            n = int(sl[i].split('(')[1].split(')')[0])
                            csatlist = slab.atlist[:]
                            csatlist.sort(key=lambda atom: atom.pos[2])
                            while n > 0:
                                at = csatlist.pop()
                                if at.el == plist[1]:
                                    atnums.append(at.oriN)
                                    n -= 1
                    else:
                        logger.error('PARAMETERS file: Problem with '
                                      'SITE_DEF input format')
                        raise
                newdict[sl[0]] = atnums
            rpars.SITE_DEF[plist[1]]=newdict
        elif param == 'SUPERLATTICE':
            if not 'M' in plist:
                if slab == DEFAULT:
                    logger.warning('PARAMETERS file: SUPERLATTICE '
                            'parameter appears to be in Wood notation, '
                            'but no slab was passed; cannot calculate '
                            'bulk unit cell!')
                    rpars.setHaltingLevel(2)
                else:
                    rpars.SUPERLATTICE = tl.readWoodsNotation(value,
                                                              slab.ucell)
                    rpars.superlattice_defined = True
            else:
                sublists = tl.splitSublists(llist, ',')
                if not len(sublists) == 2:
                    logger.warning('PARAMETERS file: error reading '
                                     'SUPERLATTICE matrix: number of '
                                     'lines is not equal 2.')
                    rpars.setHaltingLevel(2)
                else:
                    write=True
                    nl = []
                    for sl in sublists:
                        if len(sl) == 2:
                            try:
                                nl.append([float(s) for s in sl])
                            except ValueError:
                                logger.warning('PARAMETERS file: error '
                                    'reading SUPERLATTICE matrix: could '
                                    'not convert '+str(sl)+' to floats.')
                                rpars.setHaltingLevel(2)
                                write = False
                        else:
                            logger.warning('PARAMETERS file: error '
                                    'reading SUPERLATTICE matrix: number '
                                    'of columns is not equal 2.')
                            rpars.setHaltingLevel(2)
                            write = False
                    if write:
                        rpars.SUPERLATTICE = np.array(nl,dtype=float)
                        rpars.superlattice_defined = True
        elif param == 'SUPPRESS_EXECUTION':
            s = llist[0].lower()
            if s in ['false','f']:
                rpars.SUPPRESS_EXECUTION = False
            elif s in ['true','t']:
                rpars.SUPPRESS_EXECUTION = True
            else:
                logger.warning('PARAMETERS file: SUPPRESS_EXECUTION: '
                                'Could not parse given value. Input will '
                                'be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'SYMMETRIZE_INPUT':
            s = llist[0].lower()
            if s in ['false','f']:
                rpars.SYMMETRIZE_INPUT = False
            elif s in ['true','t']:
                rpars.SYMMETRIZE_INPUT = True
            else:
                logger.warning('PARAMETERS file: SYMMETRIZE_INPUT: Could '
                           'not parse given value. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'SYMMETRY_EPS':
            try:
                f = float(llist[0])
                if f > 0:
                    rpars.SYMMETRY_EPS = f
                    if f > 1.0:
                        logger.warning('PARAMETERS file: SYMMETRY_EPS: Given '
                            'value is greater than one Ångström. This is a '
                            'very loose constraint and might lead to '
                            'incorrect symmetry detection. Be sure to check '
                            'the output!')
                        rpars.setHaltingLevel(1)
                else:
                    logger.warning('PARAMETERS file: SYMMETRY_EPS: Input '
                                'is not positive. Input will be ignored.')
                    rpars.setHaltingLevel(1)
            except ValueError:
                logger.warning('PARAMETERS file: SYMMETRY_EPS: Could not '
                        'convert value to float. Input will be ignored.')
                rpars.setHaltingLevel(1)
            if len(llist) > 1:
                try:
                    f = float(llist[1])
                    if f > 0:
                        rpars.SYMMETRY_EPS_Z = f
                        if f > 1.0:
                            logger.warning('PARAMETERS file: SYMMETRY_EPS: '
                             'Given value is greater than one Ångström. This '
                             'is a very loose constraint and might lead to '
                             'incorrect symmetry detection. Be sure to check '
                             'the output!')
                            rpars.setHaltingLevel(1)
                    else:
                        logger.warning('PARAMETERS file: SYMMETRY_EPS: '
                           'Input is not positive. Input will be ignored.')
                        rpars.setHaltingLevel(1)
                        rpars.SYMMETRY_EPS_Z = rpars.SYMMETRY_EPS
                except ValueError:
                    rpars.SYMMETRY_EPS_Z = rpars.SYMMETRY_EPS
            else:
                rpars.SYMMETRY_EPS_Z = rpars.SYMMETRY_EPS
        elif param == 'SYMMETRY_FIND_ORI':
            s = llist[0].lower()
            if s in ['false','f']:
                rpars.SYMMETRY_FIND_ORI = False
            elif s in ['true','t']:
                rpars.SYMMETRY_FIND_ORI = True
            else:
                logger.warning('PARAMETERS file: SYMMETRY_FIND_ORI: '
                    'Could not parse given value. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'SYMMETRY_FIX':
            s = llist[0].lower()
            grouplist = ["p1","p2","pm","pg","cm","rcm","pmm","pmg","pgg",
                         "cmm","rcmm","p4","p4m","p4g","p3","p3m1","p31m",
                         "p6","p6m"]
            if s == 'true':
                pass    #same as default, determine symmetry automatically
            elif s == 'false':
                rpars.SYMMETRY_FIX = 'p1'
            elif s in grouplist:
                if not s in ["cm","pmg"]:
                    rpars.SYMMETRY_FIX = s
                else:
                    logger.warning('PARAMETERS file: SYMMETRY_FIX: For '
                                    'group '+s+', direction needs to be '
                                    'specified. Input will be ignored.')
                    rpars.setHaltingLevel(1)
            elif s[0:2] in ["pm","pg","cm"] or s[0:3] == "pmg":
                #regex to read
                rgx = re.compile(r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*'
                          +r'\[\s*(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
                m = rgx.match(value.lower())
                if not m:
                    logger.warning('PARAMETERS file: SYMMETRY_FIX: Could '
                           'not parse given value. Input will be ignored.')
                    rpars.setHaltingLevel(1)
                else:
                    i1 = i2 = -2
                    group = m.group('group')
                    try:
                        i1 = int(m.group('i1'))
                        i2 = int(m.group('i2'))
                    except ValueError:
                        logger.warning('PARAMETERS file: SYMMETRY_FIX: '
                                        'Could not parse given value. '
                                        'Input will be ignored.')
                        rpars.setHaltingLevel(1)
                    if (group in ["pm","pg","cm","rcm","pmg"]
                          and i1 in range(-1,3) and i2 in range(-1,3)):
                        rpars.SYMMETRY_FIX = m.group(0)
                    else:
                        logger.warning('PARAMETERS file: SYMMETRY_FIX: '
                                        'Could not parse given value. '
                                        'Input will be ignored.')
                        rpars.setHaltingLevel(1)
            else:
                logger.warning('PARAMETERS file: SYMMETRY_FIX: Could not '
                               'parse given value. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'TENSOR_INDEX':
            i = 0
            try:
                i = int(llist[0])
            except:
                logger.warning('PARAMETERS file: TENSOR_INDEX: could not '
                                'convert value to integer.')
            if i > 0:
                rpars.TENSOR_INDEX = i
            else:
                logger.warning('PARAMETERS file: TENSOR_INDEX should be a '
                                'positive integer.')
                rpars.setHaltingLevel(1)
        elif param == 'TENSOR_OUTPUT':
            nl = tl.recombineListElements(llist, '*')
            for s in nl:
                try:
                    v = int(s)
                    if v != 0 and v != 1:
                        logger.warning('PARAMETERS file: Problem with '
                            'TENSOR_OUTPUT input format: Found value '+str(v)
                            +', expected 0 or 1. Value will be ignored.')
                        rpars.setHaltingLevel(1)
                    else:
                        rpars.TENSOR_OUTPUT.append(v)
                except ValueError:
                    if '*' in s:
                        sl = s.split('*')
                        try:
                            r = int(sl[0])
                            v = int(sl[1])
                        except ValueError:
                            logger.warning('PARAMETERS file: Problem with '
                                'TENSOR_OUTPUT input format: could not read '
                                'value '+s+', value will be ignored.')
                            rpars.setHaltingLevel(1)
                        if v != 0 and v != 1:
                            logger.warning('PARAMETERS file: Problem with '
                                'TENSOR_OUTPUT input format: Found value '
                                +str(v)+', expected 0 or 1. Value will be '
                                'ignored.')
                            rpars.setHaltingLevel(1)
                        else:
                            for i in range(0,r):
                                rpars.TENSOR_OUTPUT.append(v)
                    else:
                        logger.warning('PARAMETERS file: Problem with '
                            'TENSOR_OUTPUT input format: could not read '
                            'value '+s+', value will be ignored.')
                        rpars.setHaltingLevel(1)
        elif param == 'THEO_ENERGIES':
            if len(llist) == 1:
                # single value input - only one energy requested
                try:
                    f = float(llist[0])
                except ValueError:
                    logger.warning('PARAMETERS file: Failed to convert '
                       'THEO_ENERGIES input to floats Input will be ignored')
                    rpars.setHaltingLevel(1)
                else:
                    if f > 0:
                        rpars.THEO_ENERGIES = [f, f, 1]
                    else:
                        logger.warning('PARAMETERS file: Unexpected '
                            'input for THEO_ENERGIES. Input will be ignored.')
                    rpars.setHaltingLevel(1)
            elif len(llist) == 3:
                fl = []
                defined = 0
                for s in llist:
                    if s == "_":
                        fl.append(-1)
                    else:
                        try:
                            f = float(s)
                            if f > 0:
                                fl.append(f)
                                defined += 1
                            else:
                                logger.warning('PARAMETERS file: '
                                    'THEO_ENERGIES values have to be '
                                    'positive. Input will be ignored.')
                                rpars.setHaltingLevel(1)
                        except ValueError:
                            logger.warning('PARAMETERS file: Failed to '
                                 'convert THEO_ENERGIES input to floats. '
                                 'Input will be ignored.')
                            rpars.setHaltingLevel(1)
                if len(fl) == 3:
                    if defined < 3:
                        rpars.THEO_ENERGIES = fl
                    elif (fl[0]>0 and fl[1]>fl[0] and fl[2]>0):
                        if (fl[1]-fl[0])%fl[2] != 0:
                            #if the max is not hit by the steps exactly,
                            # correct max up to make it so
                            fl[0] -= fl[2]-(fl[1]-fl[0])%fl[2]
                            if fl[0] <= 0:
                                fl[0] = fl[0] % fl[2]
                                if fl[0] == 0:
                                    fl[0] += fl[2]
                            logger.info('THEO_ENERGIES parameter: '
                                    '(Eto-Efrom)%Estep != 0, Efrom was '
                                    'corrected to '+str(fl[0]))
                        rpars.THEO_ENERGIES = fl
                    else:
                        logger.warning('PARAMETERS file: Unexpected '
                                        'input for THEO_ENERGIES. Input '
                                        'will be ignored.')
                        rpars.setHaltingLevel(1)
                else:
                    logger.warning('PARAMETERS file: THEO_ENERGIES '
                                    'appears to have no value.')
                    rpars.setHaltingLevel(1)
            else:
                logger.warning('PARAMETERS file: Unexpected input for '
                                'THEO_ENERGIES. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'T_DEBYE':
            try:
                rpars.T_DEBYE = float(llist[0])
            except ValueError:
                logger.warning('PARAMETERS file: T_DEBYE: '
                                'Could not convert value to float. Input '
                                'will be ignored.')
                rpars.setHaltingLevel(1)
            else:
                if rpars.T_DEBYE < 0:
                    rpars.T_DEBYE = None
                    logger.warning('PARAMETERS file: T_DEBYE must '
                                    'be positive. Input will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'T_EXPERIMENT':
            try:
                rpars.T_EXPERIMENT = float(llist[0])
            except ValueError:
                logger.warning('PARAMETERS file: T_EXPERIMENT: '
                                'Could not convert value to float. Input '
                                'will be ignored.')
                rpars.setHaltingLevel(1)
            else:
                if rpars.T_EXPERIMENT < 0:
                    rpars.T_EXPERIMENT = None
                    logger.warning('PARAMETERS file: T_EXPERIMENT must '
                                    'be positive. Input will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'V0_IMAG':
            try:
                f = float(llist[0])
            except:
                logger.warning('PARAMETERS file: Failed to convert '
                        'V0_IMAG input to float. Input will be ignored.')
            else:
                if f < 0:
                    logger.warning('PARAMETERS file: V0_IMAG should be '
                        'positive real value. Input '+llist[0]+' not '
                        'accepted. Input will be ignored.')
                else:
                    rpars.V0_IMAG = f
        elif param == 'V0_REAL':
            if llist[0].lower() == 'rundgren':
                try:
                    c = []
                    for i in range(0,4):
                        c.append(float(llist[i+1]))
                    setTo = ("workfn-max("+str(round(c[0],3))
                        +", (("+str(round(c[1],3))+")+("
                        +str(round(c[2],3))+")/sqrt(EEV+workfn+("
                        +str(round(c[3],3))+"))))")
                except:
                    logger.warning("PARAMETERS file: V0_REAL parameter: "
                        "could not parse constants for Rundgren-type "
                        "function. Input will be ignored.")
                    rpars.setHaltingLevel(1)
            else:
                setTo = re.sub("(?i)EE","EEV+workfn",value)
            setTo = setTo.rstrip()
            rpars.V0_REAL = setTo
        elif param == 'V0_Z_ONSET':
            try:
                rpars.V0_Z_ONSET = float(llist[0])
            except ValueError:
                logger.warning('PARAMETERS file: V0_Z_ONSET: Could not '
                            'convert value to float. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'VIBR_AMP_SCALE':
            rpars.VIBR_AMP_SCALE.extend(value.split(","))
                                # taken at face value, interpreted later
        else:
            logger.warning('PARAMETERS file: Parameter '+param+' not '
                            'recognized.')
            rpars.setHaltingLevel(1)

    rf.close()
    logger.setLevel(loglevel)
    return(rpars)

def readDISPLACEMENTS(rp, filename="DISPLACEMENTS"):
    """Reads the DISPLACEMENTS file to rp.disp_blocks without interpreting 
    it."""
    rp.disp_blocks = []
    copyblock = None
    loopStarts = []
    loopLine = False
    try:
        with open(filename, 'r') as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.error("DISPLACEMENTS not found.")
        raise
    except:
        logger.error("Error reading DISPLACEMENTS file.")
        raise
    for (i, line) in enumerate(lines):
        if "!" in line:
            line = line.split("!")[0]
        line = line.strip()
        if re.search(r'<\s*/?\s*loop\s*>', line.lower()):
            if "=" in line:
                logger.warning("DISPLACEMENTS file: One line cannot contain "
                                "both an '=' sign and a <loop> flag.")
                rp.setHaltingLevel(2)
                continue
            loopLine = True
            ltags = re.findall(r'<\s*/?\s*loop\s*>', line.lower())
            n_blocks_real = len(rp.disp_blocks)
            if copyblock:
                n_blocks_real += 1
            for t in ltags:
                if re.match(r'<\s*loop\s*>', t):
                    loopStarts.append(n_blocks_real)
                else:
                    if len(loopStarts) == 0:
                        logger.warning("DISPLACEMENTS file: Unmatched "
                                        "</loop> flag found, skipping.")
                        rp.setHaltingLevel(2)
                        continue
                    if loopStarts[-1] == n_blocks_real:
                        logger.warning("DISPLACEMENTS file: Found empty "
                                        "loop.")
                        loopStarts.pop()
                        rp.setHaltingLevel(1)
                        continue
                    rp.disp_loops.append((loopStarts.pop(), n_blocks_real - 1))
            continue
        if not '=' in line:
            continue
        if line.startswith('=='):
            if re.match(r'==\s*s', line.lower()):  # search block
                loopLine = False
                try:
                    name = (line[re.match(r'==\s*s(earch)?\s+', 
                                          line.lower()).span()[1]:].strip())
                except:
                    name = ""
                names = [n for (_, n) in rp.disp_blocks]
                if not name or name in names:  # get unique name
                    if not name:
                        i = 1
                    else:
                        i = 2
                    nn = (name + " {}".format(i)).strip()
                    while nn in names:
                        i += 1
                        nn = (name + " {}".format(i)).strip()
                    name = nn
                if not rp.disp_blocks or len(rp.disp_blocks[-1][0]) != 0:
                    if copyblock:
                        rp.disp_blocks.append(copyblock)
                        copyblock = None
                    rp.disp_blocks.append(([], name))
                else:
                    rp.disp_blocks[-1] = ([], name)
                continue
        elif loopLine and not len(rp.disp_blocks) == 0:
            logger.error("DISPLACEMENTS file line {}: Line with <loop> or "
                          "</loop> flags must be followed by a new search "
                          "block or end of file!".format(i+1))
            return("Syntax error")
        loopLine = False
        if len(rp.disp_blocks) == 0:
            rp.disp_blocks.append(([], ""))
        p = line.split("=")[0]
        if ("xy" in p or "ab" in p) and not "[" in p:  # shorthand for 2 blocks
            if not copyblock:
                if rp.disp_blocks[-1][1]:
                    names = (rp.disp_blocks[-1][1]+" [1 0]", 
                             rp.disp_blocks[-1][1]+" [0 1]")
                else:
                    names = ("", "")
                rp.disp_blocks[-1] = (rp.disp_blocks[-1][0], names[0])
                copyblock = (rp.disp_blocks[-1][0][:], names[1])
            if "xy" in p:
                s = "xy"
            else:
                s = "ab"
            copyblock[0].append(line.replace(s, s+"[0 1]", 1))
            line = line.replace(s, s+"[1 0]", 1)

        elif copyblock:
            copyblock[0].append(line)
        rp.disp_blocks[-1][0].append(line)
    if copyblock:
        rp.disp_blocks.append(copyblock)
    if len(rp.disp_blocks[-1][0]) == 0 and len(rp.disp_blocks) > 1:
        rp.disp_blocks = rp.disp_blocks[:-1]
    if len(loopStarts) != 0:
        logger.warning("DISPLACEMENTS file: Unmatched <loop> flags found, "
                        "loops are still open at end of file.")
        rp.setHaltingLevel(2)
    return 0

def readDISPLACEMENTS_block(rp, sl, dispblock):
    """Reads a block from the DISPLACEMENTS file and adds the information to 
    all atoms in the slab."""
    abst = np.transpose(sl.ucell[0:2,0:2])
    uCellState = sl.uCellMod # if the unit cell gets modified by SYM_DELTA,
                             #   restore it afterwards
    (lines, name) = dispblock
    #read DISPLACEMENTS:
    mode = 0  #0: not reading; 1: reading geo, 2: reading vib, 3: reading occ
    regex = False   #read regular expressions as-is or not
    constraints = []  # collect constraints here first, interpret after loop
    for line in lines:
        if "!" in line:
            line = line.split("!")[0]
        line = line.strip()
        if '=' in line:
            if line[0] == '=':
                llist = line[1:].split()
                if llist[0][0].lower() == 'g':
                    mode = 1  # geometry
                elif llist[0][0].lower() == 'v':
                    mode = 2  # vibration
                elif llist[0][0].lower() == 'o':
                    mode = 3  # occupation
                elif llist[0][0].lower() == 'c':
                    mode = 4  # constraint
                elif llist[0][0].lower() == 'r':
                    regex = True
                    if len(llist) >= 2:
                        if llist[1].lower() == 'off':
                            regex = False
                else:
                    logger.warning("DISPLACEMENTS: Found line starting with "
                                    "'=', but didn't recognize flag.")
                    rp.setHaltingLevel(1)
                continue
            elif mode == 0:
                continue
            else:
                pside = line.split('=')[0].strip()
                if pside:
                     #get rid of spaces and check the leftmost entry.
                    plist = tl.linelist(pside)
                    if plist:
                        param = plist[0]
                        if param[0] == '!':
                            continue
                else:
                    continue
        else:
            continue
        value = line.split('=')[1].strip()
        try:
            llist = tl.linelist(value)
        except IndexError:
            logger.warning("DISPLACEMENTS file: " + param + " appears to "
                            "have no value")
            rp.setHaltingLevel(1)
            continue
        if not llist:
            logger.warning("DISPLACEMENTS file: " + param + " appears to "
                            "have no value")
            rp.setHaltingLevel(1)
            continue
        if param.lower() == "sym_delta":
            s = llist[0].lower()
            grouplist = ["p1","p2","pm","pg","cm","rcm","pmm","pmg","pgg",
                         "cmm","rcmm","p4","p4m","p4g","p3","p3m1","p31m",
                         "p6","p6m"]
            targetsym = ""
            if s[0] == "t":
                # True - go to highest symmetry
                targetsym = "found"
            elif s[0] == "f":
                # False - set p1
                targetsym = "p1"
            elif s in grouplist:
                if not s in ["cm","pmg"]:
                    targetsym = s
                else:
                    logger.warning('DISPLACEMENTS file: SYM_DELTA: For '
                                    'group '+s+', direction needs to be '
                                    'specified. Input will be ignored.')
                    rp.setHaltingLevel(2)
            elif s[0:2] in ["pm","pg","cm"] or s[0:3] == "pmg":
                #regex to read
                rgx = re.compile(r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*'
                          +r'\[\s*(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
                m = rgx.match(value.lower())
                if not m:
                    logger.warning('DISPLACEMENTS file: SYM_DELTA: Could '
                           'not parse given value. Input will be ignored.')
                    rp.setHaltingLevel(2)
                else:
                    i1 = i2 = -2
                    group = m.group('group')
                    try:
                        i1 = int(m.group('i1'))
                        i2 = int(m.group('i2'))
                    except ValueError:
                        logger.warning('DISPLACEMENTS file: SYM_DELTA: '
                                        'Could not parse given value. '
                                        'Input will be ignored.')
                        rp.setHaltingLevel(1)
                    if (group in ["pm","pg","cm","rcm","pmg"]
                          and i1 in range(-1,3) and i2 in range(-1,3)):
                        targetsym = m.group(0)
                    else:
                        logger.warning('DISPLACEMENTS file: SYM_DELTA: '
                                        'Could not parse given value. '
                                        'Input will be ignored.')
                        rp.setHaltingLevel(2)
            else:
                logger.warning('DISPLACEMENTS file: SYM_DELTA: Could not '
                               'parse given value. Input will be ignored.')
                rp.setHaltingLevel(2)
            if targetsym != "":
                sl.revertUnitCell(uCellState)
                sl.setSymmetry(rp, targetsym)
                sl.enforceSymmetry(rp, movement=False, rotcell=False)
            continue
        # if we are still here, then it's a "normal" line
        ats = []
        if mode == 1:
            rgx = re.compile(r'\s*(?P<nums>(([-:0-9]|((l|L)\([-:0-9\s]+\))'
                    r')\s*)*)(?P<dir>(z|((ab|xy)\s*)?\[[-0-9]+\s+[-0-9]+\]|'
                    r'(azi|r)\(?((ab|xy)\s*)?\[[-0-9\.]+\s+[-0-9\.]+\]\)?)?'
                    r'(\s*offset)?)')
        else:
            rgx = re.compile(r'(?P<nums>(([-:0-9]|((l|L)\([-:0-9\s]+\)))\s*)*)'
                             r'(\s*offset)?')
        if mode == 4:
            if param.lower() == "geo":
                ctype = 1
            elif param.lower() == "vib":
                ctype = 2
            elif param.lower() == "occ":
                ctype = 3
            else:
                logger.warning("DISPLACEMENTS file: CONSTRAINT flag "+param+
                        " not recognized. Input will be ignored.")
                rp.setHaltingLevel(2)
                continue
            pspec = pside.split(maxsplit=1)[1]
        else:
            pspec = pside
        _break = True
        for s in pspec.split(","):
            s = s.strip()
            sname = s.split()[0]
            spec = s.split(sname,1)[1].strip().lower()
            m = rgx.match(spec)
            if not m or (mode == 1 and not m.group("dir")):
                logger.warning('DISPLACEMENTS file: could not parse as '
                    'numbers or direction, skipping line: '+pside)
                rp.setHaltingLevel(2)
                break
            ats.append((sname, m.group("nums")))
            if mode == 1:
                dr = m.group("dir")
                if not dr:
                    logger.warning('DISPLACEMENTS file: could '
                        'not parse direction, skipping line: '+pside)
                    rp.setHaltingLevel(2)
                    break
                # quick consistency check on direction:
                if "[" in dr:
                    if not "azi" in dr and not "r" in dr:
                        try:
                            dirints = (int(dr.split("[")[1].split(" ")[0]),
                                       int(dr.split(" ")[1].split("]")[0]))
                        except ValueError:
                            logger.warning('DISPLACEMENTS file: could '
                                'not parse direction "'+dr+'", skipping '
                                'line: '+pside)
                            rp.setHaltingLevel(2)
                            break
                    else:
                        try:
                            dirfloats = (float(dr.split("[")[1]
                                               .split(" ")[0]),
                                     float(dr.split(" ")[1].split("]")[0]))
                        except ValueError:
                            logger.warning('DISPLACEMENTS file: could '
                                'not parse direction "'+dr+'", skipping '
                                'line: '+pside)
                            rp.setHaltingLevel(2)
                            break
        else:
            _break = False
        if _break:
            rp.setHaltingLevel(2)
            continue
        # get list of atoms to manipulate
        targetAtEls = []
        _break = True
        for (sname, nums) in ats:
            numlist = tl.linelist(nums)
            if numlist:
                # get proper numlist as integers
                # first need to recombine expressions like "L(2 3)"
                nl = []
                partstring = ""
                bracketcount = 0
                for s in numlist:
                    partstring += s
                    bracketcount += s.count("(") - s.count(")")
                    if bracketcount == 0:
                        nl.append(partstring)
                        partstring = ""
                    else:
                        partstring += " "
                if bracketcount != 0:
                    logger.warning("DISPLACEMENTS file: Bracket mismatch in "
                                    "line " + pside)
                    if len(partstring) > 0:
                        nl.append(partstring)
                numlist = []
                for s in nl:
                    if not "l(" in s:
                        l = tl.readIntRange(s)
                        if len(l) > 0:
                            numlist.extend(l)
                        else:
                            logger.warning('DISPLACEMENTS file: could not '
                                'parse integer range, skipping line: '+pside)
                            rp.setHaltingLevel(2)
                            break
                    else:
                        m = re.match(r'l\((?P<laynum>[-:0-9\s]+)\)', s)
                        if not m:
                            logger.warning('DISPLACEMENTS file: could not '
                                    'parse layer expression, skipping line: '
                                    +pside)
                            rp.setHaltingLevel(2)
                            break
                        l = tl.readIntRange(m.group("laynum"))
                        if len(l) == 0:
                            logger.warning('DISPLACEMENTS file: could '
                                    'not parse layer expression, skipping '
                                    'line: '+pside)
                            rp.setHaltingLevel(2)
                            break
                        for ln in l:
                            if ln > len(sl.layers):
                                logger.warning('DISPLACEMENTS file: '
                                    'layer number out of bounds, skipping '
                                    'line: '+pside)
                                rp.setHaltingLevel(2)
                                break
                            numlist.extend([at.oriN for at in sl.atlist
                                        if at.layer == sl.layers[ln-1]])
                else:  # loop finished without beak
                    _break = False
                if _break:
                    break
                _break = True
            # interpret label:
            targetsites = []
            # first try with POSCAR element:
            for site in sl.sitelist:
                s = sname
                if not "_" in s:
                    if regex:
                        s += ".*"
                    else:
                        s += "*"
                if not regex:
                    s = re.escape(s)   #double-slash non-literal characters
                    #if regular expressions are not enabled, we want to
                    #  still interpret * as "any number of any characters"
                    s = s.replace('\\*','.*')
                m = re.match(s, site.label)
                if m:
                    #if the length of the matched text == the site label,
                    #  it's a real match
                    if m.end(0) == len(site.label):
                        targetsites.append(site)
            if len(targetsites) > 0:
                targetel = ""   # use ALL elements of site
            elif mode == 3:
                # if we're reading occupation, only POSCAR elements are
                #   allowed -> error
                logger.warning('DISPLACEMENTS file: could not parse '
                    +sname+' as POSCAR element or site label, skipping '
                    'line.')
                rp.setHaltingLevel(2)
                break
            else:
                # check ELEMENT_MIX elements:
                el = sname.split("_")[0]
                s = el
                if not regex:
                    s = re.escape(s)   #double-slash
                    s = s.replace('\\*','.*')
                cels = [cel for cel in sl.chemelem if re.match(s, cel)]
                if el not in sl.elements and len(cels) > 0:
                    pel = ""
                    for k in rp.ELEMENT_MIX:
                        if any(e in rp.ELEMENT_MIX[k] for e in cels):
                            pel = k
                            targetel = [e for e in cels 
                                        if e in rp.ELEMENT_MIX[k]][0]
                            break
                    s = sname.replace(sname.split("_")[0], pel, 1)
                    if s in sl.elements:  # is element, all sites
                        targetsites.extend([site for site in sl.sitelist if 
                                            site.el == s])
                    else:
                        for site in sl.sitelist:
                            if not regex:
                                s = re.escape(s)   #double-slash
                                s = s.replace('\\*','.*')
                            m = re.match(s, site.label)
                            if m:
                                if m.end(0) == len(site.label):
                                    targetsites.append(site)
                else:
                    logger.warning('DISPLACEMENTS file: could not parse '
                        +sname+' as element or site label, skipping '
                        'line.')
                    rp.setHaltingLevel(2)
                    break
                if len(targetsites) == 0:
                    logger.warning('DISPLACEMENTS file: could not parse '
                        +sname+' as element or site label, skipping '
                        'line.')
                    rp.setHaltingLevel(2)
                    break
            for at in sl.atlist:
                if ((at.oriN in numlist or len(numlist) == 0)
                      and at.site in targetsites):
                    targetAtEls.append((at, targetel))
        else:  # loop finished without break
            _break = False
        if _break:
            continue  # error message has already happened
        if len(targetAtEls) == 0:
            logger.warning('DISPLACEMENTS file: no atoms found for '
                'line, line will be skipped: '+pside)
            rp.setHaltingLevel(1)
            continue
        if mode == 2 or (mode == 1 and not "offset" in dr):
            # geometrical or vibrational displacement, get range:
            try:
                fl = [float(s) for s in llist]
            except:
                if ("offset" in pside and llist[0].lower() 
                                            in ["clear", "original"]):
                    fl = [0.]
                elif "offset" in pside and llist[0].lower().startswith("cont"):
                    continue    # instruction to use value from previous search
                                # -> default, ignore
                else:
                    fl = []
                    logger.warning('DISPLACEMENTS file: could not parse '
                        +value+' as list of floats, skipping line.')
                    rp.setHaltingLevel(1)
                    continue
            if len(fl) < 3:
                if mode == 2 and len(fl) == 1:  # interpret as static offset
                    fl = [fl[0], fl[0], 1]
                else:
                    logger.warning('DISPLACEMENTS file: too few values '
                        'found, skipping line: '+pside)
                    rp.setHaltingLevel(1)
                    continue
            steps = abs(int(round((fl[1]-fl[0]) / fl[2])))+1
            mid = (fl[1]+fl[0]) / 2
            if steps % 2 == 0:      # even number of steps, extend range
                steps += 1
            drange = np.arange(mid - ((steps-1)/2*fl[2]),
                               mid + ((steps-1)/2*fl[2])+1e-6, fl[2])
            if fl[1] < fl[0]:
                drange = drange[::-1]  #reverse
        if mode == 1 and "offset" in dr:
            # geo offset, get value
            try:
                offval = float(llist[0])
            except:
                if llist[0].lower() in ["clear", "original"]:
                    offval = 0.
                elif llist[0].lower().startswith("cont"):
                    continue
                else:
                    logger.warning('DISPLACEMENTS file: could not parse '
                                    +value+' as float, skipping line.')
                    rp.setHaltingLevel(1)
                    continue
            drange = [offval]
        if mode == 1:
            if dr.strip() == "offset":
                if offval != 0.:
                    logger.warning('DISPLACEMENTS file: cannot assign '
                            'geo offset : no direction given: '+pside)
                    rp.setHaltingLevel(1)
                    continue
                else:
                    for (at, targetel) in targetAtEls:
                        at.clearOffset(1, targetel)
                    continue
            # geometrical displacement, get direction vector:
            drvec = np.array([0.,0.,0.])
            if not "azi" in dr and not "r" in dr:
                if "z" in dr:
                    drvec[2] = -1. #invert direction: positive = away from bulk
                elif "xy" in dr:
                    v = np.array(dirints)
                    v = v / np.linalg.norm(v)
                    drvec[:2] = v
                else:
                    v = dirints[0]*abst[0] + dirints[1]*abst[1]
                    v = v / np.linalg.norm(v)
                    drvec[:2] = v
                disprange = [f*drvec for f in drange]
                for (at, targetel) in targetAtEls:
                    # check whether displacement is allowed
                    if "z" in dr:
                        allowed = True
                    elif type(at.freedir) == int:
                        if at.freedir == 0:
                            allowed = False
                        else:
                            allowed = True
                    else:
                        freev = at.freedir[0]*abst[0] + at.freedir[1]*abst[1]
                        freev = freev / np.linalg.norm(freev)
                        if abs(abs(np.dot(drvec[:2],freev))-1) < 1e-5:
                            allowed = True
                        else:
                            allowed = False
                    if allowed:
                        if not "offset" in dr:
                            at.assignDisp(mode, disprange, targetel)
                        else:
                            at.clearOffset(1, targetel)
                            at.assignDisp(4, disprange, targetel)
                    else:
                        logger.warning("In-plane displacement assignment for "
                            "atom {} is forbidden by symmetry and will be "
                            "skipped. See 'FreeDir' in POSCAR file. To "
                            "apply this displacement, use either the "
                            "SYMMETRY_FIX parameter to lower the symmetry, "
                            "or use SYM_DELTA in the DISPLACEMENTS file to "
                            "allow symmetry breaking for this atom."
                            .format(at.oriN))
            else:
                if "xy" in dr:
                    c = np.array(dirfloats)
                else:
                    c = dirfloats[0]*abst[0] + dirfloats[1]*abst[1]
                for (at, targetel) in targetAtEls:
                    if "r" in dr:
                        v = at.cartpos[:2] - c
                        v = v / np.linalg.norm(v)
                        allowed = True
                        if type(at.freedir) == int and at.freedir == 0:
                            allowed = False
                        elif type(at.freedir) != int:
                            freev = (at.freedir[0]*abst[0] + 
                                     at.freedir[1]*abst[1])
                            freev = freev / np.linalg.norm(freev)
                            if abs(abs(np.dot(v,freev))-1) > 1e-5:
                                allowed = False
                        if allowed:
                            drvec = np.array([0.,0.,0.])
                            drvec[:2] = v
                            disprange = [f*drvec for f in drange]
                            if not "offset" in dr:
                                at.assignDisp(mode, disprange, targetel)
                            else:
                                at.assignDisp(4, disprange, targetel)
                    else:  # azi
                        # allowed only for completely free atoms
                        if type(at.freedir) == int and at.freedir == 1:
                            r = at.cartpos[:2] - c
                            disprange = []
                            for d in drange:
                                # disprange based on rotating vector around c
                                ar = d / np.linalg.norm(r)  # rotation angle
                                tm = np.array([[np.cos(ar)-1, -np.sin(ar)],
                                               [np.sin(ar), np.cos(ar)-1]])
                                disprange.append(np.append(np.dot(tm, r), 0.))
                            if not "offset" in dr:
                                at.assignDisp(mode, disprange, targetel)
                            else:
                                at.clearOffset(1, targetel)
                                at.assignDisp(4, disprange, targetel)
                        else:
                            logger.warning("In-plane azimuthal displacement "
                                "assignment for atom {} is forbidden by "
                                "symmetry and will be skipped. See 'FreeDir' "
                                "in POSCAR file. To apply this displacement, "
                                "use either the SYMMETRY_FIX parameter to "
                                "lower the symmetry, or use SYM_DELTA in the "
                                "DISPLACEMENTS file to allow symmetry "
                                "breaking for this atom.".format(at.oriN))
        elif mode == 2:
            # vibrational displacement, apply:
            for (at, targetel) in targetAtEls:
                if "offset" in pside:
                    at.clearOffset(mode, targetel)
                if not (len(drange) == 1 and drange[0] == 0.):
                    at.assignDisp(mode, drange, targetel)
        elif mode == 3:
            # occupations, get ranges:
            sublists = tl.splitSublists(llist, ',')
            maxsteps = 0
            _break = True
            for subl in sublists:
                # get range
                try:
                    fl = [float(s) for s in subl[-3:]]
                except:
                    fl = []
                    if len(subl) == 2:
                        try:
                            fl = [float(subl[-1])]
                        except:
                            pass
                    if len(fl) == 1:
                        fl = [fl[0], fl[0], 1]
                    if len(fl) < 3:
                        logger.warning('DISPLACEMENTS file: could not parse '
                                   +value+' as list of floats, skipping line.')
                        rp.setHaltingLevel(1)
                        break
                if len(fl) < 3:
                    if len(fl) == 1:  # interpret as static offset
                        fl = [fl[0], fl[0], 1]
                    else:
                        logger.warning('DISPLACEMENTS file: too few values '
                                        'found, skipping line: '+pside)
                        rp.setHaltingLevel(1)
                        break
                steps = abs(int(round((fl[1]-fl[0]) / fl[2])))+1
                if steps > maxsteps:
                    maxsteps = steps
            else:
                _break = False
            if _break:
                continue
            for subl in sublists:
                try:
                    fl = [float(s) for s in subl[-3:]]
                except:
                    f = float(subl[-1])  # checked before, will work
                    fl = [f, f, 1]
                steps = abs(int(round((fl[1]-fl[0]) / fl[2])))+1
                if steps > maxsteps:
                    steps = maxsteps
                    logger.warning('DISPLACEMENTS file: inconsistent step '
                        'numbers for occupancies of '+pside+'. Decreasing '
                        'spacing for '+subl[0]+' to make number of steps '
                        'equal.')
                mid = (fl[1]+fl[0]) / 2
                drange = np.arange(mid-((steps-1)/2*fl[2]),
                               mid+((steps-1)/2*fl[2])+1e-5, fl[2])
                if fl[1] < fl[0]:
                    drange = drange[::-1]  #reverse
                # get element
                if not "+" in "".join(subl):
                    if subl[0] in sl.chemelem:
                        targetel = subl[0]
                    else:
                        logger.warning('DISPLACEMENTS file: '+subl[0]+'not '
                                    'found in element list. No assignment '
                                    'will be made.')
                        rp.setHaltingLevel(1)
                        continue
                    for (at, _) in targetAtEls:
                        if targetel in at.disp_occ:
                            at.assignDisp(mode, drange, targetel)
                        else:
                            logger.warning('DISPLACEMENTS file: '
                                'trying to address atom number '
                                +str(at.oriN)+' with wrong element. '
                                'Atom will be skipped.')
                            rp.setHaltingLevel(1)
                else:
                    # !!! currently unused, not in wiki - delete?
                    elparts = "".join(subl[:-3]).split("+")
                    elweights = []
                    elsum = 0.0
                    rgx = re.compile(r'\s*(?P<number>[0-9\.]+)\s*'
                                     +r'(?P<name>[a-zA-Z]+)')
                    for part in elparts:
                        m = rgx.match(part)
                        if not m:
                            logger.warning('DISPLACEMENTS file: could '
                                'not parse '+part+' as number and '
                                'element. No assignment will be made.')
                            rp.setHaltingLevel(1)
                            break
                        try:
                            f = float(m.group("number"))
                        except ValueError:
                            logger.warning('DISPLACEMENTS file: '
                             'could not parse '+m.group("number")
                             +' as float. No assignment will be made.')
                            rp.setHaltingLevel(1)
                            break
                        el = m.group("name")
                        if not el in sl.chemelem:
                            logger.warning('DISPLACEMENTS file: '+el
                                +' not found in element list. No '
                                'assignment will be made.')
                            rp.setHaltingLevel(1)
                            break
                        else:
                            elweights.append((el, f))
                            elsum += f
                    else:  # loop finished without break
                        for (el, f) in elweights:
                            ndr = [v*f/elsum for v in drange]
                            for (at, targetel) in targetAtEls:
                                if el in at.disp_occ:
                                    at.assignDisp(mode, ndr, el)
                                else:
                                    logger.warning('DISPLACEMENTS file: '
                                        'trying to address atom number '
                                        +str(at.oriN)+' with wrong element. '
                                        'Atom will be skipped.')
                                    rp.setHaltingLevel(1)
        elif mode == 4:
            constraints.append((targetAtEls, ctype, value))
    sl.revertUnitCell(uCellState)  # if modified by SYM_DELTA, go back
    # now read constraints
    for (targetAtEls, ctype, value) in constraints:
        if value.lower().startswith("link"):
            linkto = targetAtEls[0]
            for (at, el) in targetAtEls:  # [1:]:
                at.assignConstraint(ctype, targetel = el, linkAtEl = linkto)
        elif value.lower().startswith("ind("):
            try:
                ind = int(value.lower().split("ind(")[1].split(")")[0])
            except:
                logger.warning("DISPLACEMENTS file: Could not convert "+value
                                +" to index, skipping line.")
                continue
            for (at, el) in targetAtEls:
                at.assignConstraint(ctype, targetel = el, index = ind)
        else:
            if ctype == 3:
                occEl = value.split()[0]
                value = value.split(maxsplit=1)[-1]
            try:
                f = float(value.split()[0])
            except:
                logger.warning("DISPLACEMENTS file: Could not convert "+value
                                +" to float, skipping line.")
                continue
            for (at, el) in targetAtEls:
                if ctype == 3:
                    el = occEl
                at.assignConstraint(ctype, targetel = el, value = f)
    # Reading done, consistency check...
    o = ""
    if name:
        o = " "+name
    logger.debug("DISPLACEMENTS block{} was read successfully".format(o))
    for at in sl.atlist:
        if len(at.disp_occ) > 5:
            logger.error('DISPLACEMENTS file: Atom '+str(at.oriN)+' has '
                'occupations defined for more than 5 elements, which is not '
                'supported by TensErLEED.')
            return 1
        occlists = []
        error = False
        for k in at.disp_occ:
            occlists.append(at.disp_occ[k])
        for i in range(0,len(occlists[0])):
            totalocc = 0.
            for l in occlists:
                if len(l) <= i:
                    error = True
                    break
                else:
                    totalocc += l[i]
            if error:
                logger.error('DISPLACEMENTS file: Lengths of occupation '
                    'lists per element differ for atom '+at.oriN)
                return 1
            if totalocc > 1 + 1e-4:    # some tolerance
                logger.error('DISPLACEMENTS file: Occupations for atom '
                    +str(at.oriN)+' sum to more than 1 for step '+str(i+1)+'.')
                return 1
            elif totalocc < 1 - 1e-4 and len(at.disp_occ) > 4:
                logger.error('DISPLACEMENTS file: Occupations for atom '
                    +str(at.oriN)+' sum to less than 1 for step '+str(i+1)+', '
                    'but a vacancy cannot be added since there are already '
                    '5 elements.')
                return 1
    return 0

def checkVIBROCC(rp, slab, generate=False):
    """Fills default values and does consistency check of site vibrational
    amplitudes and occupations. If vibrational amplitudes are automatically
    calculated here, returns True, else returns False."""
    vibAmpGenerated = False
    for site in slab.sitelist:
        #check if the vibrational amplitude is defined for the main element(s)
        if not site.el in site.vibamp:
            if not site.el in rp.ELEMENT_MIX:
                if generate:
                    if site.getVibAmp(rp, site.el) == 0:
                        vibAmpGenerated = True
                    else:
                        logger.warning('VIBROCC file: Failed to get '
                            +site.el+' vibrational amplitude for site '
                            +site.label)
                        rp.setHaltingLevel(2)
                else:
                    logger.warning('VIBROCC file: No '+site.el+
                        ' vibrational amplitude defined for site '+site.label)
                    rp.setHaltingLevel(2)
            else:
                v = -1
                for subel in rp.ELEMENT_MIX[site.el]:
                    if subel in site.vibamp:
                        v = site.vibamp[subel]
                    elif generate:
                        if site.getVibAmp(rp, subel) == 0:
                            vibAmpGenerated = True
                            v = site.vibamp[subel]
                if v == -1:
                    logger.warning('VIBROCC file: No vibrational amplitude '
                                'defined for any of the main elements in site '
                                +site.el+'_'+site.name)
                    rp.setHaltingLevel(2)
                else:
                    # vibrational amplitudes were defined for some of the
                    #  main elements but not all. Fill the other values
                    for subel in rp.ELEMENT_MIX[site.el]:
                        if not subel in site.vibamp:
                            logger.warning('VIBROCC file: No '+subel
                                    +' vibrational amplitude defined for site '
                                    +site.el+'_'+site.name+'. Using '
                                    'vibrational amplitude from other element '
                                    'in site.')
                            rp.setHaltingLevel(1)
                            site.vibamp[subel] = v
        if site.el in rp.ELEMENT_MIX:  #just use first element in ELEMENT_MIX
            mainva = site.vibamp[rp.ELEMENT_MIX[site.el][0]]
        else:
            mainva = site.vibamp[site.el]
        #for other elements, fill vibrational elements with that of the main
        # element (probably not present)
        for el in slab.chemelem:
            if not el in site.vibamp:
                site.vibamp[el] = mainva
        # check and fill occupations:
        o = 0.
        for el in slab.chemelem:
            if el in site.occ:
                o += site.occ[el]
        if 'Vac' in site.occ:
            o += site.occ['Vac']
        if not site.el in site.occ:
            if site.el in rp.ELEMENT_MIX:
                if not any([(el in site.occ) 
                             for el in rp.ELEMENT_MIX[site.el]]):
                    for el in rp.ELEMENT_MIX[site.el]:
                        site.occ[el] = (1. - o) / len(rp.ELEMENT_MIX[site.el])
            else:
                site.occ[site.el] = 1. - o
            o = 1.
        for el in slab.chemelem:
            if not el in site.occ:
                site.occ[el] = 0.
        if o < 0.99:
            logger.debug('VIBROCC file: Site '+site.label+' has a total '
                'occupation less than one. Interpreting remaining {:.2f} '
                'as vacancies.'.format(1-o))
        elif o > 1.0:
            if o > 1.01:
                logger.warning('VIBROCC file: Site '+site.label+' has a '
                        +'total occupation greater than one ({:.2f}). '
                        'Occupations will be re-scaled to 1.'.format(o))
                rp.setHaltingLevel(2)
            for el in site.occ:
                site.occ[el] *= 1/o
        # finally, check whether we have any vibrational amplitudes or
        #  occupations assigned to non-existant elements, and if so, drop
        #  them and warn the user about it:
        dl = []
        for el in site.vibamp:
            if not el in slab.chemelem:
                logger.warning('VIBROCC file: Site '+site.el+'_'+site.name
                        +' has a vibrational amplitude defined for an unknown '
                        'element, which will be dropped ('+el+').')
                rp.setHaltingLevel(1)
                dl.append[el]
        for el in dl:
            site.vibamp.pop(el, None)
        dl = []
        for el in site.occ:
            if not el in slab.chemelem and el != "Vac":
                logger.warning('VIBROCC file: Site '+site.el+'_'+site.name
                        +' has an occupation defined for an unknown element, '
                        'which will be dropped ('+el+').')
                rp.setHaltingLevel(1)
                dl.append[el]
        for el in dl:
            site.occ.pop(el, None)
    logger.debug("VIBROCC value consistency check finished")
    return vibAmpGenerated

def readVIBROCC(rp, slab, filename='VIBROCC'):
    """Reads VIBROCC and adds the information to all sites in the slab.
    If vibrational amplitudes are automatically calculated here, returns True,
    else returns False."""
    if not (rp.T_EXPERIMENT is None or rp.T_DEBYE is None):
        generate = True
    else:
        generate = False
    vibAmpGenerated = False
    # open input file
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        if generate:
            logger.info("No VIBROCC file found, generating vibrational "
                         "amplitudes from temperature...")
            for site in slab.sitelist:
                if not site.el in rp.ELEMENT_MIX:
                    if site.getVibAmp(rp, site.el) != 0:
                        logger.error("Failed to generate vibrational "
                            "amplitude for "+site.el+" in site "+site.label)
                        raise
                else:
                    for subel in rp.ELEMENT_MIX[site.el]:
                        if site.getVibAmp(rp, subel) != 0:
                            logger.error("Failed to generate vibrational "
                                "amplitude for " + site.el + " in site "
                                + site.label)
                            raise
            try:
                checkVIBROCC(rp, slab, generate=True)
            except:
                raise
            return True
        else:
            logger.error("VIBROCC not found.")
            raise
    #read VIBROCC:
    mode = 0  #0: not reading; 1: vib.amp., 2: occ., 3: search offsets
    regex = False   #read regular expressions as-is or not
    for line in rf:
        line = line.lstrip()
        if "!" in line:
            line = line.split("!")[0]
        if not '=' in line:
            continue
        else:
            if line[0] == '=':
                llist = line[1:].split()
                if llist[0][0].lower() == 'v':
                    mode = 1
                elif llist[0][0].lower() == 'o':
                    mode = 2
                elif llist[0][0].lower() == 's':
                    mode = 3
                elif llist[0][0].lower() == 'r':
                    regex = True
                    if len(llist) >= 2:
                        if llist[1].lower() == 'off':
                            regex = False
                    continue
                else:
                    logger.warning("VIBROCC: Found line starting with '=', "
                                    "but didn't recognize vibrational "
                                    "amplitude or occupation")
                    rp.setHaltingLevel(1)
                continue
            elif mode == 0:
                continue
            else:
                param = line.split('=')[0]
                if not param:
                    continue
                else:
                    plist = param.split()
                    if plist:
                        param = plist[0]
                    else:
                        continue
        if mode in [1,2]:
            try:
                llist = tl.linelist(line.split('=')[1])
            except IndexError:
                logger.warning('VIBROCC file: ' + param + ' appears to have '
                                'no value')
                rp.setHaltingLevel(1)
                continue
            if not llist:
                logger.warning('VIBROCC file: ' + param + ' appears to have '
                                'no value')
                rp.setHaltingLevel(1)
                continue
            # first get values on the right
            sublists = tl.splitSublists(tl.readToExc(llist), ',')
            # read parameter on the left
            targetsites = []
            for site in slab.sitelist:
                if regex:         #if regular expressions are enabled, take the
                    prep = param     #  parameter input at face value
                else:
                    prep = re.escape(param)   #double-slash non-literal characters
                    #if regular expressions are not enabled, we want to still
                    #  interpret * as "any number of any characters":
                    prep = prep.replace('\\*','.*')
                m = re.match(prep, site.label)
                if m:
                    #if the length of the matched text == the site label,
                    #  it's a real match
                    if m.end(0) == len(site.label):
                        targetsites.append(site)
            if len(targetsites) == 0:
                logger.warning('VIBROCC file: No sites matching '+param
                                +' found, line will be skipped.')
                rp.setHaltingLevel(1)
                continue
            else:
                for site in targetsites:
                    if mode == 1:   #decide which dictionary to write to
                        td = site.vibamp
                    else:
                        td = site.occ
                    for sl in sublists:
                        if len(sl) == 1:
                            #if there's only one value, it should be a float
                            #  for the main site element
                            try:
                                if site.el in rp.ELEMENT_MIX:
                                    for subel in rp.ELEMENT_MIX[site.el]:
                                        td[subel] = float(sl[0])
                                else:
                                    td[site.el] = float(sl[0])
                            except:
                                logger.error('VIBROCC file: Error reading '
                                        'value '+sl[0]+' at parameter '+param)
                                raise
                        else:
                            el = sl[0].capitalize()
                            try:
                                td[el] = float(sl[1])
                            except KeyError:
                                logger.error('VIBROCC file: Element '+el+
                                        'not recognized at parameter '+param)
                            except:
                                logger.error('VIBROCC file: Error reading '
                                        'value '+sl[1]+' at parameter '+param)
                                raise
        if mode == 3:
            try:
                ind = int(plist[1])
            except:
                logger.error('VIBROCC file: Error converting to atom number: '
                              +plist[1])
                continue
            else:
                targetatlist = [at for at in slab.atlist if at.oriN == ind]
                if len(targetatlist) == 1:
                    targetat = targetatlist[0]
                else:
                    logger.error('VIBROCC file: {} atoms found with number {}'
                                  .format(len(targetatlist), ind))
                    continue
            if param.lower() in ["pos", "geo"]:
                om = 1
                targetdict = targetat.offset_geo
            elif param.lower() == "vib":
                om = 2
                targetdict = targetat.offset_vib
            elif param.lower() == "occ":
                om = 3
                targetdict = targetat.offset_occ
            else:
                logger.error('VIBROCC file: Flag not recognized: '+param)
                continue
            s = line.split("=")[1]
            subls = s.split(",")
            for l in subls:
                ll = l.split()
                if (len(ll) != 4 and om == 1) or (len(ll) != 2 and om != 1):
                    logger.error('VIBROCC file: Wrong number of values in '
                                  'sublist '+l)
                    continue
                else:
                    el = ll[0]
                    if not el in slab.chemelem:
                        logger.error('VIBROCC file: Element '+el+' not '
                                      'recognized')
                        continue
                    values = []
                    try:
                        for s in ll[1:]:
                            values.append(float(s))
                    except:
                        logger.error('VIBROCC file: Could not convert values '
                                      'to floats: '+l)
                        continue
            if om == 1:
                value = np.array(values)
                value[2] *= -1 # invert z
            else:
                value = values[0]
            targetdict[el] = value
    logger.debug("VIBROCC file was read successfully")
    # now fill up default values & do consistency checks:
    vibAmpGenerated = checkVIBROCC(rp, slab, generate=generate)
    if vibAmpGenerated:
        return True
    else:
        if rp.T_EXPERIMENT is not None:
            logger.warning("Parameter T_EXPERIMENT is defined but unused.")
        if rp.T_DEBYE is not None:
            logger.warning("Parameter T_DEBYE is defined but unused.")
    return False

def readIVBEAMS(filename='IVBEAMS'):
    """Reads an IVBEAMS file and returns a list of beams (using Beam class)"""
    # open input file
    try:
        with open(filename, 'r') as rf:
            ivbeamlines = rf.readlines()
    except:
        raise
    #read IVBEAMS:
    linenum = 1		# iterates the current line being read
    hklist = []
    for line in ivbeamlines:
        # ignore brackets and vbars, except as spacers
        line = line.replace("(", " ")
        line = line.replace(")", " ")
        line = line.replace("|", " ")
        llist = tl.linelist(line)
        if len(llist) == 1 and linenum != 1:
            logger.warning('A line with only one element was found in '
                            'IVBEAMS and will be skipped: '+line)
        elif len(llist) >= 2:
            f = [None,None]
            for i in range(0,2):
                try:
                    f[i] = float(llist[i])
                except ValueError:
                    if '/' in llist[i]:
                        try:
                            f[i] = (float(llist[i].split('/')[0])
                                    / float(llist[i].split('/')[1]))
                        except ValueError:
                            if linenum != 1:
                                logger.error('Error reading IVBEAMS line: '
                                              +line)
                                raise
                    else:
                        if linenum != 1:
                            logger.error('Error reading IVBEAMS line: '+line)
                            raise
            if linenum != 1:
                if not (f[0],f[1]) in hklist and not None in f:
                    hklist.append((f[0],f[1]))
            else:
                # check whether there is data in first line by mistake
                if not None in f:  # data was read
                    logger.warning('It looks like the first line in the '
                        'IVBEAMS file may contain data. Note that the first '
                        'line should be a header line, so the data will not '
                        'be read.')
        linenum += 1
    beams = []
    for hk in hklist:
        beams.append(tl.Beam(hk))
    logger.debug("IVBEAMS file was read successfully")
    return beams

def sortIVBEAMS(sl, rp):
    """Sorts the beams in IVBEAMS such that they appear in the same order as
    in the _BEAMLIST. Returns the sorted list."""
    # read BEAMLIST
    if rp.beamlist == []:
        logger.warning("sortIVBEAMS routine: no beamlist passed, "
            "attempting to read _BEAMLIST directly.")
        try:
            with open('_BEAMLIST', 'r') as rf:
                rp.beamlist = rf.readlines()
        except FileNotFoundError:
            logger.error("_BEAMLIST not found.")
            raise
    err = 1e-3           #since beams are saved as floats, give error tolerance
    symeq = tl.getSymEqBeams(sl, rp)
    # first, get beamlist as floats
    blfs = []
    for line in rp.beamlist:
        llist = tl.linelist(line)
        if len(llist) > 1:
            fl = []
            for i in range(0,2):
                try:
                    fl.append(float(llist[i]))
                except:
                    pass
            if len(fl) == 2:
                blfs.append(fl)
    # now figure out if there are beams in IVBEAMS without a partner
    for (ind, ib) in enumerate(rp.ivbeams):
        found = any([ib.isEqual_hk(lb, eps=err) for lb in blfs])
        if found:
            continue
        eqbl = []   # hk of equivalent beams
        for (hk, i) in symeq:
            if ib.isEqual_hk(hk, eps=err):
                eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
        for eqb in eqbl:
            for lb in blfs:
                if abs(eqb[0]-lb[0]) < err and abs(eqb[1]-lb[1]) < err:
                    # rename the beam
                    if any([ib2.isEqual_hk(eqb, eps=err) for ib2 
                                                          in rp.ivbeams]):
                        logger.debug("Beam " + ib.label + " is not in "
                        "BEAMLIST, but an equivalent beam already is. "
                        "Beam will be dropped.")
                    else:
                        b = tl.Beam(eqb)
                        logger.debug("Beam " + ib.label + " is not in "
                            "BEAMLIST, renaming to equivalent beam "
                             + b.label + ".")
                        rp.ivbeams[ind] = b
                    found = True
                    break
            if found:
                break
        if not found:
            logger.warning('IVBEAMS contains beam ' + ib.label + ', which '
                'was not found in the _BEAMLIST file. Beam will be dropped.')
    # now sort
    ivsorted = []
    for lb in blfs:
        ivsorted.extend([ib for ib in rp.ivbeams if ib.isEqual_hk(lb, 
                                                                  eps=err)])
    return ivsorted

def readPHASESHIFTS(sl, rp, readfile='_PHASESHIFTS', check=True):
    """"Reads from a _PHASESHIFTS file, returns the data as a list of tuples
    (E, enps), where enps is a list of lists, containing one list of values
    (for different L) for each element. Therefore, len(phaseshifts) is the
    number of energies found, len(phaseshifts[0][1]) should match the number
    of elements, and len(phaseshifts[0][1][0]) is the number of different
    values of L in the phaseshift file. The "check" parameters controls
    whether the phaseshifts that were found should be checked against the
    parameters / slab. If it is set to False, then passing "None" as sl and rp
    will work."""
    rf74x10 = ff.FortranRecordReader('10F7.4')
    ri3 = ff.FortranRecordReader('I3')

    try:
        rf = open(readfile, 'r')
    except FileNotFoundError:
        logger.error("_PHASESHIFTS file not found.")
        raise

    filelines = []
    for line in rf:
        filelines.append(line[:-1])
    rf.close()

    try:
        nel = ri3.read(filelines[0])[0]
    except:
        logger.error("Exception while trying to read _PHASESHIFTS: could not "
                      "find number of blocks in first line.")
        raise
    phaseshifts = []

    firstline = filelines[0]
    readline = 1
    linesperblock = 0
    while readline < len(filelines):
        if linesperblock:
            en = rf74x10.read(filelines[readline])[0]
            enps = []
            for i in range(0,nel):
                elps = []
                for j in range(0,linesperblock):
                    llist = rf74x10.read(filelines[readline
                                                   +(i*linesperblock)+j+1])
                    llist = [f for f in llist if f is not None]
                    elps.extend(llist)
                enps.append(elps)
            phaseshifts.append((en,enps))
            readline += linesperblock*nel+1
        else:
            #first check how many lines until the next energy:
            lineit = 1
            llist = rf74x10.read(filelines[readline+lineit])
            llist = [f for f in llist if f is not None]
            longestline = len(llist)
            shortestline = longestline
            lastlen = longestline
            cont = True
            while cont:
                lineit += 1
                llist = rf74x10.read(filelines[readline+lineit])
                llist = [f for f in llist if f is not None]
                if len(llist) == 1:
                    if lastlen == 1 or (shortestline > 1
                                        and shortestline < longestline):
                        cont = False  #found next energy
                    else:
                        shortestline = 1
                        #blocklines += 1
                elif len(llist) != longestline:
                    shortestline = len(llist)
                lastlen = len(llist)
            linesperblock = int((lineit-1)/nel)
            if not linesperblock or (((lineit-1)/nel) - linesperblock != 0.0):
                logger.warning("Error while trying to read _PHASESHIFTS: "
                    "Could not parse file: The number of blocks may not match "
                    "the number given in the first line. A new _PHASESHIFTS "
                    "file will be generated.")
                rp.setHaltingLevel(1)
                return ("", [], True, True)
            # don't increase readline -> read the same block again afterwards

    if not check:
        newpsGen, newpsWrite = False, False
    else:
        # check whether the phaseshifts that were found fit the data:
        newpsGen, newpsWrite = True, True
                    # recommend that new values should be generated / written
        psblocks = 0
        for el in sl.elements:
            if el in rp.ELEMENT_MIX:
                n = len(rp.ELEMENT_MIX[el])
            else:
                n = 1
            psblocks += n*len([s for s in sl.sitelist if s.el == el])
        # check for MUFTIN parameters:
        muftin = True
        llist = tl.linelist(firstline)
        if len(llist) >= 6:
            for i in range(1,5):
                try:
                    float(llist[i])
                except ValueError:
                    muftin = False
        else:
            muftin = False
        if rp.V0_REAL == "default" and not muftin:
            logger.warning("Could not convert first line of "
                "_PHASESHIFTS file to MUFTIN parameters. A new "
                "_PHASESHIFTS file will be generated.")
            rp.setHaltingLevel(1)
        elif len(phaseshifts[0][1]) == psblocks:
            logger.debug("Found "+str(psblocks)+" blocks in _PHASESHIFTS "
                          "file, which is consistent with PARAMETERS.")
            newpsGen, newpsWrite = False, False
        elif len(phaseshifts[0][1]) == sl.nelem:
            logger.warning("Found fewer blocks than expected in the "
                "_PHASESHIFTS file. However, the number of blocks matches "
                "the number of chemical elements. A new _PHASESHIFTS file "
                "will be generated, assuming that each block in the old "
                "file should be used for all atoms of one element.")
            rp.setHaltingLevel(1)
            oldps = phaseshifts[:]
            phaseshifts = []
            for (en,oldenps) in oldps:
                enps = []
                j = 0   #block index in old enps
                for el in sl.elements:
                    if el in rp.ELEMENT_MIX:
                        m = len(rp.ELEMENT_MIX[el])
                    else:
                        m = 1
                    n = len([s for s in sl.sitelist if s.el == el])
                    for i in range(0,m):    # repeat for chemical elements
                        for k in range(0, n):   # repeat for sites
                            enps.append(oldenps[j])
                        j += 1  # count up the block in old enps
                phaseshifts.append((en,enps))
            newpsGen = False
            firstline = str(len(phaseshifts[0][1])).rjust(3) + firstline[3:]
        else:
            logger.warning("_PHASESHIFTS file was read but is "
                "inconsistent with PARAMETERS. A new _PHASESHIFTS file "
                "will be generated.")
            rp.setHaltingLevel(1)

        # check whether energy range is large enough:
        checkfail = False
        er = np.arange(rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+1e-4, 
                       rp.THEO_ENERGIES[2])
        psmin = round(phaseshifts[0][0]*27.2116, 2)
        psmax = round(phaseshifts[-1][0]*27.2116, 2)
        if rp.V0_REAL == "default":
            llist = tl.linelist(firstline)
            c = []
            try:
                for i in range(0,4):
                    c.append(float(llist[i+1]))
            except:
                checkfail = True
            else:
                er_inner = [e + (rp.FILAMENT_WF - max(c[0], 
                                        c[1] + (c[2]/np.sqrt(e + c[3] 
                                                          + rp.FILAMENT_WF))))
                            for e in er] # energies at which scattering occurs
        else:
            try:
                v0r = float(rp.V0_REAL)
            except:
                checkfail = True
            else:
                er_inner = [e + v0r for e in er]
        if not checkfail:
            if (psmin > min(er_inner) or psmax < max(er_inner)):
                if (psmin > min(er_inner) and psmin <= 20. 
                                          and psmax >= max(er_inner)):
                    # can lead to re-calculation of phaseshifts every run if 
                    #  V0r as calculated by EEASiSSS differs from 'real' V0r.
                    #  Don't automatically correct.
                    logger.warning("Lowest value in _PHASESHIFTS file ({:.1f} "
                        "eV) is larger than the lowest predicted scattering "
                        "energy ({:.1f} eV). If this causes problems in the "
                        "reference calculation, try deleting the _PHASESHIFTS "
                        "file to generate a new one, or increase the starting "
                        "energy in the THEO_ENERGIES parameter."
                        .format(psmin, min(er_inner)))
                else:
                    logger.warning("The energy range found in the _PHASESHIFTS"
                        " file is smaller than the energy range requested for "
                        "theoretical beams. A new _PHASESHIFTS file will be "
                        "generated.")
                    newpsGen, newpsWrite = True, True
        else:
            logger.warning("Could not check energy range in _PHASESHIFTS "
                "file. If energy range is insufficient, try deleting the "
                "_PHASESHIFTS file to generate a new one.")
    return (firstline, phaseshifts, newpsGen, newpsWrite)

def readFdOut(readfile="fd.out"):
    """Reads the fd.out file produced by the refcalc and returns a list of
    Beam objects."""
    try:
        with open(readfile, 'r') as rf:
            filelines = [l[:-1] for l in rf.readlines()]
            rf.seek(0)
            fdout = rf.read()
    except FileNotFoundError:
        logger.error("Error in readFdOut: file "+str(readfile)+" not found.")
        raise
    if filelines[0].startswith(" SOME ERROR OCCURED WHILE READING"):
        logger.error("File "+readfile+" reports error: "+filelines[0])
        raise Exception("File "+readfile+" reports error.")
    if ".  CORRECT TERMINATION" in fdout:   # happens for superpos output...
        fdout = fdout.split(".  CORRECT TERMINATION")[0]
    theobeams = []
    i = 1   # number lines as in text editor - be careful about indexing!
    nbeams = 1e10   # some large number, just to enter the while loop
    while i < nbeams+3:
        llist = tl.linelist(filelines[i-1])
        if i == 1:
            pass #header - skip
        elif i == 2:
            nbeams = int(llist[0])
        else:
            theobeams.append(tl.Beam((float(llist[1]),float(llist[2]))))
        i += 1

    # re-label the beams to get the correct number of characters and formatting
    mw = max([beam.lwidth for beam in theobeams])
    for beam in theobeams:
        beam.lwidth = mw
        beam.getLabel()

    # From now on, all the lines correspond to blocks of beam intensities.
    # Collect them, skipping empty lines
    blocks = []
    for line in filelines[i-1:]:
        llist = tl.linelist(line)
        if len(llist) > 0:  # skip empty lines
            try:
                float(llist[0])
            except:
                break  # end of data
            blocks.extend(llist)

    # Each block contains exactly nbeams+2 entries
    # reshape blocks so that each line corresponds to one block
    blocks = np.reshape(blocks, (-1, nbeams+2))

    # and parse the data
    for block in blocks:
        en = float(block[0])
        values = [float(s) for s in block[2:]]
        for (j, beam) in enumerate(theobeams):
            beam.intens[en] = values[j]
    return theobeams, fdout

def checkDelta(filename, at, el, rp):
    """Checks whether a given delta file corresponds to the requested 
    displacements of a given atom. Returns True or False."""
    eps = 1e-4
    fgeo = []                               # found geo disp
    fvib = []                               # found vib disp
    at.mergeDisp(el)
    if el == "vac":
        dgeo = [np.array([0.0, 0.0, 0.0])]  # requested geo disp
        dvib = [0.0]                        # requested vib disp
    else:
        if el in at.disp_geo:
            dgeo = at.disp_geo[el]
        else:
            dgeo = at.disp_geo["all"]
        if el in at.disp_vib:
            dvib = at.disp_vib[el]
        else:
            dvib = at.disp_vib["all"]
    try:
        with open(filename, "r") as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.error("Error reading file "+filename)
        raise
    try:
        nbeams = int(lines[1][0:3])  # number of beams
        nvar = int(lines[1][6:9])    # number of variations (geo*vib)
    except:
        logger.error("Error parsing file " + filename)
        raise
    if nbeams != len(rp.ivbeams):
        return False
    if nvar != len(dgeo)*len(dvib):
        return False
    beams = []      # read beams from delta file
    atline = 2
    for i in range(atline, len(lines)):
        try:
            fl = [float(s) for s in lines[i].split()]
        except:
            logger.error("Error parsing file "+filename)
            raise
        for j in range(0, int(len(fl)/2)):
            beams.append((fl[2*j], fl[(2*j) + 1]))
        if len(beams) == nbeams:
            atline = i+1
            break
    for (i, hk) in enumerate(beams):   # check beams
        if not rp.ivbeams[i].isEqual_hk(hk, eps=eps):
            return False
    atline += 1   # skip the line after beams
    # geo displacements start here
    rf74x10 = ff.FortranRecordReader("10F7.4")
    parselist = []
    repeats = False
    endgeo = False
    entrycount = 0
    for i in range(atline, len(lines)):
        try:
            fl = [f for f in rf74x10.read(lines[i]) if f is not None]
        except:
            logger.error("Error parsing file "+filename)
            raise
        if len(fl) < 10:
            atline = i+1
            endgeo = True  # short line -> end of geo block
        parselist = parselist + fl
        while len(parselist) >= 3:
            v = parselist[:3]
            new = np.array([v[1], v[2], v[0]])
            parselist = parselist[3:]
            entrycount += 1
            if not repeats:
                append = True
                if fgeo and np.linalg.norm(new - fgeo[0]) < eps:
                    repeats = True
                    append = False
                if append:
                    fgeo.append(new)
                if len(fgeo) == nvar:
                    atline = i+1
                    endgeo = True
                    break
            elif not endgeo:
                if (entrycount-1) % len(fgeo) == 0:  # should repeat here
                    if np.linalg.norm(new - fgeo[0]) > eps:
                        atline = i
                        endgeo = True
                        break
        if endgeo:
            break
    ngeo = len(fgeo)
    # check geometry
    if ngeo != len(dgeo):
        return False
    for (i, gd) in enumerate(fgeo):
        if np.linalg.norm(gd - dgeo[i]) > eps:
            return False
    # vib displacement starts here
    nvib = nvar / ngeo
    if int(nvib) - nvib > 1e-4:
        logger.error("Error reading file "+filename+": number of geometry "
                      "variations found does not match header.")
        return False
    nvib = int(nvib)
    for i in range(atline, len(lines)):
        try:
            fl = [f for f in rf74x10.read(lines[i]) if f is not None]
        except:
            logger.error("Error parsing file "+filename)
            raise
        parselist = parselist + fl
        while len(parselist) >= ngeo:
            fvib.append(parselist[0])
            if any([f != 0. for f in parselist[1:ngeo]]):
                logger.warning("File "+filename+": Found unexpected entries "
                        "in list of vibrational displacements.")
            parselist = parselist[ngeo:]
        if len(fvib) >= nvib:
            break
    # check vibrations:
    if el.lower() == "vac":
        voff = 0.
    else:
        voff = at.site.vibamp[el]
    for (i, f) in enumerate(fvib):
        bv = round(round(dvib[i] + voff, 4) / 0.529, 4)  # in bohr radii
                    # rounding twice to account for 1. writing to delta-input,
                    # 2. reading from DELTA file. Precision 0.529 taken from 
                    # TensErLEED GLOBAL
        if abs(f - bv) >= 1e-4:
            return False
    return True

def readOUTBEAMS(filename="EXPBEAMS.csv", sep=";", enrange=[]):
    """Reads beams from an EXPBEAMS.csv or THEOBEAMS.csv file. Returns a list
    of Beam objects. The 'sep' parameter defines the separator. If an energy
    range 'enrange' is passed, beams that contain no data within that range
    will be filtered out before returning. If a Slab and Rparams object are 
    passed, will also check whether any of the experimental beams should be 
    equivalent, and if they are, warn, discard one and raise the halting 
    level."""
    beams = []
    try:
        with open(filename, 'r') as rf:
            lines = [l[:-1] for l in rf.readlines()]
    except FileNotFoundError:
        logger.error("Error reading "+filename)
        raise
    firstline = True
    rgx = re.compile(r'[\*\(\s]*(?P<h>[-0-9/]+)\s*\|\s*(?P<k>[-0-9/]+)')
    for line in lines:
        if firstline and sep not in line:   # try some other separators
            for sep2 in [s for s in [";", ","] if s != sep]:
                if sep2 in line and len(line.split(sep2)) > 2:
                    logger.info("Found separator '"+sep2+"' in "+filename
                                 +", expected '"+sep+"'. Attempting to read "
                                 "with '"+sep2+"'.")
                    sep = sep2
                    break
        llist = line.split(sep)
        if llist[-1].strip() == "":
            llist = llist[:-1]
        if firstline:
            firstline = False
            for label in llist[1:]:
                m = rgx.match(label)
                if m == None:
                    logger.error("readOUTBEAMS: Could not parse h/k in "
                                  "label: "+label)
                    return []
                sh = m.group("h")   #string h
                sk = m.group("k")   #string k
                try:
                    if "/" in sh:
                        h = int(sh.split("/")[0]) / int(sh.split("/")[1])
                    else:
                        h = int(sh)
                    if "/" in sk:
                        k = int(sk.split("/")[0]) / int(sk.split("/")[1])
                    else:
                        k = int(sk)
                except:
                    logger.error("readOUTBEAMS: Could not parse h/k in "
                                    "label: "+label)
                    return []
                beams.append(tl.Beam((h,k)))
        elif len(line) > 1:
            try:
                en = float(llist[0])
            except:
                logger.error("readOUTBEAMS: Could not parse "+llist[0]+"as "
                              "an energy")
                return []
            for i in range(0,len(llist)):
                try:
                    f = float(llist[i+1])
                    if not np.isnan(f):
                        beams[i].intens[en] = f
                except:
                    f = None
    if len(enrange) == 2:
        remlist = []
        for b in beams:
            if ((enrange[0] > 0 and max(b.intens) < enrange[0]) or
                (enrange[1] > 0 and min(b.intens) > enrange[1])):
                # beam has no data in given interval; remove
                remlist.append(b)
                logger.warning("Experimental beam "+b.label+" contains no "
                                "data in the given energy range. The beam "
                                "will be ignored.")
        for b in remlist:
            beams.remove(b)
    # if sl and rp:           # check for equivalent beams                       # !!! CLEANUP
    #     remlist = []
    #     symeq = tl.getSymEqBeams(sl, rp)
    #     for (bi, b) in enumerate(beams):
    #         if b in remlist:
    #             continue
    #         eqbl = []   # hk of equivalent beams
    #         for (hk, i) in symeq:
    #             if b.isEqual_hk(hk):
    #                 eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
    #                 break
    #         for b2 in beams[bi+1:]:
    #             for hk in eqbl:
    #                 if b2.isEqual_hk(hk):
    #                     remlist.append(b2)
    #                     rp.setHaltingLevel(2)
    #                     w = ("Experimental beam "+b2.label+" is "
    #                         "equivalent to experimental beam "+b.label+". ")
    #                     if rp.HALTING > 2:
    #                         w += "Beam "+b2.label+" will be discarded."
    #                     else:
    #                         w += ("Remove one of them, or average them as "
    #                               "appropriate.")
    #                     logger.warning(w)
    #                     break
    #     for b in remlist:
    #         beams.remove(b)
    totalrange = 0.
    minmax = enrange[:]
    if len(minmax) < 2:
        minmax = [0, 1e5]
    if minmax[0] < 0:
        minmax[0] = 0
    if minmax[1] <= 0:
        minmax[1] = 1e5
    for b in beams:
        totalrange += (min(max(b.intens), minmax[1]) 
                       - max(min(b.intens), minmax[0]))
    logger.info("Loaded "+filename+" file containing {} beams (total energy "
                "range: {:.2g} eV).".format(len(beams), totalrange))
    return beams

def checkEXPBEAMS(sl, rp):
    remlist = []
    symeq = tl.getSymEqBeams(sl, rp)
    for (bi, b) in enumerate(rp.expbeams):
        if b in remlist:
            continue
        eqbl = []   # hk of equivalent beams
        for (hk, i) in symeq:
            if b.isEqual_hk(hk):
                eqbl.extend([hk2 for (hk2, j) in symeq if i == j])
                break
        for b2 in rp.expbeams[bi+1:]:
            for hk in eqbl:
                if b2.isEqual_hk(hk):
                    remlist.append(b2)
                    w = ("Experimental beam "+b2.label+" is "
                        "equivalent to experimental beam "+b.label+". ")
                    if rp.HALTING > 2:
                        w += "Beam "+b2.label+" will be discarded."
                    else:
                        w += ("Remove one of them, or average them as "
                              "appropriate.")
                    logger.warning(w)
                    rp.setHaltingLevel(2)
                    break
    for b in remlist:
        rp.expbeams.remove(b)

def readAUXEXPBEAMS(filename="AUXEXPBEAMS", interactive=False):
    """Reads beams from an AUXEXPBEAMS file, which already has the formatting
    required by TensErLEED. Returns a list of Beam objects. Only works if the
    comment lines are of format *( h k )*, including the stars, with h and k
    separated by whitespace. If 'interactive' is True, then the user will be 
    asked to clarify labels that cannot be read."""
    expbeams = []
    try:
        with open(filename, 'r') as rf:
            lines = [l[:-1] for l in rf.readlines()]
    except FileNotFoundError:
        logger.error("Error reading AUXEXPBEAMS.")
        raise
    read = False
    rf62x12 = ff.FortranRecordReader('12F6.2')
    rgx = re.compile(r'[\*\(\s]*(?P<h>[-0-9/]+)\s+(?P<k>[-0-9/]+)')
    for line in lines:
        if "*" in line:
            read = True
            topline = True  #next line contains number of beams and scaling
            failedToRead = False
            m = rgx.match(line)
            if m != None:
                sh = m.group("h")   #string h
                sk = m.group("k")   #string k
                try:
                    h = tl.parseMathSqrt(sh)
                    k = tl.parseMathSqrt(sk)
                except:
                    failedToRead = True
                else:
                    newbeam = tl.Beam((h, k))
                    expbeams.append(newbeam)
            else:
                failedToRead = True
            if failedToRead:
                if not interactive:
                    logger.error("readAUXEXPBEAMS: Could not parse h/k in "
                                  "line:\n"+line)
                    return []
                else:
                    print("readAUXEXPBEAMS: Could not parse h/k in "
                                  "line:\n"+line)
                    while True:
                        try:
                            hks = input("Please manually input the beam "
                                        "indices (space-separated): ")
                        except:
                            return []
                        if hks.lower() in ["quit", "exit"]:
                            return []
                        if hks and len(hks.split()) > 1:
                            try:
                                h = tl.parseMathSqrt(hks.split()[0])
                                k = tl.parseMathSqrt(hks.split()[1])
                            except:
                                print("Could not parse h/k")
                            else:
                                newbeam = tl.Beam((h, k))
                                expbeams.append(newbeam)
                                break
        elif read and topline:
            topline = False
            llist = tl.linelist(line)
            try:
                scaling = float(llist[1])
            except:
                logger.error("readAUXEXPBEAMS: Could not parse number of "
                              "beams or scaling factor in line:\n"+line)
                return []
        elif read:
            newvals = rf62x12.read(line)
            i = 0
            while i+1 < len(newvals):
                if newvals[i] == None:
                    break
                newbeam.intens[newvals[i]] = newvals[i+1] / scaling
                i += 2
    return expbeams

def readSDTL_next(filename="SD.TL", offset=0):
    """Reads SDTL from offset to end, returns new offset and new content in
    between as string."""
    try:
        with open(filename, "r") as rf:
            rf.seek(offset)
            content = rf.read()
            newoffset = rf.tell()
        return(newoffset, content)
    except:
        logger.error("Error reading SD.TL file")
        return (offset, "")     # return old offset, no content

def readSDTL_blocks(content, whichR = 0):
    """Attempts to interpret a given string as one or more blocks of an SD.TL
    file. Returns a list with one entry per block: (gen, rfacs, configs),
    where gen is the generation number and rfacs is a list of R-factors in
    that block, and configs is a list of lists of parameter values. Optional
    argument whichR specifies whether to use integer or fractional beams
    instead of average."""
    returnList = []
    blocklist = content.split("CCCCCCCCCCC    GENERATION")[1:]
    for block in blocklist:
        gen = 0
        try:
            gen = int(block[:13])
        except:
            logger.warning("While reading SD.TL, could not interpret "
                "generation number "+block[:13])
        lines = block.split("\n")
        rfacs = []
        configs = []
        for line in lines:
            if "|" in line and line[:3] != "IND":
                # this is a line with data in it
                try:
                    rav = float(line.split("|")[2 + whichR]) #average R-factor
                    rfacs.append(rav)
                    valstring = line.split("|")[-1].rstrip()
                    pars = []
                    while len(valstring) > 0:
                        pars.append(int(valstring[:4]))
                        valstring = valstring[4:]
                    configs.append(pars)
                except:
                    logger.error("Could not read line in SD.TL:\n"+line)
        if gen != 0 and len(rfacs) > 0 and len(configs) > 0:
            returnList.append((gen, rfacs, configs))
        else:
            logger.warning("A block in SD.TL was read but not understood.")
    return returnList

def readSDTL_end(filename="SD.TL"):
    """Reads the last generation block from the SD.TL file."""
    # get the last block from SD.TL:
    bwr = tl.BackwardsReader(filename)
    lines = [""]
    while not "CCCCCCCCCCC    GENERATION" in lines[-1] and len(bwr.data) > 0:
        lines.append(bwr.readline())
    lines.reverse()
    bwr.close()
    return lines

def readROUT_end(filename="ROUT"):      # UNUSED
    """Reads the last line of ROUT, returns the average R-factor as float."""
    bwr = tl.BackwardsReader(filename)
    line = ""
    while not "AVERAGE R-FACTOR =" in line and len(bwr.data) > 0:
        line = bwr.readline()
    bwr.close()
    rfac = 0
    v0rshift = 0
    rgx = re.compile(r'.+SHIFT\s*(?P<shift>[-0-9.]+).+=\s*(?P<rfac>[-0-9.]+)')
    m = rgx.match(line)
    if m:
        try:
            rfac = float(m.group("rfac"))
        except:
            logger.error("Could not read R-factor from "+filename)
        try:
            v0rshift = float(m.group("shift"))
        except:
            logger.error("Could not read inner potential shift from "
                          + filename)
    return rfac, v0rshift

def readROUT(filename="ROUT"):
    """
    Reads the ROUT file.

    Parameters
    ----------
    filename : string, optional
        DESCRIPTION. Pass if you want to read from a file other than the 
        default 'ROUT'

    Returns
    -------
    tuple (r, r_int, r_frac), float v0rshift, list rfaclist
        r, r_int, r_frac: tuple of floats:
            average R-factor for all beams, integer-order beams and 
            fractional order beams.
        v0rshift: inner potential shift corresponding to r, r_int, r_frac
        rfaclist: list of floats, r-factors per beam in order of experimental
            beams

    """
    with open(filename, 'r') as rf:
        lines = rf.readlines()
    line = ""
    i = 0
    while not "AVERAGE R-FACTOR =" in line and i+1 < len(lines):
        i += 1
        line = lines[-i]
    if line == "":
        return 0, 0, []
    rfac = 0
    rfac_int = -1
    rfac_frac = -1
    v0rshift = 0
    rgx = re.compile(r'.+SHIFT\s*(?P<shift>[-0-9.]+).+=\s*(?P<rfac>[-0-9.]+)')
    m = rgx.match(line)
    if m:
        try:
            rfac = float(m.group("rfac"))
        except:
            logger.error("Could not read R-factor from "+filename)
        try:
            v0rshift = float(m.group("shift"))
        except:
            logger.error("Could not read inner potential shift from "
                          + filename)
    else:
        return (0,0,0), 0, []
    # now read the R-factors per beam at v0rshift
    rfaclist = []
    for line in [l for l in lines if len(l) >= 70]:
        values = line[18:69].split()  # limits to 999 beams;
                                      #  use [17:69] if more are required
        try:
            index = int(values[0])
            v0r = float(values[2])
            rav = float(values[-1])
        except:
            pass    # ignore line
        else:
            if v0r == v0rshift and index > 0:
                if index != len(rfaclist)+1:
                    logger.warning("Unexpected index mismatch in readROUT. "
                                    "Reading R-factors per beam will fail.")
                else:
                    rfaclist.append(rav)
            elif v0r == v0rshift and index == -1:
                if line.startswith("AV.-INT"):
                    rfac_int = rav
                elif line.startswith("AV.-FRAC"):
                    rfac_frac = rav
    return (rfac, rfac_int, rfac_frac), v0rshift, rfaclist

def getTensorOriStates(sl, path):
    """Reads POSCAR, PARAMETERS and VIBROCC from the target path, gets the 
    original state of the atoms and sites, and stores them in the given 
    slab's atom/site oriState variables."""
    for fn in ["POSCAR", "PARAMETERS", "VIBROCC"]:
        if not os.path.isfile(os.path.join(path, fn)):
            logger.error("File "+fn+" is missing in "+path)
            return("Could not check Tensors: File missing")
    # loglevel = logging.getLogger("tleedm").level
    loglevel = logger.level
    logger.setLevel(logging.ERROR)
    dn = os.path.basename(path)
    try:
        tsl = tl.readPOSCAR(os.path.join(path, "POSCAR"))
        trp = tl.readPARAMETERS(slab = tsl, filename = 
                                os.path.join(path, "PARAMETERS"))
        tsl.fullUpdate(trp)
        tl.readVIBROCC(trp, tsl, filename = os.path.join(path, "VIBROCC"))
        tsl.fullUpdate(trp)
    except:
        logger.error("Error checking Tensors: Error while reading "
                      "input files in "+dn)
        return("Could not check Tensors: Error loading old input "
               "files")
    finally:
        # logging.getLogger("tleedm").setLevel(loglevel)
        logger.setLevel(loglevel)
    if len(tsl.atlist) != len(sl.atlist):
        logger.error("POSCAR from "+dn+" is incompatible with "
                      "current POSCAR.")
        return("Tensors file incompatible")
    for at in sl.atlist:
        tal = [tat for tat in tsl.atlist if at.oriN == tat.oriN]
        if len(tal) != 1:
            logger.error("POSCAR from "+dn+" is incompatible with "
                          "current POSCAR.")
            return("Tensors file incompatible")
        at.copyOriState(tal[0])
    if len(tsl.sitelist) != len(sl.sitelist):
        logger.error("Sites from "+dn+" input differ from current "
                      "input.")
        return("Tensors file incompatible")
    for site in sl.sitelist:
        tsitel = [s for s in tsl.sitelist if site.label == s.label]
        if len(tsitel) != 1:
            logger.error("Sites from "+dn+" input differ from "
                          "current input.")
            return("Tensors file incompatible")
        site.oriState = copy.deepcopy(tsitel[0])
    return 0