# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer

Functions for reading from and writing to the PARAMETERS file
"""

import logging
import numpy as np
import re
import shutil
import os

import tleedmlib as tl

logger = logging.getLogger("tleedm.files.parameters")

# list of allowed parameters
knownParams = ['ATTENUATION_EPS', 'BEAM_INCIDENCE', 'BULKDOUBLING_EPS',
    'BULKDOUBLING_MAX', 'BULK_REPEAT', 'DOMAIN', 'ELEMENT_MIX', 
    'ELEMENT_RENAME', 'FILAMENT_WF', 'FORTRAN_COMP', 'HALTING', 
    'IV_SHIFT_RANGE', 'LAYER_CUTS', 'LAYER_STACK_VERTICAL', 'LMAX', 
    'LOG_DEBUG', 'LOG_SEARCH', 'N_BULK_LAYERS', 'N_CORES', 'PHASESHIFT_EPS', 
    'PLOT_COLORS_RFACTOR', 'RUN', 'R_FACTOR_SMOOTH', 'R_FACTOR_TYPE', 
    'SCREEN_APERTURE', 'SEARCH_BEAMS', 'SEARCH_CONVERGENCE', 'SEARCH_CULL', 
    'SEARCH_MAX_GEN', 'SEARCH_POPULATION', 'SEARCH_START', 'SITE_DEF', 
    'SUPERLATTICE', 'SUPPRESS_EXECUTION', 'SYMMETRIZE_INPUT', 'SYMMETRY_EPS', 
    'SYMMETRY_FIND_ORI', 'SYMMETRY_FIX', 'TENSOR_INDEX', 'TENSOR_OUTPUT', 
    'THEO_ENERGIES', 'T_DEBYE', 'T_EXPERIMENT', 'V0_IMAG', 'V0_REAL', 
    'V0_Z_ONSET', 'VIBR_AMP_SCALE']
paramAlias = {}
for p in knownParams:
    paramAlias[p.lower().replace("_","")] = p

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
        if (param not in knownParams and 
                param.lower().replace("_","") in paramAlias):
            param = paramAlias[param.lower().replace("_","")]
        if param not in knownParams:
            continue
        try:
            value = line.split('=', maxsplit=1)[1].rstrip()
            llist = value.split()  #read the stuff to the right of "="
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


def readPARAMETERS(filename='PARAMETERS'):
    """Reads a PARAMETERS file and returns an Rparams object with the
    raw input, without interpretation"""
    # open input file
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
        line = line.lstrip()
        if line.lstrip().startswith("SEARCH_KILL"):
            if not re.match(r"\s*SEARCH_KILL\s*=\s*[Ff](alse)?", line):
                logger.warning('PARAMETERS file: SEARCH_KILL is set at start '
                        'of program. This means the search will be stopped as '
                        'soon as it starts. Delete SEARCH_KILL from '
                        'PARAMETERS to avoid this.')
        if not "=" in line:
            continue     #ignore all lines that don't have an "=" sign at all
        param = line.split('=')[0]      #parameter is defined left of "="
        if not param:
            continue
        #get rid of spaces and check the leftmost entry.
        plist = param.split()
        param = plist[0]
        if (param not in knownParams and 
                param.lower().replace("_","") in paramAlias):
            param = paramAlias(param.lower().replace("_",""))
        if param not in knownParams:
            logger.warning('PARAMETERS file: Parameter '+param+' not '
                            'recognized.')
            rpars.setHaltingLevel(1)
            continue
        value = line.split('=', maxsplit=1)[1].rstrip()
        if not value:
            logger.warning('PARAMETERS file: ' + param + ' appears to '
                            'have no value')
            rpars.setHaltingLevel(1)
            continue
        if not param in rpars.readParams:
            rpars.readParams[param] = []
        rpars.readParams[param].append((plist, value))
    rf.close()
    return rpars
    
def interpretPARAMETERS(rpars, slab=None, silent=False):
    """Interprets the string values in an Rparams object, read previously by 
    readPARAMETERS, to fill the parameter variables. Returns 0 on success."""
    loglevel = logger.level
    if silent:
        logger.setLevel(logging.ERROR)
    # track some parameters while reading
    searchConvRead = False
    # define order that parameters should be read in
    orderedParams = ["LOG_DEBUG", "RUN"]
    checkParams = [p for p in orderedParams if p in rpars.readParams]
    checkParams.extend([p for p in knownParams if (p in rpars.readParams and 
                                                     not p in checkParams)])
    domainsIgnoreParams = ['BULK_REPEAT', 'ELEMENT_MIX', 'ELEMENT_RENAME',
        'LAYER_CUTS', 'N_BULK_LAYERS', 'SITE_DEF', 'TENSOR_INDEX', 
        'TENSOR_OUTPUT']
    for (param, plist, value) in [(param, lel[0], lel[1]) 
                                  for param in checkParams
                                  for lel in rpars.readParams[param]]:
        if rpars.hasDomains and param in domainsIgnoreParams:
            continue
        llist = value.split()
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
                    sublists = tl.base.splitSublists(llist, ',')
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
                vec = tl.base.readVector(s, slab.ucell)
                if vec is None:
                    logger.warning('PARAMETERS file: BULK_REPEAT: '
                        'Could not parse input expression. Input will '
                        'be ignored.')
                else:
                    rpars.BULK_REPEAT = vec
        elif param == 'DOMAIN':
            # check name
            if len(plist) > 1:
                name = plist[1]
            else:
                name = ""
            names = [n for (n, _) in rpars.DOMAINS]
            if name in names:
                logger.warning('PARAMETERS file: Multiple sources defined '
                    'for DOMAIN {}. Last entry will be ignored.'.format(name))
                rpars.setHaltingLevel(2)
                continue
            if not name:  # get unique name
                i = 1
                while str(i) in names:
                    i += 1
                name = str(i)
            # check path
            if os.path.exists(value):
                path = value
            elif os.path.isfile(value + ".zip"):
                path = value + ".zip"
            else:
                logger.warning('PARAMETERS file: Value for DOMAIN {} could '
                        'not be interpreted as either a path or a .zip file.'
                        .format(name))
                rpars.setHaltingLevel(2)
                continue
            rpars.DOMAINS.append((name, path))
        elif param == 'ELEMENT_MIX':
            ptl = [el.lower() for el in tl.leedbase.periodic_table]
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
            ptl = [el.lower() for el in tl.leedbase.periodic_table]
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
                if not silent:
                    logger.setLevel(logging.DEBUG)
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
                l = tl.base.readIntRange(s)
                if len(l) > 0:
                    rl.extend(l)
                else:
                    logger.warning('PARAMETERS file: RUN: Could not '
                            'interpret value '+s+', skipping value...')
                    rpars.setHaltingLevel(2)
            if len(rl) > 0:
                if 4 in rl:
                    logger.info('Found domain search.')
                    rpars.hasDomains = True
                i = 0
                while i < len(rl):
                    if rl[i] not in [0,1,2,3,4,11,31]:
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
        elif param == 'SCREEN_APERTURE':
            try:
                f = float(llist[0])
            except:
                logger.warning('PARAMETERS file: SCREEN_APERTURE: could not '
                               'convern value to float. Input will be '
                               'ignored.')
                rpars.setHaltingLevel(1)
            else:
                if 0 <= f <= 180:
                    rpars.SCREEN_APERTURE = f
                else:
                    logger.warning('PARAMETERS file: SCREEN_APERTURE value '
                        'must be between 0 and 180 degrees. Input {} will be '
                        'ignored.'.format(llist[0]))
                    rpars.setHaltingLevel(1)
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
                    if not searchConvRead:  # clear default values
                        rpars.SEARCH_MAX_DGEN = {"all": 0, "best": 0, "dec": 0}
                        searchConvRead = True
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
            sublists = tl.base.splitSublists(llist, ',')
            for sl in sublists:
                atnums = []
                for i in range(1, len(sl)):
                    l = tl.base.readIntRange(sl[i])
                    if len(l) > 0:
                        atnums.extend(l)
                    elif "top(" in sl[i]:
                        if slab is None:
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
                if slab is None:
                    logger.warning('PARAMETERS file: SUPERLATTICE '
                            'parameter appears to be in Wood notation, '
                            'but no slab was passed; cannot calculate '
                            'bulk unit cell!')
                    rpars.setHaltingLevel(2)
                else:
                    rpars.SUPERLATTICE = tl.leedbase.readWoodsNotation(value,
                                                              slab.ucell)
                    rpars.superlattice_defined = True
            else:
                sublists = tl.base.splitSublists(llist, ',')
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
            nl = tl.base.recombineListElements(llist, '*')
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
    logger.setLevel(loglevel)
    return 0

def modifyPARAMETERS(rp, modpar, new="", comment=""):
    """Looks for 'modpar' in the PARAMETERS file, comments that line out, and 
    replaces it by the string specified by 'new'"""
    oriname = "PARAMETERS_ori_"+rp.timestamp
    if not oriname in rp.manifest:
        try:
            shutil.copy2("PARAMETERS", oriname)
        except:
            logger.error("modifyPARAMETERS: Could not copy PARAMETERS file "
                "to PARAMETERS_ori. Proceeding, original file will be lost.")
        rp.manifest.append(oriname)
    if not "PARAMETERS" in rp.manifest:
        rp.manifest.append("PARAMETERS")
    output = ""
    headerPrinted = False

    try:
        with open("PARAMETERS", "r") as rf:
            plines = rf.readlines()
    except:
        logger.error("Error reading PARAMETERS file.")
        raise
    found = False
    for line in plines:
        if "! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #" in line:
            headerPrinted = True
        valid = False
        param = ""
        if "=" in line:    #ignore all lines that don't have an "=" sign at all
            param = line.split('=')[0]        #parameter is defined left of "="
            if param: 
                valid = True
                plist = param.split()  
                if plist: 
                    param = plist[0]
                if param[0] == '!': 
                    valid = False
        if valid and param == modpar:
            found = True
            if new:
                if comment == "":
                    comment = "line automatically changed to:"
                output += "!"+line[:-1] + " ! " + comment + "\n"
                output += new + "\n"
            else:
                if comment == "":
                    comment = "line commented out automatically"
                output += ("!"+line.rstrip()).ljust(35)+ " ! " + comment + "\n"
        else:
            output += line
    if not found:
        if not headerPrinted:
            output += """

! ######################################################
! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #
! ######################################################
"""
            headerPrinted = True
        output += "\n" + modpar + " = " + new + " ! " + comment
    try:
        with open("PARAMETERS", "w") as wf:
            wf.write(output)
    except:
        logger.error("modifyPARAMETERS: Failed to write PARAMETERS file.")
        raise
    return