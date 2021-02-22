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
knownParams = [
    'ATTENUATION_EPS', 'BEAM_INCIDENCE', 'BULKDOUBLING_EPS',
    'BULKDOUBLING_MAX', 'BULK_REPEAT', 'DOMAIN', 'DOMAIN_STEP', 'ELEMENT_MIX',
    'ELEMENT_RENAME', 'FILAMENT_WF', 'FORTRAN_COMP', 'HALTING',
    'IV_SHIFT_RANGE', 'LAYER_CUTS', 'LAYER_STACK_VERTICAL', 'LMAX',
    'LOG_DEBUG', 'LOG_SEARCH', 'N_BULK_LAYERS', 'N_CORES', 'PARABOLA_FIT',
    'PHASESHIFT_EPS', 'PLOT_COLORS_RFACTOR', 'RUN', 'R_FACTOR_SMOOTH',
    'R_FACTOR_TYPE', 'SCREEN_APERTURE', 'SEARCH_BEAMS', 'SEARCH_CONVERGENCE',
    'SEARCH_CULL', 'SEARCH_MAX_GEN', 'SEARCH_POPULATION', 'SEARCH_START',
    'SITE_DEF', 'SUPERLATTICE', 'SUPPRESS_EXECUTION', 'SYMMETRIZE_INPUT',
    'SYMMETRY_CELL_TRANSFORM', 'SYMMETRY_EPS', 'SYMMETRY_FIND_ORI',
    'SYMMETRY_FIX', 'TENSOR_INDEX', 'TENSOR_OUTPUT', 'THEO_ENERGIES',
    'TL_VERSION', 'T_DEBYE', 'T_EXPERIMENT', 'V0_IMAG', 'V0_REAL',
    'V0_Z_ONSET', 'VIBR_AMP_SCALE']
paramAlias = {}  # keys should be all lowercase, with no underscores
for p in knownParams:
    paramAlias[p.lower().replace("_", "")] = p


def updatePARAMETERS(rp, filename='PARAMETERS'):
    """
    Reads PARAMETERS file again, but ignores everything not concerning the
    search or STOP. Updates the given Rparams object accordingly.

    Parameters
    ----------
    rp : Rparams
        Parameters for current run, as defined previously. Will be updated if
        parameters have changed.
    filename : string, optional
        The file to be read. The default is 'PARAMETERS'.

    Returns
    -------
    None.

    """
    try:
        with open(filename, 'r') as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.warning("updatePARAMETERS routine: PARAMETERS file not found.")
        return
    for line in lines:
        if "!" in line:
            line = line.split("!")[0].rstrip()
        for param in ["SEARCH_KILL", "STOP"]:  # SEARCH_KILL is legacy name
            if line.lstrip().upper().startswith(param):
                if not re.match(r"\s*"+param+r"\s*=\s*[Ff](alse)?", line):
                    rp.STOP = True
        if "=" not in line:
            continue  # ignore all lines that don't have an "=" sign at all
        param = line.split('=')[0]        # parameter is defined left of "="
        if param:
            # get rid of spaces and check the leftmost entry.
            plist = param.split()
            if plist:
                param = plist[0]
        if (param not in knownParams and
                param.lower().replace("_", "") in paramAlias):
            param = paramAlias[param.lower().replace("_", "")]
        if param not in knownParams:
            continue
        try:
            value = line.split('=', maxsplit=1)[1].rstrip()
            llist = value.split()  # read the stuff to the right of "="
        except IndexError:
            llist = []
        if not llist:
            continue
        if param == 'SEARCH_CONVERGENCE':
            flags = plist[1:]
            if not flags or flags[0].lower() not in ['dgen', 'gaussian']:
                continue
            fl = [None, None]
            for (i, s) in enumerate(llist[:2]):
                try:
                    fl[i] = float(s)
                except Exception:
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
    """
    Reads a PARAMETERS file and returns an Rparams object with the
    raw input, without interpretation.

    Parameters
    ----------
    filename : str, optional
        The file to be read. The default is 'PARAMETERS'.

    Returns
    -------
    rpars : Rparams
        Object storing parameters for current run. Will contain parameters
        read in this function, and some 'global' parameters defined at runtime.
    """
    try:
        rf = open(filename, 'r')
    except FileNotFoundError:
        logger.error("PARAMETERS not found.")
        raise
    # read PARAMETERS:
    rpars = tl.Rparams()
    for line in rf:
        if "!" in line:
            line = line.split("!")[0].rstrip()
        line = line.lstrip()
        for param in ["STOP", "SEARCH_KILL"]:
            if line.lstrip().upper().startswith(param) and not re.match(
                    r"\s*" + param + r"\s*=\s*[Ff](alse)?", line):
                logger.warning(
                    'PARAMETERS file: {0} was set at start of '
                    'program. Modifying PARAMETERS to disable {0}; re-insert '
                    'it if you actually want to stop.'.format(param))
                modifyPARAMETERS(rpars, param, comment='Disabled at program '
                                 'start', path=os.path.dirname(filename),
                                 suppress_ori=True)
        if "=" not in line:
            continue     # ignore all lines that don't have an "=" sign at all
        param = line.split('=')[0]      # parameter is defined left of "="
        if not param:
            continue
        # get rid of spaces and check the leftmost entry.
        plist = param.split()
        param = plist[0]
        if (param not in knownParams and
                param.lower().replace("_", "") in paramAlias):
            param = paramAlias(param.lower().replace("_", ""))
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
        if param not in rpars.readParams:
            rpars.readParams[param] = []
        rpars.readParams[param].append((plist, value))
    rf.close()
    return rpars


def interpretPARAMETERS(rpars, slab=None, silent=False):
    """
    Interprets the string values in an Rparams object, read previously by
    readPARAMETERS, to fill the parameter variables.

    Parameters
    ----------
    rpars : Rparams
        Object storing parameters for current run. Created previously by
        readPARAMETERS, and should already contain raw string data.
    slab : Slab, optional
        Slab object with elements and atomic position data. If not passed, some
        parameters will not be interpreted.
    silent : bool, optional
        If True, less output will be printed. The default is False.

    Raises
    ------
    ValueError
        Raised if the value for some parameter cannot be used.

    Returns
    -------
    None

    """

    def setBoolParameter(rp, param, value, varname=None,
                         addAllowedValues={False: [], True: []},
                         haltingOnFail=1):
        """
        Generic function for setting a given parameter to a boolean value.

        Parameters
        ----------
        rp : Rparams
            The Rparams object for which the parameter should be set
        param : str
            The name of the parameter in PARAMETERS
        value : str
            The value as found in the PARAMETERS file
        varname : str, optional
            The variable name in the Rparams class, if it differs from 'param'
        addAllowedValues : dict, optional
            Additional string values which should be interpreted as False or
            True. By default, 'f' and 'false' are False, 't' and 'true' are
            True.
        haltingOnFail : int, optional
            What halting level should be set if the 'fail' event is triggered.

        Returns
        ----------
        int
            0 if value was set, 1 otherwise

        """
        s = value.lower()
        allowedValues = {False: ['false', 'f'], True: ['true', 't']}
        for key in allowedValues:
            allowedValues[key].extend(addAllowedValues[key])
            if s in allowedValues[key]:
                v = key
                break
        else:
            logger.warning('PARAMETERS file: {}: Could not interpret given '
                           'value. Input will be ignored.'.format(param))
            rpars.setHaltingLevel(haltingOnFail)
            return 1
        if not varname:
            varname = param
        setattr(rp, varname, v)
        return 0

    def setNumericalParameter(rp, param, value, varname=None, type_=float,
                              range_=(None, None), haltingOnFail=1,
                              range_exclude=(False, False),
                              outOfRangeEvent=('fail', 'fail')):
        """
        Generic function for setting a given parameter to an integer or float
        value.

        Parameters
        ----------
        rp : Rparams
            The Rparams object for which the parameter should be set
        param : str
            The name of the parameter in PARAMETERS
        value : str
            The value as found in the PARAMETERS file
        varname : str, optional
            The variable name in the Rparams class, if it differs from 'param'
        type_ : int or float, optional
            Which type the variable must take
        range_ : tuple, values can be int or float or None, optional
            Defines upper and lower limits for the value. If either is set to
            None, only the other one will be interpreted. The default is None
            for both.
        range_exclude : tuple (bool, bool), optional
            Whether the upper and lower limits themselves are allowed values
            or not. The default is False for both (range is inclusive).
        outOfRangeEvent : str 'fail' or 'set', optional
            What to do if the given value lies outside the range. For 'fail',
            the value will be ignored. For 'modulo', the value will be
            brought within the range by '%' operation. For 'set', the value
            will be set to the lowest allowed value. The default is
            ('fail', 'fail').
        haltingOnFail : int, optional
            What halting level should be set if the 'fail' event is triggered.

        Returns
        ----------
        int
            0 if value was set, 1 otherwise

        """

        try:
            v = type_(value)
        except ValueError:
            logger.warning('PARAMETERS file: {}: Could not convert value to '
                           '{}. Input will be ignored.'
                           .format(param, type_.__name__))
            rpars.setHaltingLevel(haltingOnFail)
            return 1
        outOfRange = [False, False]
        if range_[0] is not None and (v < range_[0] or
                                      (range_exclude[0] and v == range_[0])):
            outOfRange[0] = True
        if range_[1] is not None and (v > range_[1] or
                                      (range_exclude[1] and v == range_[1])):
            outOfRange[1] = True
        rangeChars = ({False: "[", True: "]"}, {False: "]", True: "["})
        if all([v is not None for v in range_]):
            outOfRangeStr = 'not in range {}{}, {}{}'.format(
                                        rangeChars[0][range_exclude[0]],
                                        range_[0], range_[1],
                                        rangeChars[1][range_exclude[1]])
        elif range_[0] is not None:
            if range_exclude[0]:
                outOfRangeStr = 'less than or equal to {}'.format(range_[0])
            else:
                'less than {}'.format(range_[0])
        else:
            if range_exclude[1]:
                outOfRangeStr = 'larger than or equal to {}'.format(range_[1])
            else:
                'larger than {}'.format(range_[1])
        if any([outOfRange[i] and outOfRangeEvent[i] == 'fail'
                for i in range(0, 2)]):
            logger.warning('PARAMETERS file: {}: Value {} is {}. Input will '
                           'be ignored.'.format(param, v, outOfRangeStr))
            rpars.setHaltingLevel(haltingOnFail)
            return 1
        else:
            for i in range(0, 2):
                if outOfRange[i]:
                    if outOfRangeEvent == 'modulo':
                        if not range_[1]:
                            raise ValueError(
                                'Cannot use outOfRangeEvent modulo if upper '
                                'limit is undefined.')
                        if not range_[0]:
                            range_[0] = 0
                        setTo = (((v - range_[0]) % (range_[1] - range_[0]))
                                 + range_[0])
                    elif range_exclude[i]:
                        mult = 1
                        if i == 1:
                            mult = -1
                        if type_ == float:
                            setTo = range_[i] + mult*1e-4
                        else:
                            setTo = range_[i] + mult
                    else:
                        setTo = range_[i]
                    logger.warning('PARAMETERS file: {}: Value {} is {}. '
                                   'Value will be set to {}.'
                                   .format(param, v, outOfRangeStr, setTo))
                    v = setTo
                    break
        if not varname:
            varname = param
        setattr(rp, varname, v)
        return 0

    loglevel = logger.level
    if silent:
        logger.setLevel(logging.ERROR)
    # track some parameters while reading
    searchConvRead = False
    # define order that parameters should be read in
    orderedParams = ["LOG_DEBUG", "RUN"]
    checkParams = [p for p in orderedParams if p in rpars.readParams]
    checkParams.extend([p for p in knownParams if (p in rpars.readParams and
                                                   p not in checkParams)])
    domainsIgnoreParams = [
        'BULK_REPEAT', 'ELEMENT_MIX', 'ELEMENT_RENAME',
        'LAYER_CUTS', 'N_BULK_LAYERS', 'SITE_DEF', 'SUPERLATTICE',
        'SYMMETRY_CELL_TRANSFORM', 'TENSOR_INDEX', 'TENSOR_OUTPUT']
    for (param, plist, value) in [(param, lel[0], lel[1])
                                  for param in checkParams
                                  for lel in rpars.readParams[param]]:
        if (4 in rpars.RUN or rpars.domainParams) and (param in
                                                       domainsIgnoreParams):
            continue
        llist = value.split()
        # parameters not interpreted right now
        if param == 'VIBR_AMP_SCALE':
            rpars.VIBR_AMP_SCALE.extend(value.split(","))
        # simple bool parameters
        elif param in ['LOG_DEBUG', 'LOG_SEARCH', 'SUPPRESS_EXECUTION',
                       'SYMMETRIZE_INPUT', 'SYMMETRY_FIND_ORI']:
            setBoolParameter(rpars, param, llist[0])
        # slightly more complicated bools
        elif param == 'LAYER_STACK_VERTICAL':
            setBoolParameter(rpars, param, llist[0],
                             addAllowedValues={False: ['c'], True: ['z']})
        # positive-only integers
        elif param in ['BULKDOUBLING_MAX', 'N_CORES', 'SEARCH_MAX_GEN',
                       'TENSOR_INDEX']:
            setNumericalParameter(rpars, param, llist[0], type_=int,
                                  range_=(1, None))
        # positive-only floats
        elif param in ['T_DEBYE', 'T_EXPERIMENT', 'V0_IMAG', 'TL_VERSION']:
            setNumericalParameter(rpars, param, llist[0], range_=(0, None))
        # simple numericals
        elif param == 'V0_Z_ONSET':
            setNumericalParameter(rpars, param, llist[0])
        elif param == 'ATTENUATION_EPS':
            setNumericalParameter(rpars, param, llist[0], range_=(0.0001, 1),
                                  range_exclude=(False, True))
        elif param == 'BULKDOUBLING_EPS':
            setNumericalParameter(rpars, param, llist[0],
                                  range_=(0.0001, None),
                                  outOfRangeEvent=('set', 'fail'))

        elif param == 'HALTING':
            setNumericalParameter(rpars, param, llist[0], type_=int,
                                  range_=(1, 3))
        elif param == 'LMAX':
            r = setNumericalParameter(rpars, param, llist[0], type_=int,
                                      range_=(1, 15),
                                      outOfRangeEvent=('set', 'set'))
            if r == 0 and rpars.PHASESHIFT_EPS != 0:
                logger.warning(
                    'PARAMETERS file: Both LMAX and '
                    'PHASESHIFT_EPS are being defined. PHASESHIFT_EPS '
                    'will be ignored.')
        elif param == 'N_BULK_LAYERS':
            setNumericalParameter(rpars, param, llist[0], type_=int,
                                  range_=(1, 2), haltingOnFail=2)
        elif param == 'R_FACTOR_SMOOTH':
            setNumericalParameter(rpars, param, llist[0], type_=int,
                                  range_=(0, 999))
        elif param == 'R_FACTOR_TYPE':
            setNumericalParameter(rpars, param, llist[0], type_=int,
                                  range_=(1, 2))
        elif param == 'SCREEN_APERTURE':
            setNumericalParameter(rpars, param, llist[0], range_=(0, 180))
        elif param == 'SEARCH_POPULATION':
            r = setNumericalParameter(rpars, param, llist[0], type_=int,
                                      range_=(1, None))
            if r == 0 and rpars.SEARCH_POPULATION < 16:
                logger.warning('SEARCH_POPULATION is very small. A '
                               'minimum value of 16 is recommended.')
        elif param == 'SYMMETRY_EPS':
            r = setNumericalParameter(rpars, param, llist[0],
                                      range_=(1e-10, None))
            if r == 0 and rpars.SYMMETRY_EPS > 1.0:
                logger.warning(
                    'PARAMETERS file: SYMMETRY_EPS: Given '
                    'value is greater than one Ångström. This is a '
                    'very loose constraint and might lead to '
                    'incorrect symmetry detection. Be sure to check '
                    'the output!')
                rpars.setHaltingLevel(1)
            if len(llist) > 1:
                r = setNumericalParameter(rpars, param, llist[0],
                                          range_=(1e-10, None),
                                          varname='SYMMETRY_EPS_Z')
                if r == 0 and rpars.SYMMETRY_EPS_Z > 1.0:
                    logger.warning(
                        'PARAMETERS file: SYMMETRY_EPS: Given '
                        'value for z is greater than one Ångström. This is a '
                        'very loose constraint and might lead to '
                        'incorrect symmetry detection. Be sure to check '
                        'the output!')
                    rpars.setHaltingLevel(1)
                if r != 0:
                    rpars.SYMMETRY_EPS_Z = rpars.SYMMETRY_EPS
            else:
                rpars.SYMMETRY_EPS_Z = rpars.SYMMETRY_EPS
        # non-trivial parameters
        elif param == 'BEAM_INCIDENCE':
            range_ = {'THETA': (-90, 90), 'PHI': (0, 360)}
            outOfRangeEvent = {'THETA': ('fail', 'fail'),
                               'PHI': ('modulo', 'modulo')}
            if ',' in value:
                sublists = tl.base.splitSublists(llist, ',')
                for sl in sublists:
                    for name in ['THETA', 'PHI']:
                        if sl[0].upper() == name:
                            setNumericalParameter(
                                rpars, 'BEAM_INCIDENCE {}'.format(name),
                                sl[1], range_=range_[name],
                                outOfRangeEvent=outOfRangeEvent[name],
                                varname=name)
                            break
                    else:
                        logger.warning(
                            'PARAMETERS file: BEAM_INCIDENCE: Unknown flag '
                            'found. Input will be ignored.')
                        rpars.setHaltingLevel(1)
            else:
                for ind, name in enumerate(['THETA', 'PHI']):
                    setNumericalParameter(
                        rpars, 'BEAM_INCIDENCE {}'.format(name),
                        llist[ind], range_=range_[name], varname=name,
                        outOfRangeEvent=outOfRangeEvent[name])
        elif param == 'BULK_REPEAT':
            s = value.lower()
            if "[" not in s:
                if "(" not in s:
                    try:
                        rpars.BULK_REPEAT = abs(float(llist[0]))
                    except ValueError:
                        logger.warning(
                            'PARAMETERS file: BULK_REPEAT: Could not convert '
                            'value to float. Input will be ignored.')
                        rpars.setHaltingLevel(1)
                else:
                    m = re.match(r'\s*(c|z)\(\s*(?P<val>[0-9.]+)\s*\)', s)
                    if not m:
                        logger.warning(
                            'PARAMETERS file: BULK_REPEAT: Could not parse '
                            'input expression. Input will be ignored.')
                    else:
                        try:
                            v = abs(float(m.group("val")))
                        except Exception:
                            logger.warning(
                                'PARAMETERS file: BULK_REPEAT: Could not '
                                'convert value to float. Input will be '
                                'ignored.')
                            rpars.setHaltingLevel(1)
                        else:
                            if "z" in s:
                                rpars.BULK_REPEAT = v
                            else:  # c
                                rpars.BULK_REPEAT = slab.ucell[2, 2] * v
            else:  # vector
                vec = tl.base.readVector(s, slab.ucell)
                if vec is None:
                    logger.warning(
                        'PARAMETERS file: BULK_REPEAT: Could not parse input '
                        'expression. Input will be ignored.')
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
                logger.warning(
                    'PARAMETERS file: Multiple sources defined '
                    'for DOMAIN {}. Last entry will be ignored.'.format(name))
                rpars.setHaltingLevel(2)
                continue
            if not name:  # get unique name
                i = 1
                while str(i) in names:
                    i += 1
                name = str(i)
            # check path
            value = value.strip()
            if os.path.exists(value):
                path = value
            elif os.path.isfile(value + ".zip"):
                path = value + ".zip"
            else:
                logger.warning(
                    'PARAMETERS file: Value for DOMAIN {} could not be '
                    'interpreted as either a path or a .zip file.'
                    .format(name))
                rpars.setHaltingLevel(2)
                continue
            rpars.DOMAINS.append((name, path))
        elif param == 'DOMAIN_STEP':
            try:
                i = int(llist[0])
            except ValueError:
                logger.warning(
                    'PARAMETERS file: DOMAIN_STEP: Could not '
                    'convert value to integer. Input will be ignored.')
                rpars.setHaltingLevel(1)
            if not (1 <= i <= 100):
                logger.warning('PARAMETERS file: DOMAIN_STEP: Invalid '
                               'value given. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            if 100 % i != 0:
                j = i-1
                while 100 % j != 0:
                    j -= 1
                logger.warning(
                    'PARAMETERS file: DOMAIN_STEP: 100 is not '
                    'divisible by given value {}. DOMAIN_STEP will be set to '
                    '{} instead.'.format(i, j))
                rpars.setHaltingLevel(1)
                rpars.DOMAIN_STEP = j
            else:
                rpars.DOMAIN_STEP = i
        elif param == 'ELEMENT_MIX':
            ptl = [el.lower() for el in tl.leedbase.periodic_table]
            found = False
            for el in llist:
                if el.lower() not in ptl:
                    logger.warning(
                        'PARAMETERS file: ELEMENT_MIX for {0}: {1} not found '
                        'in periodic table. ELEMENT_MIX for {0} will be '
                        'ignored.'.format(plist[1], el))
                    rpars.setHaltingLevel(1)
                    found = True
            if not found:
                rpars.ELEMENT_MIX[plist[1]] = [el.capitalize()
                                               for el in llist]
        elif param == 'ELEMENT_RENAME':
            ptl = [el.lower() for el in tl.leedbase.periodic_table]
            if llist[0].lower() not in ptl:
                logger.warning(
                    'PARAMETERS file: ELEMENT_RENAME for {0}: {1} not found '
                    'in periodic table. ELEMENT_RENAME for {0} will be '
                    'ignored.'.format(plist[1], llist[0]))
                rpars.setHaltingLevel(1)
            else:
                rpars.ELEMENT_RENAME[plist[1]] = llist[0].capitalize()
        elif param == 'FILAMENT_WF':
            if llist[0].lower() == 'w':
                rpars.FILAMENT_WF = 4.5
            elif llist[0].lower() == 'lab6':
                rpars.FILAMENT_WF = 2.65
            else:
                setNumericalParameter(rpars, param, llist[0])
        elif param == 'FORTRAN_COMP':
            if len(plist) <= 1 and llist[0].lower() in ["ifort", "gfortran"]:
                rpars.getFortranComp(comp=llist[0].lower())
            elif (len(plist) > 1 and plist[1].lower() == "mpi"
                  and llist[0].lower() in ["mpifort", "mpiifort"]):
                rpars.getFortranMpiComp(comp=llist[0].lower())
            else:
                delim = llist[0][0]     # should be quotation marks
                if delim not in ["'", '"']:
                    logger.warning(
                        'PARAMETERS file: FORTRAN_COMP '
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
                    logger.warning(
                        'PARAMETERS file: FORTRAN_COMP parameter: '
                        'Could not interpret flags. Value will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'IV_SHIFT_RANGE':
            if len(llist) in (2, 3):
                fl = []
                try:
                    fl = [float(s) for s in llist]
                except ValueError:
                    logger.warning(
                        'PARAMETERS file: Failed to convert IV_SHIFT_RANGE '
                        'input to floats Input will be ignored')
                    rpars.setHaltingLevel(1)
                else:
                    if fl[1] >= fl[0]:
                        for i in range(0, 2):
                            rpars.IV_SHIFT_RANGE[i] = fl[i]
                    else:
                        logger.warning(
                            'PARAMETERS file: IV_SHIFT_RANGE '
                            'end energy has to be greater than or equal '
                            'to start energy. Input will be ignored.')
                        rpars.setHaltingLevel(1)
                    if len(fl) == 3:
                        if fl[2] > 0:
                            rpars.IV_SHIFT_RANGE[2] = fl[2]
                        else:
                            logger.warning(
                                'PARAMETERS file: IV_SHIFT_RANGE step has to '
                                'be positive. Input will be ignored.')
                            rpars.setHaltingLevel(1)
            else:
                logger.warning('PARAMETERS file: Unexpected number of '
                               'values for IV_SHIFT_RANGE. Input will be '
                               'ignored.')
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
            for (i, s) in enumerate(llist):
                if "dz" in s.lower() or "dc" in s.lower():
                    m = rgx.match(value.lower())
                    if m:
                        try:
                            float(m.group('cutoff'))
                            llist[i] = m.group(0)
                        except Exception:
                            logger.warning(
                                'PARAMETERS file: LAYER_CUTS parameter: Could '
                                'not parse function ' + s + '. Input will be '
                                'ignored.')
                            rpars.setHaltingLevel(1)
                            continue
                elif not (s == "<" or s == ">"):
                    try:
                        float(s)
                    except Exception:
                        logger.warning('PARAMETERS file: LAYER_CUTS '
                                       'parameter: Error parsing values. '
                                       'Input will be ignored.')
                        rpars.setHaltingLevel(1)
                        continue
            rpars.LAYER_CUTS = llist
        elif param == 'PARABOLA_FIT':
            sublists = tl.base.splitSublists(llist, ',')
            for sl in sublists:
                flag = sl[0].lower()
                if flag.lower() == 'localise':
                    flag = 'localize'
                value_error = ('PARAMETERS file: PARABOLA_FIT: Value {} is '
                               'not valid for flag {}. Value will be ignored.'
                               .format(sl[1], sl[0]))
                if flag == 'type':
                    if sl[1].lower() in ('linear', 'linearregression', 'lasso',
                                         'ridge', 'elasticnet'):
                        rpars.PARABOLA_FIT['type'] = sl[1]
                    else:
                        logger.warning(value_error)
                        rpars.setHaltingLevel(1)
                elif flag in ('alpha', 'mincurv', 'localize'):
                    try:
                        f = float(sl[1])
                    except ValueError:
                        f = -1
                    if f >= 0:
                        rpars.PARABOLA_FIT[flag] = f
                    else:
                        logger.warning(value_error)
                        rpars.setHaltingLevel(1)
        elif param == 'PHASESHIFT_EPS':
            if rpars.LMAX != 0:
                logger.warning('PARAMETERS file: Both LMAX and '
                               'PHASESHIFT_EPS are being defined. '
                               'PHASESHIFT_EPS will be ignored.')
            else:
                try:
                    f = float(llist[0])
                except ValueError:
                    s = llist[0].lower()[0]
                    if s == 'r':    # rough
                        f = 0.1
                    elif s == 'n':  # normal
                        f = 0.05
                    elif s == 'f':  # fine
                        f = 0.01
                    elif s == 'e':  # extrafine
                        f = 0.001
                    else:
                        logger.warning('PARAMETERS file: PHASESHIFT_EPS: '
                                       'Could not convert value to float. '
                                       'Input will be ignored.')
                        rpars.setHaltingLevel(1)
                if f > 0 and f < 1:
                    rpars.PHASESHIFT_EPS = f
                else:
                    logger.warning(
                        'PARAMETERS file: PHASESHIFT_EPS: Unexpected value '
                        '(should be between 0 and 1). Input will be ignored.')
                    rpars.setHaltingLevel(1)
        elif param == 'PLOT_COLORS_RFACTOR':
            if len(llist) >= 2:
                if len(llist) > 2:
                    logger.warning(
                        'PARAMETERS file: PLOT_COLORS_RFACTOR '
                        'parameter: Expected two values, found {}. First two '
                        'values will be used.'.format(len(llist)))
                    rpars.setHaltingLevel(1)
                rpars.PLOT_COLORS_RFACTOR = (llist[0], llist[1])
            else:
                logger.warning(
                    'PARAMETERS file: PLOT_COLORS_RFACTOR '
                    'parameter: Expected two values, found {}. Input will be '
                    'ignored.'.format(len(llist)))
                rpars.setHaltingLevel(1)
        elif param == 'RUN':
            rl = []
            for s in llist:
                ir = tl.base.readIntRange(s)
                if len(ir) > 0:
                    rl.extend(ir)
                else:
                    logger.warning(
                        'PARAMETERS file: RUN: Could not interpret value '
                        + s + ', skipping value...')
                    rpars.setHaltingLevel(2)
            if len(rl) > 0:
                if 4 in rl:
                    logger.info('Found domain search.')
                i = 0
                while i < len(rl):
                    if rl[i] not in (0, 1, 2, 3, 4, 11, 12, 31):
                        logger.warning(
                            'PARAMETERS file: RUN: Value {} does not '
                            'correspond to a segment and will be skipped.'
                            .format(rl[i]))
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
        elif param == 'SEARCH_BEAMS':
            if llist[0][0].lower() in ["0", "a"]:
                rpars.SEARCH_BEAMS = 0
            elif llist[0][0].lower() in ["1", "i"]:
                rpars.SEARCH_BEAMS = 1
            elif llist[0][0].lower() in ["2", "f"]:
                rpars.SEARCH_BEAMS = 2
            else:
                logger.warning('PARAMETERS file: SEARCH_BEAMS: value not '
                               'recognized. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'SEARCH_CONVERGENCE':
            if len(plist) == 1:
                if value.lower().strip() == 'off':
                    rpars.GAUSSIAN_WIDTH_SCALING = 1.
                    continue
                else:
                    logger.warning(
                        'PARAMETERS file: SEARCH_CONVERGENCE: '
                        'no flag given, value ' + value + ' not recognized.')
                    rpars.setHaltingLevel(1)
                    continue
            flags = plist[1:]
            if flags[0].lower() not in ['dgen', 'gaussian']:
                logger.warning(
                    'PARAMETERS file: SEARCH_CONVERGENCE: flag "'
                    + flags[0] + '" not recognized.')
                rpars.setHaltingLevel(1)
                continue
            fl = [None, None]
            for (i, s) in enumerate(llist[:2]):
                try:
                    fl[i] = float(s)
                except ValueError:
                    logger.warning(
                        'PARAMETERS file: SEARCH_CONVERGENCE gaussian: could '
                        'not convert value to float.')
                    rpars.setHaltingLevel(1)
            if flags[0].lower() == 'gaussian':
                if fl[0] is not None and fl[0] > 0:
                    rpars.GAUSSIAN_WIDTH = fl[0]
                else:
                    logger.warning(
                        'PARAMETERS file: SEARCH_CONVERGENCE '
                        'gaussian should be a positive number.')
                    rpars.setHaltingLevel(1)
                if fl[1] is not None:
                    if 0 < fl[1] <= 1:
                        rpars.GAUSSIAN_WIDTH_SCALING = fl[1]
                    else:
                        logger.warning(
                            'PARAMETERS file: SEARCH_CONVERGENCE gaussian: '
                            'scaling value should be in range ]0, 1[')
                        rpars.setHaltingLevel(1)
            elif flags[0].lower() == 'dgen':
                if len(flags) == 1:
                    target = 'dec'
                elif flags[1].lower() in ['dec', 'best', 'all']:
                    target = flags[1].lower()
                else:
                    logger.warning(
                        'PARAMETERS file: SEARCH CONVERGENCE '
                        'dgen: flag "' + flags[1] + '" not recognized.')
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
                        logger.warning(
                            'PARAMETERS file: SEARCH_CONVERGENCE dgen '+target
                            + ': scaling value cannot be smaller than 1.')
                        rpars.setHaltingLevel(1)
        elif param == 'SEARCH_CULL':
            try:
                f = float(llist[0])
            except ValueError:
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
        elif param == 'SEARCH_START':
            s = llist[0].lower()
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
                    ir = tl.base.readIntRange(sl[i])
                    if len(ir) > 0:
                        atnums.extend(ir)
                    elif "top(" in sl[i]:
                        if slab is None:
                            logger.warning(
                                'PARAMETERS file: SITE_DEF '
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
            rpars.SITE_DEF[plist[1]] = newdict
        elif param in ['SUPERLATTICE', 'SYMMETRY_CELL_TRANSFORM']:
            if 'M' not in plist:
                if slab is None:
                    logger.warning(
                        'PARAMETERS file: {} parameter appears to be in Wood '
                        'notation, but no slab was passed; cannot calculate '
                        'bulk unit cell!'.format(param))
                    rpars.setHaltingLevel(2)
                else:
                    setattr(rpars, param, tl.leedbase.readWoodsNotation(
                        value, slab.ucell))
                    if param == 'SUPERLATTICE':
                        rpars.superlattice_defined = True
            else:
                sublists = tl.base.splitSublists(llist, ',')
                if not len(sublists) == 2:
                    logger.warning('PARAMETERS file: error reading {} '
                                   'matrix: number of lines is not equal 2.'
                                   .format(param))
                    rpars.setHaltingLevel(2)
                else:
                    write = True
                    nl = []
                    for sl in sublists:
                        if len(sl) == 2:
                            try:
                                nl.append([float(s) for s in sl])
                            except ValueError:
                                logger.warning(
                                    'PARAMETERS file: error '
                                    'reading {} matrix: could not convert {} '
                                    'to floats.'.format(param, sl))
                                rpars.setHaltingLevel(2)
                                write = False
                        else:
                            logger.warning(
                                'PARAMETERS file: error reading {} matrix: '
                                'number of columns is not equal 2.'
                                .format(param))
                            rpars.setHaltingLevel(2)
                            write = False
                    if write:
                        setattr(rpars, param, np.array(nl, dtype=float))
                        if param == 'SUPERLATTICE':
                            rpars.superlattice_defined = True
        elif param == 'SYMMETRY_FIX':
            s = llist[0].lower()
            grouplist = [
                "p1", "p2", "pm", "pg", "cm", "rcm", "pmm", "pmg", "pgg",
                "cmm", "rcmm", "p4", "p4m", "p4g", "p3", "p3m1", "p31m", "p6",
                "p6m"]
            if s == 'true':
                pass    # same as default, determine symmetry automatically
            elif s == 'false':
                rpars.SYMMETRY_FIX = 'p1'
            elif s in grouplist:
                if s not in ("cm", "pmg"):
                    rpars.SYMMETRY_FIX = s
                else:
                    logger.warning('PARAMETERS file: SYMMETRY_FIX: For '
                                   'group '+s+', direction needs to be '
                                   'specified. Input will be ignored.')
                    rpars.setHaltingLevel(1)
            elif s[0:2] in ["pm", "pg", "cm"] or s[0:3] == "pmg":
                # regex to read
                rgx = re.compile(
                    r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*'
                    + r'\[\s*(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
                m = rgx.match(value.lower())
                if not m:
                    logger.warning(
                        'PARAMETERS file: SYMMETRY_FIX: Could '
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
                    if (group in ["pm", "pg", "cm", "rcm", "pmg"]
                            and i1 in range(-1, 3) and i2 in range(-1, 3)):
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
        elif param == 'TENSOR_OUTPUT':
            nl = tl.base.recombineListElements(llist, '*')
            for s in nl:
                try:
                    v = int(s)
                    if v != 0 and v != 1:
                        logger.warning(
                            'PARAMETERS file: Problem with '
                            'TENSOR_OUTPUT input format: Found value '+str(v)
                            + ', expected 0 or 1. Value will be ignored.')
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
                            logger.warning(
                                'PARAMETERS file: Problem with '
                                'TENSOR_OUTPUT input format: could not read '
                                'value '+s+', value will be ignored.')
                            rpars.setHaltingLevel(1)
                        if v != 0 and v != 1:
                            logger.warning(
                                'PARAMETERS file: Problem with '
                                'TENSOR_OUTPUT input format: Found value '
                                + str(v) + ', expected 0 or 1. Value will be '
                                'ignored.')
                            rpars.setHaltingLevel(1)
                        else:
                            for i in range(0, r):
                                rpars.TENSOR_OUTPUT.append(v)
                    else:
                        logger.warning(
                            'PARAMETERS file: Problem with '
                            'TENSOR_OUTPUT input format: could not read '
                            'value '+s+', value will be ignored.')
                        rpars.setHaltingLevel(1)
        elif param == 'THEO_ENERGIES':
            if len(llist) == 1:
                # single value input - only one energy requested
                try:
                    f = float(llist[0])
                except ValueError:
                    logger.warning(
                        'PARAMETERS file: Failed to convert '
                        'THEO_ENERGIES input to floats Input will be ignored')
                    rpars.setHaltingLevel(1)
                else:
                    if f > 0:
                        rpars.THEO_ENERGIES = [f, f, 1]
                    else:
                        logger.warning(
                            'PARAMETERS file: Unexpected input for '
                            'THEO_ENERGIES. Input will be ignored.')
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
                                logger.warning(
                                    'PARAMETERS file: '
                                    'THEO_ENERGIES values have to be '
                                    'positive. Input will be ignored.')
                                rpars.setHaltingLevel(1)
                        except ValueError:
                            logger.warning(
                                'PARAMETERS file: Failed to convert '
                                'THEO_ENERGIES input to floats. '
                                'Input will be ignored.')
                            rpars.setHaltingLevel(1)
                if len(fl) == 3:
                    if defined < 3:
                        rpars.THEO_ENERGIES = fl
                    elif (fl[0] > 0 and fl[1] > fl[0] and fl[2] > 0):
                        if (fl[1] - fl[0]) % fl[2] != 0:
                            # if the max is not hit by the steps exactly,
                            #   correct max up to make it so
                            fl[0] -= fl[2] - (fl[1] - fl[0]) % fl[2]
                            if fl[0] <= 0:
                                fl[0] = fl[0] % fl[2]
                                if fl[0] == 0:
                                    fl[0] += fl[2]
                            logger.info(
                                'THEO_ENERGIES parameter: '
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
        elif param == 'V0_REAL':
            if llist[0].lower() == 'rundgren':
                try:
                    c = []
                    for i in range(0, 4):
                        c.append(float(llist[i+1]))
                    setTo = (
                        "workfn-max("+str(round(c[0], 3))
                        + ", (("+str(round(c[1], 3)) + ")+("
                        + str(round(c[2], 3)) + ")/sqrt(EEV+workfn+("
                        + str(round(c[3], 3)) + "))))")
                except ValueError:
                    logger.warning(
                        "PARAMETERS file: V0_REAL parameter: "
                        "could not parse constants for Rundgren-type "
                        "function. Input will be ignored.")
                    rpars.setHaltingLevel(1)
            else:
                setTo = re.sub("(?i)EE", "EEV+workfn", value)
            setTo = setTo.rstrip()
            rpars.V0_REAL = setTo
    logger.setLevel(loglevel)
    return


def modifyPARAMETERS(rp, modpar, new="", comment="", path="",
                     suppress_ori=False, include_left=False):
    """
    Looks for 'modpar' in the PARAMETERS file, comments that line out, and
    replaces it by the string specified by 'new'

    Parameters
    ----------
    rp : Rparams
        The run parameters object.
    modpar : str
        The parameter that should be modified.
    new : str, optional
        The new value for the parameter. If not passed (default), the parameter
        will be commented out without replacement.
    comment : str, optional
        A comment to be added in the new line in PARAMETERS.
    path : str, optional
        Where to find the PARAMETERS file that should be modified.
    suppress_ori : bool, optional
        If True, no 'PARAMETERS_original' file will be created.
    include_left : bool, optional
        If False (default), 'new' will be interpreted as only the string that
        should go on the right-hand side of the equal sign. If True, the entire
        line will be replace by 'new'.

    Returns
    -------
    None.

    """

    file = os.path.join(path, "PARAMETERS")
    oriname = "PARAMETERS_ori_"+rp.timestamp
    ori = os.path.join(path, oriname)
    if oriname not in rp.manifest and not suppress_ori:
        try:
            shutil.copy2(file, ori)
        except Exception:
            logger.error(
                "modifyPARAMETERS: Could not copy PARAMETERS file "
                "to PARAMETERS_ori. Proceeding, original file will be lost.")
        rp.manifest.append(oriname)
    if "PARAMETERS" not in rp.manifest and not path:
        rp.manifest.append("PARAMETERS")
    output = ""
    headerPrinted = False

    try:
        with open(file, "r") as rf:
            plines = rf.readlines()
    except Exception:
        logger.error("Error reading PARAMETERS file.")
        raise
    found = False
    for line in plines:
        if "! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #" in line:
            headerPrinted = True
        valid = False
        param = ""
        if "=" in line:  # ignore all lines that don't have an "=" sign at all
            param = line.split('=')[0]      # parameter is defined left of "="
            if param:
                valid = True
                plist = param.split()
                if plist:
                    param = plist[0]
                if param[0] == '!':
                    valid = False
        else:
            for p in ["STOP", "SEARCH_KILL"]:
                if line.lstrip().upper().startswith(p):
                    valid = True
                    param = p
        if valid and param == modpar:
            found = True
            if new:
                if (modpar+" = "+new.strip()) == line.split("!")[0].strip():
                    output += line
                else:
                    if comment == "":
                        comment = "line automatically changed to:"
                    output += "!"+line[:-1] + " ! " + comment + "\n"
                    if not include_left:
                        output += modpar + " = " + new + "\n"
                    else:
                        output += new + "\n"
            else:
                if comment == "":
                    comment = "line commented out automatically"
                output += (("!" + line.rstrip()).ljust(35)
                           + " ! " + comment + "\n")
        else:
            output += line
    if new and not found:
        if not headerPrinted:
            output += """

! ######################################################
! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #
! ######################################################
"""
            headerPrinted = True
        output += "\n" + modpar + " = " + new
        if comment:
            output += " ! " + comment
    try:
        with open(file, "w") as wf:
            wf.write(output)
    except Exception:
        logger.error("modifyPARAMETERS: Failed to write PARAMETERS file.")
        raise
    return
