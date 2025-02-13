"""Functions for reading and writing the VIBROCC file."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging
import re

import numpy as np

from viperleed.calc.lib.base import readToExc
from viperleed.calc.lib.base import splitSublists


logger = logging.getLogger(__name__)


def readVIBROCC(rp, slab, filename='VIBROCC', silent=False):
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
                if site.el not in rp.ELEMENT_MIX:
                    if site.getVibAmp(rp, site.el) is not None:
                        logger.error("Failed to generate vibrational "
                                     "amplitude for " + site.el + " in site "
                                     + site.label)
                        raise
                else:
                    for subel in rp.ELEMENT_MIX[site.el]:
                        if site.getVibAmp(rp, subel) is not None:
                            logger.error(
                                "Failed to generate vibrational amplitude for "
                                + site.el + " in site " + site.label)
                            raise
            try:
                checkVIBROCC(rp, slab, generate=True)
            except Exception:
                raise
            return True
        else:
            logger.error("VIBROCC not found.")
            raise
    # read VIBROCC:
    mode = 0  # 0: not reading; 1: vib.amp., 2: occ., 3: search offsets
    regex = False   # read regular expressions as-is or not
    for line in rf:
        line = line.lstrip()
        for c in ["!", "#", "%"]:    # start of comment
            line = line.split(c)[0].rstrip()
        if '=' not in line:
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
        if mode in [1, 2]:
            try:
                llist = line.split('=')[1].split()
            except IndexError:
                logger.warning(f'VIBROCC file: {param} appears to have '
                               'no value')
                rp.setHaltingLevel(1)
                continue
            if not llist:
                logger.warning(f'VIBROCC file: {param} appears to have '
                               'no value')
                rp.setHaltingLevel(1)
                continue
            # first get values on the right
            sublists = splitSublists(readToExc(llist), ',')
            # read parameter on the left
            targetsites = []
            for site in slab.sitelist:
                if regex:        # if regular expressions are enabled, take the
                    prep = param     # parameter input at face value
                else:
                    # double-slash non-literal characters
                    prep = re.escape(param)
                    # if regular expressions are not enabled, we want to still
                    #  interpret * as "any number of any characters":
                    prep = prep.replace('\\*', '.*')
                m = re.match(prep, site.label)
                if m:
                    # if the length of the matched text == the site label,
                    #  it's a real match
                    if m.end(0) == len(site.label):
                        targetsites.append(site)
            if len(targetsites) == 0:
                logger.warning('VIBROCC file: No sites matching '+param
                               + ' found, line will be skipped.')
                rp.setHaltingLevel(1)
                continue
            else:
                for site in targetsites:
                    if mode == 1:   # decide which dictionary to write to
                        td = site.vibamp
                    else:
                        td = site.occ
                    for sl in sublists:
                        if len(sl) == 1:
                            # if there's only one value, it should be a float
                            #  for the main site element
                            try:
                                if site.el in rp.ELEMENT_MIX:
                                    for subel in rp.ELEMENT_MIX[site.el]:
                                        td[subel] = float(sl[0])
                                else:
                                    td[site.el] = float(sl[0])
                            except Exception:
                                logger.error(
                                    'VIBROCC file: Error reading '
                                    'value '+sl[0]+' at parameter '+param)
                                raise
                        else:
                            el = sl[0].capitalize()
                            try:
                                td[el] = float(sl[1])
                            except KeyError:
                                logger.error(
                                    'VIBROCC file: Element ' + el + 'not '
                                    'recognized at parameter ' + param)
                            except Exception:
                                logger.error(
                                    'VIBROCC file: Error reading value '
                                    + sl[1] + ' at parameter ' + param)
                                raise
        # Search Offsets
        if mode == 3:
            try:
                ind = int(plist[1])
            except Exception:
                logger.error('VIBROCC file: Error converting to atom number: '
                             + plist[1])
                continue
            else:
                targetatlist = [at for at in slab if at.num == ind]
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
            sub_lists = s.split(",")
            for li in sub_lists:
                ll = li.split()
                if (len(ll) != 4 and om == 1) or (len(ll) != 2 and om != 1):
                    logger.error('VIBROCC file: Wrong number of values in '
                                 'sublist '+li)
                    continue
                else:
                    el = ll[0].capitalize()
                    if el not in slab.chemelem:
                        logger.error('VIBROCC file: Element '+el+' not '
                                     'recognized')
                        continue
                    values = []
                    try:
                        for s in ll[1:]:
                            values.append(float(s))
                    except ValueError:
                        logger.error('VIBROCC file: Could not convert values '
                                     'to floats: '+li)
                        continue
            if om == 1:
                value = np.array(values)
                value[2] *= -1  # invert z
            else:
                value = values[0]
            targetdict[el] = value
    rf.close()
    if not silent:
        logger.debug("VIBROCC file was read successfully")
    # now fill up default values & do consistency checks:
    vibAmpGenerated = checkVIBROCC(rp, slab, generate=generate, silent=silent)
    if vibAmpGenerated:
        return True
    else:
        if rp.T_EXPERIMENT is not None:
            logger.warning("Parameter T_EXPERIMENT is defined but unused.")
        if rp.T_DEBYE is not None:
            logger.warning("Parameter T_DEBYE is defined but unused.")
    return False


def checkVIBROCC(rp, slab, generate=False, silent=False):
    """Fills default values and does consistency check of site vibrational
    amplitudes and occupations. If vibrational amplitudes are automatically
    calculated here, returns True, else returns False."""
    vibAmpGenerated = False
    for site in slab.sitelist:
        # check if the vibrational amplitude is defined for the main element(s)
        if site.el not in site.vibamp:
            if site.el not in rp.ELEMENT_MIX:
                if generate:
                    if site.getVibAmp(rp, site.el) == 0:
                        vibAmpGenerated = True
                    else:
                        logger.warning(
                            'VIBROCC file: Failed to get ' + site.el +
                            ' vibrational amplitude for site ' + site.label)
                        rp.setHaltingLevel(2)
                else:
                    logger.warning(
                        'VIBROCC file: No ' + site.el + ' vibrational '
                        'amplitude defined for site '+site.label)
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
                    logger.error(
                        'VIBROCC file: No vibrational amplitude defined for '
                        'any of the main elements in site ' + site.label)
                    rp.setHaltingLevel(3)
                    # !!! raise exception
                else:
                    # vibrational amplitudes were defined for some of the
                    #  main elements but not all. Fill the other values
                    for subel in rp.ELEMENT_MIX[site.el]:
                        if subel not in site.vibamp:
                            logger.warning(
                                'VIBROCC file: No ' + subel + ' vibrational '
                                'amplitude defined for site ' + site.label +
                                '. Using vibrational amplitude from other '
                                'element in site.')
                            rp.setHaltingLevel(1)
                            site.vibamp[subel] = v
        if site.el in rp.ELEMENT_MIX:  # just use first element in ELEMENT_MIX
            mainva = site.vibamp[rp.ELEMENT_MIX[site.el][0]]
        else:
            mainva = site.vibamp[site.el]
        # for other elements, fill vibrational elements with that of the main
        # element (probably not present)
        for el in slab.chemelem:
            if el not in site.vibamp:
                site.vibamp[el] = mainva
        # check and fill occupations:
        o = 0.
        for el in slab.chemelem:
            if el in site.occ:
                o += site.occ[el]
        if 'Vac' in site.occ:
            o += site.occ['Vac']
        if site.el not in site.occ:
            if site.el in rp.ELEMENT_MIX:
                if not any([(el in site.occ)
                            for el in rp.ELEMENT_MIX[site.el]]):
                    for el in rp.ELEMENT_MIX[site.el]:
                        site.occ[el] = (1. - o) / len(rp.ELEMENT_MIX[site.el])
            else:
                site.occ[site.el] = 1. - o
            o = 1.
        for el in slab.chemelem:
            if el not in site.occ:
                site.occ[el] = 0.
        if o < 0.99:
            logger.debug(
                'VIBROCC file: Site '+site.label+' has a total '
                'occupation less than one. Interpreting remaining {:.2f} '
                'as vacancies.'.format(1-o))
        elif o > 1.0:
            if o > 1.01:
                logger.warning(
                    'VIBROCC file: Site '+site.label+' has a total occupation '
                    'greater than one ({:.2f}). Occupations will be re-scaled '
                    'to 1.'.format(o))
                rp.setHaltingLevel(2)
            for el in site.occ:
                site.occ[el] *= 1/o
        # finally, check whether we have any vibrational amplitudes or
        #  occupations assigned to non-existant elements, and if so, drop
        #  them and warn the user about it:
        dl = []
        for el in site.vibamp:
            if el not in slab.chemelem:
                logger.warning(
                    'VIBROCC file: Site ' + site.label + ' has a vibrational '
                    'amplitude defined for an unknown element, which will be '
                    'dropped (' + el + ').')
                rp.setHaltingLevel(1)
                dl.append(el)
        for el in dl:
            site.vibamp.pop(el, None)
        dl = []
        for el in site.occ:
            if el not in slab.chemelem and el != "Vac":
                logger.warning(
                    'VIBROCC file: Site ' + site.label + ' has an occupation '
                    'defined for an unknown element, which will be dropped ('
                    + el + ').')
                rp.setHaltingLevel(1)
                dl.append(el)
        for el in dl:
            site.occ.pop(el, None)
    if not silent:
        logger.debug("VIBROCC value consistency check finished")
    return vibAmpGenerated


def writeVIBROCC(sl, filename='VIBROCC', silent=False):
    """Writes a new VIBROCC file with the optimized parameters obtained after
    the search. The new file will not follow the ordering of the input VIBROCC
    file."""
    output = "= Vibrational Amplitudes\n"
    for site in sl.sitelist:
        ol = site.label + " = "
        targetels = [el for el in site.vibamp if (site.occ[el] != 0
                                                  or el in site.mixedEls)]
        if len(targetels) == 1:
            ol += "{:.4g}".format(site.vibamp[targetels[0]])
        else:
            for i, el in enumerate(targetels):
                ol += el + " {:.4g}".format(site.vibamp[el])
                if i < len(targetels)-1:
                    ol += ", "
        output += ol + "\n"

    output += "\n= Occupations\n"
    for site in sl.sitelist:
        write = True
        ol = site.label + " = "
        targetels = [el for el in site.occ if site.occ[el] != 0]
        if (len(targetels) == 1 and site.occ[targetels[0]] == 1
                and site.el == targetels[0] and not site.mixedEls):
            write = False
        elif len(targetels) == 1 and not site.mixedEls:
            ol += "{:.4g}".format(site.occ[targetels[0]])
        else:
            for i, el in enumerate(targetels):
                ol += el + " {:.4g}".format(site.occ[el])
                if i < len(targetels)-1:
                    ol += ", "
        if write:
            output += ol + "\n"

    output += "\n= Search offsets\n"
    # figure out which atoms to write, in which order
    offsetList = []    # metric per atom
    for at in [at for at in sl if not at.is_bulk]:
        to = 0
        # !!! TODO: Think about weights for the three
        for el in at.offset_occ:
            to += abs(at.offset_occ[el])
        for el in at.offset_vib:
            to += abs(at.offset_vib[el])
        for el in at.offset_geo:
            to += np.linalg.norm(at.offset_geo[el])
        if to >= 1e-4:  # !!! TODO: enough?
            offsetList.append((to, at))
    offsetList.sort(key=lambda tup: -tup[0])
    for (_, at) in offsetList:
        # POSITION OFFSET
        write = False
        ol = "POS {} = ".format(at.num)
        for el in at.offset_geo:
            if np.linalg.norm(at.offset_geo[el]) >= 1e-4:
                write = True
                ol += el + " {:g} {:g} {:g}, ".format(at.offset_geo[el][0],
                                                      at.offset_geo[el][1],
                                                      -at.offset_geo[el][2])
        if write:
            output += ol[:-2]+"\n"
        # VIBRATION OFFSET
        write = False
        ol = "VIB {} = ".format(at.num)
        for el in at.offset_vib:
            if abs(at.offset_vib[el]) >= 1e-3:
                write = True
                ol += el + " {:g}, ".format(at.offset_vib[el])
        if write:
            output += ol[:-2]+"\n"
        # OCCUPATION OFFSET
        write = False
        ol = "OCC {} = ".format(at.num)
        for el in at.offset_occ:
            if abs(at.offset_occ[el]) >= 1e-3:
                write = True
                ol += el + " {:g}, ".format(at.offset_occ[el])
        if write:
            output += ol[:-2]+"\n"
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write "+filename)
        raise
    if not silent:
        logger.info("Wrote to "+filename+" successfully")
    return
