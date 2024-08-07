"""Functions for writing files relevant to the superpos calculation."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging

import fortranformat as ff
import numpy as np

logger = logging.getLogger(__name__)

from viperleed.calc.lib.version import Version


def writeSuperposInput(sl, rp, config, param_name="PARAM",
                       contrin_name="superpos-CONTRIN",
                       for_error=False, only_vary=None):
    """
    Writes a PARAM and a CONTRIN file as input for the Superpos calculation.
    Returns the contents of CONTRIN, as this will also be needed as an input
    stream for the superpos calculation.

    Parameters
    ----------
    sl : Slab
        The Slab object to be modified.
    rp : Rparams
        The run parameters.
    config : list of int
        Parameter indices to be applied.
    param_name : str, optional
        Alternative name for the PARAM file.
    contrin_name : str, optional
        Alternative name for the superpos-CONTRIN file.
    for_error : bool, optional
        If True, the 'config' parameter will be ignored, and superpos will
        instead be calculated for all variations.
    only_vary : list of SearchPar, optional
        If running for error, pass a list of search parameters that should be
        varied. All other search parameters will be kept at their "centered"
        position.

    Returns
    -------
    str
        Contents of CONTRIN.

    """
    if only_vary is None:
        only_vary = rp.searchpars
    mnfiles = 0
    occupations = []    # occupation per delta file
    deltanames = []     # names of delta files
    surfaceChecks = []  # is that atom at the surface or not?
    indices = []        # vib and geo indices from search, cast to 1D array
    n_var = 1
    n_conc = 1
    if for_error:
        n_var = max([sp.steps for sp in rp.searchpars
                     if sp.mode in ("geo", "vib") and sp in only_vary])
        n_conc = max([sp.steps for sp in rp.searchpars
                      if sp.mode == "occ" and sp in only_vary])
    # determine which atoms are "at the surface"
    surfats = sl.getSurfaceAtoms(rp)
    # make lists to print
    which_atoms = [at for at in rp.search_atlist
                   if [sp for sp in rp.searchpars
                       if sp.atom == at and (sp.non_zero or sp in only_vary)]]
    for at in which_atoms:
        sps = [sp for sp in rp.searchpars if sp.atom == at]
        occpar = [sp for sp in sps if sp.mode == "occ"][0]  # exactly one
        totalocc = np.array([0.] * n_conc)
        for el in at.disp_occ:
            if at in surfats:
                surfaceChecks.append(1)
            else:
                surfaceChecks.append(0)
            if not for_error:
                occind = rp.searchpars.index(occpar)
                o = [at.disp_occ[el][config[occind]-1]]
            elif occpar not in only_vary:
                o = [at.disp_occ[el][occpar.center-1]] * n_conc
            else:
                try:
                    o = [at.disp_occ[el][j] for j in range(n_conc)]
                except IndexError:
                    logger.error("Superpos for error: Inconsistent number of "
                                 "variation steps. Found maximum {} steps,"
                                 " but only {} for {}, element {}."
                                 .format(n_conc, len(at.disp_occ[el]), at, el))
                    rp.setHaltingLevel(2)
                    return ""
            occupations.append(o)
            totalocc += np.array(o)
            pl = [sp for sp in sps if sp.el == el]
            if not pl:
                logger.error('No search parameters found '
                             f'for {at}. Aborting...')
                rp.setHaltingLevel(2)
                return ""
            deltanames.append(pl[0].deltaname)
            mnfiles += 1
            if not for_error or any(sp not in only_vary for sp in pl):
                if el in at.disp_geo:
                    k = el
                else:
                    k = "all"
                ngeo = len(at.disp_geo[k])
                vibpar = [sp for sp in pl if sp.mode == "vib"]
                if len(vibpar) > 0:
                    if not for_error:
                        vibind = config[rp.searchpars.index(vibpar[0])] - 1
                    else:
                        vibind = vibpar[0].center - 1
                else:
                    vibind = 0
                geopar = [sp for sp in pl if sp.mode == "geo"]
                if len(geopar) > 0:
                    if not for_error:
                        geoind = config[rp.searchpars.index(geopar[0])] - 1
                    else:
                        geoind = geopar[0].center - 1
                else:
                    geoind = 0
                # cast vib and geo indices from 2D to 1D:
                indices.append([(ngeo * vibind) + geoind + 1] * n_var)
            else:
                whichpar = max(pl, key=lambda x: x.steps)
                if whichpar.steps == n_var:
                    indices.append(list(range(1, n_var + 1)))
                else:
                    logger.error("Superpos for error: Inconsistent number of "
                                 "variation steps. Found maximum {} steps,"
                                 " but only {} for {}."
                                 .format(n_var, whichpar.steps, at))
                    rp.setHaltingLevel(2)
                    return ""
        vp = [sp for sp in sps if sp.el == "vac"]
        if len(vp) != 0 and np.any(totalocc < 1.0):
            occupations.append(list(1 - totalocc))
            deltanames.append(vp[0].deltaname)
            mnfiles += 1
            indices.append([1] * n_var)
            if at in surfats:
                surfaceChecks.append(1)
            else:
                surfaceChecks.append(0)

    # write PARAM
    param = """C  DIMENSIONS MAY BE CHANGED IN PARAMETER-STATEMENT
C
C  MNFILES: (maximum) number of different files containing amplitude changes
C           to be used as input files for the I(E) calculation
C  MNCONCS: (maximum) number of different combinations of occupation probabilities
C           for the lattice sites in the surface - i.e., weight factors multiplied
C           to the delta amplitudes from each file
C  MNT0:    (maximum) number of beams for which I(E) spectra are calculated
C  MNCSTEP: (maximum) number of parameter combinations (geometry x vibrations)
C           per energy tabulated in an individual delta amplitude file
C  MNATOMS: currently inactive - set to 1
C
"""
    param += "      PARAMETER (MNFILES={})\n".format(mnfiles)
    param += "      PARAMETER (MNCONCS={})\n".format(n_conc)
    param += ("      PARAMETER (MNT0={}, MNCSTEP={}, MNATOMS=1)\n"
              .format(len(rp.ivbeams), rp.mncstep))
    try:
        with open(param_name, "w") as wf:
            wf.write(param)
    except Exception:
        logger.error("Failed to write to PARAM file for superpos")
        raise

    # collect contrin output
    if rp.TL_VERSION < Version('1.7.0'):
        formatter = {'int': ff.FortranRecordWriter('I3'),
                     'occ': ff.FortranRecordWriter('10F7.4'),
                     'var': ff.FortranRecordWriter('I3'),
                     'ndelta': ff.FortranRecordWriter('100I3'),
                     }
    else:
        formatter = {'int': ff.FortranRecordWriter('I4'),
                     'occ': ff.FortranRecordWriter('10F7.4'),
                     'var': ff.FortranRecordWriter('I3'),
                     'ndelta': ff.FortranRecordWriter('1000I5'),
                     }
    contrin = ((formatter['int'].write([mnfiles])
                + formatter['int'].write([0])).ljust(16)
               + "no. of files, VarAll: all possible parameter combinations? "
               "(1/0)\n")
    contrin += (formatter['int'].write([n_conc]).ljust(16)
                + "number of concentration steps\n")
    occ_trans = map(list, zip(*occupations))  # transpose occupations
    for occ in occ_trans:
        contrin += formatter['occ'].write(occ) + "\n"
    for i in range(0, len(deltanames)):
        contrin += deltanames[i].ljust(15) + formatter['var'].write([n_var])
        contrin += (formatter['var'].write([surfaceChecks[i]])
                    + formatter['var'].write([1]) + "  FILENAME,VARIATIONS,"
                    "SURFACE?,FORMATTED? (A15,3I3)\n")
        contrin += formatter['ndelta'].write(indices[i]) + "\n"
    try:
        with open(contrin_name, "w") as wf:
            wf.write(contrin)
    except Exception:
        logger.error("Failed to write " + contrin_name + " file. Execution "
                     "will continue...")

    return contrin
