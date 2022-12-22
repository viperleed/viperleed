# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:24:16 2020

@author: Florian Kraushofer, Alexander Imre
"""

import os
import logging
import subprocess
import numpy as np
from pathlib import Path

from viperleed import fortranformat as ff

from viperleed.leedbase import (HARTREE_TO_EV,
                                BOHR_TO_ANGSTROM,
                                ANGSTROM_TO_BOHR)
from viperleed.guilib.base import get_equivalent_beams, BeamIndex

logger = logging.getLogger("tleedm.beamgen")


def runBeamGen(sl, rp, beamgensource='', domains=False):
    """Writes necessary input for the beamgen3 code, then runs it. The relevant
    output file will be renamed to BEAMLIST."""
    if rp.TL_VERSION < 1.7:
        beamgensource = os.path.join('tensorleed', 'beamgen3.out')
    else:
        beamgensource = os.path.join('tensorleed', 'beamgen.v1.7')
    beamgensource = os.path.join(rp.sourcedir, beamgensource)
    output = ''
    if rp.TL_VERSION < 1.7:
        formatter = {'uc': ff.FortranRecordWriter('2F7.4'),
                     'latmat': ff.FortranRecordWriter('2I3'),
                     'emax': ff.FortranRecordWriter('F7.1'),
                     'dmin': ff.FortranRecordWriter('F4.1'),
                     'tst': ff.FortranRecordWriter('F11.6'),
                     }
    else:
        formatter = {'uc': ff.FortranRecordWriter('2F9.4'),
                     'latmat': ff.FortranRecordWriter('2I4'),
                     'emax': ff.FortranRecordWriter('F6.1'),
                     'dmin': ff.FortranRecordWriter('F7.4'),
                     'tst': ff.FortranRecordWriter('F8.6'),
                     }
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    ucbulk = np.transpose(sl.bulkslab.ucell[:2, :2])
    output += formatter['uc'].write(ucbulk[0]).ljust(36) + 'ARA1\n'
    output += formatter['uc'].write(ucbulk[1]).ljust(36) + 'ARA2\n'
    ol = formatter['latmat'].write([int(round(f)) for f in rp.SUPERLATTICE[0]])
    output += ol.ljust(36) + 'LATMAT - overlayer\n'
    ol = formatter['latmat'].write([int(round(f)) for f in rp.SUPERLATTICE[1]])
    output += ol.ljust(36) + 'LATMAT -  matrix\n'
    if rp.TL_VERSION < 1.7:
        output += ('  1                                 SSYM - symmetry code'
                   ' - cf. van Hove / Tong 1979, always 1\n')
    if not domains:
        dmin = sl.getMinLayerSpacing() * 0.7
    else:
        dmin = min([dp.sl.getMinLayerSpacing() for dp in rp.domainParams])*0.7
    ol = (formatter['emax'].write([rp.THEO_ENERGIES[1]])
          + formatter['dmin'].write([dmin])).ljust(36)
    output += (ol + 'EMAX,DMIN - max. energy, min. interlayer distance for '
               'layer doubling\n')
    output += (formatter['tst'].write([rp.ATTENUATION_EPS]).ljust(36)
               + 'TST - convergence criterion for fd. reference calculation\n')
    output += ('99999                               KNBMAX - max. number of '
               'beams to be written')
    try:
        with open('DATA', 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write DATA for BEAMLIST generation.")
        raise
    # if os.name == 'nt':
    #     logger.error("Beamlist generation is currently not supported on "
    #                  "Windows. Use a linux shell to run beamlist generation "
    #                  "script.")
    #     raise EnvironmentError("Beamlist generation is currently not "
    #                            "supported on Windows.")
    # else:
    try:
        subprocess.call(beamgensource)
    except Exception:
        logger.error("Failed to execute beamgen script.")
        raise
    # clean up folder, rename files
    try:
        os.remove('BELIST')
        os.remove('PROT')
    except Exception:
        logger.warning("BEAMLIST generation: Failed to remove BELIST or "
                       "PROT file")
    try:
        os.rename('DATA', 'beamgen3-input')
    except Exception:
        logger.warning("Failed to rename beamlist generation input file "
                       "DATA to beamgen3-input")
    try:
        os.rename('NBLIST', 'BEAMLIST')
    except Exception:
        logger.error("Failed to rename beamlist generation output file "
                     "NBLIST to BEAMLIST")
        raise
    logger.debug("Wrote to BEAMLIST successfully.")
    return

def beamgen(sl, rp, domains=False, beamlist_name="BEAMLIST"):
    """Calculates and writes the contents for the file BEAMLIST.

    BEAMLIST contains a list of all diffraction beams that will be used 
    for internal calculations (as opposed to IVBEAMS, which contains the
    beams to be output). The file list the beams with indices as float
    values and a lower cutoff energy below which the beam is evanescent
    (i.e. does not leave the surface and has intensity 0). The used 
    format is defined in make_beamlist_string and dictated by TensErLEED
    and the legacy beamgen scripts that were used before. Note that
    BEAMLIST is not read directly by refcalc, but instead all input
    files for the refcalc are combined into one string (by collectFIN in
    iorefcalc.py) and then piped in.

    NB: The energies calculated here are slightly higher than the values
    from beamgenv3 (and beamgen.v1.7) because we use *more accurate*
    values for unit conversions.The legacy code used these rounded
    values:
    HARTREE_TO_EV = 27.21
    BOHR_TO_ANGSTROM = 0.529
    Similarly, the list of included beams may be different for the same 
    energy range, as the legacy code used rounded values and typecast a
    cutoff from float to int.

    In any case, this version should give more accurate energy values 
    and be more generous in how many beams are considered.

    Parameters
    ----------
    sl : Slab
        Slab object.
    rp : Rparams.
        Run parameters.
    domains : bool, optional
        Flag to indicate if performing a domain calculation,
        by default False.
    beamlist_name : str, optional
        Filename to be written, by default "BEAMLIST".
    """
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)

    e_max = rp.THEO_ENERGIES[1]
    surf_ucell = sl.surface_vectors
    inv_surf_vectors = sl.reciprocal_vectors

    if not domains:
        d_min = sl.getMinLayerSpacing() * 0.7
    else:
        d_min = min([dp.sl.getMinLayerSpacing() for dp in rp.domainParams])*0.7

    conv_crit = rp.ATTENUATION_EPS  # convergence criterion for refcalc
    # effective cutoff energy to use (scale to correct units)
    e_max_eff = (2*e_max +
                 (np.log(conv_crit)/(d_min*ANGSTROM_TO_BOHR))**2
                 *HARTREE_TO_EV)

    # use guilib to generate list of beams
    leedParameters = {
        'eMax': e_max_eff,
        'surfBasis': surf_ucell,
        'SUPERLATTICE': rp.SUPERLATTICE,
        'surfGroup': sl.foundplanegroup,
        'bulkGroup': sl.bulkslab.foundplanegroup,
    }
    equivalent_beams = get_equivalent_beams(leedParameters)
    # convert to float array
    indices = np.array(list(BeamIndex(beam[0]) for beam in equivalent_beams),
                       dtype="float64")
    # calculate cutoff energy for each beam
    energies = (np.sum(np.dot(indices, inv_surf_vectors)**2, axis=1)
                /2 *HARTREE_TO_EV *BOHR_TO_ANGSTROM**2)  # scale to correct units
    logger.debug(f"Highest energy considered in BEAMLIST={np.max(energies)}")

    # write to file
    beamlist_contents = make_beamlist_string(indices, energies)
    write_file_path = Path(beamlist_name)
    try:
        with open(write_file_path, 'w') as file:
            file.write(beamlist_contents)
    except Exception:
        logger.error(f"Unable to write file {beamlist_name}")
        raise

    logger.debug("Wrote to BEAMLIST successfully.")
    return


def make_beamlist_string(indices, energies):
    """Creates contents for file BEAMLIST in the format used be the 
    legacy beamgen scripts by U. Loeffler and R. Doell.

    Parameters
    ----------
    indices : np.ndarray, shape=(n_beams, 2)
        Indices (diffraction orders) of the beams.
    energies : np.ndarray, shape=(n_beams,)
        Lower cutoff energies for the beams.

    Returns
    -------
    str
        String representation of the contents of the BEAMLIST file.

    Raises
    ------
    ValueError
        If indices and energies have incompatible shapes.
    """
    n_beams = indices.shape[0]
    if not energies.shape == (n_beams,) or not indices.shape == (n_beams, 2):
        raise ValueError(
            f"Incompatible size of indices (shape={indices.shape})"
            f"and energies (shape={energies.shape}).")

    # set up Fortran format as was used by beamgen
    beamlist_format = ff.FortranRecordWriter(
        "2F10.5,2I3,10X,'E =  ',F10.4,2X,'NR.',I5"  # beamgen v1.7 had I5, beamgen v3 had I4 for some reason
        )
    # first line contains number of beams
    content = ff.FortranRecordWriter('10I3').write([n_beams]) + '\n'

    # iterate over all beams and format lines
    for nr, ((beam_h, beam_k), energy) in enumerate(zip(indices, energies)):
        # nr+1 because of Fortran indices starting at 1
        line = beamlist_format.write([beam_h, beam_k, 1, 1, energy,nr+1])
        content += line + '\n'

    return content
