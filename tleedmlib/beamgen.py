# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:24:16 2020

@author: Florian Kraushofer
"""

import os
import logging
import subprocess
import numpy as np
from viperleed import fortranformat as ff

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

def beamgen(sl, rp, domain=False):
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    
    e_max = rp.THEO_ENERGIES[1]
    surf_ucell = sl.ucell[:2, :2].T

    leedParameters = {
        'eMax': e_max,
        'surfBasis': surf_ucell,
        'SUPERLATTICE': rp.SUPERLATTICE,
        'surfGroup': sl.foundplanegroup,
        'bulkGroup': sl.bulkslab.foundplanegroup,
        'screenAperture': rp.SCREEN_APERTURE,
    }
    # TODO: domains!!
    equivalent_beams = get_equivalent_beams(leedParameters)
    indices = np.array(list(BeamIndex(beam[0]) for beam in equivalent_beams), dtype="float64")
    energies = np.sum(np.dot(indices, surf_ucell)**2, axis=1)
    # forbidden beams have negative group nr
    g_max_squared = 2*e_max + (np.log(tst)/d_min)**2
    
def evanescent_energy(beam_index, ucell):
    return np.sum(np.dot(beam_index, ucell)**2)

def make_beamlist_string(indices, energies):
    n_beams = indices.shape[0]
    if not energies.shape[0] == n_beams:
        raise ValueError(
            f"Incompatible size of indices (shape={indices.shape})"
            f"and energies (shape={energies.shape}).")

    # set up Fortran format as was used by beamgen
    beamlist_format = ff.FortranRecordWriter(
        "2F10.5,2I3,10X,'E =  ',F10.4,' EV',2X,'NR.',I4"
        )
    # first line contains number of beams
    content = ff.FortranRecordWriter('10I3').write([n_beams]) + '\n'
    
    # iterate over all beams and format lines
    for nr, ((beam_h, beam_k), energy) in enumerate(zip(indices, energies)):
        # nr+1 because of Fortran indices starting at 1
        line = beamlist_format.write([beam_h, beam_k, 1, 1, energy,nr+1])
        content += line + '\n'

    return content
