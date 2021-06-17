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

logger = logging.getLogger("tleedm.beamgen")


def runBeamGen(sl, rp, beamgensource=os.path.join('tensorleed',
                                                  'beamgen3.out'),
               domains=False):
    """Writes necessary input for the beamgen3 code, the runs it. The relevant
    output file will be renamed to BEAMLIST."""
    beamgensource = os.path.join(rp.sourcedir, beamgensource)
    output = ''
    f74x2 = ff.FortranRecordWriter('2F7.4')
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    ucbulk = np.transpose(sl.bulkslab.ucell[:2, :2])
    ol = f74x2.write(ucbulk[0])
    ol = ol.ljust(36)
    output += ol + 'ARA1\n'
    ol = f74x2.write(ucbulk[1])
    ol = ol.ljust(36)
    output += ol + 'ARA2\n'
    i3x2 = ff.FortranRecordWriter('2I3')
    ol = i3x2.write([int(round(f)) for f in rp.SUPERLATTICE[0]])
    ol = ol.ljust(36)
    output += ol + 'LATMAT - overlayer\n'
    ol = i3x2.write([int(round(f)) for f in rp.SUPERLATTICE[1]])
    ol = ol.ljust(36)
    output += ol + 'LATMAT -  matrix\n'
    output += ('  1                                 SSYM - symmetry code - cf.'
               ' van Hove / Tong 1979, always 1\n')
    if not domains:
        dmin = sl.getMinLayerSpacing() * 0.7
    else:
        dmin = min([dp.sl.getMinLayerSpacing() for dp in rp.domainParams])*0.7
    f71 = ff.FortranRecordWriter('F7.1')
    f41 = ff.FortranRecordWriter('F4.1')
    ol = f71.write([rp.THEO_ENERGIES[1]]) + f41.write([dmin])
    ol = ol.ljust(36)
    output += (ol + 'EMAX,DMIN - max. energy, min. interlayer distance for '
               'layer doubling\n')
    output += ('   {:.4f}                           TST - convergence '
               'criterion for fd. reference calculation\n'
               .format(rp.ATTENUATION_EPS))
    output += ('9999                                KNBMAX - max. number of '
               'beams to be written (may be a format problem!)')
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
