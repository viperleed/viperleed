# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 10:23:21 2020

@author: Florian Kraushofer

Functions for reading and writing PatternInfo.tlm
"""

import logging

logger = logging.getLogger("tleedm.files.patterninfo")


def writePatternInfo(sl, rp, filename="PatternInfo.tlm"):
    """Writes a PatternInfo file that can be used by the TLEEDMAP GUI utility
    to display the expected LEED pattern and show beam labelling."""
    output = "eMax = {:.2f}\n".format(rp.THEO_ENERGIES[1])
    mstring = "[[{}, {}], [{}, {}]]".format(sl.ucell[0, 0], sl.ucell[1, 0],
                                            sl.ucell[0, 1], sl.ucell[1, 1])
    output += "surfBasis = "+mstring+"\n"
    mstring = ("[[{:.0f}, {:.0f}], [{:.0f}, {:.0f}]]"
               .format(rp.SUPERLATTICE[0, 0], rp.SUPERLATTICE[0, 1],
                       rp.SUPERLATTICE[1, 0], rp.SUPERLATTICE[1, 1]))
    output += "superlattice = "+mstring+"\n"
    if sl.planegroup in ["pm", "pg", "cm", "rcm", "pmg"]:
        pgstring = sl.planegroup+str(sl.orisymplane.par)
    else:
        pgstring = sl.planegroup
    output += "surfGroup = "+pgstring+"\n"
    if sl.bulkslab is None:
        logger.error("PatternInfo.tlm: bulk slab has not been initialized.")
        raise RuntimeError("writePatternInfo called without bulk slab.")
    output += "bulkGroup = "+sl.bulkslab.foundplanegroup+"\n"
    output += "bulk3Dsym = "+sl.bulkslab.getBulk3Dstr()
    # write output
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return
