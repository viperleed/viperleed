"""Functions for reading and writing PatternInfo.tlm."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging


logger = logging.getLogger(__name__)


def writePatternInfo(sl, rp, filename="PatternInfo.tlm"):
    """Writes a PatternInfo file that can be used by the TLEEDMAP GUI utility
    to display the expected LEED pattern and show beam labelling."""
    output = f"eMax = {rp.THEO_ENERGIES.max:.2f}\n"
    mstring = "[[{}, {}], [{}, {}]]".format(*sl.ab_cell.T.ravel())
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
    output += "bulk3Dsym = "+sl.bulkslab.get_bulk_3d_str()
    # write output
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return
