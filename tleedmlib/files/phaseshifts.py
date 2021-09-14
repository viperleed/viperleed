# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 17:20:43 2020

@author: Florian Kraushofer

Functions for reading and writing the PHASESHIFTS file
"""

import logging
import numpy as np
import os

from viperleed import fortranformat as ff

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
except Exception:
    plotting = False
else:
    plotting = True

logger = logging.getLogger("tleedm.files.phaseshifts")


def readPHASESHIFTS(sl, rp, readfile='PHASESHIFTS', check=True,
                    ignoreEnRange=False):
    """Reads from a PHASESHIFTS file, returns the data as a list of tuples
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

    # legacy - allow "_PHASESHIFTS"
    if (readfile == 'PHASESHIFTS' and not os.path.isfile('PHASESHIFTS')
            and os.path.isfile('_PHASESHIFTS')):
        logger.info("Found no PHASESHIFTS file, but found legacy file named "
                    "_PHASESHIFTS. Renaming _PHASESHIFTS to PHASESHIFTS.")
        os.rename('_PHASESHIFTS', 'PHASESHIFTS')

    try:
        rf = open(readfile, 'r')
    except FileNotFoundError:
        logger.error("PHASESHIFTS file not found.")
        raise

    filelines = []
    for line in rf:
        filelines.append(line[:-1])
    rf.close()

    try:
        nel = ri3.read(filelines[0])[0]
    except Exception:
        logger.error("Exception while trying to read PHASESHIFTS: could not "
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
            for i in range(0, nel):
                elps = []
                for j in range(0, linesperblock):
                    llist = rf74x10.read(filelines[readline
                                                   + (i*linesperblock)+j+1])
                    llist = [f for f in llist if f is not None]
                    elps.extend(llist)
                enps.append(elps)
            phaseshifts.append((en, enps))
            readline += linesperblock*nel+1
        else:
            # first check how many lines until the next energy:
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
                        cont = False  # found next energy
                    else:
                        shortestline = 1
                elif len(llist) != longestline:
                    shortestline = len(llist)
                lastlen = len(llist)
            linesperblock = int((lineit-1)/nel)
            if not linesperblock or (((lineit-1)/nel) - linesperblock != 0.0):
                logger.warning(
                    "Error while trying to read PHASESHIFTS: "
                    "Could not parse file: The number of blocks may not match "
                    "the number given in the first line. A new PHASESHIFTS "
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
        llist = firstline.split()
        if len(llist) >= 6:
            for i in range(1, 5):
                try:
                    float(llist[i])
                except ValueError:
                    muftin = False
        else:
            muftin = False
        if rp.V0_REAL == "default" and not muftin:
            logger.warning(
                "Could not convert first line of PHASESHIFTS file to MUFTIN "
                "parameters. A new PHASESHIFTS file will be generated.")
            rp.setHaltingLevel(1)
        elif len(phaseshifts[0][1]) == psblocks:
            logger.debug("Found "+str(psblocks)+" blocks in PHASESHIFTS "
                         "file, which is consistent with PARAMETERS.")
            newpsGen, newpsWrite = False, False
        elif len(phaseshifts[0][1]) == len(sl.chemelem):
            logger.warning(
                "Found fewer blocks than expected in the "
                "PHASESHIFTS file. However, the number of blocks matches "
                "the number of chemical elements. A new PHASESHIFTS file "
                "will be generated, assuming that each block in the old "
                "file should be used for all atoms of one element.")
            rp.setHaltingLevel(1)
            oldps = phaseshifts[:]
            phaseshifts = []
            for (en, oldenps) in oldps:
                enps = []
                j = 0   # block index in old enps
                for el in sl.elements:
                    if el in rp.ELEMENT_MIX:
                        m = len(rp.ELEMENT_MIX[el])
                    else:
                        m = 1
                    n = len([s for s in sl.sitelist if s.el == el])
                    for i in range(0, m):    # repeat for chemical elements
                        for k in range(0, n):   # repeat for sites
                            enps.append(oldenps[j])
                        j += 1  # count up the block in old enps
                phaseshifts.append((en, enps))
            newpsGen = False
            firstline = str(len(phaseshifts[0][1])).rjust(3) + firstline[3:]
        else:
            logger.warning(
                "PHASESHIFTS file was read but is inconsistent with "
                "PARAMETERS. A new PHASESHIFTS file will be generated.")
            rp.setHaltingLevel(1)

    if check and not ignoreEnRange:
        # check whether energy range is large enough:
        checkfail = False
        er = np.arange(rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+1e-4,
                       rp.THEO_ENERGIES[2])
        psmin = round(phaseshifts[0][0]*27.211396, 2)
        psmax = round(phaseshifts[-1][0]*27.211396, 2)
        if rp.V0_REAL == "default":
            llist = firstline.split()
            c = []
            try:
                for i in range(0, 4):
                    c.append(float(llist[i+1]))
            except (ValueError, IndexError):
                checkfail = True
            else:
                er_inner = [e + (rp.FILAMENT_WF - max(c[0],
                                 c[1] + (c[2]/np.sqrt(e + c[3]
                                                      + rp.FILAMENT_WF))))
                            for e in er]  # energies at which scattering occurs
        else:
            try:
                v0r = float(rp.V0_REAL)
            except ValueError:
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
                    logger.warning(
                        "Lowest value in PHASESHIFTS file ({:.1f} "
                        "eV) is larger than the lowest predicted scattering "
                        "energy ({:.1f} eV). If this causes problems in the "
                        "reference calculation, try deleting the PHASESHIFTS "
                        "file to generate a new one, or increase the starting "
                        "energy in the THEO_ENERGIES parameter."
                        .format(psmin, min(er_inner)))
                else:
                    logger.warning(
                        "The energy range found in the PHASESHIFTS"
                        " file is smaller than the energy range requested for "
                        "theoretical beams. A new PHASESHIFTS file will be "
                        "generated.")
                    newpsGen, newpsWrite = True, True
        else:
            logger.warning(
                "Could not check energy range in PHASESHIFTS "
                "file. If energy range is insufficient, try deleting the "
                "PHASESHIFTS file to generate a new one.")
    return (firstline, phaseshifts, newpsGen, newpsWrite)


def writePHASESHIFTS(firstline, phaseshifts, filename='PHASESHIFTS'):
    """Takes phaseshift data and writes it to a PHASESHIFTS file."""
    output = firstline
    if output[-1] != "\n":
        output += "\n"
    f74x10 = ff.FortranRecordWriter('10F7.4')
    f74 = ff.FortranRecordWriter('F7.4')
    for (en, enps) in phaseshifts:
        output += f74.write([en])+"\n"
        for block in enps:
            output += f74x10.write(block)+"\n"
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
        logger.debug("Wrote to "+filename+" successfully.")
    except Exception:
        logger.error("Exception while writing PHASESHIFTS file: ",
                     exc_info=True)
        raise
    return


def plot_phaseshifts(sl, rp, filename="Phaseshifts_plots.pdf"):
    """
    Outputs plots of the phaseshifts in rp.phaseshifts in a pdf file.

    Parameters
    ----------
    sl : Slab
        Slab object, used to determine elements and sites represented in
        the phaseshifts.
    rp : Rparams
        Run parameters. Phaseshifts are read from rp.phaseshifts.
    filename : str, optional
        Name of the output file. The default is "Phaseshifts_plots.pdf".

    Returns
    -------
    None.

    """
    global plotting
    if not plotting:
        logger.debug("Necessary modules for plotting not found. Skipping "
                     "Phaseshift plotting.")
        return
    ps_labels = []
    for el in sl.elements:
        # this reproduces the order of blocks contained in PHASESHIFTS:
        if el in rp.ELEMENT_MIX:
            chemelList = rp.ELEMENT_MIX[el]
        elif el in rp.ELEMENT_RENAME:
            chemelList = [rp.ELEMENT_RENAME[el]]
        else:
            chemelList = [el]
        ps_labels.extend([cel + " in " + s.label + " site"
                          for cel in chemelList
                          for s in sl.sitelist if s.el == el])
    energies = np.array([ps[0]*27.211396 for ps in rp.phaseshifts])
    ps_vals = np.array([ps[1] for ps in rp.phaseshifts])

    figsize = (7, 4)
    linewidth = 1
    nlplot = min(np.shape(ps_vals)[-1], rp.LMAX[1]+1)

    figs = []
    # colors = ["#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
    #           "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
    #           "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d"]
    #   colorblind safe 16 - not great
    colors = ["#000000", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]  # colorblind safe 8
    styles = ["-", "--", ":", "-."]  # for when we run out of colors
    for i, label in enumerate(ps_labels):
        fig, ax = plt.subplots(figsize=figsize)
        figs.append(fig)
        fig.subplots_adjust(left=0.1, right=0.98,
                            bottom=0.11, top=0.92)
        ax.set_title(label)
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Phaseshifts (rad)")
        ax.grid(True, linewidth=0.1*linewidth)
        [sp.set_linewidth(0.7*linewidth) for sp in ax.spines.values()]
        ax.tick_params(bottom=True, top=True, left=True, right=True,
                       axis='both', direction='in', width=0.7*linewidth)
        ax.set_xlim((np.min(energies), np.max(energies)))
        for j in range(0, nlplot):
            ax.plot(energies, ps_vals[:, i, j],
                    linewidth=linewidth, label="L = {}".format(j),
                    c=colors[j % len(colors)],
                    ls=styles[(j // len(colors)) % len(styles)])
        legend = ax.legend(ncol=(nlplot // 8 + 1))
        legend.get_frame().set_linewidth(linewidth*0.7)

    try:
        pdf = PdfPages(filename)
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
        pdf.close()
    except PermissionError:
        logger.error("plot_phaseshifts: Cannot open file {}. Aborting."
                     .format(filename))
    return
