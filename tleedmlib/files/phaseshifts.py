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

from viperleed.tleedmlib.leedbase import (get_atomic_number,
                                          get_element_symbol)

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
_HARTREE_TO_EV = 27.211396

def readPHASESHIFTS(sl, rp, readfile='PHASESHIFTS', check=True,
                    ignoreEnRange=False):
    """Reads from a PHASESHIFTS file.
    
    Parameters
    ----------
    sl: Slab
    rp: RunParameters
    readfile: str, optional
        filename to be read. Default is 'PHASESHIFTS'
    check: bool, optional
        Wether to check for consistence agains sl and rp. Default is True. If False,
        sl and rp can be None.
    ignoreEnRange: bool, optional
        Check wether the energy range in readfile is sufficient to cover the
        energy range requested in rp. Default is False.
        
    Returns
    -------
    firstline: str
        Header line of the readfile, containing Rundgren parameters.
    phaseshifts: list of tuple
        Each element is (energy, pahseshifts at energy) where 
        phaseshifts at energy is [[el0_L0, el0_L1, ...], [el1_L0, ...], ...]
        where eli_Lj is the phaseshift for element i and angular momentum j.
        Therefore, len(phaseshifts) is the number of energies found, 
        len(phaseshifts[0][1]) should match the number of elements, and
        len(phaseshifts[0][1][0]) is the number of different
        values of L in the phaseshift file.
    newpsGen: bool
        Wether the inconsitency found requires a full recalculation of the phaseshifts.
        This happens if: 
            1) an error occurs while parsing,
            2) if check, rp.V0_real was not given and could not interpret firstline,
            3) LMAX of readfile is smaller than required in rp,
            4) if number of blocks inconsitent with elements in sl.
    newpsWrite: bool
        Wether the inconsitency found requires writing a new PHASESHIFTS file.
        This is distinct from newpsGen as we can at times generate a new file 
        from the information read.
        This is the case if:
            5) There are fewer blocks than expected, but the number of blocks
                matches the number of chemical elements in sl. This means however,
                that all sites with the same chemical element will use the same
                phaseshift.
"""
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
    else:
        filelines = rf.readlines()
    finally:
        rf.close()

    try:
        nel = ri3.read(filelines[0])[0]
    except Exception:
        logger.error("Exception while trying to read PHASESHIFTS: could not "
                     "find number of blocks in first line.")
        raise
    phaseshifts = []

    firstline = filelines[0]
    line_idx = 1
    try:
        lines_per_block = find_block_size(filelines[1:])
    except ValueError as err:
        logger.warning(
            f"Error while trying to read PHASESHIFTS: {err} "
            "A new PHASESHIFTS file will be generated."
        )
        rp.setHaltingLevel(1)
        return "", [], True, True
    
    # check block length is consitent with number of species
    if (lines_per_block-1)%nel:
        logger.warning(
            "Error while trying to read PHASESHIFTS: "
            "Could not parse file: The number of blocks may not match "
            "the number given in the first line. A new PHASESHIFTS "
            "file will be generated."
        )
        rp.setHaltingLevel(1)
        return "", [], True, True

    lines_per_element = (lines_per_block-1)//nel

    while line_idx < len(filelines):
        current_line = filelines[line_idx].strip()
        if not current_line:
            line_idx += 1
            continue
        energy = rf74x10.read(filelines[line_idx])[0] # read energy
        enps = []
        for i in range(nel):
            elps = []
            for j in range(lines_per_element):
                llist = rf74x10.read(filelines[line_idx
                                                + (i*lines_per_element)+j+1])
                llist = [f for f in llist if f is not None]
                elps.extend(llist)
            enps.append(elps)
        phaseshifts.append((energy, enps))
        line_idx += lines_per_element*nel+1

    if not check:
        return firstline, phaseshifts, False, False

    # try extracting muffin tin parameters:
    muftin = []
    llist = firstline.split()
    if len(llist) >= 6:
        for i in range(1, 5):
            try:
                muftin.append(float(llist[i]))
            except (ValueError, IndexError):
                muftin = []
                break
    else:
        muftin = []

    (phaseshifts,
     firstline,
     newpsGen,
     newpsWrite) = __check_consistency_rp_elements(sl, rp, phaseshifts, firstline, muftin)

    if not ignoreEnRange:
       newpsGen, newpsWrite = __check_consistency_energy_range(rp, phaseshifts, muftin, newpsGen, newpsWrite)

    # check consitency of phaseshift values
    __check_consitency_element_order(sl, phaseshifts, eps=rp.PHASESHIFT_EPS)

    return firstline, phaseshifts, newpsGen, newpsWrite


def __check_consistency_rp_elements(sl, rp, phaseshifts, firstline, muftin):
    """Check whether the phaseshifts that were found fit the number 
    of elements and LMAX expected.

    Parameters
    ----------
    sl : Slab
        Surface Slab object.
    rp : Rparams
        Run paramters.
    phaseshifts : list
        Nested list of phaseshifts.
    firstline : str
        First line of the phaseshifts file.
    muftin : ??
        TODO

    Returns
    -------
    ??
        TODO
    """
    newpsGen, newpsWrite = True, True # should new values should be generated / written to file?

    n_l_values = len(phaseshifts[0][1][0])
    n_el_and_sites = len(phaseshifts[0][1])

    n_el_and_sites_expected = 0
    for el in sl.elements:
        if el in rp.ELEMENT_MIX:
            n = len(rp.ELEMENT_MIX[el])
        else:
            n = 1
        n_el_and_sites_expected += n*len([s for s in sl.sitelist if s.el == el])

    if rp.V0_REAL == "default" and not muftin:
        logger.warning(
            "Could not convert first line of PHASESHIFTS file to MUFTIN "
            "parameters. A new PHASESHIFTS file will be generated.")
        rp.setHaltingLevel(1)

    elif n_el_and_sites == n_el_and_sites_expected:
        logger.debug("Found "+str(n_el_and_sites_expected)+" blocks in PHASESHIFTS "
                        "file, which is consistent with PARAMETERS.")
        newpsGen, newpsWrite = False, False

    # Check that the phaseshifts read in have sufficient lmax
    elif n_l_values < rp.LMAX[1] + 1:
        logger.warning(
            "Maximum angular momentum LMAX in PHASESHIFTS "
            "file is lower than required by PARAMETERS. A "
            "new PHASESHIFTS file will be generated."
            )
        rp.setHaltineLevel(1)

    elif n_el_and_sites == len(sl.chemelem):
        logger.warning(
            "Found fewer blocks than expected in the "
            "PHASESHIFTS file. However, the number of blocks matches "
            "the number of chemical elements. A new PHASESHIFTS file "
            "will be generated, assuming that each block in the old "
            "file should be used for all atoms of one element.")
        rp.setHaltingLevel(1)
        oldps = phaseshifts[:]
        phaseshifts = []
        for (energy, oldenps) in oldps:
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
            phaseshifts.append((energy, enps))
        newpsGen = False
        firstline = str(len(phaseshifts[0][1])).rjust(3) + firstline[3:]

    else:
        logger.warning(
            "PHASESHIFTS file was read but is inconsistent with "
            "PARAMETERS. A new PHASESHIFTS file will be generated.")
        rp.setHaltingLevel(1)

    return phaseshifts, firstline, newpsGen, newpsWrite


def __check_consistency_energy_range(rp, phaseshifts, muftin, newpsGen, newpsWrite):
    # check whether energy range is large enough:
    checkfail = False
    er = np.arange(rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+1e-4,
                    rp.THEO_ENERGIES[2])
    psmin = round(phaseshifts[0][0]*_HARTREE_TO_EV, 2)
    psmax = round(phaseshifts[-1][0]*_HARTREE_TO_EV, 2)
    if rp.V0_REAL == "default" or isinstance(rp.V0_REAL, list):
        if isinstance(rp.V0_REAL, list):
            c = rp.V0_REAL
        else:
            c = muftin
            checkfail = not bool(muftin)
        if c and not checkfail:
            # energies at which scattering occurs
            er_inner = er + rp.FILAMENT_WF
            er_inner -= np.maximum(
                c[0],
                c[1] + (c[2]/np.sqrt(er_inner + c[3]))
                )
    else:
        try:
            er_inner = er + float(rp.V0_REAL)
        except (ValueError, TypeError):
            checkfail = True

    if checkfail:
        logger.warning(
            "Could not check energy range in PHASESHIFTS "
            "file. If energy range is insufficient, try deleting the "
            "PHASESHIFTS file to generate a new one."
            )
        return newpsGen, newpsWrite

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
                "The energy range found in the PHASESHIFTS "
                "file is smaller than the energy range requested for "
                "theoretical beams. A new PHASESHIFTS file will be "
                "generated.")
            newpsGen, newpsWrite = True, True

    return newpsGen, newpsWrite


def __check_consitency_element_order(sl, phaseshifts, eps):
    """Tries to determine if sites/elements may have been assigned the
    wrong phaseshift.

    In general at high energies and at high LMAX, heavier elements 
    (higher atomic number) have larger elastic scattering phaseshifts 
    than lighter atoms. If this is not the case in the read in 
    PHASESHIFTS file, the ordering in the file may be wrong.
    This is a very common user mistake, that often occurs when using a
    different POSCAR with pre-existing PHASESHIFTS.

    Parameters
    ----------
    sl : Slab
        Surface Slab object. Contains site types to check against.
    phaseshifts : list
        Phaseshifts list. Can be take from rp.phaseshifts.
    eps : float
        Epsilon to use when comparing phaseshift values.
        Can be taken from rp.PHASESHIFT_EPS.

    Returns
    -------
    set(str)
        Set of elements which may be assigned the wrong phaseshifts.
        If no inconsitency is detected, the set is empty.
    """
    n_ps_energies = len(phaseshifts)
    n_l_values = len(phaseshifts[0][1][0])
    n_sites = len(phaseshifts[0][1])

    # get atomic number of all site types in the order they appear in sl.sitelist
    atomic_numbers = tuple(get_atomic_number(site.el) for site in sl.sitelist)

    # get phaseshifts at highest LMAX that are larger than eps
    en_idx = -1 # look at highest energy only, TODO: enough?
    for l_value in reversed(range(4, n_l_values)): # TODO: may want a higher cutoff than 4?
        ps_sites = np.array(list(
            phaseshifts[en_idx][1][site_idx][l_value]
                    for site_idx in range(n_sites)
                    ))
        if all(abs(ps_sites) > eps):
            break
        # could not find any energy where all phaseshifts higher than eps
        # try lower angular momentum

    affected_elements = set()
    ps_pairs = list(zip(atomic_numbers, ps_sites))
    at_number = lambda item : item[0]
    pair_ps = lambda item : item[1]
    ps_pairs.sort(key=at_number) # sort by atomic number
    # go through all pairs and check that heavier elements
    # have larger phaseshifts
    prev_pair = ps_pairs[0]
    for pair in ps_pairs:
        # same element
        if at_number(pair) == at_number(prev_pair):
            continue
        # different element, where
        # at_number(pair) > at_number(prev_pair)
        # because the list was sorted
        elif abs(pair_ps(pair)) >= abs(pair_ps(prev_pair)) - eps:
            continue
        affected_elements.update(at_number(prev_pair), at_number(pair))
    affected_elements = set(get_element_symbol(el) for el in affected_elements)

    if affected_elements:
        elements_str = ", ".join(affected_elements)
        logger.warning(
            "Detected inconsitency in PHASESHIFTS: "
            "Found larger phaseshifts for lighter elements "
            f"than for heavier elements at LMAX={l_value}. "
            "This may mean that the blocks in PHASESHIFTS "
            "are in the wrong order! "
            f"Affected elements: {elements_str}."
        )
    return affected_elements


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
    energies = np.array([ps[0]*_HARTREE_TO_EV for ps in rp.phaseshifts])
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


def find_block_size(filelines):
    """Returns the periodicity of lengths found in filelines.
    
    """
    
    lens_as_chars = "".join(chr(len(l.rstrip())) for l in filelines)
    
    # allow empty lines at the end
    pass
    while lens_as_chars.endswith(chr(0)):
        lens_as_chars = lens_as_chars[:-1]
    pass
    if chr(0) in lens_as_chars:
        raise ValueError("Empty line found.")
    
    # see stackoverflow.com/questions/29481088
    period = (lens_as_chars+lens_as_chars).find(lens_as_chars, 1, -1)
    pass
    if period == -1:
        raise ValueError("Could not identify block.")
    
    return period
