"""Functions for reading and writing the PHASESHIFTS file."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-18'
__license__ = 'GPLv3+'

from itertools import combinations
import logging
import os
from pathlib import Path

import fortranformat as ff
import numpy as np

from viperleed.calc.lib.leedbase import HARTREE_TO_EV
from viperleed.calc.lib.matplotlib_utils import CAN_PLOT
from viperleed.calc.lib.matplotlib_utils import log_without_matplotlib
from viperleed.calc.lib.matplotlib_utils import prepare_matplotlib_for_calc
from viperleed.calc.lib.periodic_table import (get_atomic_number,
                                               get_element_symbol)

if CAN_PLOT:
    prepare_matplotlib_for_calc()
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages


logger = logging.getLogger(__name__)


def readPHASESHIFTS(sl, rp, readfile='PHASESHIFTS', check=True,
                    ignoreEnRange=False):
    """Read from a PHASESHIFTS file.

    Parameters
    ----------
    sl: Slab
    rp: RunParameters
    readfile: str, optional
        filename to be read. Default is 'PHASESHIFTS'.
    check: bool, optional
        Wether to check for consistence agains sl and rp.
        If False, sl and rp can be None. Default is True.
    ignoreEnRange: bool, optional
        Check wether the energy range in readfile is sufficient to
        cover the energy range requested in rp. Default is False.

    Returns
    -------
    firstline: str
        Header line of the readfile, containing Rundgren parameters.
    phaseshifts: list of tuple
        Each element is (energy, phaseshifts_at_energy) with
        phaseshifts_at_energy = [[el0_L0, el0_L1, ...], [el1_L0, ...], ...]
        where eli_Lj is the phaseshift for element i and angular momentum j.
        Therefore, len(phaseshifts) is the number of energies found,
        len(phaseshifts[0][1]) should match the number of sites, and
        len(phaseshifts[0][1][0]) is the number of different
        values of L in the phaseshift file.
    newpsGen: bool
        Whether the inconsitency found requires a full recalculation
        of the phaseshifts. This happens if:
            1) an error occurs while parsing,
            2) if check, rp.V0_real was not given and could not
               interpret firstline,
            3) LMAX of readfile is smaller than required in rp,
            4) if number of blocks inconsitent with elements in sl.
    newpsWrite: bool
        Whether the inconsitency found requires writing a new
        PHASESHIFTS file. This is distinct from newpsGen as we
        can at times generate a new file from the information read.
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
    except (ValueError, IndexError):
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
    if (lines_per_block - 1) % nel:
        logger.warning(
            "Error while trying to read PHASESHIFTS: "
            "Could not parse file: The number of blocks may not match "
            "the number given in the first line. A new PHASESHIFTS "
            "file will be generated."
        )
        rp.setHaltingLevel(1)
        return "", [], True, True

    lines_per_element = (lines_per_block - 1) // nel

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
     newpsWrite) = __check_consistency_rp_elements(sl, rp, phaseshifts,
                                                   firstline, muftin)

    if not ignoreEnRange and not newpsGen:
        newpsGen = __check_consistency_energy_range(rp, phaseshifts,
                                                    muftin, newpsGen)
        newpsWrite |= newpsGen

    # Check consitency of phaseshift values (unless we
    # already know we have to make a new file anyway)
    if not newpsGen:
        __check_consistency_element_order(rp, sl, phaseshifts)

    return firstline, phaseshifts, newpsGen, newpsWrite


def __check_consistency_rp_elements(sl, rp, phaseshifts, firstline, muftin):
    """Check that phaseshifts fit the expected number of elements and LMAX.

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
    muftin : list
        Rundgren coefficients interpreted from the firstline.

    Returns
    -------
    phaseshifts: list of tuple
        Each element is (energy, phaseshifts_at_energy) with
        phaseshifts_at_energy = [[el0_L0, el0_L1, ...], [el1_L0, ...], ...]
        where eli_Lj is the phaseshift for element i and angular momentum j.
        Therefore, len(phaseshifts) is the number of energies found,
        len(phaseshifts[0][1]) should match the number of sites, and
        len(phaseshifts[0][1][0]) is the number of different
        values of L in the phaseshift file.
    firstline: str
        Header line of the file read, containing Rundgren parameters,
        modified if newpsWrite is True.
    newpsGen: bool
        Wether the inconsitency found requires a full recalculation
        of the phaseshifts. This happens if:
            1) rp.V0_real was not given and could not interpret firstline,
            2) LMAX in phaseshifts is smaller than required in rp,
            3) if number of blocks inconsitent with elements in sl.
    newpsWrite: bool
        Wether the inconsitency found requires writing a new
        PHASESHIFTS file. This is distinct from newpsGen as we
        can at times generate a new file from the information read.
        This is the case if:
            4) There are fewer blocks than expected, but the number of
               blocks matches the number of chemical elements in sl.
               This means however, that all sites with the same element
               will use the same phaseshift.
    """
    newpsGen, newpsWrite = True, True

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
            "parameters. A new PHASESHIFTS file will be generated."
            )
        rp.setHaltingLevel(1)

    elif n_el_and_sites == n_el_and_sites_expected:
        logger.debug(f"Found {n_el_and_sites_expected} blocks in PHASESHIFTS "
                     "file, which is consistent with PARAMETERS.")
        newpsGen, newpsWrite = False, False

    # Check that the phaseshifts read in have sufficient lmax
    elif n_l_values < rp.LMAX.max + 1:  # +1 because of L=0
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
            "file should be used for all atoms of one element."
            )
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
            "PARAMETERS. A new PHASESHIFTS file will be generated."
            )
        rp.setHaltingLevel(1)

    return phaseshifts, firstline, newpsGen, newpsWrite


def __check_consistency_energy_range(rp, phaseshifts, muftin, newpsGen):
    """Check that the energy range of phaseshifts is large enough for rp."""
    checkfail = False

    er = np.arange(rp.THEO_ENERGIES.start,
                   rp.THEO_ENERGIES.stop + 1e-4,
                   rp.THEO_ENERGIES.step)
    psmin = round(phaseshifts[0][0] * HARTREE_TO_EV, 2)
    psmax = round(phaseshifts[-1][0] * HARTREE_TO_EV, 2)

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
        return newpsGen

    if (psmin > min(er_inner) or psmax < max(er_inner)):
        if (psmin > min(er_inner) and psmin <= 20.
                and psmax >= max(er_inner)):
            # can lead to re-calculation of phaseshifts every run if
            # V0r as calculated by EEASiSSS differs from 'real' V0r.
            # Don't automatically correct.
            logger.warning(
                f"Lowest value in PHASESHIFTS file ({psmin:.1f} eV) is "
                "larger than the lowest predicted scattering energy "
                f"({min(er_inner):.1f} eV). If this causes problems in the "
                "reference calculation, try deleting the PHASESHIFTS "
                "file to generate a new one, or increase the starting "
                "energy in the THEO_ENERGIES parameter."
                )
        else:
            logger.warning(
                "The energy range found in the PHASESHIFTS "
                "file is smaller than the energy range requested for "
                "theoretical beams. A new PHASESHIFTS file will be "
                "generated."
                )
            return True

    return newpsGen


def __check_consistency_element_order(rp, sl, phaseshifts,
                                     eps=None, l_max_cutoff=3):
    """Determine if elements may have been assigned wrong phaseshifts.

    In general at high energies and at high LMAX, heavier elements
    (higher atomic number) have larger elastic scattering phaseshifts
    than lighter atoms. If this is not the case in the read in
    PHASESHIFTS file, the ordering in the file may be wrong.
    This is a very common user mistake, that often occurs when
    using a different POSCAR with pre-existing PHASESHIFTS.

    Warnings will be issued if an inconsistency is found for
    elements whose atomic numbers differ by at least 2. This
    means, e.g., that Cr--Ni will not issue warnings for Fe,
    but Ti and Cu will.

    Parameters
    ----------
    rp : Rparams
        Run parameters.
    sl : Slab
        Surface Slab object. Contains site types to check against.
    phaseshifts : list
        Phaseshifts list, as read in from a PHASESHIFTS file.
    eps : float, optional
        Tolerance to use when comparing phaseshift values. If None
        or not given, eps == rp.PHASESHIFT_EPS. Default is None.
    l_max_cutoff : int, optional
        Lower cutoff for angular momentum quantum number to be used for
        comparison of phaseshifts. Heavier elements should always have
        higher phaseshifts in the limit of high energy and high angular
        momentum. Behaviour for low energy/low angular momentum is not
        as clear cut, as phaseshifts may cross zero and are pi periodic.
        Default is 3.

    Returns
    -------
    may_have_wrong_phaseshifts : set
        Set of chemical elements (str) which may be assigned the wrong
        phaseshifts. If no inconsitency is detected, the set is empty.
    """
    if not eps:
        eps = rp.PHASESHIFT_EPS
    n_l_values = len(phaseshifts[0][1][0])
    n_sites_in_ps = len(phaseshifts[0][1])  # No. species in PHASESHIFTS

    # Get atomic number of all site types in the order
    # they appear in sl.sitelist, taking into account
    # that some sites may have mixed occupation
    element_mix = []
    for site in sl.sitelist:
        element_mix.extend(
            site.mixedEls  # Mixed occupation
            or [site.el]   # Only one element
            )
    # get atomic numbers from element symbols
    real_el = lambda el: (rp.ELEMENT_RENAME[el] if el in rp.ELEMENT_RENAME
                          else el)
    atomic_numbers = tuple(get_atomic_number(real_el(el))
                           for el in element_mix)

    # Get phaseshifts at highest LMAX that are larger than eps
    en_idx = -1  # look at highest energy only
    for l_value in reversed(range(l_max_cutoff, n_l_values)):
        ps_sites = np.array(list(
            phaseshifts[en_idx][1][site_idx][l_value]
            for site_idx in range(n_sites_in_ps)
            ))
        if all(abs(ps_sites) > eps):
            break
    else:
        # At least one of the phaseshifts are smaller than eps at this energy
        logger.warning(
            "Could not check consistency of PHASESHIFTS file: "
            f"PHASESHIFTS for some sites are smaller than {eps} at the largest "
            "energy for LMAX >= {l_max_cutoff}. This may happen if you are "
            "using very light scatterers (e.g. hydrogen)."
            )
        rp.setHaltingLevel(1)
        return set()

    affected_elements = set()
    ps_pairs = list(zip(atomic_numbers, ps_sites))
    at_number = lambda item : item[0]
    pair_ps = lambda item : item[1]
    ps_pairs.sort(key=at_number) # sort by atomic number

    # Go through all pairs and check that heavier
    # elements have larger phaseshifts
    for pair_1, pair_2 in combinations(ps_pairs, 2):
        at_num_change = at_number(pair_1) - at_number(pair_2)
        if at_num_change <= 0:  # Same element or redundant iteration
            continue
        # Different element. Check that phaseshifts are also
        # ordered, and complain only if elements are not
        # close neighbors
        if (abs(pair_ps(pair_1)) >= abs(pair_ps(pair_2)) - eps
                and at_num_change > 2):
            continue
        affected_elements.update((at_number(pair_1), at_number(pair_2)))
    # get element symbols
    may_have_wrong_phaseshifts = set(get_element_symbol(el)
                                     for el in affected_elements)

    if may_have_wrong_phaseshifts:
        elements_str = ", ".join(affected_elements)
        logger.warning(
            "Detected inconsitency in PHASESHIFTS: "
            "Found larger phaseshifts for lighter elements "
            f"than for heavier elements at LMAX={l_value}. "
            "This may mean that the blocks in PHASESHIFTS "
            "are in the wrong order! "
            f"Affected elements: {elements_str}. "
            "You may want to regenerate the PHASESHIFTS file or "
            "reorder the blocks (see utility ModifyPhaseshifts). "
        )
        rp.setHaltingLevel(1)
    return may_have_wrong_phaseshifts


def writePHASESHIFTS(firstline, phaseshifts, file_path=Path()/'PHASESHIFTS'):
    """Takes phaseshift data and writes it to a PHASESHIFTS file."""
    _file_path = Path(file_path)
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
        with open(_file_path, 'w') as wf:
            wf.write(output)
        logger.debug(f"Wrote to {_file_path} successfully.")
    except Exception:
        logger.error("Exception while writing PHASESHIFTS file: ",
                     exc_info=True)
        raise
    return


@log_without_matplotlib(logger, msg='Skipping Phaseshift plotting.')
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
    energies = np.array([ps[0]*HARTREE_TO_EV for ps in rp.phaseshifts])
    ps_vals = np.array([ps[1] for ps in rp.phaseshifts])

    figsize = (7, 4)
    linewidth = 1
    nlplot = min(np.shape(ps_vals)[-1], rp.LMAX.max+1)

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
                    linewidth=linewidth, label=fr"$\ell$ = {j}",
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
    """Returns the periodicity of lengths found in filelines."""
    lens_as_chars = "".join(chr(len(l.rstrip())) for l in filelines)

    # Allow empty lines at the end
    while lens_as_chars.endswith(chr(0)):
        lens_as_chars = lens_as_chars[:-1]
    if chr(0) in lens_as_chars:
        raise ValueError("Empty line found.")

    # See stackoverflow.com/questions/29481088
    period = (lens_as_chars * 2).find(lens_as_chars, 1, -1)
    if period == -1:
        raise ValueError("Could not identify block.")

    return period
