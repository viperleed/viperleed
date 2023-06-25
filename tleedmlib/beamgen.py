# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:24:16 2020

@author: Florian Kraushofer
@author: Alexander Imre
"""

import os
import logging
import subprocess
from pathlib import Path

import fortranformat as ff
import numpy as np


from viperleed.tleedmlib.leedbase import (HARTREE_TO_EV,
                                BOHR_TO_ANGSTROM,
                                ANGSTROM_TO_BOHR)
from viperleed.tleedmlib import symmetry
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


def generate_beamlist(sl, rp, domains=False, beamlist_name="BEAMLIST"):
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
        symmetry.findSymmetry(sl.bulkslab, rp)

    e_max = rp.THEO_ENERGIES[1]
    surf_ucell = sl.surface_vectors
    inv_bulk_surf_vectors = sl.bulkslab.reciprocal_vectors

    if not domains:
        d_min = sl.getMinLayerSpacing() * 0.7
    else:
        d_min = min([dp.sl.getMinLayerSpacing() for dp in rp.domainParams])*0.7

    conv_crit = rp.ATTENUATION_EPS  # convergence criterion for refcalc
    # effective cutoff energy to use (scale to correct units)
    e_max_eff = (e_max +
                 (np.log(conv_crit)/(d_min*ANGSTROM_TO_BOHR))**2
                 / 2 *HARTREE_TO_EV)

    # use guilib to generate list of beams
    leedParameters = {
        'eMax': e_max_eff,
        'surfBasis': surf_ucell,
        'SUPERLATTICE': rp.SUPERLATTICE,
        'surfGroup': sl.foundplanegroup,
        'bulkGroup': sl.bulkslab.foundplanegroup,
        'screenAperture': 180,  # all beams, since this is for internal calculation
    }
    # use **only** beams from domain specified in rp.SUPERLATTICE
    equivalent_beams = get_equivalent_beams(leedParameters, domains=0)
    # strip away symmetry group information
    beam_indices_raw = list(BeamIndex(beam[0]) for beam in equivalent_beams)
    subset_classes, reduced_indices = get_beam_scattering_subsets(beam_indices_raw) # TODO: create test case to check that len(subset_classes) == np.linalg.det(rp.SUPERLATTICE)

    # sort beams into scattering subsets
    beam_subsets = [[] for set in range(len(subset_classes))]
    for index, red_index in zip(beam_indices_raw, reduced_indices):
        applicable_subset = subset_classes.index(red_index)
        beam_subsets[applicable_subset].append(index)

    all_energies = []
    all_indices_arr = []
    beamlist_content = ""
    # for every subset calculate energies, sort and generate partial string
    for beam_indices in beam_subsets:
        # convert to float array
        indices_arr = np.array(beam_indices, dtype="float64")
        # calculate cutoff energy for each beam
        energies = (np.sum(np.dot(indices_arr, inv_bulk_surf_vectors)**2, axis=1)
                    /2 *HARTREE_TO_EV *BOHR_TO_ANGSTROM**2)  # scale to correct units

        # sort beams by energy (same as sorting by |G|)
        sorting_order = np.argsort(energies)
        energies, indices_arr = energies[sorting_order], indices_arr[sorting_order]

        # generate file contents for beam subset
        all_indices_arr.append(indices_arr)
        all_energies.append(energies)
    beamlist_content = make_beamlist_string(all_indices_arr, all_energies)
    max_energy = max((np.max(energies) for energies in all_energies))
    logger.debug(f"Highest energy considered in BEAMLIST: {max_energy:.2f}eV")

    # write to file
    write_file_path = Path(beamlist_name)
    try:
        with open(write_file_path, 'w') as file:
            file.write(beamlist_content)
    except Exception:
        logger.error(f"Unable to write file {beamlist_name}")
        raise

    logger.debug("Wrote to BEAMLIST successfully.")
    return


def make_beamlist_string(all_indices, all_energies):
    """Creates contents for file BEAMLIST for each beamset in the format
    used be the legacy beamgen scripts by U. Loeffler and R. Doell.

    Parameters
    ----------
    all_indices : list(np.ndarray, shape=(n_beams_subset, 2))
        Indices (diffraction orders) of the beams. n_beams_subset is
        the number of beams in each subset.
    all_energies : np.ndarray, shape=(n_beams_subset,)
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
    # set up Fortran format as was used by beamgen
    beamlist_format = ff.FortranRecordWriter(
        "2F10.5,2I3,10X,'E =  ',F10.4,2X,'NR.',I5"
        # beamgen v1.7 had I5, beamgen v3 had I4 for some reason
        )
    beam_nr = 1
    content = ""
    for indices, energies in zip(all_indices, all_energies):
        n_beams = indices.shape[0]
        if not energies.shape == (n_beams,) or not indices.shape == (n_beams,
                                                                     2):
            raise ValueError(
                f"Incompatible size of indices (shape={indices.shape})"
                f"and energies (shape={energies.shape}).")

        # first line contains number of beams
        content += ff.FortranRecordWriter('10I3').write([n_beams]) + '\n'
        # TODO: why limit to 999 beams?

        # iterate over all beams and format lines
        for (beam_h, beam_k), energy in zip(indices, energies):
            line = beamlist_format.write([beam_h, beam_k, 1, 1, energy,
                                          beam_nr])
            content += line + '\n'
            beam_nr += 1

    return content


def get_beam_scattering_subsets(beam_indices_raw):
    """Takes a list of beam_indices and returns the beam scattering
    subsets and reduced indices for sorting.

    LEED diffraction beams are grouped into subsets for the computation
    of reflection/transmission matrices. In the full dynamic scattering
    calculation (refcalc), one needs to consider that one beam can be
    scattered into another. However, this is only possible, if the
    beam wave-vectors are related by the *bulk* unit cell. I.e. a beam
    (1/2,0) can be scattered into (3/2,0), but not into (1,0).
    This property can be used in the refcalc, to simplify calculations
    by making the reflection/transmission matrices block-diagonal.
    To enable this, we need to group the beams accordingly in BEAMLIST.

    This function takes a list of beams (as BeamIndex objects) and
    calculates reduced indices via h_red = h%1, k_red = k%1 (wrapping
    the beams back in the first Brillouin zone). It then takes a set
    of the reduced indices, to generate unique identifiers of the
    subsets. The first (and possibly only) subset contains, by
    definition, the integer beams, starting with (0|0).

    Parameters
    ----------
    beam_indices_raw : [BeamIndex]
        List of beam indices.

    Returns
    -------
    tuple
        A tuple containing the unique first Brillouin zone beam indices
        for the beam subsets. The length gives the number of subsets.
    list
        A list of corresponding reduced indices for beam indices_raw.
    """

    reduced_indices = [(h%1,k%1) for (h,k) in beam_indices_raw]
    subset_classes = set(reduced_indices)

    # sort order of subsets by |(h_red, k_red)|^2
    by_h_k_red = lambda index: index[0]**2 + index[1]**2
    subset_classes = tuple(sorted(subset_classes, key= by_h_k_red))

    return subset_classes, reduced_indices