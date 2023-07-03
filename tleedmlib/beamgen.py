# -*- coding: utf-8 -*-
"""
Module beamgen of viperleed.tleedmlib.

Creates the BEAMLIST file for TensErLEED calculations.

@author: Alexander Imre
@author: Florian Kraushofer

Original version by Florian Kraushofer (2020) was a wrapper for the fortran
beamgen script by Uli Löffler et al. Complete refactro by Alexander Imre (2023)
that removes the dependency on the fortran script and instead accomplishes the
same in python. This is more flexible and allows for more accurate calculations.
"""

import logging
from pathlib import Path

import numpy as np
import fortranformat as ff

from viperleed.guilib.base import get_equivalent_beams, BeamIndex
from viperleed.tleedmlib import symmetry
from viperleed.tleedmlib.leedbase import (HARTREE_TO_EV,
                                          BOHR_TO_ANGSTROM,
                                          ANGSTROM_TO_BOHR)

logger = logging.getLogger('tleedm.beamgen')


def calc_and_write_beamlist(sl, rp, domains=False, beamlist_name='BEAMLIST'):
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
    rp : Rparams
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
        d_min = min((dp.sl.getMinLayerSpacing() for dp in rp.domainParams))*0.7

    # convergence criterion for refcalc;
    # beams are propagated between layers if relative change is larger than this
    conv_crit = rp.ATTENUATION_EPS
    # effective cutoff energy to use (scale to correct units)
    e_max_eff = (e_max +
                 (np.log(conv_crit)/(d_min*ANGSTROM_TO_BOHR))**2
                 / 2 *HARTREE_TO_EV)

    # use guilib to generate list of beams
    leed_parameters = {
        'eMax': e_max_eff,
        'surfBasis': surf_ucell,
        'SUPERLATTICE': rp.SUPERLATTICE,
        'surfGroup': sl.foundplanegroup,
        'bulkGroup': sl.bulkslab.foundplanegroup,
        'screenAperture': 180,  # all beams, because internal calculation
    }
    # use **only** beams from domain specified in rp.SUPERLATTICE
    # beams come pre-sorted from get_equivalent_beams()
    equivalent_beams = get_equivalent_beams(leed_parameters, domains=0)
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
    beamlist_content = ''
    # for every subset calculate energies, sort and generate partial string
    for beam_indices in beam_subsets:
        # convert to float array
        indices_arr = np.array(beam_indices, dtype='float64')
        # calculate cutoff energy for each beam and scale to correct units
        energies = (np.sum(np.dot(indices_arr, inv_bulk_surf_vectors)**2,
                           axis=1)
                    /2 *HARTREE_TO_EV *BOHR_TO_ANGSTROM**2) 

        # generate file contents for beam subset
        all_indices_arr.append(indices_arr)
        all_energies.append(energies)
    beamlist_content = make_beamlist_string(all_indices_arr,
                                            all_energies,
                                            rp.TL_VERSION)
    max_energy = np.max(all_energies)
    logger.debug(f'Highest energy considered in BEAMLIST: {max_energy:.2f}eV')

    # write to file
    write_file_path = Path(beamlist_name)
    try:
        with open(write_file_path, 'w', encoding='utf-8') as file:
            file.write(beamlist_content)
    except Exception:
        logger.error(f'Unable to write file {beamlist_name}')
        raise

    logger.debug('Wrote to BEAMLIST successfully.')


def make_beamlist_string(all_indices, all_energies, tl_version=1.7):
    """Creates contents for file BEAMLIST for each beamset in the format
    used be the legacy beamgen scripts by U. Loeffler and R. Doell.

    Parameters
    ----------
    all_indices : list(np.ndarray, shape=(n_beams_subset, 2))
        Indices (diffraction orders) of the beams. n_beams_subset is
        the number of beams in each subset.
    all_energies : np.ndarray, shape=(n_beams_subset,)
        Lower cutoff energies for the beams.
    tl_version : float, optional
        Version of TensErLEED, by default 1.7. To be taken from
        Rparams.TL_VERSION. This values decides the format of the output string.

    Returns
    -------
    str
        String representation of the contents of the BEAMLIST file.

    Raises
    ------
    ValueError
        If indices and energies have incompatible shapes.
    """
    # Set up Fortran format as was used by beamgen.
    # TensErLEED v1.7 and higher used beamgen v1.7; earlier versions used v3
    # This matters because the format changed slightly:
    # beamgen v1.7 had I5, beamgen v3 had I4 for some reason
    if tl_version >= 1.7:
        beamlist_format = ff.FortranRecordWriter(
            "2F10.5,2I3,10X,'E =  ',F10.4,2X,'NR.',I4"
            )
    else:
        beamlist_format = ff.FortranRecordWriter(
            "2F10.5,2I3,10X,'E =  ',F10.4,2X,'NR.',I5"
            )
    beam_nr = 1
    content = ''
    for indices, energies in zip(all_indices, all_energies):
        n_beams = indices.shape[0]
        if not energies.shape == (n_beams,) or not indices.shape == (n_beams,
                                                                     2):
            raise ValueError(
                f'Incompatible size of indices (shape={indices.shape})'
                f'and energies (shape={energies.shape}).')

        # first line contains number of beams
        content += ff.FortranRecordWriter('10I3').write([n_beams]) + '\n'
        # TODO: why limit to 999 beams?

        # iterate over all beams and format lines
        for beam_hk, energy in zip(indices, energies):
            line = beamlist_format.write([*beam_hk, 1, 1, energy,
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
    beam_indices_raw : list(BeamIndex)
        List of beam indices.

    Returns
    -------
    subset_classes : list
        A tuple containing the unique first Brillouin zone beam indices
        for the beam subsets. The length gives the number of subsets.
    reduced_indices : list
        A list of corresponding reduced indices for beam indices_raw.
    """

    reduced_indices = [(h%1,k%1) for (h,k) in beam_indices_raw]
    subset_classes = set(reduced_indices)

    # sort order of subsets by |(h_red, k_red)|^2
    subset_classes = sorted(subset_classes,
                            key= lambda index: index[0]**2 + index[1]**2)

    return subset_classes, reduced_indices
