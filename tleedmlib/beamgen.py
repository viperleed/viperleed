# -*- coding: utf-8 -*-
"""Module beamgen of viperleed.tleedmlib.

@author: Alexander Imre
@author: Florian Kraushofer

Creates the BEAMLIST file for TensErLEED calculations.

Original version by Florian Kraushofer (2020) was a wrapper for
the Fortran beamgen script by Uli Löffler et al. Complete refactor
by Alexander Imre (2023) that removes the dependency on the Fortran
script and instead accomplishes the same in python. This is more
flexible and allows for more accurate calculations.
"""

import logging
from pathlib import Path

import fortranformat as ff
import numpy as np

from viperleed.guilib.base import get_equivalent_beams, BeamIndex
from viperleed.tleedmlib import symmetry
from viperleed.tleedmlib.leedbase import BOHR_TO_ANGSTROM, HARTREE_TO_EV
from viperleed.tleedmlib.leedbase import getLEEDdict

H_BAR_SQ_OVER_2M = 0.5 * HARTREE_TO_EV * BOHR_TO_ANGSTROM**2  # h**2/2m
logger = logging.getLogger('tleedm.beamgen')


def calc_and_write_beamlist(slab, rpars, domains=False,
                            beamlist_name='BEAMLIST'):
    """Calculates and writes the contents for the file BEAMLIST.

    BEAMLIST contains a list of all diffraction beams that will be
    used for internal calculations (as opposed to IVBEAMS, which
    contains the beams to be output). The file lists the beams with
    indices as float values and a lower cutoff energy below which
    the beam is evanescent (i.e., does not leave the surface and
    has intensity 0). The used format is defined in make_beamlist_lines
    and dictated by TensErLEED and the legacy beamgen scripts that were
    used before. Note that BEAMLIST is not read directly by refcalc,
    but instead all input files for the refcalc are combined into one
    string (by collectFIN in iorefcalc.py) and then piped in.

    NB: The energies calculated here are slightly higher than those
    from beamgenv3 (and beamgen.v1.7) because we use *more accurate*
    values for unit conversions. The legacy code used these rounded
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
    slab : Slab
        Slab object. A .bulkslab will be added if it was not
        yet created for this slab.
    rpars : Rparams
        Run parameters.
    domains : bool, optional
        Flag to indicate if performing a domain calculation,
        by default False.
    beamlist_name : str or pathlike, optional
        Filename to be written, by default "BEAMLIST".
    """
    if slab.bulkslab is None:
        slab.bulkslab = slab.makeBulkSlab(rpars)
        symmetry.findSymmetry(slab.bulkslab, rpars)

    # Use guilib to generate list of beams
    leed_parameters = getLEEDdict(slab, rpars)
    leed_parameters['screenAperture'] = 180  # All possible beams
    leed_parameters['eMax'] = _get_emax_for_evanescent_beams(slab, rpars,
                                                             domains)

    # use **only** beams from domain specified in rpars.SUPERLATTICE
    # beams come pre-sorted from get_equivalent_beams()
    equivalent_beams = get_equivalent_beams(leed_parameters, domains=0)

    # Log beam groups for debugging
    _log_beamgroups(equivalent_beams)

    # Strip away symmetry information, and sort into scattering subsets
    beam_subsets = _split_into_subsets([BeamIndex(hk)  # str->Fraction
                                        for hk, *_ in equivalent_beams])

    inv_bulk_surf_vectors = slab.bulkslab.reciprocal_vectors
    all_energies = []
    all_indices_arr = []
    for beam_indices in beam_subsets.values():
        indices_arr = np.array(beam_indices, dtype='float64')  # from Fraction
        # Calculate cutoff energy for each beam
        energies = (                                                            # TODO: we could probably remove the energies from BEAMLIST completely. It seems they are not used in TensErLEED (see subroutine READIN in lib.tleed.f). Would need to remove it from here, and readBEAMLIST in beams.py.
            H_BAR_SQ_OVER_2M
            * np.sum(np.dot(indices_arr, inv_bulk_surf_vectors)**2, axis=1)
            )

        # generate file contents for beam subset
        all_indices_arr.append(indices_arr)
        all_energies.append(energies)
    beamlist_content = make_beamlist_string(all_indices_arr,
                                            all_energies,
                                            rpars.TL_VERSION)
    # get highest energy considered; groups may have different shapes
    max_energy = max(np.max(group) for group in all_energies)
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


def _get_emax_for_evanescent_beams(slab, rpars, domains):
    """Return an energy cut-off that will generate also evanescent beams.

    The full-dynamic calculation considers scattering of propagating
    beams as well as evanescent ones. Among the latter, only those
    that survive attenuation when propagating between two layers are
    considered. Beams are considered surviving if their amplitude is
    attenuated by a factor less than ATTENUATION_EPS. Here we increase
    the highest energy by the corresponding factor.

    Parameters
    ----------
    slab : Slab
        The structure for which beams are to be calculated.
    rpars : Rparams
        Run parameters.
    domains : bool
        Flag to indicate if performing a domain calculation.

    Returns
    -------
    e_max : float
        The energy (in electronvolts) that generates also the
        correct evanescent beams.
    """
    if not domains:
        d_min = slab.getMinLayerSpacing()
    else:
        d_min = min(dp.sl.getMinLayerSpacing() for dp in rpars.domainParams)
    d_min *= 0.7                                                                # TODO: may want to complain if this is small as it will give a huge load of beams (and may mean different LAYER_CUTS should be used).

    e_max = rpars.THEO_ENERGIES[1]
    e_max += H_BAR_SQ_OVER_2M * (np.log(rpars.ATTENUATION_EPS) / d_min)**2
    return e_max


def _log_beamgroups(equivalent_beams):
    """Log information on beam groups if desired."""
    if logger.getEffectiveLevel() > 5:
        return

    full_log_msg = ['Equivalent beams:',
                    '(   h     |   k     ),group,']
    for beam_hk, beam_group, *_ in equivalent_beams:
        index = BeamIndex(beam_hk)
        line = f'{format(index, "(4,4)s")}, {beam_group:4},\n'
        full_log_msg.append(line)

    # Split log message into two parts
    # Thee first 15 lines are intended for log-level verbose...
    logger.log(level=5, msg='\n'.join(full_log_msg[:15]))
    # ...the rest only at very verbose
    logger.log(level=1, msg='\n'.join(full_log_msg[15:]))


def make_beamlist_lines(all_indices, all_energies, tl_version):
    """Creates contents for file BEAMLIST for each beamset in the format
    used be the legacy beamgen scripts by U. Loeffler and R. Doell.

    Beam sets that are not related by a bulk translation are
    separated into different "blocks". The format is the same
    as used be the legacy beamgen scripts by U. Loeffler and R.
    Doell. See also help(get_beam_scattering_subsets).

    Parameters
    ----------
    all_indices : list
        Elements are numpy arrays with shape (n_beams_subset, 2).
        Each element is one subset, containing the hk-indices
        (diffraction orders) of the beams. n_beams_subset is the
        number of beams in each subset.
    all_energies : list
        Lower cut-off energies for the beams. Each element is one
        subset. Elements are numpy arrays with shape (n_beams_subset,).
    tl_version : float
        Version of TensErLEED. To be taken from Rparams.TL_VERSION.
        This values decides the format of the output string.

    Returns
    -------
    str
        String representation of the contents of the BEAMLIST file.

    Raises
    ------
    ValueError
        If indices and energies have incompatible shapes.
    ValueError
        If one of the beam sets exceeds the maximum size
        currently supported in the Fortran code.
    """
    # Set up Fortran format as was used by beamgen.
    # TensErLEED v1.7 and higher used beamgen v1.7; earlier versions
    # used v3. This matters because the format changed slightly:
    # beamgen v1.7 had I5, beamgen v3 had I4. Also the maximum
    # number of beams per beamset has changed since.
    if tl_version < 1.7:
        fmt = {
            'beamset': ff.FortranRecordWriter(
                "2F10.5,2I3,10X,'E =  ',F10.4,2X,'NR.',I4"
                ),
            'n_beams': ff.FortranRecordWriter('15I4')
            }
        n_beams_max = 9999
    else:
        fmt = {
            'beamset': ff.FortranRecordWriter(
                "2F10.5,2I3,10X,'E =  ',F10.4,2X,'NR.',I5"
                ),
            'n_beams': ff.FortranRecordWriter('15I5')
            }
        n_beams_max = 99999

    beam_nr = 1
    content = ''
    for indices, energies in zip(all_indices, all_energies):
        n_beams = indices.shape[0]
        if energies.shape != (n_beams,) or indices.shape != (n_beams, 2):
            raise ValueError(
                f'Incompatible size of indices (shape={indices.shape})'
                f'and energies (shape={energies.shape}).'
                )

        # first line contains number of beams
        if n_beams > n_beams_max:
            raise ValueError(
                f'Too many beams in beam set ({n_beams}). TensErLEED v'
                f'{tl_version:.2f} code supports at most {n_beams_max}.'
                )
        content += fmt['n_beams'].write([n_beams]) + '\n'

        # iterate over all beams and format lines
        for beam_hk, energy in zip(indices, energies):
            line = fmt['beamset'].write([*beam_hk, 1, 1, energy, beam_nr])
            content += line + '\n'
            beam_nr += 1

    return content


def get_beam_scattering_subsets(beam_indices_raw):
    """Return beam scattering subsets and reduced indices for a list of beams.

    LEED diffraction beams are grouped into subsets for the computation
    of reflection/transmission matrices. In the full-dynamic scattering
    calculation (refcalc), one needs to consider that one beam can be
    scattered into another. However, this is only possible, if the
    beam wave-vectors are related by the *bulk* unit cell. I.e., beam
    (1/2, 0) can be scattered into (3/2, 0), but not into (1, 0).
    This property can be used in the refcalc to simplify calculations
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
    beam_indices_raw : Sequence
        Beam indices in (h, k) form.

    Returns
    -------
    subset_classes : list
        The unique first-Brillouin-zone beam indices for the beam
        subsets, sorted by their length. The length gives the number
        of subsets.
    reduced_indices : list
        Reduced version of the indices in beam_indices_raw.
    """
    reduced_indices = [(h%1, k%1) for (h, k) in beam_indices_raw]
    subset_classes = set(reduced_indices)

    # sort order of subsets by |(h_red, k_red)|^2
    subset_classes = sorted(subset_classes, key=lambda hk: hk[0]**2 + hk[1]**2)

    return subset_classes, reduced_indices


def _split_into_subsets(beam_indices_raw):
    """Split raw beam indices into bulk-independent subsets.

    See also help(get_beam_scattering_subsets).

    Parameters
    ----------
    beam_indices_raw : Sequence
        Beam indices in (h, k) form.

    Returns
    -------
    beam_subsets : dict
        Keys are the unique first-Brillouin-zone beam indices for
        the beam subsets. They are insertion-ordered by their norm.
        Values are all the indices in beam_indices_raw that belong
        to the specific subset.
    """
    (subset_classes,
     reduced_indices) = get_beam_scattering_subsets(beam_indices_raw)           # TODO: create test case to check that len(subset_classes) == np.linalg.det(rp.SUPERLATTICE)
    beam_subsets = {red_index: [] for red_index in subset_classes}
    for raw_index, red_index in zip(beam_indices_raw, reduced_indices):
        beam_subsets[red_index].append(raw_index)
    return beam_subsets
