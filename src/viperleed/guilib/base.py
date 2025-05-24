"""
======================================
  ViPErLEED Graphical User Interface
======================================
      *** module guilib.base ***

Created: 2020-01-11
Author: Michele Riva

This module contains base functions and classes shared by the whole GUI
"""

import ast
import copy
import sys
import re
from collections import defaultdict

import numpy as np

from viperleed.guilib.classes.beamindex import BeamIndex
from viperleed.guilib.classes.planegroup import PlaneGroup
from viperleed.guilib.leedsim.classes.oldleedpatterns import LEEDPattern


###############################################################################
#                                  FUNCTIONS                                  #
###############################################################################


def get_equivalent_beams(leed_parameters, *other_leed_parameters, domains=None):
    """
    Generates a sorted list of LEED beams, including their fractional
    indices and their symmetry-equivalence

    Parameters
    ----------
    leedParameters: dictionary
      The following keys are needed
      - 'eMax': float
                maximum primary beam energy used in the TensErLEED calculation

      - 'surfBasis': 2x2 numpy array of floats
                     unit vectors in Cartesian coordinates defining the basis
                     vectors aS and bS of the surface unit cell.
                     aS = basis[0], bS = basis[1]

      - 'SUPERLATTICE': 2x2 numpy array of ints
                        Superlattice matrix that defines the relation between
                        the bulk unit vectors (aB, bB) and the surface ones.
                        m = [[m11, m12],[m21, m22]]
                        aS = m11*aB + m12*bB
                        bS = m12*aB + m22*bB

      - 'surfGroup': string
                     plane group of the surface

      - 'bulkGroup': string
                     plane group of the bulk

      The following keys are optional
      - 'bulk3Dsym': string or array-like
            This parameter is used to describe the isomorphic part (i.e.,
            neglecting translations) of screw axes and glide planes orthogonal
            to the surface.
            -- when passing a string, one the following formats is
               required (with or without white spaces):
                   * "r(#, #, ...), m([#, #], ...)"
                   * "m([#, #], ...), r(#, #, ...)"
                   * "r(#, #, ...)"
                   * "m([#, #], ...)"
               The quantities above are:
                   * "r(...)" -- list of rotation orders for screw axes
                                (acceptable: 2, 3, 4, 6)
                   * "m(...)" -- list of in-plane directions lying on the glide
                                 plane, expressed in 'fractional coordinates' of
                                 the bulk unit vectors, i.e., "[i,j]" represents
                                 the vector i*a + j*b. Acceptable: [1,0], [0,1],
                                 [1,1], [1,-1], [1,2], [2,1]
            -- when passing an array-like, it should be a 'list' of 2x2
               'matrices' with integer entries (floats will be rounded and cast
               to int).
            NB: by 2020-06-24, no check is performed on whether the operations
                above are actually compatible with the shape of the unit cell!
      - 'screenAperture': float
            This parameter can be used to define the aperture of the solid angle
            captured by the LEED screen in degrees. Acceptable values are
            between zero and 180

    *other_leed_parameters: dictionaries, optional
      unpacked list of LEED parameters for additional reconstructions
      that may be present on the surface as an incoherent superposition with the
      mandatory first argument. See the leed_parameters argument for details of
      the keys. The function will return the list of equivalent superposed
      beams, i.e., superposed beams are symmetry-equivalent only if they are the
      sum of beams equivalent in each of the single structures. When multiple
      parameters are passed:
      - the maximum energy used is the largest among the eMax values among the
        parameters
      - the same screenAperture will be used for all, corresponding to the
        largest among the values passed
      - consistency of the bulk unit cells, symmetry groups and bulk3Dsym is
        checked

    domains: iterable, int or None, default=None
             Domain indices for which the beams get exported
             The indices (zero-based) follow the following convention:
             - dom = 0 -> domain whose SUPERLATTICE is passed as a key of
                          leedParameters
             - dom = 1 ... -> domains generated from SUPERLATTICE as a result
                              of mirror operations.
                              Glide planes count as mirrors;
                              Only one mirror/glide is used among those
                              parallel to each other
             - dom = ... -> domains generated from SUPERLATTICE as a result of
                            rotation operations.
                            Only one rotation axis is used among the those
                            equivalent to each other (e.g., only C4 at origin,
                            not C4 at cell center)
             Since the ordering of the domains is not necessarily the same as
             the order of operations, one needs to rely on the user making the
             right choice when asking for domains.

             If None, all domains are used. If a single integer, only that
             domain is considered.

             If multiple LEED parameter dictionaries are passed, domains should
             be a list of lists with as many entries as there are LEED
             parameters, None or a single integer.

    Returns
    -------
    list of tuples [(name_0, id_0), (name_1, id_1), ...]

        name_i: str
                fractional indices of the beam in the form
                '(num_h/den_h, num_k/den_k)'

        id_i: int
              representing the grouping index for symmetry-equivalent beams.
              All symmetry-equivalent beams share the same abs(id).
                id < 0 for glide-extinct beams,
                id = 0 for spot (0, 0),
                id > 0 for non-extinct beams
              For spots originating from the superposition of different
              domains, id is positive even if some of the domains contribute
              with glide-extinct beams.
              ids are generated based on the (h, k) values, so that beams
              closer to the (0, 0) spot have smaller abs(id)

        The entries in the list are sorted in this order:
            abs(id)  [i.e., energy], -h-k, -h
        The sorting choice for the (h, k) indices makes 'nicer' beams
        appear earlier (e.g., 1|0 before -1|0).
    """
    # process arguments to allow the function to accept the LEED parameters
    # being passed as a single "list" of dictionaries, or as an unpacked "list"
    # of dictionaries
    if isinstance(leed_parameters, dict):
        leed_parameters = (leed_parameters, *other_leed_parameters)

    # check that the domains keyword argument is consistent with the number of
    # parameters passed
    if domains is None:
        domains = [None]*len(leed_parameters)
    elif isinstance(domains, int):
        domains = [[domains]]*len(leed_parameters)
    elif not hasattr(domains, '__len__'):
        raise TypeError("The keyword argument 'domains' should be either an "
                        f"iterable, an integer or None. Found {type(domains)} "
                        "instead.")
    elif len(domains) != len(leed_parameters):
        raise ValueError("Not as many domains as LEED parameters passed."
                         f"Expected {len(leed_parameters)}, found "
                         f"len(domains)")

    leed_parameters, all_leed = check_multi_leed_params(leed_parameters)

    for i, leed in enumerate(all_leed):
        if domains[i] is None:
            domains[i] = range(leed.n_domains)

    # Now use a logic similar to the one in LEEDPattern.get_equivalentSpots
    # for figuring out which superimposed beams are indeed equivalent (i.e.,
    # they are superposition of equivalent beams), starting from all the beams
    # of each structure. Also, keep track of which beams are extinct, as this
    # information needs to pass on to the end.
    all_beams = []
    all_extinct = []
    for doms, leed in zip(domains, all_leed):
        fract, groups, *_ = zip(*leed.get_equivalentSpots(domains=doms))

        # use the group indices to create dictionaries of index: beams
        beams = defaultdict(set)
        extinct = []
        for beam, group in zip(fract, groups):
            if group < 0:
                extinct.append(beam)
            beams[group].add(beam)

        # now re-index the dictionaries such that there are as many keys as
        # beams and each entry is
        #    beam: list of equivalent beams (including beam itself)
        all_beams.append({beam: eq_beams
                          for eq_beams in beams.values()
                          for beam in eq_beams})
        all_extinct.append(extinct)

    all_beams_cp = copy.deepcopy(all_beams)

    # Flatten the list of all beams, keeping only uniques
    flat_beams = set(beam for beams in all_beams for beam in beams)

    # and iterate through each with the same logics as in
    # LEEDPattern.get_equivalentSpots
    # Example on square bulk:
    #   p(2x2)-pmm + c(2x2)-pm[1 0]
    #   overlapping spots are {1 | 0}, {1/2 | 1/2}, and {1 | 1}
    #
    #   p(2x2): ( 1 | 1), (-1 |  1), (-1 | -1), (1 | -1) equivalent
    #   c(2x2): ( 1 | 1), ( 1 | -1) equivalent
    #           (-1 | 1), (-1 | -1) equivalent
    #   -> ( 1 | 1), ( 1 | -1) equivalent
    #      (-1 | 1), (-1 | -1) equivalent, but not equivalent to the others
    #
    #
    #   THIS IS THE SAME CODE AS IN LEEDPattern, probably will consolidate the
    #   two later on
    #
    eq_beams = []
    for beam in flat_beams:
        # from each structure, if there is a beam <beam>, take the list of all
        # those equivalent to it.
        # Example: beam = (1 | 1)
        #   -> take [[( 1 | 1), (-1 |  1), (-1 | -1), (1 | -1)],
        #            [( 1 | 1), ( 1 | -1)]]
        beam_lists = [beams_dict[beam]
                      for beams_dict in all_beams_cp
                      if beam in beams_dict]
        if beam_lists:
            # if there are beams to process in the list, take the set
            # intersection of the elements of the list, i.e., all those in
            # common to all structures
            common_beams = set.intersection(*beam_lists)

            # now "mark as processed" in all the structures the beams coming
            # from the intersection by removing them
            for beams_dict in all_beams_cp:
                for processed_beam in common_beams:
                    beams_dict.pop(processed_beam, None)
            eq_beams.append(common_beams)

    # sort within each equivalence group
    llst = [sorted(list(beams), key=all_leed[0].beamsSortCriterion)
            for beams in eq_beams]
    # and by energy
    sorted_beams = sorted(llst, key=all_leed[0].sortEnergy)

    # and fix the indices, also accounting for extinct beams
    beams_with_indices = []
    for i, beams in enumerate(sorted_beams):
        for beam in beams:
            # find which structures overlap
            overlapping_structs = [s + 1 for s in range(len(all_leed))
                                   if beam in all_beams[s]]
            group_idx = i

            # figure out whether the beam is extinct in all the structures, in
            # which case the index goes negative
            n_extinct = len([1 for s in overlapping_structs
                             if beam in all_extinct[s-1]])
            if n_extinct == len(overlapping_structs) and n_extinct > 0:
                group_idx *= -1
            beams_with_indices.append((beam, group_idx))
    return beams_with_indices


def project_to_first_domain(beam_list, leed_parameters, *other_leed_parameters,
                            domains=None):
    """
    Given a list of beams and a one or more dictionaries of parameters defining
    the LEED geometry, returns a list of beams that 'projects' the input onto
    the first domain. I.e., the output contains only beams that belong to the
    first domain and are equivalent to those given as an input. At present
    [2020-05-04] this function does not account for non-normal beam incidence.

    Parameters
    ----------
    leed_parameters:  dictionary
        See docstring of 'get_equivalent_beams' function for details
    beam_list:  list of tuples, list of 2-element lists, or Nx2 array
        Each element in the list is an (h, k) pair of fractional beam indices
        with type(h) == type(k) == type(fractions.Fraction)

    Returns
    -------
        list of tuples. Each tuple is a beam 'projected' to the first domain.
        Duplicates are removed, and only one of the symmetry-equivalent beams is
        retained (according to the symmetry group of the first domain)
    """
    # process arguments to allow the function to accept the LEED parameters
    # being passed as a single "list" of dictionaries, or as an unpacked "list"
    # of dictionaries
    if isinstance(leed_parameters, dict):
        leed_parameters = (leed_parameters, *other_leed_parameters)

    leed_parameters, (leed, *_) = check_multi_leed_params(leed_parameters)

    ops = leed.reciprocal_lattices['bulk'].group.operations(include_3d=True)
    superlattices = [params['SUPERLATTICE'] for params in leed_parameters]

    def is_in_first_domain(beam):
        """
        Returns True if beam belongs to the first domain, False otherwise
        """
        # only those beams that, processed through one of the superlattice
        # matrices give integer indices belong to the first domain
        for m in superlattices:
            indices = np.dot(beam, m.T)
            if all(i.denominator == 1 for i in indices):
                return True
        return False

    def group_beams(beams_in):
        """
        Given a list of beams in the form ('indices', group_index) generates a
        dictionary with key=group_index and value=[beams with group_index]
        """
        # Using a defaultdict(list) makes each key an empty list when created
        # and allows to extend it with the beams. Each beam is a BeamIndex
        ddict = defaultdict(list)
        for *beams, group in beams_in:
            ddict[group].extend(BeamIndex(beam) for beam in beams)
        return ddict

    def to_inequivalent(beams_in):
        """
        Given a list of beams of the first domain, returns another list
        containing only the inequivalent beams
        """
        # get the list of all equivalent beams for the first domain only
        # It's a list of tuples of the form ('index', group)
        eq_beams = get_equivalent_beams(leed_parameters, domains=0)

        # construct a dictionary of the form beam_group: [*beams].
        group_dict = group_beams(eq_beams)

        # Now create a new dictionary with as many keys as beams, and for each
        # beam, use the first beam in its beam group as value. This allows later
        # to replace each beam in the original list with only one of the beams
        # equivalent to it.
        beams_dict = {beam: beams[0]
                      for group, beams in group_dict.items()
                      for beam in beams}

        # (1) replace each beam in beams_in with the first one of its beam
        #     group, as per the beams_dict dictionary created above
        # (2) remove duplicates by transforming to a set
        # (3) typecast each element (a BeamIndex) to a tuple
        try:
            ineq_beams = set(beams_dict[beam] for beam in beams_in)
        except KeyError:
            # catching the KeyError here should point out to the user that
            # either something is wrong in the list of the experimental beams
            # given, or something is wrong with the SUPERLATTICE matrix
            for beam in beams_in:
                try:
                    beams_dict[beam]
                except KeyError:
                    err = (f"Beam {beam} is incompatible with all the current "
                           "SUPERLATTICE matrices:\n"
                           + "\n".join(str(m) for m in superlattices))

                    # Check if the reason why the beam was not found is that
                    # it would lie outside the LEED screen
                    b = leed.bulk_basis
                    g = np.linalg.norm(np.dot(beam, b)) * 1e10  # 1/m
                    el_m = 9.109e-31    # kg
                    el_q = 1.60218e-19  # C
                    hbar = 1.05457e-34  # J*s

                    # calculate the exit angle
                    s_angle = np.sqrt(hbar**2 * g**2
                                      /(2 * el_m * el_q * leed.max_energy))
                    ang = 2*np.degrees(np.arcsin(s_angle))
                    aperture = leed_parameters[0].get('screenAperture', 110)
                    if ang > aperture:
                        err += ("\nThe beam would need a minimum LEED screen "
                                f"aperture of {ang:.1f} deg at Emax="
                                f"{leed.max_energy} eV, while you are using "
                                f"{aperture} deg. You may have swapped the "
                                "in-plane unit vectors, or you may want to "
                                "set a larger aperture with the "
                                "SCREEN_APERTURE PARAMETER")
                    raise ValueError(err)
        return [tuple(beam) for beam in ineq_beams]

    projected = []
    for op in ops:
        proj = np.dot(beam_list, np.asarray(op).T)
        # keep only the beams that are in the first domain
        first_dom = np.asarray([is_in_first_domain(beam) for beam in proj])
        projected.extend(tuple(beam) for beam in proj[first_dom])
    # remove duplicates
    projected = list(set(projected))

    inequivalent = to_inequivalent(projected)

    # now decide if one should keep all the beams in projected, or only the
    # inequivalent ones -> might be nicer to do with a generator to avoid
    # creating the whole list!
    #
    # We will keep all of them if there is a beam group of the whole system that
    # has a number of beams not integer multiple of the number of beams in the
    # inequivalent list belonging to the same group. This prevents beam
    # averaging errors when, e.g., a spot is composed of 3 beams but only 2
    # would result from reducing to inequivalent
    all_groups = group_beams(
        get_equivalent_beams(leed_parameters, domains=domains)
        )
    for beams in all_groups.values():
        n = len(beams)
        n_ineq = sum(beam in beams for beam in inequivalent)
        if n_ineq == 0:
            continue
        if n % n_ineq != 0:
            return projected
    return inequivalent


def check_type(value, typ):
    """Basic type check for value.

    Parameters
    ----------
    value: object
           The object whose type needs to be checked.
    typ: str
         One of 'str', 'int', 'float, 'number',  'list', 'tuple', 'dict',
         'ndarray', 'arraylike'

    Returns
    -------
    True if type is correct, raises a TypeError otherwise
    """
    is_type = False

    if typ == 'str':
        is_type = isinstance(value, str)
    elif typ == 'int':
        is_type = isinstance(value, int)
    elif typ == 'list':
        is_type = isinstance(value, list)
    elif typ == 'tuple':
        is_type = isinstance(value, tuple)
    elif typ == 'dict':
        is_type = isinstance(value, dict)
    elif typ == 'float':
        is_type = isinstance(value, float)
    elif typ == 'number':
        is_type = check_type(value, 'float') or check_type(value, 'int')
    elif typ == 'ndarray':
        is_type = isinstance(value, np.ndarray)
    elif typ == 'arraylike':
        is_type = hasattr(value, '__iter__') \
                  and not any(isinstance(value, typ) for typ in [str, dict])
    else:
        raise ValueError(f"{typ} is an unknown type")
    return is_type


def check_leed_params(leed_parameters):                                         # TODO: remove: Check is done with LEEDParameters
    """
    Check if all parameters in the leed_parameters dictionary are acceptable
    """
    needed_keys = ('eMax', 'surfBasis', 'SUPERLATTICE',
                   'surfGroup', 'bulkGroup')

    if any(key not in leed_parameters for key in needed_keys):
        missing = ', '.join(key for key in needed_keys
                            if key not in leed_parameters)
        raise NameError("Not enough keys in LEED parameters dictionary. "
                        f"Missing: {missing}")

    # Check types first
    if not check_type(leed_parameters['eMax'], 'number'):
        raise TypeError("Maximum LEED energy should be a number. "
                        f"Found {type(leed_parameters['eMax'])} instead")
    if not check_type(leed_parameters['surfBasis'], 'arraylike'):
        raise TypeError("Surface basis should be array-like. "
                        f"Found {type(leed_parameters['surfBasis'])} instead")
    if not check_type(leed_parameters['SUPERLATTICE'], 'arraylike'):
        raise TypeError("SUPERLATTICE should be array-like. "
                        f"Found {type(leed_parameters['SUPERLATTICE'])} "
                        "instead")
    for group in (leed_parameters['surfGroup'], leed_parameters['bulkGroup']):
        if not (check_type(group, 'str')
                or isinstance(group, PlaneGroup)):
            raise TypeError("Plane group should be a string. "
                            f"found {type(group)} instead")

    aperture = leed_parameters.get('screenAperture', 110.0)
    if not isinstance(aperture, (int, float)):
        raise TypeError("screenAperture should be a floating point number")

    # Then check some requirements on the values
    if leed_parameters['eMax'] < 0:
        raise ValueError("Maximum LEED energy should be positive.")
    if not np.shape(leed_parameters['surfBasis']) == (2, 2):
        raise ValueError("Lattice basis needs to have a (2, 2) shape. "
                         f"Found {np.shape(leed_parameters['surfBasis'])} "
                         "instead")
    if not np.shape(leed_parameters['SUPERLATTICE']) == (2, 2):
        raise ValueError("SUPERLATTICE needs to have a (2, 2) shape. "
                         f"Found {np.shape(leed_parameters['SUPERLATTICE'])} "
                         "instead")
    if aperture < 0 or aperture > 180:
        raise ValueError("screenAperture should be between 0 and 180. "
                         f"Found {aperture} instead.")
    # type and format checking for the plane groups is done in the PlaneGroup
    # instance constructor directly

    return True


def check_multi_leed_params(leed_parameters):                                   # TODO: remove. Checks can be done with LEEDParametersList
    """
    Checks consistency of the LEED parameters passed as a list of dictionaries,
    and returns a consistent version. The returned version has the same eMax and
    screenAperture for all (both the max among the values). Aside from
    checking acceptable values for the parameters (as in check_leed_params, it
    also checks that the bulk lattices are the same for all. Raises errors
    otherwise.

    Parameters
    ----------
    leed_parameters: list of dictionaries


    Returns
    -------
    (consistent_leed_parameters, leed_patterns)

    - consistent_leed_parameters: list of dictionaries with eMax and
                                  ScreenApertures equal for all
    - leed_patterns: list of LEEDPattern
                     patterns generated from consistent_leed_parameters
    """

    # check that each single parameter passed has all the necessary contents,
    # and:
    # - take the largest energy as eMax for all
    # - take the largest screenAperture as the screenAperture for all
    emax = 0
    aperture = 0
    for params in leed_parameters:
        if not check_leed_params(params):
            return None
        emax = max(params['eMax'], emax)
        aperture = max(params.get('screenAperture', 0), aperture)

    leed_patterns = []
    for params in leed_parameters:
        # update eMax and screenAperture
        params['eMax'] = emax
        if aperture > 0:
            params['screenAperture'] = aperture
        leed_patterns.append(LEEDPattern(params))

    # check consistency of bulk
    bulk = leed_patterns[0].reciprocal_lattices['bulk']
    b_ops = set(bulk.group.operations(include_3d=True))
    for leed in leed_patterns[1:]:
        this_bulk = leed.reciprocal_lattices['bulk']
        # bulk lattices should be the same
        if not np.allclose(bulk.basis, this_bulk.basis):
            raise ValueError("Inconsistent bulk bases found in the input "
                             "parameters")
        # bulk groups also should be the same. Not sure how to handle the
        # case in which some of the parameters does not contain the bulk3Dsym
        # but others do. Right now it raises errors.
        if (not bulk.group == this_bulk.group
            or not b_ops == set(this_bulk.group.operations(include_3d=True))):
            raise ValueError("Inconsistent symmetry operations of bulk lattices"
                             " in the input parameters")
    return (leed_parameters, leed_patterns)


def catch_gui_crash():
    """
    Function that allows to catch exceptions that cause the GUI to crash and
    to print them to terminal
    """
    sys._excepthook = sys.excepthook
    def exception_hook(exctype, value, traceback):
        print("########## Caught an exception! ##########")
        # TODO: here one would like to rather open a parent-less
        # QMessageBox reporting the exception. Perhaps even log
        # the event to disk (or have the option in the message).
        sys._excepthook(exctype, value, traceback)
        sys.exit(1)
    sys.excepthook = exception_hook


def check_py_version(version_to_check, check_what='earlier'):
    """Check the interpreter version against a given input.

    Parameters
    ----------
    version_to_check : str
        A string of the form "x.y.z" with major.minor.patch.
        Patch and minor can also be left out.
    check_what : {'earlier', 'later', 'same'}, default='earlier'
        'earlier' checks interpreter_version < version_to_check
        'later' checks interpreter_version > version_to_check
        'same'  checks interpreter_version == version_to_check

    Returns
    -------
    bool
        True if the condition is verified
    """
    if check_what not in ('earlier', 'later', 'same'):
        raise ValueError("check_py_version: Invalid check_what. Should be "
                         "'earlier', 'later', or 'same'.")

    major = sys.version_info.major
    minor = sys.version_info.minor
    patch = sys.version_info.micro

    version_to_check = version_to_check.split('.')
    if len(version_to_check) == 3:
        # major.minor.patch
        py_version = 10000*major + 100*minor + patch
        version_to_check = (10000*int(version_to_check[0])
                            + 100*int(version_to_check[1])
                            + int(version_to_check[2]))
    elif len(version_to_check) == 2:
        # major.minor
        py_version = 100*major + minor
        version_to_check = (100*int(version_to_check[0])
                            + int(version_to_check[1]))
    elif len(version_to_check) == 1:
        # major
        py_version = major
        version_to_check = int(version_to_check[0])
    else:
        raise ValueError("check_py_version: invalid version")

    if check_what == 'earlier':
        return py_version < version_to_check
    if check_what == 'later':
        return py_version > version_to_check
    return py_version == version_to_check


def parallel(v1, v2):
    """
    Check whether vectors v1 and v2 (of arbitrary number of components) are
    parallel to one another

    Parameters
    ----------
    v1, v2 : iterables of real numbers, or None
        If either vector is None or has only zero entries the comparison
        returns False

    Returns
    -------
    v1 parallel v2 : bool

    Raises
    ------
    ValueError if the two vectors have different length
    """
    if v1 is None or v2 is None:
        return False
    if len(v1) != len(v2):
        raise ValueError("parallel: can compare only vectors of equal length")
    # check null entries
    v12 = np.asarray((v1, v2))
    null_entries = np.abs(v12) < 1e-8
    if any(sum(null) == len(v) for null, v in zip(null_entries, v12)):
        # one of the vectors is identically zero
        return False
    if any(null_entries[0] != null_entries[1]):
        # the vectors have zeros at different positions
        return False
    nonnul_1 = v12[0, ~null_entries[0]]
    nonnul_2 = v12[1, ~null_entries[1]]
    ratios = nonnul_1/nonnul_2
    return all(np.abs(ratios - ratios[0]) < 1e-8)


def orientation(vector, zero_pi=True, precision=4):
    """
    Returns the angle in degrees between the vector and the x axis

    Parameters
    ----------
    vector : numpy.ndarray
        vector = [x, y] is the vector for which the orientation is computed
    zero_pi : bool, default=True
        If True, returns angles only in the [0, 180] range, i.e., for
        angle < 0, return 180 - angle. Otherwise angle is in [-180, 180]
    precision : int
        number of decimal places for rounding angles in degrees

    Returns
    -------
    angle : float
        angle in degrees
    """
    angle = round(np.degrees(np.arctan2(vector[1], vector[0])), precision)
    # if angle >= 0 or not zero_pi:
        # return angle
    if not zero_pi:
        return angle
    return angle % 180
