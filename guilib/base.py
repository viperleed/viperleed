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
# from fractions import Fraction
from quicktions import Fraction  # faster version of Fraction (~ factor of 2)
from collections import defaultdict
from warnings import warn as warning   # eventually will replace with logging

import numpy as np

# from viperleed.guilib.leedsim.classes import LEEDPattern
# from guilib.leedsim.classes import LEEDPattern
# from guilib import profile_lines, profile_calls
from viperleed import guilib as gl


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
                           "\n".join(str(m) for m in superlattices))

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
    for group, beams in all_groups.items():
        n = len(beams)
        n_ineq = sum(beam in beams for beam in inequivalent)
        if n_ineq == 0:
            continue
        if n % n_ineq != 0:
            return projected
    return inequivalent


def check_type(value, typ):
    """
    Basic type check for value

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


def check_leed_params(leed_parameters):
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


def check_multi_leed_params(leed_parameters):
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
    - leed_patterns: list of gl.LEEDPattern
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
        leed_patterns.append(gl.LEEDPattern(params))
    
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
        print(exctype, value, traceback)
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


def string_matrix_to_numpy(str_matrix, dtype=float, needs_shape=tuple()):
    """
    Takes a string representing a matrix with float values and parses it into a
    numpy array of the right shape. It can also check if the string represents
    a matrix of a specific shape.

    Parameters
    ----------
    - str_matrix: str
                  String to be parsed to extract a matrix. It needs to be
                  formatted as an array-like input acceptable for numpy
    - dtype: data type; default float
             Data type that the numpy array returned will have
    - needs_shape: array-like; default = tuple()
                   if this keyword argument is given, the function checks that
                   the shape of the parsed matrix matches the input

    Returns
    -------
    np.array or None

    If needs_shape is given, the function returns None in case the shapes do not
    match. Otherwise always returns a numpy array.
    """
    if not check_type(str_matrix, 'str'):
        raise TypeError("Argument 0 of string_matrix should be string. "
                        f"Found {type(str_matrix)} instead")
    if not check_type(needs_shape, 'arraylike'):
        raise TypeError("Argument 1 of string_matrix should be array-like. "
                        f"Found {type(needs_shape)} instead")
    try:
        matrix = np.asarray(ast.literal_eval(str_matrix), dtype=float)
    except SyntaxError as expt:
        raise RuntimeError("Could not extract a matrix "
                           f"from string {str_matrix}. "
                           "Most likely a bracket or a comma is missing")

    # if dtype=int, round first, then typecast
    if dtype == int:
        matrix = np.round(matrix).astype(int)

    if len(needs_shape) > 0 and matrix.shape != tuple(needs_shape):
        return None
    return matrix


def format_floats(format_specs, *numbers):
    """
    Formats floats to have them all aligned to the point, and according to
    format_spec. format_spec is in the form

    [[fill]align][sign][#][0][minimumwidth][.precision][type]

    - fill and align will be used to pad the overall resulting string
    - sign, #, 0 are discarded
    - minimumwidth is the minimum number of characters used to represent the
      integer parts. Defaults to the minimum number necessary to have all
      numbers aligned to the decimal point
    - precision is the number of decimal digits. Defaults to 5.
    """
    int_len = integer_part_length(*numbers)
            
    # standard pattern for format_specs
    pattern = (r"^(?P<align>[<>=^]\d+)?"
               r"[\+\- ]?"
               r"[#]?"
               r"[0]?"
               r"(?P<minwidth>\d*)"
               r"([.](?P<prec>\d+))?"
               r"[bdcoxXneEfFgG\%]?$")
    m = re.match(pattern, format_specs)
    if m is None:
        raise TypeError("Incorrect format specifier for type f.")
    # process the format specifications:
    # number of digits:
    prec = m.group('prec')
    if prec:
        prec = int(prec)
    else:
        prec = 5
    # minimum width of integer part
    # (will left-pad with spaces if needed):
    min_w = m.group('minwidth')
    if min_w:
        min_w = int(min_w)
    else:
        min_w = int_len
    tot_w = max(min_w + prec, int_len + prec) + 1  # +1 for decimal point
    
    format = f">{tot_w}.{prec}f"
    raw = ','.join(f"{float(number):{format}}" for number in numbers)
    align = ''
    if m.group('align'):
        align = f"{m.group('align')}"
    
    return f"{raw:{align}}"


def integer_part_length(*numbers):
    return max(len(f"{int(number)}") for number in numbers)


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
    print(null_entries, ~null_entries, v12, v12[~null_entries])
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


def screen_radius(energy, aperture):                                          # WILL BE MOVED, renamed, and docstring needs work
    """
    Returns the radius of the projection of a LEED screen with a given
    aperture that contains LEED spots with 'energy'

    Parameters
    ----------
    energy : float
        primary energy of the electrons
    aperture: float, default=110.0
              degrees of aperture of the solid angle captured by the
              LEED screen. The current (2020-08-17) default value is taken
              from the dimensions of the screen of the ErLEED optics. The
              MCP SpectaLEED from Omicron appears to have an equivalent
              aperture of ~80°
    """
    electron_mass = 9.109e-31  # kg
    electron_charge = 1.60218e-19   # C
    hbar = 1.05457e-34  # J*s
    # The next one is the prefactor in AA^-1 eV^(-1/2)
    rt2me_hbar = np.sqrt(2*electron_mass*electron_charge) / hbar*1e-10

    return rt2me_hbar*np.sqrt(energy)*np.sin(np.radians(aperture)/2)

################################################################################
#                                   CLASSES                                    #
################################################################################


class BeamIndex(tuple):
    """
    Convenience class to store a 2-element tuple that represents a Miller index
    for a LEED beam. Each index is a Fraction.
    """
    separators = (',', '|')

    def __new__(cls, *indices, denominator=1, from_numerators=False):
        """
        Parameters
        ----------
        indices : str, or iterable of str or number
            indices of the beam. Can be passed as a single argument or as two
            arguments.
            When a single argument, it should be either a string of the form
            'idx1, idx' or 'idx1 | idx2' (spaces don't count), or a 2-element
            iterable with indices.
            When two arguments, both should be either a string or a number,
            i.e., mixed input of the form ('idx1', idx2) is not processed
            correctly
        denominator : int (default=1)
            this is used for speeding up instantiation when indices are given
            as numbers rather than strings, and it is not used at all for string
            inputs. It should be the largest common denominator between the
            indices.
        from_numerators : bool (default=False)
            Use True when passing only the numerators as indices. The
            denominator is taken from the denominator optional parameter. When
            passing numerators only, it is most efficient to give the indices
            as ints rather than floats
        """
        if not isinstance(denominator, int):
            raise TypeError("BeamIndex: denominator must be an int")
        indices = cls.__process_indices(indices)
        if isinstance(indices[0], str):
            beam_indices = [Fraction(index) for index in indices]
        else:
            if from_numerators:
                # the indices passed are already the numerators
                numerators = indices
            else:
                # the indices passed are fractional, get the numerators
                numerators = (round(indices[0]*denominator),
                              round(indices[1]*denominator))
            beam_indices = [Fraction(num, denominator) for num in numerators]
        return super().__new__(cls, beam_indices)

    @staticmethod
    def __process_indices(indices):
        """
        Checks and processes the indices as needed

        Returns
        -------
        indices : iterable, 2 elements

        Raises
        ------
        TypeError
        ValueError
        """
        if len(indices) > 2:
            raise ValueError("BeamIndex accepts at most two indices, "
                             f"{len(indices)} given instead.")
        if len(indices) == 1:
            temp = indices[0]
            if isinstance(temp, str):
                for separator in BeamIndex.separators:
                    indices = temp.split(separator)
                    if len(indices) == 2:
                        # found an acceptable separator
                        break
            else:
                if not hasattr(temp, '__len__'):
                    raise TypeError("BeamIndex: when one argument given, "
                                    "it should be a string or a 2-element "
                                    "array-like.")
                indices = temp
        if len(indices) != 2:
            raise ValueError("BeamIndex: too many/few indices. "
                             "Exactly 2 indices should be given. "
                             f"Found {len(indices)} instead.")
        return indices

    def __str__(self):
        return f"{', '.join(str(index) for index in self)}"

    def __repr__(self):
        return f"BeamIndex({', '.join(str(index) for index in self)})"
    
    def __format__(self, format_specs):
        """
        Customized formatting of BeamIndex. format_specs is a standard format
        specifier in the form:

        [[fill]align][sign][#][0][minimumwidth][.precision][type]

        The following formats apply:
        - type == "s": 
            returns "(num/den|num/den)" where the width of each of the h,k
            fields is dictated by the minimumwidth specifier in format_specs. In
            this case, minimumwidth should be of the form "(w_num,w_den)", so
            that the indices can be aligned at the slash. If this is omitted,
            the width of both fields is equal, and equal to the longest among
            the two.
            If h or k are integers, their "/den" is omitted, and replaced with
            white spaces if the other is fractional.
            Negative signs are printed, positive ones are replaced with spaces
            if the other index is negative
        - type == "f":
            returns f"({float(h)},{float(k)})". If .precision is not given,
            uses 5 digits. If minimumwidth is given, it is treated as the
            minimum width of the integer part only. The two indices are aligned
            on the decimal point.
        - all others return f"({str(self)}:{format_specs}})"
        """
        if format_specs == '':
            return str(self)
        if format_specs[-1] not in ('s', 'f'):
            # basic format
            return f"({str(self):{format_specs}})"
        if format_specs[-1] == 'f':
            return f"({format_floats(format_specs, *self)})"

        # Case 's':
        num_min_len, den_min_len = self.get_format_lengths('s')

        # now search the specifier for something like "(\d+,\d+)" to be
        # interpreted as the minimum widths of the two fields, to update
        # the minimum lengths of the fields
        m = re.search(r"((?P<num>\d+),(?P<den>\d+))", format_specs)
        if m is not None:
            num_min_len = max(num_min_len, int(m.group('num')))
            den_min_len = max(den_min_len, int(m.group('den')))
            format_specs = format_specs.replace(
                f"({m.group('num')},{m.group('den')})", '')

        raws = ['', '']
        for i, hk in enumerate(self):
            # numerator is right-justified in its field
            raws[i] = f"{hk.numerator:>{num_min_len}}"

            # denominator a bit more complicated, as it depends on whether
            # it is == 1 or not
            if hk.denominator == 1:
                n_white = den_min_len
                n_white += 1 if den_min_len else 0  # slash if needed
                raws[i] += ' '*(n_white)
            else:
                raws[i] += f'/{hk.denominator:<{den_min_len}}'
        return f"{f'({raws[0]}|{raws[1]})':{format_specs}}"

    @property
    def numerators(self):
        return tuple(index.numerator for index in self)
    
    def get_format_lengths(self, str_or_float):
        """
        Returns the minimum number of characters that can represent the
        BeamIndex as a fraction or as a float. The two cases behave differently:
        - str_or_float = 's':
            returns a 2-tuple with the minimum lengths of numerator and
            denominator, such that both h and k are aligned at the slash (if
            fractional) or at the lowest significant digit of the numerator
        - str_or_float = 'f':
            returns a 1-tuple with the minimum length of the string
            representation of the integer parts that can represent both
        """
        if str_or_float not in ('s', 'f'):
            raise ValueError("Invalid format specifier. Should be 's' or 'f'.")
        if str_or_float == 'f':
            return (integer_part_length(*self),)
        # case 's'
        num_min_len = max(len(str(hk.numerator)) for hk in self)
        dens = [hk.denominator for hk in self]
        if all(den == 1 for den in dens):
            den_min_len = 0
        else:
            den_min_len = max(len(str(den)) for den in dens)
        return num_min_len, den_min_len


class PlaneGroup():
    """
    PlaneGroup(group='p1')

    Class representing a planar 2D group

    Parameters
    ----------
    group: str, default 'p1'
           Hermann-Mauguin notation for the 2D plane group (with some
           extensions). Acceptable values:
           'p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]', 'pg[0 1]', 'cm[1 0]',
           'cm[0 1]', 'cm[1 1]', 'cm[1 -1]', 'cm[1 2]', 'cm[2 1]', 'rcm[1 0]',
           'rcm[0 1]', 'pmm', 'pmg[1 0]', 'pmg[0 1]', 'pgg', 'cmm', 'rcmm',
           'p4', 'p4m', 'p4g', 'p3', 'p3m1', 'p31m', 'p6', 'p6m'

           Please refer to viperleed/docs/plane_groups.pdf for more information

    Attributes
    ----------
    group: str
           same as input parameter

    Methods
    -------
    operations() Returns tuple with group operations

    get_subgroups() Returns set of strings with the subgroups of group

    The point group operations are represented by 2x2 matrices (np ndarrays),
    and correspond to applying the operation in 'fractional' coordinates.
    The naming follows the following convention:
    - identity: E
    - rotations: Cn, Cmn
              * Cn -> 2pi/n counterclockwise rotation
              * Cmn -> -2pi/n counterclockwise rotation
    - mirrors: Mx, My, M45, Mm45 or Mij
              * Orthogonal bases:
                Mx: mirror across line through first unit vector
                My: mirror across line through second unit vector
                M45: mirror across line at +45deg with respect to first
                     unit vector
                Mm45: mirror across line at -45deg with respect to first
                      unit vector
              * Other bases:
                Mij: mirror across a line through vector v = i*a + j*b,
                     where a and b are the basis unit vectors
    """
    # These two are good for all cells
    E = (1, 0), (0, 1)
    C2 = (-1, 0), (0, -1)    # = -E

    # These two are good for rectangular and square cells
    Mx = (1, 0), (0, -1)
    My = (-1, 0), (0, 1)     # = -Mx

    # These are good for square cells
    C4 = (0, -1), (1, 0)
    Cm4 = (0, 1), (-1, 0)    # = -C4
    M45 = (0, 1), (1, 0)
    Mm45 = (0, -1), (-1, 0)  # = -M45

    # These are good for rhombic and hex (both obtuse)
    M11 = (0, 1), (1, 0)     # = M45
    M1m1 = (0, -1), (-1, 0)  # = -M11 = Mm45
    M01 = (-1, -1), (0, 1)
    M10 = (1, 0), (-1, -1)

    # And these are good for hex only (obtuse)
    C6 = (1, 1), (-1, 0)
    Cm6 = (0, -1), (1, 1)
    C3 = (0, 1), (-1, -1)
    Cm3 = (-1, -1), (1, 0)
    M21 = (1, 1), (0, -1)
    M12 = (-1, 0), (1, 1)

    allGroups = {'p1': (E,),
                 'p2': (E, C2),
                 'pm[1 0]': (E, Mx), 'pm[0 1]': (E, My),
                 'pg[1 0]': (E, Mx), 'pg[0 1]': (E, My),
                 'cm[1 0]': (E, M10), 'cm[0 1]': (E, M01),
                 'cm[1 1]': (E, M11), 'cm[1 -1]': (E, M1m1),
                 'cm[1 2]': (E, M12), 'cm[2 1]': (E, M21),
                 'rcm[1 0]': (E, Mx), 'rcm[0 1]': (E, My),
                 'pmm': (E, Mx, My, C2),
                 'pmg[1 0]': (E, Mx, My, C2), 'pmg[0 1]': (E, Mx, My, C2),
                 'pgg': (E, C2, Mx, My),
                 'cmm': (E, C2, M11, M1m1), 'rcmm': (E, C2, Mx, My),
                 'p4': (E, C2, C4, Cm4),
                 'p4m': (E, Mx, My, M45, Mm45, C2, C4, Cm4),
                 'p4g': (E, Mx, My, M45, Mm45, C2, C4, Cm4),
                 'p3': (E, C3, Cm3),
                 'p3m1': (E, M12, M21, M1m1, C3, Cm3),
                 'p31m': (E, M10, M11, M01, C3, Cm3),
                 'p6': (E, C6, Cm6, C3, C2, Cm3),
                 'p6m': (E, M10, M01, M12, M21, M11, M1m1,
                         C6, Cm6, C3, C2, Cm3)
                 }
    # The following two dictionaries are used in screws_glides to convert
    # 1) rotation orders of screws into a tuple of the corresponding matrices
    screw_ops = {'2': (C2,),
                 '3': (C3, Cm3),
                 '4': (C4, Cm4),    # Probably need also to add C2
                 '6': (C6, Cm6)}    # Probably also need C3, Cm3 and C2
    # 2) direction contained in glide planes into the corresponding matrices
    glide_ops = {'[1,0]': M10,
                 '[0,1]': M01,
                 'x': Mx,
                 'y': My,
                 '[1,1]': M11,
                 '[1,-1]': M1m1,
                 '[1,2]': M12,
                 '[2,1]': M21}

    groupsForShape = {
        'Oblique': ('p1', 'p2'),
        'Rectangular': ('p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]',
                        'pg[0 1]', 'rcm[1 0]', 'rcm[0 1]', 'pmm', 'pmg[1 0]',
                        'pmg[0 1]', 'pgg', 'rcmm'),
        'Square': ('p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]', 'pg[0 1]',
                   'cm[1 1]', 'cm[1 -1]', 'pmm', 'pmg[1 0]', 'pmg[0 1]',
                   'pgg', 'cmm', 'p4', 'p4m', 'p4g'),
        'Rhombic': ('p1', 'p2', 'cm[1 1]', 'cm[1 -1]', 'cmm'),
        'Hexagonal': ('p1', 'p2', 'cm[0 1]', 'cm[1 0]', 'cm[1 1]',
                      'cm[1 -1]', 'cm[1 2]', 'cm[2 1]', 'cmm', 'p3', 'p3m1',
                      'p31m', 'p6', 'p6m')
        }

    subgroups = {
        'p1': {'p1'},
        'p2': {'p1', 'p2'},
        'pm[1 0]': {'p1', 'pm[1 0]'},
        'pm[0 1]': {'p1', 'pm[0 1]'},
        'pg[1 0]': {'p1', 'pg[1 0]'},
        'pg[0 1]': {'p1', 'pg[0 1]'},
        'cm[1 0]': {'p1', 'cm[1 0]'},
        'cm[0 1]': {'p1', 'cm[0 1]'},
        'cm[1 1]': {'p1', 'cm[1 1]'},
        'cm[1 -1]': {'p1', 'cm[1 -1]'},
        'cm[1 2]': {'p1', 'cm[1 2]'},
        'cm[2 1]': {'p1', 'cm[2 1]'},
        'rcm[1 0]': {'p1', 'pm[1 0]', 'pg[1 0]', 'rcm[1 0]'},
        'rcm[0 1]': {'p1', 'pm[0 1]', 'pg[0 1]', 'rcm[0 1]'},
        'pmm': {'p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pmm'},
        'pmg[1 0]': {'p1', 'p2', 'pg[1 0]', 'pm[0 1]', 'pmg[1 0]'},
        'pmg[0 1]': {'p1', 'p2', 'pg[0 1]', 'pm[1 0]', 'pmg[0 1]'},
        'pgg': {'p1', 'p2', 'pg[1 0]', 'pg[0 1]', 'pgg'},
        'cmm': {'p1', 'p2', 'cm[1 1]', 'cm[1 -1]', 'cmm'},
        'rcmm': {'p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]', 'pg[0 1]',
                 'rcm[1 0]', 'rcm[0 1]', 'pmm', 'pmg[1 0]', 'pmg[0 1]', 'pgg',
                 'rcmm'},
        'p4': {'p1', 'p2', 'p4'},
        'p4m': {'p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'cm[1 1]', 'cm[1 -1]', 'pmm',
                'cmm', 'p4', 'p4m'},
        'p4g': {'p1', 'p2', 'pg[1 0]', 'pg[0 1]', 'cm[1 1]', 'cm[1 -1]', 'pgg',
                'cmm', 'p4', 'p4g'},
        'p3': {'p1', 'p3'},
        'p3m1': {'p1', 'cm[1 -1]', 'cm[2 1]', 'cm[1 2]', 'p3', 'p3m1'},
        'p31m': {'p1', 'cm[1 0]', 'cm[0 1]', 'cm[1 1]', 'p3', 'p31m'},
        'p6': {'p1', 'p2', 'p3', 'p6'},
        'p6m': {'p1', 'p2', 'cm[1 0]', 'cm[0 1]', 'cm[1 1]', 'cm[1 -1]',
                'cm[2 1]', 'cm[1 2]', 'p3', 'p3m1', 'p31m', 'p6', 'p6m'}
        }

    allOps = (E, C2, C4, Cm4, C3, Cm3, C6, Cm6, Mx, My, M45, Mm45, M21, M12,
              M01, M10, M11, M1m1)  # this should not be needed anymore

    def __init__(self, group='p1'):
        if isinstance(group, PlaneGroup):
            bulk_3d = group.screws_glides
            group = group.group
        else:
            bulk_3d = tuple()
        group = self.check_group_name(group)
        self.group = group
        
        # The next one will be a tuple of the 2x2 matrices representing the
        # isomorphism part of screws and glide planes perpendicular to the
        # surface
        self.__ops_3d = bulk_3d

    def __repr__(self):
        """
        Representation of PlaneGroup
        """
        return f"viperleed.PlaneGroup(group={self.group!r})"

    def __str__(self):
        """
        string representation of PlaneGroup
        """
        return self.group

    def __eq__(self, other):
        """
        Equality method for PlaneGroup instances.
        
        Most likely one would like to also check the set of symmetry operations,
        as this would allow to include checks of the screws and glides as well.

        Notice that Python 3 handles correctly (and in a faster way) cases in
        which the not-equal dunder is not reimplemented
        """
        if not isinstance(other, PlaneGroup):
            # Can't compare instances of other classes
            return NotImplemented
        return self.group == other.group
    
    def same_operations(self, other, include_3d=False):
        """
        Returns whether self has the same operations as another PlaneGroup
        instance. The comparison can include or not (depending on the value of
        include_3d) also the 3D bulk operations. This does not check whether
        the Hermann-Mauguin names are the same. Use self == other to test that.
        """
        if not isinstance(other, PlaneGroup):
            return NotImplemented
        self_ops = set(self.operations(include_3d))
        other_ops = set(other.operations(include_3d))
        return self_ops == other_ops

    def check_group_name(self, group):
        """
        Fixes spaces in the group name and checks that the fixed name is an
        acceptable plane group
        """
        if not isinstance(group, str):
            raise TypeError("'group' should be a string. "
                            f"Found {type(group)} instead")
        group_re = re.compile(
            r"""
            (?P<hermann>[\w]+)      # hermann mauguin
            (?:[\[]\s*              # optional direction opening bracket
            (?P<dir1> [-]*[\d]+)    # first direction
            [\s]*                   # optional space
            (?P<dir2> [-]*[\d]+)    # second direction
            [\]])*                  # optional direction closing bracket
            """,
            re.VERBOSE)
        
        group_match = group_re.match(group)
        if group_match is None:
            raise ValueError(f"{group} is not an acceptable plane group.")
        group = group_match.group('hermann')
        
        if group_match.group('dir1') is not None:
            group = (f"{group}[{group_match.group('dir1')} "
                     f"{group_match.group('dir2')}]")

        if group not in self.allGroups.keys():
            raise ValueError(f"{group} is not an acceptable plane group.")
        
        return group
    
    @property
    def screws_glides(self):
        """
        Returns a tuple containing 2x2 tuples of integers representing the
        isomorphic part of screw axes and glide planes perpendicular to the
        surface
        """
        return self.__ops_3d
    
    @screws_glides.setter
    def screws_glides(self, input):
        """
        Set the isomorphic part of screws/glides.
        Parameters
        ----------
        input, is a 1- or 2-tuple where
        first item: str or array-like 
                 - a string of the following forms (with or without spaces):
                   "r(#, #, ...), m([#, #], [#, #], ...)"
                   "m([#, #], [#, #], ...), r(#, #, ...)"
                   "r(#, #, ...)"
                   "m([#, #], [#, #], ...)"
                   "None"
                   For screws, the list in "r()" provides the order of rotations
                   For glides, an entry "[i,j]" means the glide plane leaves
                   the in-plane direction i*a + j*b unmodified.
                 - an array-like, containing a list of 2x2 matrices of integers
        second item:
                 one of the lattice types: 'Oblique', 'Rectangular', 'Square',
                 'Hexagonal', or 'Rhombic'. This parameter is required only if
                 the first argument is a string containing definitions of glide
                 planes.
        """
        # The next if/else handles the input parameter so that, if only a
        # 1-tuple or a single string is given, the shape of the unit cell is
        # defaulted to None.
        # NB: I may change the behavior later, and rather make this a
        #     set_3d_ops method, while incorporating the getter into
        #     the operations method below
        if not input:
            self.__ops_3d = tuple()
            return
        if isinstance(input, (tuple, list, np.array)):
            if len(input) > 2 or len(input) == 0:
                raise ValueError("PlaneGroup.screws_glides: requires at most "
                                 f"a 2-tuple. {len(input)} items found.")
            elif len(input) == 2:
                shape = input[1]
                input = input[0]
            else:
                input = input[0]
                shape = None
        else:
            shape = None
        
        if not input:
            self.__ops_3d = tuple()
            return

        if not isinstance(input, (str, np.ndarray, tuple, list)):
            raise ValueError("PlaneGroup.screws_glides: Invalid input. "
                             "Need either a string or an array-like. "
                             f"Found {type(input)} instead.")

        if isinstance(input, str):
            if input.lower() == "none":
                self.__ops_3d = tuple()
                return

            ops = []  # this will contain the 2x2 tuples of the new operations

            # search the following patterns
            screw_re = r"[rR][\(](?P<screws>[\d\,\s]+)[\)]"
            glide_re = r"[mM][\(](?P<glides>[\d\[\]\,\-\s]+)[\)]"
            
            found_screws = re.search(screw_re, input)
            found_glides = re.search(glide_re, input)
            if not (found_screws or found_glides):
                raise ValueError("PlaneGroup.screws_glides: Invalid input.")
            if found_screws:  # found some screws
                screws = found_screws.group('screws').replace(' ',
                                                              '').split(',')
                if any(s not in self.screw_ops.keys() for s in screws):
                    raise ValueError("PlaneGroup.screws_glides: Invalid "
                                     "rotation order in the input. Only 2-, "
                                     "3-, 4-, and 6-fold orders allowed.")
                [ops.extend(self.screw_ops[screw]) for screw in screws]
            if found_glides:  # found some glide planes
                if shape is None:
                    raise ValueError("PlaneGroup.screws_glides: cell shape is "
                                     "required when glide planes are given.")
                # parse glide planes by removing spaces and splitting on commas
                glides = found_glides.group('glides').replace(' ','').split(',')

                # now g should contain an even number of elements:
                # each odd element is of the form "[#", each even "#]"
                for g_odd, g_even in zip(glides[::2], glides[1::2]):
                    if g_odd[0] != "[" or g_even[-1] != "]":
                        raise ValueError("PlaneGroup.screws_glides: some "
                                         "glide plane directions are not in "
                                         "the form '[i, j]'")
                    key = f"{g_odd},{g_even}"
                    # the only keys that require special attention are "[1,0]"
                    # and "[0,1]", since the matrices depend on the shape of
                    # the cell
                    if key == '[1,0]' and shape in ('Square', 'Rectangular'):
                        key = 'x'
                    if key == '[0,1]' and shape in ('Square', 'Rectangular'):
                        key = 'y'
                    try:
                        ops.append(self.glide_ops[key])
                    except KeyError:
                        raise ValueError(
                            f"PlaneGroup.screws_glides: invalid direction {key}"
                            " for glide plane. The only directions allowed "
                            f"are: {self.glide_ops.keys() - ['x', 'y']}")
            self.__ops_3d = tuple(ops)
            return

        # Otherwise, the input is an array-like, which should correspond to a
        # 1D list of 2x2 matrices with integer values

        if len(np.shape(input)) != 3 or np.shape(input)[1:] != (2, 2):
            raise ValueError("PlaneGroup.screws_glides: an array-like input "
                             "should be a 1D 'list' of 2x2 'matrices'. Found "
                             f"incompatible shape {np.shape(input)}.")

        if any(np.abs(mij % 1) > 1e-4 for mij in np.ravel(input)):
            raise ValueError("PlaneGroup.screws_glides: an array-like input "
                             "should contain only integer-valued matrices")

        self.__ops_3d = tuple(  # * make the whole list a tuple
            tuple(              # * make each array of the input into a tuple
                map(tuple,      # * map each line of each matrix to a tuple
                    np.array(mi).round().astype(int))  # after rounding to int
                )
            for mi in input
            )

    def operations(self, include_3d=False):
        """
        Returns a tuple of 2x2 tuples representing the operations of the point
        group associated with the 2D plane group
        """
        ops = list(self.allGroups[self.group])
        if include_3d:
            ops.extend(self.screws_glides)
        return tuple(ops)

    def get_subgroups(self):
        """
        Returns the 2D subgroups of PlaneGroup as a set of strings
        """
        return self.subgroups[self.group]
    
    def is_valid_group(self, group, cell_shape):
        """
        Checks whether group is a valid group for a given cell_shape. No type
        checking is done on group other than string and PlaneGroup, under the
        assumption that group will be used to create a PlaneGroup itself, which
        does the type checking in the constructor
        """
        if cell_shape not in self.groupsForShape:
            raise ValueError(f"PlaneGroup: unknown lattice shape {cell_shape}")
        valid_group = True
        if isinstance(group, str):
            group = self.check_group_name(group)
            valid_group = group in self.groupsForShape[cell_shape]
        if isinstance(group, PlaneGroup):
            valid_group = group.group in self.groupsForShape[cell_shape]
        return valid_group
    
    def transform(self, transform, inverse=None, include_3d=False):
        """
        Returns a tuple of the group operations 'projected' to a new coordinate
        system, whose coordinates are expressed by the 2x2 array-like
        transformation matrix given in the "transform" parameter. No assumption
        is made on the coordinate transform (i.e., the operation matrices
        returned may have non-integer values)

        Returns
        -------
        tuple of numpy.ndarrays
        """
        if np.shape(transform) != (2, 2):
            raise ValueError("PlaneGroup.transform requires a 2x2 array-like as"
                             " coordinate transform matrix. "
                             f" Found shape {np.shape(transform)} instead.")
        if inverse is None:
            inverse = np.linalg.inv(transform)
        elif np.shape(inverse) != (2, 2):
            raise ValueError("PlaneGroup.transform requires a 2x2 array-like as"
                             " the inverse of the coordinate transform. "
                             f" Found shape {np.shape(transform)} instead.")
        elif not np.allclose(np.dot(transform, inverse), ((1, 0), (0, 1))):
            raise ValueError("PlaneGroup.transform transformation matrix and "
                             "inverse are inconsistent.")
        return tuple(np.linalg.multi_dot((transform,
                                          op,
                                          inverse))
                     for op in self.operations(include_3d))


class Lattice():
    """
    Lattice(basis, space='real', group='p1', limit=1)

    Class for a 2D lattice. Includes basis and its shape, plane group, whether
    it is a real- or reciprocal-space lattice, and an array of lattice points.

    Parameters
    ----------
    basis : 2x2 array-like
        basis vectors a and b of the lattice, with a = basis[0], b=basis[1]
        the units are assumed to be Angstrom for real space lattices, and
        2*pi/Angstrom for reciprocal-space lattices
    space : str, default='real'
        accepts 'real' or 'reciprocal' for real and reciprocal space
        lattices, respectively
    group : str, default='p1'
        plane group in Hermann-Mauguin notation. See also the PlaneGroup
        class for acceptable values.
    limit : int, default=1
        radius used to limit the number of lattice points generated. Only
        lattice points closer to the origin than limit will be produced.

    Attributes
    ----------
    basis : 2x2 numpy.ndarray
        basis, same convention as with class instantiation

    space : str, read-only
        'real' or 'reciprocal'

    cell_shape : str, read-only
        shape of the unit cell. Can be 'Oblique', 'Rectangular',
        'Square', 'Rhombic', or 'Hexagonal'

    group : PlaneGroup
        a PlaneGroup instance representing the plane group of the lattice.
        Can be set with Lattice.group = str or PlaneGroup, or Lattice.group(str)

    lattice : numpy.ndarray, shape == (..., 2)
        array of lattice points (x, y) [for real-space] or (gx, gy) [for
        reciprocal space]

    hk : numpy.ndarray, shape == (..., 2)
        array of (h, k) indices that generate the points in lattice, i.e.,
        lattice[i] = hk[i, 0]*basis[0] + hk[i, 1]*basis[1]

    n_beams : int, read-only
        Return number of lattice points

    real_basis : 2x2 numpy.ndarray, read only
        returns a copy of the real-space basis, independently of whether
        the lattice is real or reciprocal. Use Lattice.basis for changing
        the basis

    reciprocal_basis : 2x2 numpy.ndarray, read-only
        returns the reciprocal of Lattice.basis. Use Lattice.basis for changing
        the basis

    lattice_parameters : 3-tuple of floats, read-only
        return lattice parameters as (length_a, length_b, angle), where angle
        is the angle between the basis vectors in degrees. Use Lattice.basis for
        changing the basis

    special_directions : list of (2,) numpy.ndarrays
        Returns vectors along the directions of the mirror planes, in the
        same coordinate system as Lattice.basis

    Public methods
    --------------
    get_rotated_lattice(angle)
        returns a copy of Lattice.lattice rotated by angle (degrees), without
        actually modifying Lattice

    get_rotated_basis(angle)
        returns a copy of Lattice.basis rotated by angle, without actually 
        modifying Lattice

    high_symm_transform()
        return matrix transform that gives highest symmetry basis
        (with the same lattice). Does not transform the Lattice though.
        Use Lattice.make_high_symmetry() to make high symmetry, as this also
        returns the same transform

    make_high_symmetry()
        transform Lattice so the basis has the highest possible symmetry,
        returns the transformation used

    transform(transform, as_copy=True)
        transforms the Lattice by applying transform TO THE LEFT of basis,
        i.e., the transform matrix passed is meant to act on ROW VECTORS.
        Can either modify Lattice itself (as_copy=False) or a copy of it.
        In any case, it returns the transformed Lattice instance

    equivalent_points(special_direction=None, superlattice=None)
        Returns a dictionary with (h, k): star of (h, k) entries.
        If special_direction is given, only the mirror plane containing it is
        used to compute the star of each lattice point. Use superlattice to
        return indices expressed with respect to a different basis
    
    Private methods
    ---------------
    __get_cell_shape()
        finds the shape of lattice unit cell. Use Lattice.cell_shape, unless you
        suspect that the shape is not up to date

    __generate_lattice()
        generates the lattice points
    """

    def __init__(self, basis, space='real', group='p1', limit=1):
        if not check_type(basis, 'arraylike'):
            raise TypeError(f"Lattice basis should be array-like, "
                            f"not {type(basis)}")
        if not np.shape(basis) == (2, 2):
            raise ValueError("Lattice basis needs to have a (2, 2) shape. "
                             f"Found {np.shape(basis)} instead")
        if not check_type(space, 'str'):
            raise TypeError(f"Lattice space should be a string, "
                            f"not {type(space)}")
        if space not in ('real', 'reciprocal'):
            raise ValueError(f"Lattice space {space} unknown")
        if not check_type(limit, 'number'):
            raise TypeError(f"Lattice limit should be a scalar, "
                            f"not {type(limit)}")

        self._basis = np.asarray(basis)
        self._space = space
        self._shape = self.__get_cell_shape()  # __get_cell_shape
        
        # check if the plane group given is consistent with the cell shape
        if not PlaneGroup().is_valid_group(group, self.cell_shape):
            raise ValueError(f"Lattice: invalid group {group} for lattice "
                             f"shape {self.cell_shape}")
        self._group = PlaneGroup(group)
        self._limit = limit
        self.lattice, self.hk = self.__generate_lattice()
    
    def __repr__(self):
        txt = (f"{self.cell_shape} "
               + f"viperleed.Lattice({self.basis}, ".replace('\n', '')
               + f"space={self.space}, group={self.group}, "
               + f"limit={self._limit})")
        return txt
    
    @property
    def basis(self):
        """
        Returns the lattice basis as a 2x2 numpy.ndarray, with a = basis[0] and
        b = basis[1]
        """
        return self._basis

    @basis.setter
    def basis(self, basis):
        """
        Sets the lattice basis to basis and updates the other attributes
        
        Parameters
        ----------
        basis : 2x2 array-like
        """
        if not check_type(basis, 'arraylike'):
            raise TypeError("basis must be array-like. "
                            f"Found {type(basis)} instead.")
        if not np.shape(basis) == (2, 2):
            raise ValueError("Lattice basis needs to have a (2, 2) shape. "
                             f"Found {np.shape(basis)} instead.")
        self._basis = np.array(basis)
        self._shape = self.__get_cell_shape()
        self.lattice, self.hk = self.__generate_lattice()

    @property
    def lattice_parameters(self):
        """
        Returns the length of the lattice vectors and the angle alpha between
        them as a tuple (norm_a, norm_b, alpha)
        """
        norm_a, norm_b = np.linalg.norm(self.basis, axis=-1)
        alpha = np.arccos(np.dot(self.basis[0], self.basis[1])/(norm_a*norm_b))
        return norm_a, norm_b, np.degrees(alpha)

    @property
    def reciprocal_basis(self):
        """
        Returns the reciprocal of the lattice basis. Notice that if
        self.space == 'real', this returns the reciprocal space basis.
        If self.space == 'reciprocal', it returns the reciprocal of
        the reciprocal, i.e., the real-lattice basis.

        This is different from the return value of real_basis
        """        
        # Working out the cross products
        #   a* = 2pi/area [b x z][0:2] = 2pi/area [b12, -b11]
        #   b* = 2pi/area [z x a][0:2] = 2pi/area [-a12, a11]
        # where
        #   area = det(basis).
        # Thus
        #  | a* |                    | b12  -a12| T
        #  |    | = 2pi/det(basis) * |          |   = 2pi * basis^(-T)
        #  | b* |                    |-b11   a11|
        
        return 2*np.pi*np.linalg.inv(self.basis).T

    @property
    def real_basis(self):
        """
        Always returns the real-space basis of the lattice as a 2x2 np.array
        """
        if self.space == 'reciprocal':
            return self.reciprocal_basis
        # Not sure why this was like this, but seems useless. Probably was
        # just to yield a copy.
        # return np.dot(self.basis, [[1, 0], [0, 1]])
        return self.basis.copy()

    @property
    def cell_shape(self):
        """
        Returns the shape of the unit cell as a string. The value returned is
        'Oblique', 'Rectangular', 'Square', 'Rhombic', or 'Hexagonal'.
        """
        return self._shape

    @property
    def space(self):
        """
        Returns the space of the lattice as a string. The value returned is
        either 'real' or 'reciprocal'. This property cannot be changed after
        initialization.
        """
        return self._space

    @property
    def group(self):
        """
        Returns the plane group of the lattice as a PlaneGroup instance
        """
        return self._group

    @group.setter
    def group(self, group):
        """
        Set the plane group to a PlaneGroup instance with Hermann-Mauguin
        symbol equal to group. See PlaneGroup for a list of acceptable input
        parameters.
        """
        self._group = PlaneGroup(group)

    @property
    def n_beams(self):
        """
        Returns
        -------
        int
        the number of lattice points (for real-space) or LEED beams (for
        reciprocal-space)
        """
        return len(self.hk)
    
    @property
    def special_directions(self):
        """
        Returns vectors along the directions of the mirror planes, in the
        same coordinate system as self.real_basis. The list returned is as long
        as there are operations in self.group, with None instead of a direction
        vector in correspondence of rotations

        Returns
        -------
        list of numpy.ndarrays
        """
        # Work on the real-basis coordinates, as directions are the
        # same in real and reciprocal space. Using a reciprocal-space basis
        # would require to transpose-invert the symmetry operations.
        # Also notice that we're not including the 3d operations, as this
        # makes sense only for a 'surface' lattice, and the role of 3d
        # operations is merely that of determining how many domains may
        # occur in LEED << IS THIS CORRECT??
        transform = np.linalg.inv(self.real_basis)
        ops = self.group.transform(transform)

        # The special directions are those parallel to the eigenvectors
        # of the mirrors with eigenvalue == 1 (use 1e-4 tolerance).
        # In fact, this is the direction that is unchanged when applying 
        # the mirror operation. Rotations do not give any special direction.
        directions = []
        for op in ops:
            if np.linalg.det(op) > 0:  # Rotation
                directions.append(None)
                continue
            eigval, eigvecs = np.linalg.eig(op)
            # eigvecs has the eigenvectors as columns, ordered in the same
            # way as the eigenvalues in eigvals. Notice the use of .extend()
            # rather than .append(). This is because eigvecs is a matrix,
            # and its 1-element slice is a (1, 2) array. Using .extend() adds
            # a (2,)-shaped array.
            directions.extend(eigvecs.T[np.abs(eigval - 1) < 1e-4])
        return directions

    def __get_cell_shape(self):
        """
        Determine the shape of the lattice ('Oblique', 'Rhombic', 'Hexagonal',
        'Rectangular', or 'Square'), also changing the real-space basis in an
        handedness-conserving fashion so that the angle between the basis
        vectors in real space is obtuse (for hexagonal and rhombic only)
        
        Returns
        -------
        string
        """

        basis = self._basis

        # Relative tolerance factor within which things are assumed to be equal
        eps = 1e-3

        # cosine of angle between vectors
        cosine = np.dot(basis[0], basis[1])/(np.linalg.norm(basis[0])
                                             * np.linalg.norm(basis[1]))

        # Mismatch between the length of the two vectors
        delta = np.linalg.norm(basis[0])/np.linalg.norm(basis[1]) - 1
        if abs(cosine) < eps:  # angle is 90°
            if np.abs(delta) < eps:
                return 'Square'
            return 'Rectangular'
        if np.abs(delta) < eps:  # rhombic or hex
            if self.space == 'real':
                if cosine > eps:  # angle is acute -> make it obtuse
                    print('Warning: the input might be corrupted,'
                          'since the real-space basis is not obtuse!')
                    transform = (0, -1), (1, 0)  # this keeps the handedness
                else:
                    transform = (1, 0), (0, 1)
                self._basis = np.dot(transform, self._basis)
                basis = self._basis
                cosine = np.dot(basis[0],
                                basis[1])/(np.linalg.norm(basis[0])
                                           * np.linalg.norm(basis[1]))
            else:
                if cosine > eps:
                    # The reciprocal lattice will be acute, since
                    # the real one is obtuse.
                    cosine *= -1

            if abs(cosine + 0.5) < eps:  # angle is 120deg
                return 'Hexagonal'
            return 'Rhombic'
        return 'Oblique'

    def __generate_lattice(self):
        """
        Generates a list of lattice points given a basis.

        Parameters
        ----------
        limit : scalar
            Determines which portion of the lattice is generated.
            For space == 'real', the lattice is generated up to a radius
            of 1.5*limit.
                One should then plot from -limit to +limit. This should
                cover any post-rotation of the lattice that the user
                might later request.
            For space == 'reciprocal', the lattice is generated up to a
            radius of limit

        basis : 2x2 array-like, default = self.basis
            Contains the unit vectors as basis[0] and basis[1]

        space : string, {'real', 'reciprocal'}, default = self.space
            Determines whether a real-space or a reciprocal-space
            lattice is generated.
            This only affects the behavior of limit (see above)

        Returns
        -------
        Tuple (lat, hk)

        lat : 1D numpy.ndarray
            Lattice points generated

        hk : 1D numpy.ndarray, hk.shape == lat.shape
            Integer indices of the lattice points generated.
            The lattice points are lat = h*basis[0] + k*basis[1].
            hk[i] is a 2-element array of indices of lattice point i
        """

        limit = self._limit
        basis = self.basis
        space = self.space

        if space == 'real':
            limit *= 1.5

        # get limit for the loops that follow
        shortest = min(*np.linalg.norm(basis, axis=-1),
                       np.linalg.norm(basis[0] + basis[1])/2,
                       np.linalg.norm(basis[0] - basis[1])/2)
        h_max = int(np.ceil(limit/shortest))

        # create grid of indices
                                                                                # TODO: use BeamIndex here, rather than a bare tuple. Even better,
                                                                                #       once I have a BeamList (subclass of ndarray, see
                                                                                #       https://numpy.org/doc/stable/user/basics.subclassing.html)
                                                                                #       this can be a BeamList, that can preserve the overall type
                                                                                #       as well as the type of its elements (which will be BeamIndex)
        hk = np.mgrid[-h_max:h_max+1, -h_max:h_max+1].reshape(2,-1).T

        # Create lattice:
        # Notice that, given a row vector of indices (h, k) and a basis in
        # matrix form B = [[a1, a2],[b1, b2]], the corresponding lattice
        # point can be obtained as the row vector
        #    L = [L1, L2] = [h, k]*B
        lattice = np.dot(hk, basis)

        # Now find all those lattice points that lie within limit, and use
        # this as a mask for the output
        mask = np.linalg.norm(lattice, axis=1) <= limit

        return lattice[mask], hk[mask]

    def get_rotated_lattice(self, angle):
        """
        Returns a copy of Lattice.lattice rotated by angle

        Parameters
        ----------
        angle: float
               Rotation angle in degrees. Positive values will give a
               counterclockwise rotation

        Returns
        -------
        np.array, same shape as Lattice.lattice
        """
        angle = np.radians(angle)

        rot = np.array([[np.cos(angle), np.sin(angle)],
                        [-np.sin(angle), np.cos(angle)]])

        return np.dot(self.lattice, rot)

    def get_rotated_basis(self, angle):
        """
        Returns a rotated copy of the lattice basis

        Parameters
        ----------
        angle: float
               Rotation angle in degrees. Positive values will give a
               counterclockwise rotation

        Returns
        -------
        rotated basis: 2x2 np.array
        """
        angle = np.radians(angle)
        rot = np.array([[np.cos(angle), np.sin(angle)],
                        [-np.sin(angle), np.cos(angle)]])
        return np.dot(self.basis, rot)

    def high_symm_transform(self):
        """
        Returns a 2x2 np.array that brings the basis lattice into the highest
        possible symmetry configuration when left-multiplied to the lattice
        basis. The highest-symmetry basis generates the same lattice as the
        any other basis.

        The transform can bring an oblique lattice into square, rectangular,
        hexagonal or rhombic, or make the basis vectors as close to orthogonal
        as possible.

        This function DOES NOT transform the lattice. Call make_high_symmetry()
        for that.
        """
        # The following line should also redefine the basis so that it is
        # obtuse for real-space rhombic and hexagonal lattices
        _shape = self.cell_shape

        # Will always work on the real-space lattice for convenience,
        # then convert back to the reciprocal one in case the lattice was
        # reciprocal in the first place
        basis = self.real_basis

        # In what follows, t_elem is used to define a specific elementary
        # operation to be performed on the lattice basis. This is
        # left-multiplied to t_overall at each elementary step, so that
        # t_overall contains the overall transformation

        if _shape != 'ob':
            return (1, 0), (0, 1)

        # Lattice is oblique.
        # Transform lattice to have the shortest two vectors, with angle
        # closest to 90°.
        # This might bring it to rect, hex or rhombic.
        #
        # If neither, will anyway transform to have the closest to rect.

        # ALGORITHM for reduction to closest to rect:
        # This is a discrete version of Gram-Schmidt's algorithm to find
        # orthogonal bases
        # At each iteration:
        #   - order vectors by norm, the shortest first
        #   - determine the projection of the second on the first,
        #     and calculate the nearest integer kk
        #   - subtract from the second the projection calculated above
        #   - check whether now the second is the smallest.
        #     If yes, repeat, otherwise finished.

        # swap keeps track of whether the first and second vectors are
        # swapped at the end of this passage
        swap = (1, 0), (0, 1)
        t_overall = (1, 0), (0, 1)
        while True:
            # Swap vectors if needed to get the shortest first
            if np.linalg.norm(basis[0]) > np.linalg.norm(basis[1]):
                t_elem = (0, 1), (1, 0)
            else:
                t_elem = (1, 0), (0, 1)
            swap = np.dot(t_elem, swap)
            t_overall = np.dot(t_elem, t_overall)
            basis = np.dot(t_elem, basis)
            projection = np.dot(basis[0], basis[1])/np.dot(basis[0],
                                                           basis[0])
            projection = int(np.round(projection))
            t_elem = (1, 0), (-projection, 1)
            t_overall = np.dot(t_elem, t_overall)
            basis = np.dot(t_elem, basis)
            if np.linalg.norm(basis[0]) <= np.linalg.norm(basis[1]):
                break
        # Swap vectors back if they were overall swapped
        t_overall = np.dot(swap, t_overall)
        basis = np.dot(swap, basis)

        # END OF ALGORITHM. Now the lattice is closest to rectangular.
        # It might be still any shape (square, rect, hex, rhombic, oblique)

        # Create a dummy lattice with the new basis,
        # to check which shape it has
        _shape = Lattice(basis).cell_shape

        # If it's still oblique, try to see if it can be transformed to hex
        # or rhombic by choosing "a" not to be the shortest vector of all.
        #
        # If possible, keep the new transformation. Otherwise, stick to the
        # one that makes it closest to rectangular
        #
        # All the operations that follow are stored in a matrix t_second,
        # to be later left-multiplied to t_overall to get the full
        # transformation
        #
        if _shape != 'ob':
            t_second = (1, 0), (0, 1)
        else:
            # lattice is still oblique, even if closest to rectangular
            #
            # Re-swapping guarantees that that the matrix has on the first
            # line the shortest possible vector,
            # and on the second line the second shortest possible vector.
            #
            # The only possible combinations that can lead to a
            # rhombic/hex are a'=b+a or a'=b-a, depending
            # on whether the angle is acute or obtuse, respectively
            #
            t_second = swap
            basis = np.dot(swap, basis)

            t_elem = [[-int(np.sign(np.dot(basis[0], basis[1]))), 1],
                      [0, 1]]
            t_second = np.dot(t_elem, t_second)
            basis = np.dot(t_elem, basis)
            dummy2 = Lattice(basis)

            # The following line might change acute into obtuse
            # --> check if it was the case
            # PROBABLY IT CANNOT HAPPEN ANYWAY!! IT SHOULD RATHER BE
            # CHECKED ON DUMMY ABOVE
            _shape = dummy2.cell_shape
            sign_before = np.dot(basis[0], basis[1])
            sign_after = np.dot(dummy2.basis[0], dummy2.basis[1])
            if sign_before*sign_after < 0:
                # sign did change -> basis went from acute to obtuse
                t_second = np.dot([[0, -1], [1, 0]], t_second)

            if _shape == 'ob':
                # lattice is still oblique, no transformation is needed
                # (will keep the one closest to rect)
                t_second = (1, 0), (0, 1)
        t_overall = np.dot(t_second, t_overall)

        # Finally update the transformation matrix to account for the correct
        # space of the lattice
        if self.space == 'reciprocal':
            t_overall = np.linalg.inv(t_overall).T
        return t_overall
    
    def is_high_symmetry(self):
        """
        Checks whether the lattice is in the highest symmetry possible
        """
        return np.array_equal(((1, 0), (0, 1)), self.high_symm_transform())

    def make_high_symmetry(self):
        """
        Transform lattice basis to the one that realizes the highest possible
        symmetry and update the value of the attributes. Notice that if the
        basis can be transformed to higher symmetry, determining whether its
        plane group has higher symmetry than before is not possible.

        Returns
        -------
        np.array, shape = (2, 2)
        transformation matrix that has been left-multiplied to the basis
        """
        transform = self.high_symm_transform()

        if not np.array_equal([[1, 0], [0, 1]], transform):
            # lattice can be higher symmetry, i.e., it was oblique
            self.basis = np.dot(transform, self.basis)
            self._shape = self.__get_cell_shape()
            self.lattice, self.hk = self.__generate_lattice()
            # no need to change the group, since group was at most p2.
            # Inferring whether the group has higher symmetry is not possible

        return transform
    
    def transform(self, transform, as_copy=True):
        """
        Modifies the basis and lattice points according to transform, and
        returns a copy (if copy=True) or directly modifies the values of self.
        
        Returns
        -------
        Lattice
        Either a reference to self (copy=False) or a reference to a new Lattice.
        in both cases the lattice returned is updated with the transform
        """
        new_lattice = self
        if as_copy:
            new_lattice = copy.deepcopy(self)
        new_lattice.basis = np.dot(transform, self.basis)
        return new_lattice
    
    


