"""
======================================
  ViPErLEED Graphical User Interface
======================================
      *** module guilib.base ***

Created: 2020-01-11
Author: Michele Riva

This module contains base functions and classes shared by the whole GUI
"""

import sys
import ast
import re
from fractions import Fraction
from collections import defaultdict

import numpy as np

from guilib.leedsim.classes import LEEDPattern

###############################################################################
#                                  FUNCTIONS                                  #
###############################################################################


def get_equivalent_beams(leed_parameters, domains=None):
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

    domains: iterable or None, default=None
             List of domain indices for which the beams get exported
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
             right choice when asking for domains

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

        The entries in the list are sorted in this order: id, h, k
    """

    if not check_leed_params(leed_parameters):
        return None

    leed = LEEDPattern(leed_parameters)

    if domains is None:
        domains = range(leed.nDoms)
    elif not hasattr(domains, '__len__'):
        raise ValueError(f"The keyword argument 'domains' should be either "
                         f"an iterable or None. Found {type(domains)} "
                         "instead.")

    fract, groups, *_ = zip(*leed.get_equivalentSpots(domains=domains))
    return list(zip(fract, groups))


def project_to_first_domain(leed_parameters, beam_list, domains=None):
    """
    Given a list of beams and a dictionary of parameters defining the LEED
    geometry, returns a list of beams that 'projects' the input onto the first
    domain. I.e., the output contains only beams that belong to the first domain
    and are equivalent to those given as an input. At present [2020-05-04] this
    function does not account for non-normal beam incidence.

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
    leed = LEEDPattern(leed_parameters)
    ops = leed.bulkR.group.group_ops(include_3d=True)
    m = leed_parameters['SUPERLATTICE']

    def is_in_first_domain(beam):
        """
        Returns True if beam belongs to the first domain, False otherwise
        """
        # only those beams that, processed through the superlattice
        # matrix give integer indices belong to the first domain
        indices = np.dot(beam, m.T)
        return all(i.denominator == 1 for i in indices)
    
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
        eq_beams = get_equivalent_beams(leed_parameters, domains=[0])
        
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
                    err = (f"Beam {beam} is incompatible with the current "
                           f"SUPERLATTICE matrix \n{m}")

                    # Check if the reason why the beam was not found is that
                    # it would lie outside the LEED screen
                    leed = LEEDPattern(leed_parameters)
                    b = leed.get_BulkBasis()
                    g = np.linalg.norm(np.dot(beam, b))
                    el_m = 9.109e-31    # kg
                    el_q = 1.60218e-19  # C
                    hbar = 1.05457e-34  # J*s
                    
                    # calculate the exit angle
                    s_angle = np.sqrt(hbar**2 * g**2
                                      /(2 * el_m * el_q * leed.maxEnergy))
                    ang = 2*np.degrees(np.arcsin(s_angle))
                    aperture = leed_parameters.get('screenAperture', 110)
                    if ang > aperture:
                        err += ("\nThe beam would need a minimum LEED screen "
                                f"aperture of {ang} deg at Emax="
                                f"{leed.maxEnergy} eV, while you are using "
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
        get_equivalent_beams(leed_parameters, domains=domains))
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
        if not check_type(group, 'str'):
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

    if len(needs_shape)>0 and matrix.shape != tuple(needs_shape):
        return None
    return matrix


################################################################################
#                                   CLASSES                                    #
################################################################################


class BeamIndex(tuple):
    """
    Convenience class to store a 2-element tuple that represents a Miller index
    for a LEED beam. Each index is a fractions.Fraction
    """
    
    def __new__(cls, *indices):
        if len(indices) > 2:
            raise ValueError("BeamIndex accepts at most two indices, "
                             f"{len(indices)} given instead.")
        if len(indices) == 1:
            t = indices[0]
            if isinstance(t, str):
                indices = t.split(',')
            else:
                if not hasattr(t, '__len__'):
                    raise TypeError("BeamIndex: when one argument given, "
                                    "it should be a string or a 2-element "
                                    "array-like.")
                indices = t
            if len(indices) != 2:
                raise ValueError("BeamIndex: too many/few indices. "
                                 "Exactly 2 indices should be given. "
                                 f"Found {len(indices)} instead.")
        return super().__new__(cls, (Fraction(index) for index in indices))

    def __str__(self):
        return f"{', '.join(str(index) for index in self)}"

    def __repr__(self):
        return f"BeamIndex({', '.join(str(index) for index in self)})"


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
    group_ops() Returns tuple with group operations

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
    # EDIT 2020-06-22: all the operations below are now defined as tuples of
    #                  tuples rather than numpy arrays. This makes them
    #                  hash-able, and thus usable as set elements.
    #                  Also edited: M45 and Mm45 were exchanged.

    # These two are good for all cells
    E = ((1, 0), (0, 1))
    C2 = ((-1, 0), (0, -1))    # = -E

    # These two are good for rectangular and square cells
    Mx = ((1, 0), (0, -1))
    My = ((-1, 0), (0, 1))     # = -Mx

    # These are good for square cells
    C4 = ((0, -1), (1, 0))
    Cm4 = ((0, 1), (-1, 0))    # = -C4
    M45 = ((0, 1), (1, 0))
    Mm45 = ((0, -1), (-1, 0))  # = -M45

    # These are good for rhombic and hex (both obtuse)
    M11 = ((0, 1), (1, 0))
    M1m1 = ((0, -1), (-1, 0))  # = -M11
    M01 = ((-1, -1), (0, 1))
    M10 = ((1, 0), (-1, -1))

    # And these are good for hex only (obtuse)
    C6 = ((1, 1), (-1, 0))
    Cm6 = ((0, -1), (1, 1))
    C3 = ((0, 1), (-1, -1))
    Cm3 = ((-1, -1), (1, 0))
    M21 = ((1, 1), (0, -1))
    M12 = ((-1, 0), (1, 1))

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
                 '4': (C4, Cm4),
                 '6': (C6, Cm6)}
    # 2) direction contained in glide planes into the corresponding matrices
    glide_ops = {'[1,0]': M10,
                 '[1,1]': M11,
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
              M01, M10, M11, M1m1)

    def __init__(self, group='p1'):
        group = self.check_group_name(group)
        self.group = group
        
        # The next one will be a tuple of the 2x2 matrices representing the
        # isomorphism part of screws and glide planes perpendicular to the
        # surface
        self.__ops_3d = tuple()

    def __repr__(self):
        """
        Representation of PlaneGroup
        """
        return f"viperleed.guilib.base.PlaneGroup(group={self.group!r})"

    def __str__(self):
        """
        string representation of PlaneGroup
        """
        return self.group

    def __eq__(self, other):
        """
        Equality method for PlaneGroup instances.
        """
        if not isinstance(other, PlaneGroup):
            # Other instances are never equal
            return NotImplemented
        return self.group == other.group

    def __ne__(self, other):
        """
        Not equality method for PlaneGroup instances.
        """
        return not self == other

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
            (?:[\[]                 # optional direction opening bracket
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
        #     the group_ops method below
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
                
        if not isinstance(input, (str, np.array, tuple, list)):
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
            
            s = re.search(screw_re, input)
            g = re.search(glide_re, input)
            if not (s or g):
                raise ValueError("PlaneGroup.screws_glides: Invalid input.")
            if s:  # found some screws
                screws = s.group('screws').split(',')
                if any(s not in self.screw_ops.keys() for s in screws):
                    raise ValueError("PlaneGroup.screws_glides: Invalid "
                                     "rotation order in the input. Only 2-, "
                                     "3-, 4-, and 6-fold orders allowed.")
                [ops.extend(self.screw_ops[screw]) for screw in screws]
            if g:  # found some glide planes
                if shape is None:
                    raise ValueError("PlaneGroup.screws_glides: cell shape is "
                                     "required when glide planes are given.")
                # parse glide planes by removing spaces and splitting on commas
                g = g.group('glides').replace(' ','').split(',')

                # now g should contain an even number of elements:
                # each odd element is of the form "[#", each even "#]"
                for g_odd, g_even in zip(g[::2], g[1::2]):
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

        if len(np.shape(input)) !=3 or np.shape(input)[1,2] != (2, 2):
            raise ValueError("PlaneGroup.screws_glides: an array-like input "
                             "should be a 1D 'list' of 2x2 'matrices'. Found "
                             f"incompatible shape {np.shape(input)}.")

        if any(np.abs(mij%1) > 1e-4 for mij in np.flatten(input)):
            raise ValueError("PlaneGroup.screws_glides: an array-like input "
                             "should contain only integer-valued matrices")

        self.__ops_3d = tuple(  # * make the whole list a tuple
            tuple(              # * make each array of the input into a tuple
                map(tuple,      # * map each line of each matrix to a tuple
                    np.array(mi).round().astype(int))  # after rounding to int
                )
            for mi in input
            )

    def group_ops(self, include_3d=False):
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
    
    def transform(self, transform):
        """
        Returns a tuple of the group operations 'projected' to a new coordinate
        system, whose coordinates are expressed by the 2x2 array-like
        transformation matrix given in the "transform" parameter. No assumption
        is made on the coordinate transform (i.e., the operation matrices
        returned may have non-integer values).
        """
        if np.shape(transform) != (2, 2):
            raise ValueError("PlaneGroup.transform_to_coordinates requires a "
                             "2x2 array-like as coordinate transform matrix. "
                             f" Found shape {np.shape(transform)} instead.")
        return tuple(np.linalg.multi_dot((transform,
                                          op,
                                          np.linalg.inv(transform)))
                     for op in self.group_ops())


class Lattice():
    """
    Lattice(basis, space='real', group='p1', limit=1)

    Class for a 2D lattice. Includes basis and its shape, plane group, whether
    it is a real- or reciprocal-space lattice, and an array of lattice points.

    Parameters
    ----------
    basis: 2x2 np.array
           basis vectors a and b of the lattice, with a = basis[0], b=basis[1]
           the units are assumed to be Angstrom for real space lattices, and
           2*pi/Angstrom for reciprocal-space lattices
    space: str, default: 'real'
           accepts 'real' or 'reciprocal' for real and reciprocal space
           lattices, respectively
    group: str, default: 'p1'
           plane group in Hermann-Mauguin notation. See also the PlaneGroup
           class for acceptable values.
    limit: int, default: 1
           radius used to limit the number of lattice points generated. Only
           lattice points closer to the origin than limit will be produced.

    Properties
    ----------
    basis: 2x2 np.array
           basis, same convention as with class instantiation
    space: str
           'real' or 'reciprocal'
    type: str
          shape of the unit cell. Can be 'Oblique', 'Rectangular',
          'Square', 'Rhombic', or 'Hexagonal'
    group: PlaneGroup
           a PlaneGroup instance representing the plane group of the lattice
    lattice: np.array, shape = (..., 2)
             array of lattice points (x, y) [for real-space] (gx, gy) [for
             reciprocal space]
    hk: np.array, shape = (..., 2)
            array of (h, k) indices that generate the points in lattice, i.e.,
            lattice[i] = hk[i, 0]*basis[0] + hk[i, 1]*basis[1]

    Methods
    -------
    nbeams() Return number of lattice points

    get_lattice_type() finds type of lattice, rather use Lattice.type, unless
                       you suspect that the type is not up to date

    generate_lattice() generates the lattice points

    get_rotated_lattice(angle) returns a copy of Lattice.lattice rotated
                               by angle (degrees)

    get_rotated_basis(angle) returns a copy of Lattice.basis rotated by angle

    reciprocal_basis() returns the reciprocal of Lattice.basis

    real_basis() returns a copy of the real-space basis, independently of the
                 fact that the lattice is real or reciprocal

    lattice_parameters() return lattice parameters

    set_group(group) sets the plane group to group

    high_symm_transform() return matrix transform that gives highest symmetry
                          basis (with the same lattice).

    make_high_symmetry() transform Lattice so the basis has the highest
                         possible symmetry
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

        self._basis = np.array(basis)
        self._space = space
        self._type = self.__get_lattice_type()
        self._group = PlaneGroup(group)
        self._limit = limit
        self.lattice, self.hk = self.__generate_lattice()
    
    @property
    def basis(self):
        """
        Returns the lattice basis as a 2x2 np.array, with a = basis[0] and
        b = basis[1]
        """
        return self._basis

    @basis.setter
    def basis(self, basis):
        """
        Sets the lattice basis to basis and updates the other attributes
        
        Parameters
        ----------
        basis: 2x2 array like
        """
        if not check_type(basis, 'arraylike'):
            raise TypeError("basis must be array-like. "
                            f"Found {type(basis)} instead.")
        if not np.shape(basis) == (2, 2):
            raise ValueError("Lattice basis needs to have a (2, 2) shape. "
                             f"Found {np.shape(basis)} instead.")
        self._basis = np.array(basis)
        self._type = self.__get_lattice_type()
        self.lattice, self.hk = self.__generate_lattice()

    @property  # NB: this one does not have a setter method, so it's read-only
    def type(self):
        """
        Returns the shape of the unit cell as a string. The value returned is
        'Oblique', 'Rectangular', 'Square', 'Rhombic', or 'Hexagonal'.
        """
        return self._type

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

    def nbeams(self):
        """
        Returns
        -------
        int
        the number of lattice points (for real-space) or LEED beams (for
        reciprocal-space)
        """
        return len(self.hk)

    def __get_lattice_type(self):
        """
        Returns
        -------
        string
        """

        basis = self._basis

        # Relative tolerance factor within which things are assumed to be equal
        eps = 1e-5

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
                    print('Warning: the input file might be corrupted,'
                          'since the real-space basis is not obtuse!')
                    transform = [[0, -1], [1, 0]]  # this keeps the handedness
                else:
                    transform = [[1, 0], [0, 1]]
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

            if abs(cosine + 1/2) < eps:  # angle is 120deg
                return 'Hexagonal'
            return 'Rhombic'
        return 'Oblique'

    def __generate_lattice(self):
        """
        Generates a list of lattice points given a basis.

        Parameters
        ----------
        limit:  scalar
                Determines which portion of the lattice is generated.
                For space == 'real', the lattice is generated up to a radius
                of 1.5*limit.
                    One should then plot from -limit to +limit. This should
                    cover any post-rotation of the lattice that the user
                    might later request.
                For space == 'reciprocal', the lattice is generated up to a
                radius of limit

        basis:  2x2 array-like, default = self.basis
                Contains the unit vectors as basis[0] and basis[1]

        space:   string, {'real', 'reciprocal'}, default = self.space
                Determines whether a real-space or a reciprocal-space
                lattice is generated.
                This only affects the behavior of limit (see above)

        Returns
        -------
        Tuple (lat, hk)

        lat: 1D np.array
             Lattice points generated

        hk:  1D np.array, hk.shape == lat.shape
             Integer indices of the lattice points generated.
             The lattice points are lat = h*basis[0] + k*basis[1]
        """

        limit = self._limit
        basis = self.basis
        space = self.space

        if space == 'real':
            limit *= 1.5

        # get limit for the loops that follow
        shortest = min(np.linalg.norm(basis[0]),
                       np.linalg.norm(basis[1]),
                       np.linalg.norm(basis[0] + basis[1])/2,
                       np.linalg.norm(basis[0] - basis[1])/2)
        h_max = int(np.ceil(limit/shortest))

        # create grid of indices
        hk = np.array([(i, j)
                       for i in range(-h_max, h_max+1)
                       for j in range(-h_max, h_max+1)])

        # Create lattice:
        # Notice that, given a row vector of indices (h, k) and a basis in
        # matrix form B = [[a1, a2],[b1, b2]], the corresponding lattice
        # point can be obtained as the row vector
        #    L = [L1, L2] = [h, k]*B
        lattice = np.dot(hk, basis)

        # Now find all those lattice points that lie within limit, and use
        # this as a mask for the output
        mask = np.linalg.norm(lattice, axis=1) <= limit

        return (lattice[mask], hk[mask])

    def get_rotated_lattice(self, angle):
        """
        Returns a copy of Lattice.lattice rotates by angle

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

    def reciprocal_basis(self):
        """
        Returns the reciprocal of the lattice basis. Notice that if
        self.space == 'real', this returns the reciprocal space basis.
        If self.space == 'reciprocal', it returns the reciprocal of
        the reciprocal, i.e., the real-lattice basis.

        This is different from the return value of real_basis
        """
        # 1st and 2nd surface unit cell vectors, including rotation.
        # Shape = (3, 1)
        a_3d = np.array([*self.basis[0], 0])
        b_3d = np.array([*self.basis[1], 0])

        z_3d = np.array([0, 0, 1])

        uc_area = np.vdot(a_3d, np.cross(b_3d, z_3d))  # Area of unit cell
        a_recipr_3d = (2*np.pi/uc_area)*np.cross(b_3d, z_3d)
        b_recipr_3d = (2*np.pi/uc_area)*np.cross(z_3d, a_3d)

        return np.array([a_recipr_3d[0:2], b_recipr_3d[0:2]])

    def real_basis(self):
        """
        Always returns the real-space basis of the lattice as a 2x2 np.array
        """
        if self.space == 'reciprocal':
            return self.reciprocal_basis()
        return np.dot(self.basis, [[1, 0], [0, 1]])

    def lattice_parameters(self):
        """
        Returns the length of the lattice vectors and the angle alpha between
        them as a tuple (norm_a, norm_b, alpha)
        """
        norm_a = np.linalg.norm(self.basis[0])
        norm_b = np.linalg.norm(self.basis[1])
        alpha = np.arccos(np.dot(self.basis[0], self.basis[1])/(norm_a*norm_b))
        return (norm_a, norm_b, np.degrees(alpha))

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
        _typ = self.type

        # Will always work on the real-space lattice for convenience,
        # then convert back to the reciprocal one in case the lattice was
        # reciprocal in the first place
        basis = self.real_basis()

        # In what follows, t_elem is used to define a specific elementary
        # operation to be performed on the lattice basis. This is
        # left-multiplied to t_overall at each elementary step, so that
        # t_overall contains the overall transformation

        if _typ == 'ob':
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
            swap = [[1, 0], [0, 1]]
            t_overall = [[1, 0], [0, 1]]
            while True:
                # Swap vectors if needed to get the shortest first
                if np.linalg.norm(basis[0]) > np.linalg.norm(basis[1]):
                    t_elem = [[0, 1], [1, 0]]
                else:
                    t_elem = [[1, 0], [0, 1]]
                swap = np.dot(t_elem, swap)
                t_overall = np.dot(t_elem, t_overall)
                basis = np.dot(t_elem, basis)
                projection = np.dot(basis[0], basis[1])/np.dot(basis[0],
                                                               basis[0])
                projection = int(np.round(projection))
                t_elem = [[1, 0], [-projection, 1]]
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
            _typ = Lattice(basis).type

            # If it's still oblique, try to see if it can be transformed to hex
            # or rhombic by choosing "a" not to be the shortest vector of all.
            #
            # If possible, keep the new transformation. Otherwise, stick to the
            # one that makes it closest to rectangular
            #
            # All the operations that follow are stored in a matrix t_second,
            # to be later left-multiplied to t to get the full transformation
            #
            if _typ == 'ob':
                # lattice is still oblique, even if closest to rectangular

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
                _typ = dummy2.type
                sign_before = np.dot(basis[0], basis[1])
                sign_after = np.dot(dummy2.basis[0], dummy2.basis[1])
                if sign_before*sign_after < 0:
                    # sign did change -> basis went from acute to obtuse
                    t_second = np.dot([[0, -1], [1, 0]], t_second)

                if _typ == 'ob':
                    # lattice is still oblique, no transformation is needed
                    # (will keep the one closest to rect)
                    t_second = [[1, 0], [0, 1]]
            else:
                t_second = [[1, 0], [0, 1]]
            t_overall = np.dot(t_second, t_overall)
        else:
            t_overall = [[1, 0], [0, 1]]

        # Finally update the transformation matrix to account for the correct
        # space of the lattice
        if self.space == 'reciprocal':
            t_overall = np.linalg.inv(t_overall).transpose()
        return t_overall

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
            self.type = self.get_lattice_type()
            self.lattice, self.hk = self.generate_lattice()
            # no need to change the group, since group was at most p2.
            # Inferring whether the group has higher symmetry is not possible

        return transform
