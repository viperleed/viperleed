"""
======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes.leedpattern ***

Defines the LEEDPattern and LEEDsubpattern classes, used for displaying a
LEEDPattern. This is the NEW version of the classes, and replaces the old
implementations that is currently kept in oldleedpattern.py

Author: Michele Riva
Created: 2021-03-13
"""

import numpy as np

from viperleed import guilib as gl


class LEEDPattern:
    """
    LEEDPattern handles the LEED pattern of multiple structural domains, each of
    which may have multiple symmetry-induced domains. It also stores all the
    information needed to plot the LEED pattern.
    
    Parameters
    ----------
    It can be constructed from at least one instance of the following:
    - dict : assumed to be a LEED parameters dictionary
    - ConfigParser : assumed to be a LEED parameters dictionary
    - LEEDParameters
    - LEEDPattern
    - LEEDParametersList
    When passing multiple inputs, also an heterogeneous collection of the above
    is acceptable. However, the inputs will be checked for compatibility, i.e.,
    they should all have the same bulk lattice (both basis, group and 3D
    symmetry operations should match)
    
    duplicates : bool, optional (default=False)
        upon instantiation or when appending, determines whether or not to
        keep duplicates, i.e., input parameters that would generate identical
        LEED patterns
    """
    def __new__(cls, *args, **kwargs):
        """
        The only purpose of this is to handle the case in which a single
        LEEDPattern instance is passed. In that case the same instance is
        returned as is
        """
        if len(args) == 1 and isinstance(args[0], LEEDPattern):
            return args[0]
        return super().__new__(cls)
    
    # @gl.profile_calls(print_args=[30])
    def __init__(self, leed, *other_leeds, **kwargs):
        # Treat the special case of when a single LEEDPattern is passed:
        # __new__ already returns a reference to the LEEDPattern, so we
        # just skip everything
        if len(other_leeds) == 0 and isinstance(leed, LEEDPattern):
            return None

        # In all other cases, process the arguments one at a time,
        # collecting the necessary information. In fact, the only special
        # case to be treated is when one of the arguments is a LEEDPattern,
        # in which case we take its .parameters attribute (a LEEDParametersList)
        leeds = (leed, *other_leeds)
        for i, leed in enumerate(leeds):
            if isinstance(leed, LEEDPattern):
                leeds[i] = leeds[i].parameters

        # process acceptable keyword arguments
        duplicates = kwargs.get('duplicates', False)
        
        # and prepare the LEEDParametersList that defines the LEED pattern,
        # i.e., a list of STRUCTURAL domains. Consistency of bulk lattices
        # is delegated to LEEDParametersList
        self.__parameters = gl.LEEDParametersList(leeds, duplicates)

        # Construct the structural domains
        # TODO: self.domains = 
        

    @property
    def parameters(self):
        """
        Returns the LEEDParametersList that is used to build this LEED pattern
        """
        return self.__parameters

    @property
    def max_energy(self):
        """
        Maximum LEED energy
        """
        return self.__parameters['eMax']

    @property
    def n_domains(self):
        """
        Number of distinct domains produced by the symmetry operations of the
        bulk
        """
        return len(self.domains['operations'])

    @property
    def superlattices(self):                                                    # TODO: REMOVE?
        """
        Numpy ndarray of superlattice matrices generating the symmetry-related
        domains. The first element is the one generating the lattice whose basis
        is given in the constructor as a LEED parameter
        """
        superlattice = self.__parameters['SUPERLATTICE']
        operations = self.domains['operations']                                 # TODO: needs to change
        if operations is None:
            operations = [gl.PlaneGroup.E]  # identity matrix
        return np.einsum('ij,mjk->mik', superlattice, operations)
    
    @property
    def primary_beam_theta(self):
        """
        Polar angle of incidence with respect to the direction perpendicular
        to the surface
        """
        return self.__parameters['beamIncidence'][0]
    
    @property
    def primary_beam_phi(self):
        """
        Azimuthal angle of incidence of the primary beam with respect to the
        x axis. Positive counterclockwise when looking down on the surface.
        """
        return self.__parameters['beamIncidence'][1]

    def screen_radius(self, en):
        return gl.screen_radius(en, self.__parameters['screenAperture'])

    def rotate(self, angle, which='surf'):
        """
        This function is a mess and needs to be split into two functions.
        It currently does two different things: if which='bulk', it returns a
        copy of the bulk lattice. If which='surf' it rotates all the LEED
        subpatterns. It would probably be nicer to also have the bulk pattern
        be a special subpattern.
        """
        which = which[0:4]
        if which not in ['surf', 'bulk']:
            raise ValueError("LEEDPattern: Only 'surf', 'surface', or "
                             "'bulk' are acceptable values for rotate()")
        if hasattr(angle, '__len__'):
            raise ValueError('Angle input is not a scalar')

        if which == 'surf':
            angle = np.radians(angle)
            rot = np.array([[np.cos(angle), np.sin(angle)],
                            [-np.sin(angle), np.cos(angle)]])
            for pat in [*self.firstLEED, *self.domsLEED]:
                pat.transformCoords(rot)
            return None
        return self.reciprocal_lattices['bulk'].get_rotated_lattice(angle)


class LEEDStructuralDomains:
    pass

# TODO: LEEDsubpattern
