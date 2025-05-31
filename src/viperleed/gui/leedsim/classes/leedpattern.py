"""Module leedpattern of viperleed.gui.leedsim.classes.

Defines the LEEDPattern class, used for displaying the LEED pattern,
potentially originating from multiple structural domains.
This is the NEW version of the classes, and replaces the
old implementations that is currently kept in oldleedpattern.py
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-03-13'
__license__ = 'GPLv3+'

import numpy as np

from viperleed.gui.classes.beamindex import BeamIndex
from viperleed.gui.classes.lattice2d import Lattice2D
from viperleed.gui.leedsim.classes.leedsubpattern import (
    CachedLEEDSubpatternConstructor,
    LEEDSubpattern,
    )
from viperleed.gui.leedsim.classes.leedparameters import LEEDParametersList
from viperleed.gui.leedsim.classes.structdomains import LEEDStructuralDomains
from viperleed.gui.leedsim.utils import screen_radius
from viperleed.gui.leedsim.utils import sort_hk

from viperleed.gui import decorators as dev_


# TODO: may need to implement a way to change the radii to allow
#       non-normal incidence to work decently

class LEEDPattern:
    """Handle the LEED pattern of multiple structural domains.

    LEEDPattern handles the LEED pattern of multiple structural domains,
    each of which may have multiple symmetry-induced domains. It also
    stores all the information needed to plot the LEED pattern.

    Parameters
    ----------
    It can be constructed from at least one instance of the following:
    - dict : assumed to be a LEED parameters dictionary
    - ConfigParser : assumed to be a LEED parameters dictionary
    - LEEDParameters
    - LEEDPattern
    - LEEDParametersList
    When passing multiple inputs, also an heterogeneous collection of
    the above is acceptable. However, the inputs will be checked for
    compatibility, i.e., they should all have the same bulk lattice
    (both basis, group and 3D symmetry operations should match)

    keep_duplicates : bool, optional (default=False)
        upon instantiation or when appending, determines whether or not
        to keep duplicates, i.e., input parameters that would generate
        identical LEED patterns
    """

    def __new__(cls, *args, **_kwargs):
        """When only one LEEDPattern instance is given, return it."""
        if len(args) == 1 and isinstance(args[0], LEEDPattern):
            return args[0]
        return super().__new__(cls)

    # @dev_.profile_calls(print_args=[30])
    # @dev_.profile_lines
    @dev_.exec_time
    def __init__(self, leed, *other_leeds, **_kwargs):
        """Initialize LEEDPattern.

        Parameters
        ----------
        leed : dict, ConfigParser, LEEDParameters, LEEDPattern, or
            LEEDParametersList
            Parameters defining a LEED pattern.
        *other_leeds : same as leed, optional
            Other parameters defining additional LEED patterns.
            Combined in a LEEDParametersList together with those given
            in leed.
        keep_duplicates : bool, default=False
            Whether duplicates should be kept in the input lists, where
            parameters are considered duplicates if they would produce
            the very same LEED pattern, including symmetry-induced
            domains.
        **_kwargs
            Other keyword arguments that are not used.

        Returns
        -------
        None.

        """
        # Skip everything when a single LEEDPattern is passed, as
        # __new__ already takes care of this
        if not other_leeds and isinstance(leed, LEEDPattern):
            return

        # The only other special case to be treated is when one of the
        # arguments is a LEEDPattern
        leeds = [leed, *other_leeds]
        for i, leed_param in enumerate(leeds):
            if isinstance(leed_param, LEEDPattern):
                leeds[i] = leeds[i].parameters

        # process acceptable keyword arguments
        keep_duplicates = _kwargs.get('keep_duplicates', False)

        # and prepare the LEEDParametersList that defines the LEED
        # pattern, i.e., a list of STRUCTURAL domains. Consistency of
        # bulk lattices is delegated to LEEDParametersList
        params = LEEDParametersList(leeds, keep_duplicates)

        # Construct the structural domains
        self.domains = LEEDStructuralDomains(params)
        self.__bulk = self.bulk

        # In-plane azimuthal rotation in degrees; positive
        # counterclockwise.  Notice that this is independent from the
        # azimuthal angle of the incidence beam: it is merely a
        # graphical in-plane rotation of the whole pattern.  The
        # underlying lattices are not changed.  This rotation is only
        # applied to the sub-patterns.
        self.__rotation = 0

        # And build an few LEEDSubpattern instances (cached) that one
        # would normally like to see, i.e.:
        # - The bulk patterns, stored in a special attribute.
        self.__bulk_patterns = self.__build_bulk_subpatterns()

        # - Each structural domain separately, including all
        #   symmetry-induced domains
        # - The first domain of each structural domain separately
        # - All structural domains together, including all
        #   symmetry-induced domains.  This will be the starting
        #   configuration when a LEEDPattern is built.
        self.__subpatterns = []   # assigned in .build_subpatterns()
        for struct_id in self.domains.domain_ids:
            all_sym_domains = {struct_id: None}
            first_sym_domain = {struct_id: 0}
            print(f'\n{struct_id} ALL')
            self.build_subpatterns(domains_ids=all_sym_domains)
            print(f'\n{struct_id} FIRST')
            self.build_subpatterns(domains_ids=first_sym_domain)
        print('\nEVERYTHING')
        self.build_subpatterns()

    @property
    def bulk(self):
        """Reciprocal-space bulk Lattice2D."""
        if hasattr(self, '__bulk'):                                             # TODO: ugly
            return self.__bulk

        max_radius = self.screen_radius(self.max_energy)
        return Lattice2D(self.domains.bulk_basis, space='reciprocal',
                         group=self.parameters[0]['bulkGroup'],
                         limit=max_radius)

    @property
    def max_energy(self):
        """Maximum LEED energy."""
        return self.parameters[0]['eMax']

    @property
    def n_structures(self):
        """Return number of structural domains."""
        return self.domains.n_domains

    @property
    def n_symmetry_domains(self):
        """List the number of symmetry-induced domains for each structure."""
        return tuple(dom.n_domains for dom in self.domains)

    @property
    def parameters(self):
        """LEEDParametersList that is used to build this LEED pattern."""
        return self.domains.parameters

    @property
    def pattern(self):
        """Return the LEED pattern as a dict of LEEDSubpatterns.

        It also applies the current rotation to each of the
        sub-patterns.  One needs to call .build_subpatterns asking for
        the desired domains before being able to retrieve the correct
        pattern!

        Returns
        -------
        dict
            keys : {'bulk', 'surf'}
            values : list
                List of LEEDSubpattern objects.
        """
        self.rotate(self.__rotation)
        return {'bulk': self.__bulk_patterns, 'surf': self.__subpatterns}

    @property
    def primary_beam_angles(self):
        """Polar and azimuthal angle of incidence for the primary beam.

        Returns
        -------
        theta : float
            Polar angle of incidence in degrees, measured with
            respect to the direction perpendicular to the surface.
            Range [0,90].
        phi : float
            Azimuthal angle of incidence in degrees, measured with
            respect to the positive x axis in the Cartesian reference
            of the bulk.  Positive counterclockwise when looking down
            on the surface.  Range [0, 360)
        """
        return self.parameters[0]['beamIncidence']

    @primary_beam_angles.setter
    def primary_beam_angles(self, theta=None, phi=None):
        """Set direction angles of primary beam.

        May take some time, as it triggers a recomputation of the
        beam equivalence and superposition. If either angle is not
        given, the one in self.parameters is taken.

        Parameters
        ----------
        theta : float, optional
            Polar incidence angle in degrees, measured with
            respect to the direction perpendicular to the surface.
            Range [-90, 90].
        phi : float, optional
            Azimuthal angle of incidence in degrees, measured with
            respect to the positive x axis in the Cartesian reference
            of the bulk.  Positive counterclockwise when looking down
            on the surface.

        Returns
        -------
        None.
        """
        self.domains.set_beam_incidence(theta, phi)

    @property
    def superlattices(self):
        """Superlattice matrices for each structural domain."""
        return [dom.superlattices for dom in self.domains]

    @property
    def __equivalent_beams_dict(self):
        """Return the dictionary of beam equivalences."""
        return self.domains.eq_beams_last_config.equivalent_beams_dict

    # @dev_.profile_lines
    @dev_.exec_time
    def build_subpatterns(self, domains_ids=None):
        """Build LEED subpatterns given a domain specification.

        This function needs to be called each time one wants to plot
        a new combination of domains. It internally sets the
        information that is needed for plotting via matplotlib

        Parameters
        ----------
        domains_ids : dict, list, int or None (default=None)
            When passing a dictionary, the keys are used to select
            the structural domains based on their name; names not
            found are skipped. When passing a list it should have as
            many elements as there are structural domains. In both
            cases each value/element can be either None (selects all
            domains), an integer, or a list of integers, selecting
            which symmetry-equivalent domains to use. When a single
            integer is passed, this is used as the positional index
            for the symmetry domains in all the structural domains.

        Returns
        -------
        None.
        """
        # Rebuild the underlying LEEDEquivalentBeams object with
        # the configuration requested
        self.domains.equivalent_spots(domains_ids)

        # And create the subpatterns
        subpatts = CachedLEEDSubpatternConstructor(
            self, domains_ids=self.domains.last_domains_used
            )

        self.__subpatterns = subpatts.patterns

    @dev_.exec_time
    def beams_equivalent_to(self, beam, in_format='', out_format=''):
        """Given a beam, return beams that are equivalent to it.

        This function cannot be used before issuing once
        self.equivalent_spots at the same angular conditions
        and with the same domains.

        Parameters
        ----------
        beam : 2-element array-like of numbers
            Some reference to the beam. Which reference is passed
            should be explicitly specified with the in_format parameter
        in_format : {'fractional', 'numerator', 'g'}
            'fractional':
                when passing the full bulk fractional indices. In this
                case, beam can be a list-like with two int, or two
                floats, or a string in the form 'num1/den1, num2/den2'
                ('/den_i' can be skipped if == 1)
            'numerator':
                when passing only the numerator of the bulk fractional
                index.  In this case, beam should be a tuple of
                integers.  The denominator is taken automatically from
                the superlattice matrix
            'g':
                when passing the vector position in Cartesian
                coordinates.  In this case, beam should be expressed in
                a reference system with z orthogonal to the surface (DO
                WE NEED THIS? PROBABLY NOT AFTER I FIGURE OUT THE
                GENERAL POSITIONS ON THE SCREEN), and in the same
                in-plane reference as the bulk basis
        out_format : str, optional (default='')
            Same values as in_format. 'numerator' probably doesn't make
            much sense. When not given, the same format as in_format is
            used.

        Returns
        -------
        list
            As many elements as there are beams equivalent to beam,
            each one is a (beam, index, text) tuple, with
            beam : tuple
                formatted as requested in out_format
            index : BeamIndex
                fractional indices
            text : str
                form "i+(j)+..." where i,j,... are the indices of the
                domains that contribute to the beam, and parenthesized
                ones are extinct contributions

        Raises
        ------
        ValueError
            If beam is not one of the beams of the LEED pattern.
        """
        # Potentially useful for the hovering annotations.
        # They need the following info:
        # - beams equivalent to beam
        # - for all beams:
        #   * some type of coordinate, probably best to give out the g
        #     vectors? it will need to be processed later on to get the
        #     right rotation (depending on the view angle), and,
        #     consequently the pixel coordinates in the canvas
        #   * the fractional index + a list of the domains that
        #     contribute to that beam, perhaps already formatted as a
        #     string?
        if self.domains.fractional:
            self_format = 'fractional'
        else:
            self_format = 'numerator'

        # Convert the beam to the internal format used by self.domains
        beam_indices = self.__reformat_beams([beam], in_format=in_format,
                                             out_format=self_format)[0]

        if not self.contains_beam(beam_indices):
            raise ValueError(f"beams_equivalent_to: beam {beam} "
                             f"(transformed to {beam_indices}) "
                             "not found. It may not be an acceptable "
                             "index for the system, or may lie outside "
                             "the energy range")

        # Get the beams equivalent to beam, sorting by (h, k).
        eq_beams = sorted(self.__equivalent_beams_dict[beam_indices],
                          key=sort_hk)

        # Now get the names of the domains that overlap at each of the
        # beams, also accounting for whether the domain contributes
        # with an extinct spot
        overlapping = self.domains_contributing_to_beams(eq_beams)

        # Fractional-formatted beams
        fractional = self.__reformat_beams(eq_beams, in_format=self_format,
                                           out_format='fractional')

        # And the beams formatted as requested
        if not out_format:
            out_format = in_format
        if out_format != 'fractional':
            out_beams = self.__reformat_beams(eq_beams, in_format=self_format,
                                              out_format=out_format)
        else:
            out_beams = fractional
        return list(zip(out_beams, fractional, overlapping))

    def contains_beam(self, beam_indices):
        """Return whether a beam is part of self.

        Parameters
        ----------
        beam_indices : tuple-like
            2-element tuple-like with indices in the same format as in
            self.domains, i.e., fractional- or numerators-only.

        Returns
        -------
        bool
            True if beam_indices is a beam of the LEED pattern
        """
        last_eq = self.domains.eq_beams_last_config
        return beam_indices in last_eq.equivalent_beams_dict

    # @dev_.profile_lines
    @dev_.exec_time
    def domains_contributing_to_beams(self, beams_indices):
        """Return which domains contribute to given beams.

        Parameters
        ----------
        beams_indices : iterable
            Each element is the index of a beam, with the same
            format as in self.domains (i.e., fractional or
            numerator-only).

        Returns
        -------
        list of str
            One string per beam; each string has the form:
            "struct1:sym1+sym2+...; struct2:sym1+sym2+...".
            "struct_id: ..." is parenthesized if the domain
            contributes only in an extinct manner, i.e., its
            part will look like "(struct_id:...)".

            THIS IS NOT TRUE, BUT IT WOULD BE NICE PERHAPS:
            Similarly, those symmetry-induced
            domains that contribute in an extinct manner in a
            non-fully-extinct domain are parenthesized.
        """
        def __format(overlapped, extinct):
            """Format list of domains as per __doc__."""
            overlapped = dict(overlapped)
            extinct = dict(extinct)
            formatted = []
            for struct_id, symm_ids in overlapped.items():
                ext_ids = extinct.get(struct_id, [])

                if all(sym_id in ext_ids for sym_id in symm_ids):
                    this_format = f'{struct_id}:'
                    this_format += '+'.join(str(ext_id) for ext_id in ext_ids)
                    this_format = f'({this_format})'
                else:
                    this_format = '+'.join(str(sym_id)
                                           for sym_id in symm_ids
                                           if sym_id not in ext_ids)
                    if this_format and ext_ids:
                        this_format += '+'
                    if ext_ids:
                        this_format += "({})".format(
                            '+'.join(str(ext_id) for ext_id in ext_ids)
                            )
                    this_format = f'{struct_id}:' + this_format
                formatted.append(this_format)
            return '; '.join(formatted)

        # beams, _, domains, extinct_domains = zip(*self.domains.indexed_beams)
        beams, domains, extinct_domains = [self.domains.indexed_beams[k]
                                           for k in ('beams',
                                                     'overlap_domains',
                                                     'extinct_domains')]
        overlapping = []
        for beam in beams_indices:
            idx = beams.index(beam)
            overlapping.append(__format(domains[idx], extinct_domains[idx]))
        return overlapping

    def rotate(self, angle):
        """Rotate the whole LEED pattern by a given angle.

        This function rotates each of the sub-patterns.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees, measured with respect to the
            positive x axis in the same Cartesian reference as the bulk
            basis. Positive counterclockwise.

        Raises
        ------
        TypeError
            If angle is not a number
        """
        try:
            float(angle)
        except TypeError as err:
            raise TypeError("Rotation angle should be a scalar, not "
                            f"{type(angle).__name__!r}") from err
        self.__rotation = angle
        angle = np.radians(angle)
        sine, cosine = np.sin(angle), np.cos(angle)

        rotation = (cosine, sine), (-sine, cosine)
        for pattern in (*self.__bulk_patterns, *self.__subpatterns):
            pattern.transform_beams(rotation)

    def screen_radius(self, energy):
        """Radius of LEED screen at given energy."""
        return screen_radius(energy, self.parameters[0]['screenAperture'])

    # @dev_.exec_time
    def __build_bulk_subpatterns(self):
        """Return the subpatterns for the bulk.

        Returns
        -------
        subpatterns : list
            A list that contains at most two LEEDSubpattern
            elements.  The first one is the subpattern for
            non-extinct bulk beams and is always present.
            The second one contains extinct beams, and exists
            only if the bulk plane group is pg, pmg or p4g.
        """
        # pylint: disable=compare-to-zero

        # TODO: should one account for the beamIncidence already here,
        #       or should we do that in the .pattern property (i.e.,
        #       only when a pattern is actually requested)?
        indices, beams = self.bulk.hk, self.bulk.lattice
        group_name = self.bulk.group.group
        if 'g' not in group_name:  # no glide
            extinct_beams = []
            non_extinct_beams = beams
        else:
            # Get boolean arrays that pick the extinct spots along the
            # [1 0] and [0 1] directions
            extinct_10 = (indices[:, 1] == 0) & (indices[:, 0] % 2 == 1)
            extinct_01 = (indices[:, 0] == 0) & (indices[:, 1] % 2 == 1)

            # pg, pmg, p4g are the only options, pg and pmg may also
            # have directions
            if '[1 0]' in group_name:
                extinct_mask = extinct_10
            elif '[0 1]' in group_name:
                extinct_mask = extinct_01
            else:
                extinct_mask = extinct_10 | extinct_01
            extinct_beams = beams[extinct_mask]
            non_extinct_beams = beams[~extinct_mask]

        bulk_patterns = [LEEDSubpattern(self, non_extinct_beams, '__b')]
        if extinct_beams:
            bulk_patterns.append(LEEDSubpattern(self, extinct_beams, '__b__e'))
        return bulk_patterns

    def __reformat_beams(self, beams, **kwargs):
        """Reformat beams to fractional, numerator or g vector.

        The beams can be passed either as a reciprocal-space vector
        or as a 2-tuple-like of indices.

        Parameters
        ----------
        beams : iterable
            Each element is a 2-element array-like of numbers,
            representing a beam.  Which format is passed should
            be explicitly specified with the in_format parameter.
        in_format : {'fractional', 'numerator', 'g'}
            'fractional':
                when passing the full bulk fractional indices. In this
                case, each beam can be a list-like with two int, or two
                floats, or a string in the form 'num1/den1, num2/den2'
                ('/den_i' can be skipped if == 1).
            'numerator':
                when passing only the numerator of the bulk fractional
                index.  In this case, beam should be a tuple of
                integers.  The denominator is taken automatically
                from the superlattice matrices.  This input format
                is considered invalid if the LEEDPattern contains
                multiple structural domains with differently sized
                cells.
            'g': ### THIS MAY NOT WORK!!
                when passing the vector position in Cartesian
                coordinates.  In this case, beam should be expressed in
                a reference system with z orthogonal to the surface (DO
                WE NEED THIS? PROBABLY NOT AFTER I FIGURE OUT THE
                GENERAL POSITIONS ON THE SCREEN), and in the same
                in-plane reference as the bulk basis
        out_format : str, optional (default='')
            Same values as in_format.  When not given, the same format
            as in_format is used.
        **kwargs : dict
            Other keyword arguments are ignored.

        Returns
        -------
        reformatted_beams : numpy.ndarray or tuple of BeamIndex
            Return type depends on out_format as follows
            'out_format' == 'g' : numpy.ndarray
            'out_format' == 'fractional' : tuple of BeamIndex
            'out_format' == 'numerator' : tuple of BeamIndex

        Raises
        ------
        ValueError
            If in_format is not one of the acceptable values
        ValueError
            If either in_format or out_format is 'numerator'
            but the domains cannot be expressed as numerators
            only (i.e., there are different structures with
            differently-sized unit cells).
        """
        in_format = kwargs.get('in_format', '')
        out_format = kwargs.get('out_format', in_format)
        if in_format not in ('fractional', 'numerator', 'g'):
            raise ValueError("in_format is mandatory, and should be either "
                             "'fractional', 'numerator', or 'g'")
        # (1) Request 'numerator' but LEEDPattern needs to be expressed
        #     with fractional indices --> Error
        if 'numerator' in (in_format, out_format) and self.domains.fractional:
            raise ValueError("Invalid format 'numerator' when multiple "
                             "domains with different cells are present.")

        # (2) Request 'g' and input is already 'g' --> return as is
        if in_format == out_format == 'g':
            return np.asarray(beams)

        # (3) Request 'n' and input is already 'n' --> return as is
        if in_format == out_format == 'numerator':
            return tuple(BeamIndex(b) for b in beams)

        # (4) Input is 'g' or 'n', output is something different
        #     --> convert to fractional, so these cases can be handled
        #     the same way as 'f' -> 'f' and 'f' -> 'n'
        if in_format == 'g':
            # Since g_bulk = (h,k)_bulk @ bulk_basis, with bulk_basis
            # the reciprocal-lattice bulk basis, then:
            #       (h,k)_bulk = g_bulk @ bulk_basis^(-1).
            # The indices are the fractional ones in this case.
            beams = np.dot(beams, np.linalg.inv(self.bulk.basis))
            in_format = 'fractional'
            # A possible way to solve numerical errors when g vectors
            # are given is either: (1) store a map of (indices, g)
            # in LEEDSubPattern; this would allow to use directly the
            # 'ind' from .contains, (2) have a way in both
            # struct and sym domains to retrieve the bulk index
            # given a certain vector. May be possible using indexing
            # of .lattice and .hk from Lattice.

        # Convert everything to a list of BeamIndex, trying to limit
        # as much as possible numerical errors.  The only case in which
        # serious numerical errors will occur, and are unavoidable, is
        # when self.domains must be expressed as fractional indices and
        # we are passing in g vectors or pure floats.
        # Actually, given the current (2021-03-27) implementation of
        # BeamIndex, BeamIndex(floats) would raise errors!
        den = self.domains[0].denominator_for_bulk_beams
        if in_format == 'fractional':
            if not self.domains.fractional:
                beams = [BeamIndex(b, denominator=den) for b in beams]
            else:
                # Here numerical errors will likely occur!
                beams = [b if isinstance(b, BeamIndex) else BeamIndex(b)
                         for b in beams]
        else:  # 'numerator'
            beams = [BeamIndex(b, denominator=den, from_numerators=True)
                     for b in beams]

        if out_format == 'numerator':
            return tuple(b*den for b in beams)
        if out_format == 'fractional':
            return tuple(beams)
        return np.dot(beams, self.bulk.basis)
