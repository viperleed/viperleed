"""
======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes.symdomains ***

Defines the LEEDSymetryDomains class, used for constructing a LEED pattern

Author: Michele Riva
Created: 2021-03-13
"""

from collections.abc import Sequence
import itertools

import numpy as np

from viperleed import guilib as gl


class LEEDSymmetryDomains(Sequence):
    """
    Collection of LEED domains generated from one superlattice matrix
    and the symmetry operations of a bulk group. LEEDSymmetryDomains is an
    immutable sequence (i.e., tuple-like). Calling LEEDSymmetryDomains[i]
    returns a viperleed.Lattice of the domain
    """
    # @gl.exec_time
    def __init__(self, leed_parameters):
        self.__parameters = gl.LEEDParameters(leed_parameters)
        super().__init__()
        
        # self.__operations is a list of the bulk symmetry operations generating
        # distinct domains. Read by self.operations
        self.__operations = self.__find_domain_operations()
        
        # self.__superlattices is a numpy.ndarray of the superlattice matrices
        # that generate distinct domains. It can be accessed via the attribute
        # self.superlattices, and it is permanently set there the first time the
        # attribute is accessed
        self.__superlattices = None
        
        # The return value of self.equivalent_spots has beams expressed as
        # tuples, but including ONLY THE NUMERATOR of the beam indices with
        # respect to the bulk. Use self.denominator_for_bulk_beams to retrieve
        # the correct denominator as an integer, which accesses the private
        # property self.__denominator. self.__denominator is set once and for
        # all the first time self.denominator_for_bulk_beams is called
        self.__denominator = None

        # self.__domains is the underlying list that is accessed when issuing
        # self[index]. It's a list of gl.Lattice instances
        self.__domains = self.__build_domains()

        # self.__equiv_spots_no_superpos is a list of dictionaries of
        # dictionaries, one element per domain, ordered as in self.__domains
        # (i.e., in self).
        # The outer dictionary is labeled by 'norm', 'other', or an azimuthal
        # angle. For each of the angle settings, the inner dictionary is in the
        # form {beam_i: set(beams_i)} with beams of domain i and their
        # equivalents within the domain at the specified primary beam angle(s).
        # Used only internally, cannot be accessed
        self.__equiv_spots_no_superpos = self.__get_spot_equivalence()

        # self.__extinct_spots list of dictionaries, one element per domain,
        # same order as in self.__domains (i.e., in self). Each dictionary has
        # keys 'norm', 'other', or an azimuthal angle, and as values a numpy
        # array that lists the extinct beams of the domain at the specified
        # beam angle(s).
        # Used only internally, cannot be accessed
        self.__extinct_spots = self.__get_extinct()

        # self.__equiv_spots_no_superpos and self.__extinct_spots are used in
        # the public method self.equivalent_spots(domains, theta, phi) for
        # creating the appropriate LEEDEquivalentBeams instance (unless it's
        # already cached), whose reference is stored into self.__last_eq.
        # Used only internally, cannot be accessed
        self.__last_eq = None
    
    def __repr__(self):
        txt = (f"{self.cell_shape} "
               + f"viperleed.LEEDSymmetryDomains({self.__parameters})")
        return txt
    
    def __len__(self):
        return len(self.__domains)
    
    def __getitem__(self, el):
        # NB: if a slice is given, the result is NOT type-preserving, i.e.,
        # it will not be a LEEDSymmetryDomains but a list of Lattice(s)         # Will this be a problem?
        return self.__domains[el]

    @property
    def groups(self):
        return {'surf': self.__parameters['surfGroup'],
                'bulk': self.__parameters['bulkGroup']}

    @property
    def bulk_basis(self):
        """
        Returns the reciprocal-space basis of the bulk lattice
        """
        return np.dot(self.superlattices[0].T, self[0].basis).round(10)

    @property
    def n_domains(self):
        """
        Number of distinct domains produced by the symmetry operations of the
        bulk
        """
        return len(self)

    @property
    def operations(self):
        return self.__operations

    @property
    def superlattices(self):
        """
        Numpy ndarray of superlattice matrices generating the symmetry-related
        domains. The first element is the one generating the lattice whose basis
        is given in the constructor as a LEED parameter
        """
        if self.__superlattices is None:
            superlattice = self.__parameters['SUPERLATTICE']
            self.__superlattices = np.einsum('ij,mjk->mik',
                                             superlattice,
                                             self.operations)
        return self.__superlattices

    @property
    def g_vectors(self):
        """
        g_vectors[i] is a list of the lattice vectors of domain i in the same
        coordinate system as the bulk basis
        """
        return [dom.lattice for dom in self]

    @property
    def denominator_for_bulk_beams(self):
        if self.__denominator is None:
            self.__denominator = int(
                abs(round(np.linalg.det(self.superlattices[0])))
            )
        return self.__denominator

    def angular_offsets(self, zero_pi=False):
        """
        Returns a list of the angles between the first lattice vectors of all
        domains and those of the first domain
        """
        first_dom_angle = gl.orientation(self[0].real_basis[0], zero_pi)
        return [gl.orientation(dom.real_basis[0], zero_pi) - first_dom_angle
                for dom in self]

    # @gl.profile_lines
    def beams_equivalent_to(self, beam, **kwargs):    # this will probably move to LEEDPattern later on? As is, it's fast enough not to bother
        """
        Given a beam, returns beams that are equivalent to it. This function
        cannot be used before issuing once self.equivalent_spots at the same
        angular conditions and with the same domains.
        
        Parameters
        ----------
        beam : 2-element array-like of numbers
            Some reference to the beam. Which reference is passed should be
            explicitly specified with the in_format parameter
        in_format : str, keyword only
            Acceptable values are:
                - 'fractional':
                    when passing the full bulk fractional indices. In this case,
                    beam can be a list-like with two int, or two floats, or a
                    string in the form 'num1/den1, num2/den2' ('/den_i' can be
                    skipped if == 1)
                - 'numerator':
                    when passing only the numerator of the bulk fractional index
                    In this case, beam should be a tuple of integers. The
                    denominator is taken automatically from the superlattice
                    matrix
                - 'g':
                    when passing the vector position in Cartesian coordinates.
                    In this case, beam should be expressed in a reference system
                    with z orthogonal to the surface (DO WE NEED THIS? PROBABLY
                    NOT AFTER I FIGURE OUT THE GENERAL POSITIONS ON THE SCREEN),
                    and in the same in-plane reference as the bulk basis
        - out_format='' : str, optional (default out_format=in_format)
            Same values as in_format. 'numerator' probably doesn't make much
            sense.

        Returns
        -------
        list
            as many elements as there are beams equivalent to beam, each one is
            a (beam, index, text) tuple, with
            beam : tuple
                formatted as requested in out_format
            index : gl.BeamIndex
                fractional indices
            text : str
                form "i+(j)+..." where i,j,... are the indices of the domains
                that contribute to the beam, and parenthesized ones are extinct
                contributions
        """
        # Potentially useful for the hovering annotations.
        # They need the following info:
        # - beams equivalent to beam
        # - for all beams:
        #   * some type of coordinate, probably best to give out the g vectors?
        #     it will need to be processed later on to get the right rotation
        #     (depending on the view angle), and, consequently the pixel
        #     coordinates in the canvas
        #   * the fractional index + a list of the domains that contribute to
        #     that beam, perhaps already formatted as a string?
        if self.__last_eq is None:
            raise RuntimeError("beams_equivalent_to cannot be used before "
                               "running once self.equivalent_spots at the same "
                               "angular conditions and with the same domains")
        in_format = kwargs.get('in_format', '')
        if in_format not in ('fractional', 'numerator', 'g'):
            raise ValueError("beams_equivalent_to: in_format is mandatory, and "
                             "should be either 'fractional', 'numerator', or "
                             "'g'")
        elif in_format == 'g':
            # Since g_bulk = (h,k)_bulk @ basis_bulk, with basis_bulk the
            # reciprocal-lattice bulk basis, then
            # (h,k)_bulk = g_bulk @ basis_bulk^(-1)
            # The indices are the fractional ones in this case.
            beam_indices = np.dot(beam, np.linalg.inv(self.bulk_basis))
        else:
            beam_indices = beam

        # Now process the beam indices to have only the numerator, as this is
        # what is used in the dictionaries
        if in_format in ('fractional', 'g'):
            beam_indices = [b*self.denominator_for_bulk_beams
                            for b in beam_indices]
            # check that the indices are consistent with the bulk basis
            if any(abs(round(b) - b) > 1e-4 for b in beam_indices):
                raise ValueError("beams_equivalent_to: beam indices are "
                                 "incompatible with the bulk basis")
            beam_indices = [int(round(b)) for b in beam_indices]

        # make the indices a tuple for lookup in the dictionary of equivalent
        # beams
        beam_indices = tuple(beam_indices)
        try: 
            eq_beams = self.__last_eq.equivalent_beams_dict[beam_indices]
        except KeyError:
            raise ValueError(f"beams_equivalent_to: beam {beam} (transformed to"
                             f" {beam_indices}) not found. "
                             "It may not be an acceptable index for the "
                             "system, or may lie outside the energy range")

        # Transform eq_beams into a list, as it's a set in the dictionary, but
        # we need to keep the order for the output
        eq_beams = list(eq_beams)

        # Now get the names of the domains that overlap at each of the beams,
        # also accounting for whether the domain contributes with an extinct
        # spot
        beams, _, domains, extinct_domains = zip(*self.__last_eq.indexed_beams)
        overlapping = []
        for b in eq_beams:
            idx = beams.index(b)
            overlap = [f"({d})" if d in extinct_domains[idx] else f"{d}"
                       for d in domains[idx]]
            overlapping.append("+".join(overlap))

        # Finally process the beam indices according to the format requested
        fractional = [gl.BeamIndex(b,
                                   denominator=self.denominator_for_bulk_beams,
                                   from_numerators=True)
                      for b in eq_beams]

        out_format = kwargs.get('out_format', in_format)
        if out_format == 'fractional':
            out_beams = fractional
        elif out_format == 'numerator':
            out_beams = [b for b in eq_beams]
        elif out_format == 'g':
            out_beams = np.dot(fractional, self.bulk_basis).astype(float)
        else:
            raise ValueError("beams_equivalent_to: invalid out_format "
                             f"{out_format!r}. Expected 'fractional', "
                             "'numerator', or 'g'")
        return list(zip(out_beams, fractional, overlapping))

    # @gl.profile_lines
    def equivalent_spots(self, domains=None, theta=None, phi=None):
        """
        Returns a list of beams grouped by equivalence, given a primary beam
        direction defined by polar and azimuthal angles
        
        Parameters
        ----------
        domains : list of int, or None (default=None)
            only the beams of the domains selected by these indices will be
            output. If None, all the domains are used.
        theta, phi : number or None
            polar and azimuthal directions in degrees of the primary beam. If
            not given or None, those defined at instantiation are used
            theta should be in the [-90, 90] range
            phi is positive counterclockwise, and measured from the x axis
                in the Cartesian reference of the real-space basis of the
                first domain
        """
        # type- and value-check the input, and processing when needed
        if theta is None:
            theta = self.__parameters['beamIncidence'][0]
        if phi is None:
            phi = self.__parameters['beamIncidence'][1]
        if not all(isinstance(angle, (int, float)) for angle in (theta, phi)):
            raise TypeError("LEEDSymmetryDomains: invalid angles. "
                            "Expected a real number.")
        theta, phi = gl.conventional_angles(theta, phi)

        if domains is None:
            domains = range(self.n_domains)
        elif not all(isinstance(dom, int) for dom in domains):
            raise ValueError("Invalid domain index. "
                             "All indices should be integers")
        elif any(dom < 0 or dom >= self.n_domains for dom in domains):
            raise ValueError("Domain index out of range. "
                             "Indices should be between 0 "
                             f"and {self.n_domains}.")

        # select the bare dictionaries
        kwargs = {'domains': domains, 'theta': theta, 'phi': phi}
        domains_dicts = self.__beams_for_primary_angles('equivalent', **kwargs)
        extinct_beams = self.__beams_for_primary_angles('extinct', **kwargs)
        
        # Now, only if the beam is impinging normally, one has to rework a bit
        # the dictionary. In fact, at normal incidence, one should not only
        # consider the spot equivalence within each domain, but also the
        # equivalence among different domains. For example, the two (2x1)
        # domains on a p4xx square lattice produce spots (1/2 0) and (0 1/2),
        # respectively, that are equivalent to one another at normal incidence.
        # This means that, at normal incidence, for each domain i, all the
        # beam "lists" (values) associated with beam j (key) are the same,
        # although beam j is a different tuple, and comprise all the beams
        # for all domains that are equivalent to the j-th beam of each domain
        #
        # Will construct a list of sets, in the same order as the keys,
        # where each set is the union of the stars of the domains. This is
        # possible since all dictionaries were created with the same order
        # of insertion (actually, from the same dictionary)
        if abs(theta) < 1e-4:
            all_stars = [set.union(*stars)
                         for stars in zip(*(d.values() for d in domains_dicts))]
            # now place the unions in the domain dictionaries, preserving the
            # original keys
            domains_dicts = [dict(zip(d.keys(), all_stars))
                             for d in domains_dicts]

        # Finally defer the processing to the LEEDEquivalentBeams class, that
        # needs to know the reciprocal-space bulk basis for processing.
        # Prepare the keyword arguments before:
        kwargs = {'extinct_lists': extinct_beams, 'basis': self.bulk_basis,
                  # currently using 1-based domain indices for display purposes
                  # as it may be easier to look at from the viewpoint of the
                  # final user
                  'domain_ids': [dom + 1 for dom in domains],
                  'superlattice': self.superlattices[0],
                  'angle_key': self.__key_from_angles(theta, phi, domains)}
        self.__last_eq = gl.LEEDEquivalentBeams(domains_dicts, **kwargs)
        return self.__last_eq.indexed_beams

    def __beams_for_primary_angles(self, which_beams, theta=None, phi=None,
                                   domains=None):
        """
        Returns the correct equivalent or extinct beams given the direction
        angles of the primary beam

        Parameters
        ----------
        which_beams : str
            Can be either 'equivalent' (or any contraction down to 'eq') or
            'extinct' (or any contraction down to 'ex'). Determines what is
            returned, the dictionary of equivalent beams or the lists of
            extinct ones, respectively
        theta, phi : number or None (default=None)
            if None, the angle given as a parameter in the constructor is used.
            The azimuthal angle phi is measured with respect to the x axis
            (positive counterclockwise), in the same Cartesian coordinate frame
            as the real-space basis of the first domain. Angles are in degrees.
        domains : list of int or None (default=None)
            select which domains should be output. If None, all domains are used

        Returns
        -------
        list of dict
            if which_beams[:2] == 'eq'
            keys : tuple
                beam j of domain i
            values : set of tuple
                the star of beam j in domain i
        list
            if which_beams[:2] == 'ex'
            
        """
        # check and fix the input if needed
        which_beams = which_beams[:2]
        if which_beams not in ('eq', 'ex'):
            raise ValueError("LEEDSymmetryDomains: invalid parameter for "
                             "selecting which beams to output. Must be 'eq'"
                             "or 'ex' for equivalent or extinct")
        if theta is None:
            theta = self.__parameters['beamIncidence'][0]
        if phi is None:
            phi = self.__parameters['beamIncidence'][1]

        if domains is None:
            domains = range(self.n_domains)

        # now select which keys need to be retrieved
        key = self.__key_from_angles(theta, phi, domains)

        # and from which dictionary
        if which_beams == 'eq':
            beams = self.__equiv_spots_no_superpos
        else:
            beams = self.__extinct_spots

        # retrieve the correct dictionary keys, falling back on 'other' if the
        # angle phi is not one of the special directions
        return [beams[dom].get(key, beams[dom]['other']) for dom in domains]

    def __key_from_angles(self, theta, phi, domains):
        """
        Given the polar and azimuthal angles of incidence of the beam,
        figure out which among the keys 'norm', 'other', or any of the mirror
        direction angles needs to be retrieved
        """
        a_angle = gl.orientation(self[0].real_basis[0], zero_pi=False)
        phi = (phi - a_angle) % 180

        if theta < 1e-4:
            key = 'norm'
        elif abs(phi - round(phi)) > 1e-4:
            # will need to use 'other' for all, as all angles in the
            # dictionaries are integers: 0, 30, 45, 60, ...
            key = 'other'
        else:
            key = str(round(phi))
        if any(key in self.__equiv_spots_no_superpos[dom] for dom in domains):
            return key
        return 'other'

    def __build_domains(self):
        """
        Returns the viperleed.Lattice(s) with self.superlattices matrices
        """
        # Create dummy surface lattice just to get the reciprocal basis
        surf = gl.Lattice(self.__parameters['surfBasis'])

        # and get the maximum screen radius
        max_radius = gl.screen_radius(self.__parameters['eMax'],
                                      self.__parameters['screenAperture'])

        # Then prepare the actual reciprocal lattice of the first domain
        domains = [gl.Lattice(surf.reciprocal_basis, space='reciprocal',
                              group=self.__parameters['surfGroup'],
                              limit=max_radius)]

        # and work out those of the others.
        # The reciprocal-space basis of domain i is:
        #    Bi = Mi^(-T) Bb
        # with Bb the bulk reciprocal basis, and Mi the superlattice matrix
        # of domain i, i.e.,
        #   Mi = M0 Gi
        # where Gi is the operation (of the bulk group) that creates domain i.
        # Thus
        #   Bi = (M0 Gi)^(-T) Bb = Gi^(-T) M0^(-T) Bb = Gi^(-T) B0
        # However, when transforming the basis, Gi needs to be expressed
        # in the coordinate system of the SURFACE, i.e.,
        #   Gi = Breal0 Gi_absolute Breal0^(-1)
        #      = Breal0 Brealb^(-1) Gi_bulk Brealb Breal0^(-1)
        #      = superlattice_0 Gi_bulk superlattice_0^(-1)
        superlattice = self.superlattices[0]
        inv_superlattice = np.linalg.inv(superlattice)
        ops_t = [np.linalg.multi_dot((superlattice,
                                      op,
                                      inv_superlattice))
                 for op in self.operations]
        ops = [np.linalg.inv(op).T for op in ops_t[1:]]  # skip E (1st domain)
        domains.extend([domains[0].transform(op, as_copy=True)
                        for op in ops])
        return domains

    def __find_domain_operations(self):
        """
        Finds the symmetry operations of the bulk lattice that give distinct
        domains in LEED.
        NB: this function gives different results than what one would get from
        LEEDPat, as we're interested also in the symmetry relations between the
        intensities of the spots, while LEEDPat cares only about the presence
        or not of any spot.

        Returns
        -------
        list of operations

        -------
        """
        # The current version is based on the concept of co-sets of a group.
        # Given a group G and a subgroup H, the left co-set of H with respect
        # to the group operation g of G is
        #           gH = {g*h : h in H}.
        # Additionally, we take into account that each element of G is found in
        # exactly only one co-set (e.g., H is the identity co-set).
        #
        # The operations g_i that generate distinct co-sets are those that will
        # generate distinct lattices

        # def __array_2_tuple(arr):
            # """
            # Convenience function that is used to convert a 2x2 array to a 2x2
            # tuple
            # """
            # return tuple(map(tuple, arr))

        # 1) project the operations from the surface group to the bulk; also
        #    rounding to integers
        project_to_bulk = np.linalg.inv(self.__parameters['SUPERLATTICE'])
        surf_ops = tuple(
            op.round().astype(int)
            for op
            in self.groups['surf'].transform(project_to_bulk)
            )

        # 2) keep track of the operations of the bulk group that are already in
        #    a co-set, and of the group operations giving distinct co-sets
        coset_ops = set()
        _ops = []

        # 3) run through the bulk operations, adding the operations of new
        #    co-sets to the set above
        bulk_ops = self.groups['bulk'].operations(include_3d=True)

        for bulk_op in bulk_ops:
            if bulk_op in coset_ops:
                continue
            _ops.append(bulk_op)
            coset = set(gl.two_by_two_array_to_tuple(np.dot(bulk_op, surf_op))
                        for surf_op in surf_ops)
            coset_ops.update(coset)

        return _ops

    def __get_spot_equivalence(self):
        """
        Determines which of the LEED beams of the domains belong to the same
        equivalence class under the symmetry operations of the domain itself,
        i.e., it determines the 'star' of each of the beams.
        Returns a list of dictionary of dictionaries, one list element per each
        domain. Each key in the outer dictionary is either 'norm', the angle
        in degrees of a mirror direction (with respect to the first real-space
        lattice vector of the first domain BETTER WRT BULK??), or 'other'.
        These keys are useful for selecting equivalent spots when the primary
        beam direction is, respectively, normal to the surface, contained in a
        mirror plane, or generic.

        Returns
        -------
        list of dicts, each one
            keys : str
                'norm' and 'other' keys are always present. In addition, there
                are as many keys as self.special_directions, each one is the
                angle in degrees between the mirror direction and the first
                real-space basis vector. The angles can only be a subset of:
                '0', '30', '45', '60', '90', '120', '135', '150'
            values : dict
                keys : tuple
                    one key per each index in self.hk
                values: set of tuples
                    all beam indices symmetry-equivalent to key, including key,
                    i.e., what's known as 'the star' of the beam in key
        """
        # Will create several dictionaries, labeled depending on the beam
        # incidence directions. There always will be one for 'norm' (normal)
        # incidence, as many as there are special_directions, each of
        # which will be indexed by the angle between the special direction
        # and the first basis vector, and one for 'other'. This makes it fast
        # to later select which of the equivalence relations to choose.
        # Since the beam equivalence dictionaries are the same for all domains,
        # work on the first one only , i.e., self[0], in integer-index notation,
        # and later on translate stuff to the other domains
        all_operations = self[0].group.operations()
        a_angle = gl.orientation(self[0].real_basis[0], zero_pi=False)
        
        # Prepare the dict keys:
        labels = ['norm', 'other']
        operations = [all_operations,        # Normal incidence -> all
                      (all_operations[0],)]  # Generic -> only identity
        # Now go through the mirror directions for the other keys.
        # For a beam along the mirror, the only operations that remain are
        # identity and the mirror itself
        for op, direction in zip(all_operations, self[0].special_directions):
            if direction is None:  # the operation is a rotation
                continue
            phi = gl.orientation(direction, zero_pi=False)
            labels.append(round((phi - a_angle) % 180))
            operations.append((all_operations[0], op))
        
        # Now use the labels and the operations found above to find the star of
        # each beam of the first domain in all cases
        beams_dict_first = {}
        for label, ops in zip(labels, operations):
            # Transform all the beams with all the operations
            hk_transformed = np.einsum('ilm,mj->jil', ops, self[0].hk.T)

            # convert all beams to tuples so they can be keys for the dictionary
            hk_tuples = gl.two_by_n_array_to_tuples(self[0].hk)
            
            # and convert to tuples also the transformed beams,
            # also using a set to remove duplicates
            beams_dict_first[label] = {
                k: set(gl.two_by_n_array_to_tuples(v, axis=1))
                for k, v in zip(hk_tuples, hk_transformed)
            }
        return self.__beams_dict_to_other_domains(beams_dict_first)
    
    def __get_extinct(self):
        """
        For each domain, finds the beams extinct due to glide symmetry, and
        populates dictionaries with the in-plane directions of the primary
        beam that makes them extinct
        """
        group = self.groups['surf'].group
        # Like for __get_spot_equivalence, do the calculation on the first
        # domain, then extend to the others by relabeling
        if 'g' not in group:
            # no glide, no need to bother
            return [{'norm': [], 'other': []}]*len(self)
        
        # get boolean arrays that pick the extinct spots along the [1 0] and
        # [0 1] directions (for the first domain, but they're the same for
        # the others)
        extinct_10 = (self[0].hk[:, 1] == 0) & (self[0].hk[:, 0] % 2 == 1)
        extinct_01 = (self[0].hk[:, 0] == 0) & (self[0].hk[:, 1] % 2 == 1)
        
        beams_dict_first = {'other': []}
        # pg, pmg, p4g are the only options, and exist only for square
        # or rectangular lattices, so the angles are only 0 and/or 90
        if '[1 0]' in group:
            beams_dict_first['norm'] = self[0].hk[extinct_10]
            beams_dict_first[0] = self[0].hk[extinct_10]
        elif '[0 1]' in group:
            beams_dict_first['norm'] = self[0].hk[extinct_01]
            beams_dict_first[90] = self[0].hk[extinct_01]
        else:
            beams_dict_first['norm'] = self[0].hk[extinct_10 | extinct_01]
            beams_dict_first[0] = self[0].hk[extinct_10]
            beams_dict_first[90] = self[0].hk[extinct_01]
        return self.__beams_dict_to_other_domains(beams_dict_first)

    def __beams_dict_to_other_domains(self, beams_dict_first):
        """
        Converts a beams dictionary from the first domain into a list of beams
        dictionaries for each domain. This is done by determining angular
        offsets between domains, restricting to 0...pi, and building the offset
        dictionaries. At the same time, we also convert the integer domain
        indices into fractional bulk indices via the superlattice matrices
        
        Parameters
        ----------
        beams_dict_first : dict
            keys : str or int
                can be 'norm', 'other' or integers with the angles of some
                special directions
            values : list-like or dict
                can either be a list of integer hk values of the first domain,
                or a dictionary of the form {hk0: {hk0, hk1, ...}, ...}

        Returns
        -------
        list of dict
            one list element for each of the domains, each dict has the same
            format as the one passed as parameters, except for the fact that
            keys are all strings, and angles are adjusted to all refer to the
            first lattice vector of the first domain. The contents of the
            dictionaries are the same as passed, but converted to be the
            numerator of the fractional indices with respect to the bulk
        """
        offsets = self.angular_offsets(zero_pi=True)
        beams_dicts_all_domains = []
        mu = self.denominator_for_bulk_beams
        for offset, superlattice in zip(offsets, self.superlattices):
            ddict = {}
            for label, beams in beams_dict_first.items():
                # set up the transformation matrix that gives fractional
                # indices in the bulk basis from the integer indices of
                # each domain
                transform = np.linalg.inv(superlattice).T

                # it is more convenient (faster) to still store only the
                # numerator of the indices as tuples, as this limits numerical
                # errors. The fractional indices can be retrieved when needed
                # by dividing by self.denominator_for_bulk_beams
                transform = (mu*transform).round().astype(int)
                bulk_beam_dict = self.__beams_to_bulk(beams, transform)
                if label in ['norm', 'other']:
                    ddict[label] = bulk_beam_dict
                else:
                    ddict[str(round((label + offset) % 180))] = bulk_beam_dict
            beams_dicts_all_domains.append(ddict.copy())
        return beams_dicts_all_domains

    def __beams_to_bulk(self, beams, transform):
        """
        Take either a list [hk0, hk1, ...] or a dict of the form
        {hk0: {hk1, hk2, ...}} with hk integer indices relative to one domain,
        and transform them to a list or a dict, respectively, with fractional
        indices all over.
        """
        # In case we have a dict, extract beams from the keys (one each key)
        # as well as those from the values, i.e., beams of the star of key
        if isinstance(beams, dict):
            k, v = zip(*beams.items())
            # pack them into a flat list, as this allows to run a dot product
            # in numpy that makes calculations much faster, despite the
            # additional work to unpack/re-pack stuff
            flattened_beams = list(itertools.chain(k, *v))
        elif len(beams) == 0:
            return beams
        else:
            flattened_beams = beams

        # transform the beams. For some reason, einsum is faster than dot
        transformed_beams = np.einsum('ij,jk', flattened_beams, transform)

        # Now transform them back to a list of tuples
        transformed_beams = list(gl.two_by_n_array_to_tuples(transformed_beams,
                                                             axis=1))
        if not isinstance(beams, dict):
            return transformed_beams

        # When a dict was passed, we need to split the 
        # beams once again into keys and values
        k_transf = transformed_beams[:len(k)]
        v_transf = transformed_beams[len(k):]

        # and figure out how long each of the values used to be, so they can
        # be split again correctly among the various keys
        split_idx = list(itertools.accumulate([0, *(len(i) for i in v)]))

        # Do the splitting, also converting again to set
        v_transf = [set(v_transf[i:j])
                    for i, j in zip(split_idx, split_idx[1:])]

        # and recreate the dictionary
        return dict(zip(k_transf, v_transf))
