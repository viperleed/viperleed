"""Module symdomains of viperleed.guilib.leedsim.classes.

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


# @gl.profile_lines
# @gl.exec_time
def transform_beam_indices(indices, transform):
    """Transform indices to another basis.

    This function can be used to convert a list (or dictionary) of
    beam indices, e.g., expressed with respect to the basis of one
    domain into a list (or dictionary, respectively) of beam indices
    expressed with respect to a different basis, e.g., the one of
    the bulk.

    Parameters
    ----------
    indices : iterable or dict
        When an iterable it is assumed to have the form
        [hk0, hk1, ...]; when a dict, {hk0: [hk1, hk2, ...]}, with
        each value being again an iterable. Each of the beams should
        be a 2-element iterable. It is not safe to pass an object
        that behaves like a generator in any of the fields.
    transform : iterable
        2x2 basis-change transformation matrix that will be applied
        to each of the indices "on the right", i.e., representing an
        index as a row vector [h1, k1] the transformed index will be
        [h', k'] = [h, k] @ transform.

    Returns
    -------
    list or dict
        A list of tuples is returned when the input was an iterable.
        A dictionary of the form {hk': set{hk'}} is returned when
        a dict was passed.
    """
    if len(indices) == 0:
        return indices

    # In case we have a dict, extract indices from the keys (one
    # each key) as well as those from the values.
    if isinstance(indices, dict):
        keys, values = zip(*indices.items())
        # Pack them into a flat list, as this allows to run a dot
        # product in numpy that makes calculations much faster,
        # despite the additional work to unpack/re-pack stuff
        flattened_indices = list(itertools.chain(keys, *values))
    else:
        flattened_indices = indices

    # Transform the indices. einsum is faster than dot
    transformed_indices = np.einsum('ij,jk', flattened_indices, transform)

    # Now transform them back to a list of tuples
    transformed_indices = list(gl.two_by_n_array_to_tuples(transformed_indices,
                                                           axis=1))
    if not isinstance(indices, dict):
        return transformed_indices

    # When a dict was passed, we need to split the
    # beams once again into keys and values
    k_transf = transformed_indices[:len(keys)]
    v_transf = transformed_indices[len(keys):]

    # and figure out how long each of the values used to be, so
    # they can be split again correctly among the various keys
    split_idx = list(itertools.accumulate([0, *(len(i) for i in values)]))

    # Do the splitting, also converting again to set
    v_transf = [set(v_transf[i:j])
                for i, j in zip(split_idx, split_idx[1:])]

    # and recreate the dictionary
    return dict(zip(k_transf, v_transf))


class LEEDSymmetryDomains(Sequence):
    """Domains generated by bulk symmetry operations from a superlattice.

    Collection of LEED domains generated from one superlattice matrix
    and the symmetry operations of a bulk group.
    LEEDSymmetryDomains is an immutable sequence (i.e., tuple-like).

    Calling LEEDSymmetryDomains[i] returns a viperleed.Lattice of the
    domain
    """

    @gl.exec_time
    # @gl.profile_lines
    def __init__(self, leed_parameters):
        """Initialize LEEDSymmetryDomains instance.

        Parameters
        ----------
        leed_parameters : dict, ConfigParser or LEEDParameters
            Parameters defining the collections of domains
            related by bulk symmetry. Will be fed to
            LEEDParameters

        Returns
        -------
        None.
        """
        self.__parameters = gl.LEEDParameters(leed_parameters)

        # __const_attributes contains attributes of this instance
        # that are calculated only once (either during this
        # initialization or upon first call of one of the properties).
        # The keys are:
        # - 'bulk': a minimal bulk Lattice used to get the basis,
        #           the cell shape, and similar attributes.
        #           Read via self.bulk.
        # - 'operations': a list of the bulk symmetry operations
        #           generating distinct domains.
        #           Read via self.operations.
        # - 'superlattices': a numpy.ndarray of the superlattice
        #           matrices that generate distinct domains.
        #           Read via self.superlattices.
        # - 'denominator': the ratio of the areas of the superlattice
        #           cell with respect to the one of the bulk.  This
        #           is necessary because the return value of
        #           self.equivalent_spots() has beams expressed
        #           as tuples, and including ONLY THE NUMERATOR of
        #           the beam indices with respect to the bulk.  This
        #           makes calculations faster.
        #           Read via self.denominator_for_bulk_beams.
        self.__const_attributes = {'bulk': None,
                                   'operations': [],
                                   'superlattices': [],
                                   'denominator': None}

        # self.__domains is the underlying list that is accessed when
        # issuing self[index]. It's a list of gl.Lattice instances
        self.__domains = self.__build_domains()

        # self.__equiv_spots_no_superpos is a list of dictionaries of
        # dictionaries, one element per domain, ordered as in
        # self.__domains (i.e., in self).
        # The outer dictionary is labeled by 'norm', 'other', or an
        # azimuthal angle. For each of the angle settings, the inner
        # dictionary is in the form {beam_i: set(beams_i)} with beams
        # of domain i and their equivalents within the domain at the
        # specified primary beam angle(s).
        # Used only internally, cannot be accessed
        self.__equiv_spots_no_superpos = self.__get_spot_equivalence()

        # self.__extinct_spots is a list of dictionaries, one element
        # per domain, same order as in self.__domains (i.e., in self).
        # Each dictionary has keys 'norm', 'other', or an azimuthal
        # angle, and as values a numpy array that lists the extinct
        # beams of the domain at the specified beam angle(s).
        # Used only internally, cannot be accessed
        self.__extinct_spots = self.__get_extinct()

        # self.__equiv_spots_no_superpos and self.__extinct_spots are
        # used in the public method self.equivalent_spots(domains,
        # theta, phi) for creating the appropriate
        # gl.LEEDEquivalentBeams instance (unless it's already cached),
        # whose reference is stored into self.__last_eq.
        # Can be accessed (read-only) with self.eq_beams_last_config
        self.__last_eq = None

    def __repr__(self):
        """Return string representation of self."""
        txt = (f"{self[0].cell_shape} "
               + f"viperleed.LEEDSymmetryDomains({self.__parameters})")
        return txt

    def __len__(self):
        """Return length of self."""
        return len(self.__domains)

    def __getitem__(self, elem):
        """Item getter.

        Notice that if a slice is given, the result is NOT
        type-preserving, i.e., it will not be a LEEDSymmetryDomains
        but a list of Lattice(s)
        """
        return self.__domains[elem]

    @property
    def bulk_basis(self):
        """Reciprocal-space basis of the bulk lattice."""
        if self.__const_attributes['bulk']:
            return self.bulk.basis
        return np.dot(self.superlattices[0].T, self[0].basis).round(10)

    @property
    def bulk(self):
        """Minimal bulk Lattice. Do not use for plotting."""
        if not self.__const_attributes['bulk']:
            bulk_lattice = gl.Lattice(self.bulk_basis,
                                      group=self.groups['bulk'],
                                      space='reciprocal')
            self.__const_attributes['bulk'] = bulk_lattice
        return self.__const_attributes['bulk']

    @property
    def domain_ids(self):
        """Return an up-to-date list of domain 'ids'.

        TODO: Consider if it makes sense to return something
        else than just progressive integers. In which case
        one has to modify the way domain ids are fetched when
        constructing beam equivalences

        Returns
        -------
        list of int
        """
        return list(range(len(self)))

    @property
    def eq_beams_last_config(self):
        """Return the LEEDEquivalentBeams of the last configuration."""
        return self.__last_eq

    @property
    def extinct(self):
        """Return beams that are extinct in the last configuration.

        If the configuration (i.e., domains selected, beam-incidence
        angles) is not up to date, call .equivalent_spots() before!

        Returns
        -------
        list
            List of extinct beams as tuples of numerators

        Raises
        ------
        RuntimeError
            When this property is called before the .equivalent_spots()
            method has been invoked at least once.
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("LEEDSymmetryDomains: Need to "
                               "call .equivalent_spots() once "
                               "before .extinct is available")
        return self.eq_beams_last_config.extinct

    @property
    def groups(self):
        """Return a dictionary of the symmetry groups.

        Returns
        -------
        dict
            keys : {'surf', 'bulk'}
            values : PlaneGroup
        """
        return {'surf': self.__parameters['surfGroup'],
                'bulk': self.__parameters['bulkGroup']}

    @property
    def indexed_beams(self):
        """Return a list of indexed beams.

        If the configuration (i.e., domains selected, beam-incidence
        angles) is not up to date, call .equivalent_spots() before!

        Returns
        -------
        list
            Each entry is a tuple of the form
            (beam, group_idx, overlapping_domains, extinct_domains)

            beam : tuple
                Numerators of the indices of the beam
            group_idx : int
                Progressive index defining equivalent beams. Negative
                for fully extinct.
            overlapping_domains : list
                Domain ids for domains contributing
            extinct_domains : list
                Domain ids for domains that contribute with extinct
                beams

        Raises
        ------
        RuntimeError
            When this property is called before the equivalent_spots()
            method has been invoked at least once.
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("LEEDSymmetryDomains: Need to call "
                               ".equivalent_spots() once before "
                               ".indexed_beams is available")
        return self.eq_beams_last_config.indexed_beams

    @property
    def n_domains(self):
        """Return number of domains induced by bulk symmetry."""
        return len(self)

    @property
    def last_domains_used(self):
        """Return the domain indices used in the last configuration.

        Returns
        -------
        tuple

        Raises
        ------
        RuntimeError
            When this property is accessed before the
            equivalent_spots() method has ever been called.
            There is no "last configuration used" in this case.
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("Need to call .equivalent_spots() once "
                               "before .last_domains_used is available")
        return self.eq_beams_last_config.hash_dict['domain_ids']

    @property
    def operations(self):
        """Return a list of the bulk operations producing domains.

        Returns
        -------
        list of 2x2 tuples
        """
        if not self.__const_attributes['operations']:
            operations = self.__find_domain_operations()
            self.__const_attributes['operations'] = operations
        return self.__const_attributes['operations']

    @property
    def parameters(self):
        """Return the LEEDParameters that underlies self."""
        return self.__parameters

    @property
    def superlattices(self):
        """Return the superlattice matrices of domains.

        Returns
        -------
        numpy.ndarray
            Array of superlattice matrices generating the
            symmetry-related domains. The first element is
            the same as given in the parameters during
            initialization.
        """
        if len(self.__const_attributes['superlattices']) == 0:
            first_superlattice = self.__parameters['SUPERLATTICE']
            superlattices = np.einsum('ij,mjk->mik', first_superlattice,
                                      self.operations)
            self.__const_attributes['superlattices'] = superlattices
        return self.__const_attributes['superlattices']

    @property
    def g_vectors(self):
        """Return a list of reciprocal-space vectors for each domain.

        Returns
        -------
        list
            g_vectors[i] is a list of the lattice vectors of domain i
            in the same coordinate system as the bulk basis.
        """
        return [dom.lattice for dom in self]

    @property
    def denominator_for_bulk_beams(self):
        """Return the denominator to use for fractional indices."""
        if not self.__const_attributes['denominator']:
            den = abs(round(np.linalg.det(self.superlattices[0])))
            self.__const_attributes['denominator'] = den
        return self.__const_attributes['denominator']

    def angular_offsets(self, zero_pi=False):
        """Angular offsets between all domains and the first one.

        Returns
        -------
        list
            Angles between the first lattice vectors of
            all domains and those of the first domain
        """
        first_dom_angle = gl.orientation(self[0].real_basis[0], zero_pi)
        return [gl.orientation(dom.real_basis[0], zero_pi) - first_dom_angle
                for dom in self]

    def domains_containing_beam(self, beam):
        """Return the list of domains that contain beam."""
        # Pick one of the dicts from each domain to check
        # whether beam is in the domain. Easiest is to use
        # 'norm', since it's always there.
        beam_lists = [d['norm'].keys() for d in self.__equiv_spots_no_superpos]
        return tuple(d_id
                     for d_id, beam_list in zip(self.domain_ids, beam_lists)
                     if beam in beam_list)

    # @gl.profile_lines
    @gl.exec_time
    def equivalent_spots(self, domains=None, theta=None, phi=None):
        """Return beams grouped by equivalence, at given incidence angles.

        Parameters
        ----------
        domains : list, int, or None (default=None)
            only the beams of the domains selected by these indices
            will be output.  If None, all the domains are used. If a
            single int is passed, it is taken as the positional index
            of the only domain to be processed
        theta, phi : number or None
            polar and azimuthal directions in degrees of the primary
            beam.  If not given or None, those defined at instantiation
            are used theta should be in the [-90, 90] range.  phi is
            positive counterclockwise, and measured from the x axis in
            the Cartesian reference of the real-space basis of the
            first domain

        Returns
        -------
        list
            Each element is:
            (beam,
             beam_group_idx,
             list of domain indices of domains contributing to beam,
             list of domain indices of domains contributing as extinct)
        """
        # type- and value-check the input, and process when needed
        if theta is None:
            theta = self.__parameters['beamIncidence'][0]
        if phi is None:
            phi = self.__parameters['beamIncidence'][1]
        theta, phi = gl.conventional_angles(float(theta), float(phi))

        if domains is None:
            domains = range(self.n_domains)
        elif isinstance(domains, int):
            domains = [domains]
        elif not hasattr(domains, '__iter__'):
            raise ValueError("Invalid domain indices. Expected an "
                             "iterable, an integer, or None")
        if not all(isinstance(dom, int) for dom in domains):
            raise ValueError("Invalid domain index. "
                             "All indices should be integers")
        if any(dom < 0 or dom >= self.n_domains for dom in domains):
            raise ValueError("Domain index out of range. "
                             "Indices should be between 0 "
                             f"and {self.n_domains}.")

        # select the bare dictionaries
        kwargs = {'domains': domains, 'theta': theta, 'phi': phi}
        domains_dicts = self.__beams_for_primary_angles('equivalent', **kwargs)
        extinct_beams = self.__beams_for_primary_angles('extinct', **kwargs)

        # Now, only if the beam is impinging normally, one has to
        # rework a bit the dictionary.  In fact, at normal incidence,
        # one should not only consider the spot equivalence within each
        # domain, but also the equivalence among different domains.
        # For example, the two (2x1) domains on a p4xx square lattice
        # produce spots (1/2 0) and (0 1/2), respectively, that are
        # equivalent to one another at normal incidence.  This means
        # that, at normal incidence, for each domain i, all the beam
        # "lists" (values) associated with beam j (key) are the same,
        # although beam j is a different tuple, and comprise all the
        # beams for all domains that are equivalent to the j-th beam
        # of each domain.
        #
        # Will construct a list of sets, in the same order as the keys,
        # where each set is the union of the stars of the domains.
        # This is possible since all dictionaries were created with the
        # same order of insertion (actually, from the same dictionary)
        if theta < 1e-4:
            all_stars = [set.union(*stars)
                         for stars in zip(*(d.values()
                                            for d in domains_dicts))]
            # now place the unions in the domain dictionaries,
            # preserving the original keys
            domains_dicts = [dict(zip(d.keys(), all_stars))
                             for d in domains_dicts]

        # Finally defer the processing to the LEEDEquivalentBeams
        # class, that needs to know the reciprocal-space bulk basis for
        # processing.  Prepare the keyword arguments before:
        #
        #   TODO: Perhaps I can have
        #      'domain_ids' : [self.domain_ids[i] for i in domains]
        #   this way one can chose later to use different ids, e.g.,
        #   1-based, by simply changing self.domain_ids.
        kwargs = {'extinct_lists': extinct_beams, 'basis': self.bulk_basis,
                  'domain_ids': domains,  # 0-based indices
                  'superlattice': self.superlattices[0],
                  'angle_key': self.__key_from_angles(theta, phi, domains),
                  'caller': self}
        self.__last_eq = gl.LEEDEquivalentBeams(domains_dicts, **kwargs)
        return self.__last_eq.indexed_beams

    def set_beam_incidence(self, theta=None, phi=None):
        """Edit the beam incidence.

        Also triggers an update of the beam equivalence if
        self.equivalent_spots() was already called.

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
        if theta is None and phi is None:
            return
        if theta is None:
            theta = self.__parameters[0]['beamIncidence'][0]
        elif phi is None:
            phi = self.__parameters[0]['beamIncidence'][1]
        theta, phi = gl.conventional_angles(theta, phi)

        # Update the underlying parameters
        self.__parameters['beamIncidence'] = (theta, phi)

        # And triggers an update of the equivalent spots
        # if .equivalent_spots() was called before.
        if self.__last_eq is not None:
            self.equivalent_spots(self.last_domains_used)

    def __beams_for_primary_angles(self, which_beams, theta=None, phi=None,
                                   domains=None):
        """Return the list of beams appropriate for the angles given.

        Parameters
        ----------
        which_beams : {'eq'(uivalent), 'ex'(tinct)}
            Determines what is returned, i.e., the dictionary of
            equivalent beams or the lists of extinct ones, respectively
        theta, phi : number or None (default=None)
            if None, the angle given as a parameter in the constructor
            is used.  The azimuthal angle phi is measured with respect
            to the x axis (positive counterclockwise), in the same
            Cartesian coordinate frame as the real-space basis of the
            first domain. Angles are in degrees.
        domains : list of int or None (default=None)
            select which domains should be output. If None, all domains
            are used

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

        # retrieve the correct dictionary keys, falling back on 'other'
        # if the angle phi is not one of the special directions
        return [beams[dom].get(key, beams[dom]['other']) for dom in domains]

    def __key_from_angles(self, theta, phi, domains):
        """Pick the right dict key given the incidence angles.

        Given the polar and azimuthal angles of incidence of the beam,
        figure out which among the keys 'norm', 'other', or any of the
        mirror direction angles needs to be retrieved.

        Returns
        -------
        str : {'norm', 'other', 'azimuthal angle of mirrors/glides'}
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
        """Return the viperleed.Lattice(s) with self.superlattices matrices."""
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
        # with Bb the bulk reciprocal basis, and Mi the superlattice
        # matrix of domain i, i.e.,
        #   Mi = M0 Gi,
        # where Gi is the operation (of the bulk group) that creates
        # domain i.  Thus
        #   Bi = (M0 Gi)^(-T) Bb = Gi^(-T) M0^(-T) Bb = Gi^(-T) B0.
        # However, when transforming the basis, Gi needs to be
        # expressed in the coordinate system of the SURFACE, i.e.,
        #   Gi = Breal0 Gi_absolute Breal0^(-1)
        #      = Breal0 Brealb^(-1) Gi_bulk Brealb Breal0^(-1)
        #      = superlattice_0 Gi_bulk superlattice_0^(-1).
        superlattice = self.superlattices[0]
        inv_superlattice = np.linalg.inv(superlattice)
        ops_t = [np.linalg.multi_dot((superlattice,
                                      op,
                                      inv_superlattice))
                 for op in self.operations]
        # skip the identity, i.e., the 1st domain, as it is already
        # included above
        ops = [np.linalg.inv(op).T for op in ops_t[1:]]
        domains.extend([domains[0].transform(op, as_copy=True)
                        for op in ops])
        return domains

    def __find_domain_operations(self):
        """Find bulk symmetry operations giving distinct domains.

        Notice that this function gives different results than what
        one would get from LEEDPat, as we're interested also in the
        symmetry relations between the intensities of the spots, while
        LEEDPat cares only about the presence or not of any spot.

        Returns
        -------
        list of operations

        -------
        """
        # The current version is based on the concept of co-sets of a
        # group.  Given a group G and a subgroup H, the left co-set of
        # H with respect to the group operation g of G is
        #           gH = {g*h : h in H}.
        #   Additionally, we take into account that each element of G
        # is found in exactly only one co-set (e.g., H is the identity
        # co-set).
        #
        # The operations g_i that generate distinct co-sets are those
        # that will generate distinct lattices

        # 1) Project the operations from the surface group to the bulk;
        #    also rounding to integers
        project_to_bulk = np.linalg.inv(self.__parameters['SUPERLATTICE'])
        surf_ops = tuple(
            op.round().astype(int)
            for op
            in self.groups['surf'].transform(project_to_bulk)
            )

        # 2) Keep track of the operations of the bulk group that are
        #    already in a co-set, and of the group operations giving
        #    distinct co-sets
        coset_ops = set()
        _ops = []

        # 3) Run through the bulk operations, adding the operations of
        #    new co-sets to the set above
        bulk_ops = self.groups['bulk'].operations(include_3d=True)

        for bulk_op in bulk_ops:
            if bulk_op in coset_ops:
                continue
            _ops.append(bulk_op)
            coset = set(gl.two_by_two_array_to_tuple(np.dot(bulk_op, surf_op))
                        for surf_op in surf_ops)
            coset_ops.update(coset)

        return _ops

    # @gl.exec_time
    # @gl.profile_lines  # 37% beams_dict[label]; 63% return
    def __get_spot_equivalence(self):
        """Build equivalence classes for LEED beams.

        Determine which of the LEED beams of the domains belong to the
        same equivalence class under the symmetry operations of the
        domain itself, i.e., determine the 'star' of each of the beams.

        Returns a list of dictionary of dictionaries, one list
        element per each domain. Each key in the outer dictionary is
        either 'norm', the angle in degrees of a mirror direction (with
        respect to the first real-space lattice vector of the first
        domain), or 'other'. These keys are useful for selecting
        equivalent spots when the primary beam direction is,
        respectively, normal to the surface, contained in a mirror
        plane, or generic.

        Returns
        -------
        list of dicts, each one
            keys : str
                'norm' and 'other' keys are always present.  In
                addition, there are as many keys as there are
                self.special_directions, each one is the angle in
                degrees between the mirror direction and the first
                real-space basis vector. The angles can only be a
                subset of: {'0', '30', '45', '60', '90', '120',
                '135', '150'}
            values : dict
                keys : tuple
                    one key per each index in self.hk
                values: set of tuples
                    all beam indices symmetry-equivalent to key,
                    including key, i.e., what's known as 'the star'
                    of the beam in key
        """
        # Will create several dictionaries, labeled depending on the
        # beam incidence directions. There always will be one for
        # 'norm' (normal) incidence, as many as there are
        # special_directions, each of which will be indexed by the
        # angle between the special direction and the first basis
        # vector, and one for 'other'. This makes it fast to later
        # select which of the equivalence relations to choose.
        # Since the beam equivalence dictionaries are the same for all
        # domains, work on the first one only , i.e., self[0], in
        # integer-index notation, and later on translate stuff to the
        # other domains
        all_operations = self[0].group.operations()
        a_angle = gl.orientation(self[0].real_basis[0], zero_pi=False)

        # Prepare the dict keys:
        labels = ['norm', 'other']
        operations = [all_operations,        # Normal incidence -> all
                      (all_operations[0],)]  # Generic -> only identity
        # Now go through the mirror directions for the other keys.
        # For a beam along the mirror, the only operations that remain
        # are identity and the mirror itself
        for operation, direction in zip(all_operations,
                                        self[0].special_directions):
            if direction is None:  # the operation is a rotation
                continue
            phi = gl.orientation(direction, zero_pi=False)
            labels.append(round((phi - a_angle) % 180))
            operations.append((all_operations[0], operation))

        # Now use the labels and the operations found above to find
        # the star of each beam of the first domain in all cases
        beams_dict_first = {}
        for label, ops in zip(labels, operations):
            # Transform all the beams with all the operations
            hk_transformed = np.einsum('ilm,mj->jil', ops, self[0].hk.T)

            # convert all beams to tuples to use them as keys
            hk_tuples = gl.two_by_n_array_to_tuples(self[0].hk)

            # and convert to tuples also the transformed beams,
            # also using a set to remove duplicates
            beams_dict_first[label] = {
                k: set(gl.two_by_n_array_to_tuples(v, axis=1))
                for k, v in zip(hk_tuples, hk_transformed)
            }
        return self.__beams_dict_to_other_domains(beams_dict_first)

    def __get_extinct(self):
        """Find beams extinct due to glide symmetry.

        For each domain, find the beams extinct due to glide symmetry,
        and populate dictionaries with the in-plane directions of the
        primary beam that makes them extinct
        """
        group = self.groups['surf'].group
        # Like for __get_spot_equivalence, do the calculation on the
        # first domain, then extend to the others by relabeling
        if 'g' not in group:
            # no glide, no need to bother
            return [{'norm': [], 'other': []}]*len(self)

        # Get boolean arrays that pick the extinct spots along the
        # [1 0] and [0 1] directions (for the first domain, but they're
        # the same for the others)
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

    # @gl.profile_lines
    def __beams_dict_to_other_domains(self, beams_dict_first):
        """Translate beam dictionary of first domain to all domains.

        Convert a beams dictionary from the first domain into a list
        of beams dictionaries for each domain.  This is done by
        determining angular offsets between domains, restricting to
        [0...pi), and building the offset dictionaries.  At the same
        time, we also convert the integer domain indices into
        fractional bulk indices via the superlattice matrices

        Parameters
        ----------
        beams_dict_first : dict
            keys : str or int
                can be 'norm', 'other' or integers with the
                angles of some special directions
            values : list-like or dict
                can either be a list of integer hk values of
                the first domain, or a dictionary of the form
                {hk0: {hk0, hk1, ...}, ...}

        Returns
        -------
        list of dict
            one list element for each of the domains, each dict has the
            same format as the one passed as parameters, except for the
            fact that keys are all strings, and angles are adjusted to
            all refer to the first lattice vector of the first domain.
            The contents of the dictionaries are the same as passed,
            but converted to be the numerator of the fractional indices
            with respect to the bulk
        """
        offsets = self.angular_offsets(zero_pi=True)
        beams_dicts_all_domains = []
        den = self.denominator_for_bulk_beams
        for offset, superlattice in zip(offsets, self.superlattices):
            ddict = {}
            for label, beams in beams_dict_first.items():
                # Set up the transformation matrix that gives
                # fractional indices in the bulk basis from the
                # integer indices of each domain
                transform = np.linalg.inv(superlattice).T

                # It is more convenient (faster) to still store
                # only the numerator of the indices as tuples, as
                # this limits numerical errors.  The fractional
                # indices can be retrieved when needed by dividing
                # by self.denominator_for_bulk_beams
                transform = (den*transform).round().astype(int)
                bulk_beam_dict = transform_beam_indices(beams, transform)
                if label in ['norm', 'other']:
                    ddict[label] = bulk_beam_dict
                else:
                    ddict[str(round((label + offset) % 180))] = bulk_beam_dict
            beams_dicts_all_domains.append(ddict.copy())
        return beams_dicts_all_domains
