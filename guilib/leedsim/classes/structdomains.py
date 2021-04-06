"""Module structdomains of viperleed.guilib.leedsim.classes.

======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes.structdomains ***

Defines the LEEDStructuralDomains class, used for constructing a LEED
pattern with one or more structures, i.e., domains characterized by
superlattice matrices that cannot be reduced to one another by one of
the symmetry operations of the bulk. Those that are related by bulk
symmetry operations are handled by the LEEDSymmetryDomains class.  The
LEEDStructuralDomains class is essentially a list of LEEDSymmetryDomains
instances.

LEEDStructuralDomains used to be defined within leedpattern.py

Author: Michele Riva
Created: 2021-03-21
"""

from collections.abc import MutableSequence
import copy

from viperleed import guilib as gl


class LEEDStructuralDomains(MutableSequence):
    """
    Represent several structural domains in a LEED pattern.

    This is a container class that can be used to represent a list of
    distinct structural domains contributing to a LEED pattern. Each
    of the structural domains also has its symmetry-induced domains

    LEEDStructuralDomains(params)

    Parameters
    ----------
    params
        LEEDParametersList or anything that LEEDParametersList accepts
    """

    @gl.exec_time
    def __init__(self, params, keep_duplicates=False):
        self.__parameters = gl.LEEDParametersList(params, keep_duplicates)

        self.__list = [gl.LEEDSymmetryDomains(p) for p in self.parameters]

        # __last_eq is a reference to the last LEEDEquivalentBeams
        # instance used to combine all structural and symmetry domains.
        # __last_domains a dictionary of the last domain indices used.
        # Both are set in self.equivalent_spots()
        self.__last_eq = None
        self.__last_domains = {}

    def __delitem__(self, elem):
        """Delete item at a given index or slice."""
        del self.__list[elem]

    def __getitem__(self, elem):
        """Retrieve item at a given index or slice.

        Notice that passing a slice will return a list object, not a
        LEEDStructuralDomains object!
        """
        return self.__list[elem]

    def __setitem__(self, idx, data):
        """Set the value of element(s) at a given index or slice.

        The input is checked for consistency.

        Parameters
        ----------
        idx : int or slice
        data : anything that can instantiate a LEEDParametersList

        Raises
        ------
        ValueError
            When a slice is passed, and: the data passed contains
            duplicates, or contains data that is a duplicate of
            one of the parameters already present.
        """
        # Check if data is compatible and if there are duplicates
        new_params = copy.deepcopy(self.parameters)
        new_params[idx] = data  # raises errors if incompatible
        length_change = len(new_params) - len(self.parameters)
        if isinstance(idx, slice):
            duplicate = (length_change < len(data))
        else:
            duplicate = not length_change
        if duplicate and isinstance(idx, slice):
            raise ValueError("Data contains duplicates")
        if not duplicate:
            self.parameters[idx] = data
            if isinstance(idx, slice):
                self.__list[idx] = [gl.LEEDSymmetryDomains(d)
                                    for d in data]
            else:
                self.__list[idx] = gl.LEEDSymmetryDomains(data)

    def __len__(self):
        """Return number of structural domains."""
        return len(self.__list)

    def __repr__(self):
        """Return string representation of self."""
        txt = 'LEEDStructuralDomains([\n'
        txt += '\n'.join(repr(d) for d in self.parameters)
        txt += '\n])'
        return txt

    @property
    def bulk_basis(self):
        """Return the reciprocal-space basis of the bulk lattice."""
        return self[0].bulk_basis

    @property
    def domain_ids(self):
        """Return an up-to-date list of domain 'ids'.

        Each id is either a name, if a name is present in the
        corresponding LEEDParameters dictionary, or a string
        "S"+number, with number-1 is the index in self

        Returns
        -------
        list of str
        """
        return [param['name'] if param['name'] else f"S{i+1}"
                for i, param in enumerate(self.parameters)]

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
            List of extinct beams as either tuples or BeamIndex,
            depending on whether self.fractional is False or True

        Raises
        ------
        RuntimeError
            When this property is called before the .equivalent_spots()
            method has been invoked at least once.
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("Need to call .equivalent_spots() once "
                               "before .extinct is available")
        return self.eq_beams_last_config.extinct

    @property
    def fractional(self):
        """Return whether beam indices should be expressed fractional.

        Return whether the beam indices need to be expressed as
        fractional, or if it is enough to use the numerators. In the
        latter case the processing of equivalent beams is faster.

        Return
        ------
        bool
            True if it requires fractional indices.
        """
        denominators = [s.denominator_for_bulk_beams for s in self]
        return any(den != denominators[0] for den in denominators[1:])

    @property
    def indexed_beams(self):
        """Return a list of indexed beams.

        If the configuration (i.e., domains selected, beam-incidence
        angles) is not up to date, call .equivalent_spots() before!

        Returns
        -------
        dict
            keys : {'beams', 'group_indices',
                    'overlap_domains', 'extinct_domains'}
            values : list
                One element per each one of the beams generated
                from the parameters passed at instantiation.

        Each element in the values is:
            'beams' : tuple or BeamIndex
                Tuple with only numerators if not self.fractional.
                Full fractional index of type BeamIndex otherwise.
            'group_indices' : int
                Progressive index that uniquely identifies the beam
                equivalence class to which beam belongs.  If the beam
                is fully extinct, its group index is negative.
            'overlap_domains' : list
                List of domain identifiers specifying which domains
                contribute to this beam. Each element is a tuple
                of the form (struct_id, (symm_id1, ...)), where
                struct_id is one of self.domain_ids, and symm_id
                one of LEEDSymmetryDomains.domain_ids
            'extinct_domains' : list
                List of domains, among those that overlap, that
                only contribute to the beam in an extinct fashion.
                Each entry has the same form as in 'overlap_domains'.
                If a structure contributes only with non-extinct
                contributions (i.e., all its symmetry domains are
                non-extinct), it will not appear. If a subset of
                the symmetry domains of one structure contributes
                in an extinct fashion, only these will be listed
                in the second element of (struct_id, (sym_id1, ...)).

        Raises
        ------
        RuntimeError
            When this property is called before the equivalent_spots()
            method has been invoked at least once.
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("Need to call .equivalent_spots() once "
                               "before .indexed_beams is available")
        return self.eq_beams_last_config.indexed_beams

    @property
    def last_domains_used(self):
        """Return identifiers for the domains in the last configuration."""
        if self.eq_beams_last_config is None:
            raise RuntimeError("Need to call .equivalent_spots() once "
                               "before .last_domains_used is available")
        return self.__last_domains

    @property
    def n_domains(self):
        """Return number of distinct domains produced by bulk symmetry."""
        return len(self)

    @property
    def parameters(self):
        """Return the LEEDParametersList that underlies self."""
        return self.__parameters

    @property
    def superlattices(self):
        """Superlattice matrices of each domain.

        Returns
        -------
        Dictionary of the form {id: superlattices} for each of the
        structural domains, where superlattices is a numpy.ndarray
        of superlattice matrices
        """
        superlattices = (p.superlattices for p in self)
        return dict(zip(self.domain_ids, superlattices))

    def insert(self, idx, data, keep_duplicates=False):
        """Override list.insert allowing check for duplicates."""
        # the next raises errors if data is not compatible
        old_len = len(self.parameters)
        self.parameters.insert(idx, data, keep_duplicates)
        if len(self.parameters) > old_len:
            # do something only if something was actually inserted
            self.__list.insert(idx, data)

    def append(self, data, keep_duplicates=False):
        """Override list.append allowing check for duplicates."""
        self.insert(len(self), data, keep_duplicates)

    def extend(self, data, keep_duplicates=False):
        """Override list.extend allowing check for duplicates."""
        for datum in data:
            self.append(datum, keep_duplicates)

    @gl.exec_time
    # @gl.profile_lines
    def equivalent_spots(self, domains=None, theta=None, phi=None):
        """Return spots of selected domains that are equivalent.

        Return a list of beams grouped by equivalence, given a primary
        beam direction defined by polar and azimuthal angles, for the
        selected domains

        Parameters
        ----------
        domains : dict, list, int or None (default=None)
            When passing a dictionary, the keys are used to select
            the structural domains based on their name; names not
            found are skipped. When passing a list it should have as
            many elements as there are domains. In both cases each
            value/element can be either None (selects all domains),
            an integer, or a list of integers, selecting which
            symmetry-equivalent domains to use. When a single
            integer is passed, this is used as the positional index
            for all the domains.
        theta, phi : number or None
            polar and azimuthal directions in degrees of the primary
            beam.  If not given or None, those defined at instantiation
            are used.  theta should be in the [-90, 90] range.
            phi is positive counterclockwise, and measured from the x
            axis in the Cartesian reference of the real-space basis of
            the first domain.
        """
        # Figure out which domains we want:
        domains, structures, struct_ids = self.process_domains_input(domains)

        # Check and process the angles
        if theta is None:
            theta = self.parameters[0]['beamIncidence'][0]
        if phi is None:
            phi = self.parameters[0]['beamIncidence'][1]
        theta, phi = gl.conventional_angles(float(theta), float(phi))

        # Call the equivalent_spots() method of each of the relevant
        # LEEDSymmetryDomains instances so that the LEEDEquivalentBeams
        # instances are correctly initialized
        for symmetry_domains, structure in zip(domains, structures):
            structure.equivalent_spots(symmetry_domains, theta, phi)

        # Now each structure.eq_beams_last_config contains a
        # reference to the correct LEEDEquivalentBeams instance.
        # Combine them into another LEEDEquivalentBeams instance
        # that includes all structures.  Also build ids.
        eq_beams = [s.eq_beams_last_config for s in structures]
        symm_ids = [s.last_domains_used for s in structures]
        combined_ids = tuple(zip(struct_ids, symm_ids))
        self.__last_domains = dict(combined_ids)

        # Check if we need to pass denominators together with
        # the numerators.  This is necessary if we have multiple
        # structures with differently sized cells.
        denominators = []
        if self.fractional:
            denominators = [s.denominator_for_bulk_beams for s in structures]

        self.__last_eq = gl.LEEDEquivalentBeams(eq_beams,
                                                domain_ids=combined_ids,
                                                denominators=denominators)
        if not self.eq_beams_last_config.was_cached:
            self.__correct_beam_index()
        return self.indexed_beams

    def from_domain_id(self, domain_id):
        """Return the LEEDSymmetryDomains with a certain id.

        Parameters
        ----------
        domain_ids : str
            Idenifier of the domain. If using the positional index,
            just use self[idx].

        Returns
        -------
        LEEDSymmetryDomains
        """
        try:
            idx = self.domain_ids.index(domain_id)
        except ValueError as err:
            raise ValueError(f"{domain_id!r} not found.") from err
        return self[idx]

    def set_beam_incidence(self, theta=None, phi=None):
        """Edit the beam incidence.

        Also triggers an update of the beam equivalence
        if self.equivalent_spots() was already called.

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

        old_angles = self.parameters[0]['beamIncidence']

        # Update the underlying parameters, check if anything changed
        self.parameters.set_beam_incidence(theta, phi)

        new_angles = self.parameters[0]['beamIncidence']
        if new_angles == old_angles:
            return

        # Update the individual symmetry domains. This
        # triggers an update of the LEEDSymmetryDomains if
        # their .equivalent_spots() method was called before.
        for sym_domain in self:
            sym_domain.set_beam_incidence(*new_angles)
        # Finally update the full spot equivalence
        if self.eq_beams_last_config is not None:
            self.equivalent_spots(self.last_domains_used)

    def process_domains_input(self, domains):
        """Get correct structural domains. Used in equivalent_spots.

        Parameters
        ----------
        domains : dict, list, int or None (default=None)
            When passing a dictionary, the keys are used to select
            the structural domains based on their name; names not
            found are skipped. When passing a list it should have as
            many elements as there are domains. In both cases each
            value/element can be either None (selects all domains),
            an integer, or a list of integers, selecting which
            symmetry-equivalent domains to use. When a single
            integer is passed, it is used for all the domains.

        Returns
        -------
        domains : list
            Fixed list of domain indices
        structures : list
            Which structures were requested. Each element is a
            gl.LEEDSymmetryDomains instance
        ids : list
            Which of the ids correspond to the structures

        Raises
        ------
        ValueError
            When the number of domain indices passed is inconsistent
            (i.e., smaller or larger) than the number of structural
            domains.
        TypeError
            When the entries in the domain indices are not integers
        """
        ids = self.domain_ids
        structures = self
        if domains is None:
            domains = [None]*self.n_domains
        elif isinstance(domains, dict):
            domains_dict = domains
            ids = []
            domains = []
            structures = []
            for struct_id, symm_ids in domains_dict.items():
                try:
                    structure = self.from_domain_id(struct_id)
                except ValueError:
                    pass
                else:
                    ids.append(struct_id)
                    domains.append(symm_ids)
                    structures.append(structure)
        elif hasattr(domains, '__len__'):
            if len(domains) != self.n_domains:
                raise ValueError("equivalent_spots: Inconsistent number "
                                 "of structural domains. Expected "
                                 f"{self.n_domains}, found {len(domains)}")
        elif isinstance(domains, int):
            domains = [domains]*self.n_domains
        else:
            raise TypeError("equivalent_spots: Invalid domains "
                            f"{type(domains).__name__}")
        return domains, structures, ids

    @gl.exec_time
    # @gl.profile_lines
    def __correct_beam_index(self):
        """Correct symmetry domain indices in self.indexed_beams.

        When the LEEDEquivalentBeams combines different structural
        domains, it places in the 'overlapping list' and 'extinct list'
        the full domain ids used, e.g., [('S1', (0, 1, 2)), ...], even
        if the beam being considered only has contributions from domain
        1 of structure 'S1', as this information is not available. NOT
        TRUE! IT KNOWS IF IT WAS CREATED FROM A LIST OF LEEDEquivalentBeams
        AND IF IT WAS PASSED DENOMINATORS!!

        Here we sort out this mess by taking the LEEDSymmetryDomains
        that were used for the superposition, and selecting only
        those that contain each beam.

        This method is used only internally, and should be called
        at the end of self.equivalent_spots().
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("Cannot reindex beams if .equivalent_spots() "
                               "was never called before.")
        indexed = self.eq_beams_last_config.indexed_beams

        # Prepare the denominators to make access a bit faster
        denominators = {}
        if self.fractional:
            denominators = {struct_id: struct.denominator_for_bulk_beams
                            for struct_id, struct in zip(self.domain_ids,
                                                         self)}
        # for i, (beam, group, overlap, extinct) in enumerate(indexed):
        for i, (beam, group,
                overlap, extinct) in enumerate(zip(*indexed.values())):
            correct_overlap = []
            correct_extinct = []
            extinct = dict(extinct)
            for struct_id, symm_ids in overlap:
                dom_beam = beam
                if denominators:
                    dom_beam = dom_beam * denominators[struct_id]
                struct = self.from_domain_id(struct_id)
                corrected = struct.domains_containing_beam(dom_beam)
                correct_overlap.append((struct_id, corrected))

                extinct_ids = extinct.get(struct_id, None)
                if extinct_ids:
                    correct_extinct.append(
                        (struct_id,
                         tuple(symm_id
                               for symm_id in extinct_ids
                               if symm_id in corrected))
                        )

            # indexed[i] = (beam, group, correct_overlap, correct_extinct)
            indexed['overlap_domains'][i] = correct_overlap
            indexed['extinct_domains'][i] = correct_extinct
