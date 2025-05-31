"""Module beams_old of viperleed.gui.leedsim.classes.

Defines the LEEDSymetryDomains class, used for constructing a LEED
pattern. This module is kept for now while support for multiple
structural domains is implemented.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-03-13'
__license__ = 'GPLv3+'

import copy
from warnings import warn as warning   # eventually will replace with logging

import numpy as np


class LEEDEquivalentBeams:                                                      ## TODO: fix docstring
    """
    This class is just a processor for finding beams that are equivalent among
    the ones used for instantiation. The beams originate from the same bulk,
    whose reciprocal-space basis is also passed during instantiation.

    The class is instantiated with one or more beam dictionaries lists. Each of
    the beam dictionaries in the list represents a 'building block', e.g., the
    domains originating from symmetry or different structural domains.
    Each of the dictionaries should have the following structure:
    keys : tuple
        each key is a beam index (with respect to the basis passed) of a LEED
        beam
    values : set of tuples
        list of beams equivalent to key
    In addition a list of glide-extinct beams can also be passed upon
    instantiation. If passed, it should have the same length as the domains
    dictionaries.
    
    Alternatively, also list of LEEDEquivalentBeams instances can be passed in
    the constructor, in which case the information will be gathered from them,
    and there is no need to pass the bulk_basis or the lists of extinct beams
    
    The primary goal of the class is to combine all the dictionaries given as
    input into new data structures that contain the beam equivalences when all
    the domains given upon instantiation are combined in an incoherent fashion
    
    Since the calculations are relatively slow, and the combinations are many,
    it's a good idea to cache instances with a hash that is computed such
    that it is unique for the geometry considered. This gives speedups of more
    than a factor 100 upon re-using an already calculated geometry
    
    This class is intended only for internal use. TRUE???
    """

    hash_keys = ('basis', 'superlattice', 'domain_ids', 'angle_key')

    __cache = {}
    def __init__(self, domains, **kwargs):
        # Do some checking of the parameters passed
        # see if the domains passed are all instances of LEEDEquivalentBeams
        if all(isinstance(dom, LEEDEquivalentBeams) for dom in domains):
            # check that the elements passed are consistent with one another,
            # i.e. that the lattice bases are compatible
            if not all(np.allclose(dom.basis, domains[0].basis)
                       for dom in domains):
                raise ValueError("LEEDEquivalentBeams: domains have different"
                                 "bulk bases, and cannot be combined")
            # In this case, treat the current instance as a superposition
            # of the inputs. This allows to treat different STRUCTURES, each
            # on the same bulk, and each possibly containing several SYMMETRY-
            # related domains
            basis = domains[0].basis
            extinct_lists = [dom.extinct for dom in domains]
            domains = [dom.equivalent_beams_dict for dom in domains]

        elif not all(isinstance(dom, dict) for dom in domains):
            # The only other option is passing a list of beam-equivalence
            # dictionaries, as well as all the (now mandatory) parameter basis
            raise TypeError("LEEDEquivalentBeams: only list of dictionaries or "
                            "list of LEEDEquivalentBeams are acceptable. Not "
                            "possible to use mixed lists")
        else:
            try:
                basis = kwargs['basis']
            except KeyError:
                raise TypeError("LEEDEquivalentBeams: missing mandatory basis"
                                "parameter")
            extinct_lists = kwargs.pop('extinct_lists', []*len(domains))

        if np.shape(basis) != (2, 2):
            raise ValueError("LEEDEquivalentBeams: invalid bulk_basis "
                             + f"{basis}. ".replace('\n', '')
                             + "Expected a 2x2 array-like")
        if len(extinct_lists) != len(domains):
            raise ValueError("LEEDEquivalentBeams: incompatible number of "
                             + f"domains ({len(extinct_lists)}) found in "
                             + f"extinct_lists. Expected {len(domains)}")
        domain_ids = kwargs.get('domain_ids', range(len(domains)))
        
        self.__hash_dict = {}      # set by self.__update_hash
        
        # re-prepare the kwargs to update the hash dictionary, updating values
        # that may have changed during this initialization
        kwargs['basis'] = basis
        kwargs['domain_ids'] = domain_ids
        self.__update_hash(**kwargs)

        if self.is_cached:
            # skip all the calculations, and take the attributes from the cached
            # object
            cached = self.__cache[self.hash]
            self.__basis = cached.basis
            self.__eq_beams_dict = cached.equivalent_beams_dict
            self.__extinct = cached.extinct
            self.__indexed_beams = cached.indexed_beams
            return None

        # if never calculated before, or hash cannot be determined due to lack
        # of input, calculate again. If it can be hashed, cache self
        self.__basis = basis
        self.__eq_beams_dict = {}  # set in __build_beam_equivalence_dict
        self.__extinct = []        # set in __index_beams
        beam_groups = self.__build_beam_equivalence_dict(domains)
        beam_groups = self.__sort_beams(beam_groups)
        self.__indexed_beams = self.__index_beams(beam_groups, domains,
                                                  extinct_lists, domain_ids)
        if self.hash:
            self.__cache[self.hash] = self

    @staticmethod
    def clear_cache():
        """
        Completely clears the hash table of LEEDEquivalentBeams instances.
        This may be useful in case a completely new bulk structure is loaded,
        in case memory usage starts being an issue.
        """
        LEEDEquivalentBeams.__cache = {}

    @property
    def basis(self):
        """
        Basis vectors as a 2x2 iterable. a = basis[0], b = basis[1]
        """
        return self.__basis

    @property
    def extinct(self):
        """
        List of truly extinct beams, each is a tuple. Beams are truly extinct
        if they do not overlap, or if they originate from the superposition
        of extinct beams
        """
        return self.__extinct

    @property
    def equivalent_beams_dict(self):
        """
        Returns
        -------
        dict
            keys : tuple
                beam i
            values : set of tuples
                beams equivalent to beam i, taking the superposition of domains
                into account
        """
        return self.__eq_beams_dict

    @property
    def indexed_beams(self):
        return self.__indexed_beams

    @property
    def hash(self):
        """
        Computes the hash of the geometry defined during instantiation
        """
        if not self.__hash_dict:
            return 0
        return hash(tuple(self.__hash_dict.values()))

    @property
    def is_cached(self):
        return self.hash and (self.hash in self.__cache)

    def __update_hash(self, **kwargs):
        missing = [key for key in self.hash_keys if key not in kwargs]
        if missing:
            warning(f"LEEDEquivalentBeams: missing {missing}. No hash possible")
            return None
        self.__hash_dict = kwargs
        for k, v in self.__hash_dict.items():
            if k in ('basis', 'superlattice'):
                v = zip(*v)
            if k in ('basis', 'superlattice', 'domain_ids'):
                self.__hash_dict[k] = tuple(v)

    def __build_beam_equivalence_dict(self, beams_by_domain):
        """
        Parameters
        ----------
        beams_by_domain : list of dict
            Each dict is of the form
            key : tuple
                beam j of domain i
            value : set of tuples
                beams equivalent to beam j, possibly only considering beams
                of the domain, and certainly without accounting for
                superposition between domains

        Returns
        -------
        list of sets of tuples
            each of the elements of the list, is the collection of all beams
            that are equivalent to one another as a result of the superposition
            of the domains
            NB: this method also sets the instance attribute __eq_beams_dict,
                i.e., a dictionary of {beam: {beams equivalent to beam}}
        """
        # store a deepcopy of the dictionary, as later on we will pop some
        # of the elements, and we don't want to screw with the input
        beams_by_domain_copy = copy.deepcopy(beams_by_domain)
    
        # Flatten the list of all beams, keeping only uniques
        flat_beams = set(beam for beams in beams_by_domain for beam in beams)
        
        # Now figure out which of the beams are symmetry equivalent
        # taking superposition of spots into account:
        # spots originating from beams that are equivalent within each domain,
        # are equivalent only if all domains contribute with equivalent beams 
        # in all spots
        #
        # For example (symmetry domains on hexagonal bulk):
        #   dom1: (1 0), (-1 0), (1 -1), (-1 1) equivalent
        #   dom2: (1 0), (-1 0), (0 -1), (0 1)  equivalent
        #   dom3: (1 0), (-1 0)                 equivalent
        #   --> (1 0) and (-1 0) are equivalent
        #
        #   dom1: (1 0), (-1 0), (1 -1), (-1 1)  equivalent
        #   dom2:                (1 -1), (-1 1)  equivalent
        #   dom3: (0 1), (0 -1), (1 -1), (-1 1)  equivalent
        #   --> (1 -1) and (-1 1) are equivalent, but are
        #       not equivalent to [(1 0), (-1 0)]
        #
        # Another example (structural domains on square bulk):
        # p(2x2)-pmm + c(2x2)-pm[1 0] -- take spot (1 1)
        #   
        #   p(2x2): (1 1), (1 -1), (-1  1), (-1 -1)  equivalent
        #   c(2x2): (1 1), (1 -1)                    equivalent
        #   --> (1 1) and (1 -1) are equivalent
        #
        #   p(2x2): (1 1), (1 -1), (-1 1), (-1 -1)   equivalent
        #   c(2x2):                (-1 1), (-1 -1)   equivalent
        #   --> (-1 1) and (-1 -1) are equivalent, but are
        #       not equivalent to [(1 1), (1 -1)]
        beam_equivalence = []
        for beam in flat_beams:
            # from each structure, if there is a beam <beam>, take the list of all
            # those equivalent to it.
            # Example: beam = (1 1)
            #   -> take [[(1 1), (1 -1), (-1 -1), (1 -1)],  from p(2x2)
            #            [(1 1), (1 -1)]]                   from c(2x2)
            beam_sets = [beams_dict[beam]
                         for beams_dict in beams_by_domain_copy
                         if beam in beams_dict]
            if beam_sets:
                # if there are beams to process in the list, take the set
                # intersection, i.e., find all those in common to all structures
                common_beams = set.intersection(*beam_sets)
                
                # now "mark as processed" in all the structures the beams coming
                # from the intersection by removing them, and set the correct
                # keys of beam_equivalence
                for processed_beam in common_beams:
                    for unprocessed_beams in beams_by_domain_copy:
                        unprocessed_beams.pop(processed_beam, None)
                    if processed_beam in self.__eq_beams_dict:  # TODO: remove this after it's clear that it does not happen
                        raise KeyError("TEMP CHECK: setting multiple times the same key!")
                    self.__eq_beams_dict[processed_beam] = common_beams
                beam_equivalence.append(common_beams)
        return beam_equivalence

    def __sort_beams(self, beam_groups):
        """
        Sort beams given as a list of iterables. Each of the elements in the
        list is treated as a group of equivalent beams
        """
        # sort within each equivalence group and by energy
        return sorted(
            (sorted(beams, key=self.__sort_hk) for beams in beam_groups),
            key=self.__sort_energy
        )
    
    def __index_beams(self, beam_groups, beams_by_domain, extinct, domain_ids):
        """
        Given the (sorted) list of beam groups, assign to each of the
        groups one index. Also, determine which of the original domains
        contribute to the overlapped spots, and which of these contributions
        come from extinct beams. The domains (both in the overlapping and in
        the extinct lists) are labeled according to the domain_ids given
        NB: this method also sets up the __extinct instance attribute, that
            contains all the beams that are truly extinct, as a result of the
            superposition of extinct beams

        Parameters
        ----------
        beam_groups : iterable of iterables of tuples
            List of beam groups, already sorted
        beams_by_domain : list of iterables of tuples
            each entry is a domain, and contains all beams of the domain
        extinct : list of tuples
            Each entry in the list contains the beams that are extinct in each
            domain
        domain_ids : iterable
            contains names for the domains. Typically would be integers
        """
        
        indexed_beams = []
        for i, beams in enumerate(beam_groups):
            for beam in beams:
                # find which domains overlap
                overlapping = [d
                               for d, domain_beams in enumerate(beams_by_domain)
                               if beam in domain_beams]
                group_idx = i
                
                # figure out whether the beam is extinct in all the structures
                # that overlap, in which case the index goes negative
                extinct_domains = [d for d in overlapping if beam in extinct[d]]
                n_extinct = len(extinct_domains)
                if n_extinct == len(overlapping) and n_extinct > 0:
                    group_idx *= -1
                    self.__extinct.append(beam)
                indexed_beams.append((beam, group_idx,
                                      [domain_ids[d] for d in overlapping],
                                      [domain_ids[d] for d in extinct_domains]))
        # TODO: This format is not great. Since I'm always using it unzipped, it
        # may be worth just returning it unzipped in the first place, maybe as
        # a dictionary just to have some reasonable names, like 'beams',
        # 'groups', 'overlap domains', 'extinct domains'
        return indexed_beams

    def __sort_hk(self, index):
        sortH, sortK = index
        return -sortH-sortK, -sortH
    
    def __sort_energy(self, indices):
        # assumes that the indices passed represent beams that are equivalent to
        # one another, i.e., they all have the same length
        g = np.dot(indices[0], self.basis)
        # sort by energy, i.e., by ||g||^2:
        return g[0]**2 + g[1]**2



