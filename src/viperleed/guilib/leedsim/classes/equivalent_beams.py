"""Module beams of viperleed.guilib.leedsim.classes.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the LEEDSymetryDomains class, used for constructing a LEED pattern

Author: Michele Riva
Created: 2021-03-13
"""

import itertools
from warnings import warn as warning  # TODO: replace with logging

import numpy as np

from viperleed.guilib.classes import planegroup
from viperleed.guilib.classes.beamindex import BeamIndex
from viperleed.guilib.helpers import equal_dicts
from viperleed.guilib.helpers import two_by_two_array_to_tuple
from viperleed.guilib.leedsim.utils import sort_hk

from viperleed.guilib import decorators as dev_


@dev_.exec_time
# @dev_.profile_lines
# @dev_.profile_calls()
def beams_dict_numerators_to_fractional(beam_dict, denominator):
    """Return a dictionary of fractional indices from numerators.

    Takes a dictionary of the form {beam: set(beams)} and returns it
    as {beam/denominator: set(beams/denominator)}. The keys and the
    elements in the sets in the returned dictionary are BeamIndex
    instances. Using a (hashable) BeamIndex removes computational
    errors.

    Parameters
    ----------
    beam_dict : dict
        Dictionary containing only the numerators of the fractional
        indices of beams. Form: {beam: set(beams equivalent to beam)}
    denominator : int
        The denominator to be used to compute fractions

    Returns
    -------
    modified_beam_dict : dict
        Same as before, but each entry is a BeamIndex instance and is
        now a full fractional index, i.e., original/denominator.  Also,
        the values in the dictionary returned are a tuple, not a set!
    """
    ret = dict()
    for beam, beam_set in beam_dict.items():
        fract_beam = BeamIndex(beam, denominator=denominator,
                               from_numerators=True)
        ret[fract_beam] = beams_list_numerators_to_fractional(beam_set,
                                                              denominator)
    return ret


def beams_list_numerators_to_fractional(beam_list, denominator):
    """Return a tuple of fractional indices from numerators.

    Takes an iterable of beam indices and returns a tuple where
    each element is a fractional BeamIndex, with numerator equal
    to the beam in the original iterable, and the denominator given.

    Parameters
    ----------
    beam_list : iterable
        Iterable containing only the numerators of the fractional
        indices of beams
    denominator : int
        The denominator to be used to compute fractions

    Returns
    -------
    modified_beam_list : tuple
        Same as before, but each entry is a BeamIndex instance and
        is now a full fractional index, i.e., original/denominator.
        The return value used to be a tuple, but there is no reason
        to require ordering, as this is anyway used as a set in
        beams_dict_numerators_to_fractional, and used only for
        membership test in LEEDEquivalentBeams.__index.
    """
    return tuple(BeamIndex(beam, denominator=denominator,
                           from_numerators=True)
                 for beam in beam_list)


def same_bases(*bases):
    """Return whether the bases passed are all the same."""
    if len(bases) == 1:
        return True
    # Check for the possibility of sign changes, by
    # checking if abs(bi @ b0^-1) is the identity.
    # The multiplication is done with einsum for speed.
    inv_first = np.linalg.inv(bases[0])
    t_bases = np.einsum('ijk,kl->ijl', bases[1:], inv_first)
    return np.allclose(np.abs(t_bases) - planegroup.E, 0)


class LEEDEquivalentBeams:
    """Splits beams into equivalence classes.

    This class is a processor for finding beams that are equivalent
    among the ones used for instantiation.  The beams should originate
    from the same bulk, whose reciprocal-space basis is also passed
    during instantiation.

    The beams given are combined into a new data structure that
    contains the overall beam equivalences, i.e., the equivalence
    of the incoherent superposition of beams (if any are superposed)

    Since the calculations are relatively slow, and the combinations
    are many, instances are cached within a hash table. The hash is
    computed such that it is unique for the geometry considered.
    This gives speedups of more than a factor 100 upon re-using an
    already calculated geometry.

    This class is intended only for internal use. It does the heavy
    lifting for calculating beam eaquivalnce in the LEEDSymDomains,
    LEEDStructDomains, and LEEDPattern classes.
    """

    hash_keys = ('basis', 'superlattice', 'domain_ids', 'angle_key', 'caller')

    __dict_cache = {}  # instances created from dicts
    __self_cache = {}  # instances created from LEEDEquivalentBeams

    # TODO: pylint too-complex mccabe=11
    # @dev_.profile_lines
    @dev_.exec_time
    def __init__(self, domains, **kwargs):
        """Initialize LEEDEquivalentBeams.

        Parameters
        ----------
        domains : list
            list of LEEDEquivalentBeams or list of dict, but not
            mixed. When dictionaries are passed, each one should be
            {beam: set(beams equivalent to beam)} where each beam
            is a 2-tuple-like. When passing a single instance of
            LEEDEquivalentBeams for superposition, very little
            computation is done, but the keys 'overlap_domains'
            and 'extinct_domains' of the .indexed_beams property
            are the same as those of domains[0], i.e., they
            do not contain the domain_ids passed. The caller needs
            to fix this a fortiori.
        basis : array-like
            2x2 matrix of the bulk basis, with a, b = basis. This
            parameter is optional when domains is a list of
            LEEDEquivalentBeams, as it is extracted from domains.
            Used to compute the hash.
        extinct_lists : iterable, optional
            each element is an iterable of the beams that are fully
            extinct in each domain
        domain_ids : iterable or int, optional
            list of "names" for each of the domains passed
            (default=range(len(domains))). Used to compute the hash.
        denominators : iterable, optional
            list of denominators to be used for normalizing each of the
            beam indices passed (default=None). Passing this argument
            makes instantiation much slower the first time a new
            configuration is created. It should be used only when
            combining structures with differently sized unit cells.
            In all other cases, passing fractional indices or only
            numerators in domains and extinct_lists is faster. When
            all domains are a LEEDEquivalentBeams, denominators for
            those that .has_fractional_beams are ignored.
        superlattice : array-like, optional
            2x2 superlattice matrix or iterable of 2x2 superlattice
            matrices. Used exclusively to compute the hash. Hash cannot
            be computed without this. The value passed is ignored when
            all domains are LEEDEquivalentBeams, as it is automatically
            calculated.
        angle_key : str, optional
            can be one of {'norm', 'other', angles} where angles is one
            of the special azimuthal angles (i.e., parallel to
            mirror/glide planes). Used to compute the hash. Hash cannot
            be computed without this. The value passed is ignored when
            all domains are LEEDEquivalentBeams, as it is automatically
            calculated.
        caller : object, optional
            The instance that requested the initialization. Used for
            hashing.

        Raises
        ------
        TypeError
            * domains is a list of dictionaries but no basis is passed
            * Entries in domains are neither LEEDEquivalentBeams nor
              dict, or a mix of LEEDEquivalentBeams and dict
        ValueError
            * domains is a list of LEEDEquivalentBeams, but the bulk
              bases are incompatible
            * basis is not a 2x2 array-like
            * Number of entries in extinct_lists or denominators is
              different than those in domains
        ZeroDivisionError
            If any of the denominators is zero
        """
        domains, kwargs = self.__process_input(domains, kwargs)
        extinct_lists = kwargs.pop('extinct_lists')
        domain_ids = kwargs.get('domain_ids')

        # __hash_dict contains all the info that is used to build a
        # unique hash. It is set by self.__update_hash, and read via
        # the self.hash_dict property
        self.__hash_dict = {}
        self.__update_hash(**kwargs)

        # See if we can skip all the calculations, and
        # rather take the attributes from the cached object
        try:
            self.load_from_cache(kwargs['which_cache'])
        except RuntimeError:  # Nope, will do all calculations
            pass
        else:  # Yes, it's cached
            return

        self.__loaded_from_cache = False
        self.__basis = kwargs.get('basis')
        self.__eq_beams_dict = {}  # set in __build_beam_equivalence_dict
        self.__extinct = tuple()   # set in __index_beams

        # Before going on, if there is a 'denominators' keyword
        # argument we need to re-process domains and extinct_lists
        # to include that.  This is especially necessary to combine
        # different structures.
        # Since this slows down this part quite a bit, it is advisable
        # not to pass 'denominators' unless strictly necessary.
        denominators = kwargs.get('denominators', None)
        self.__has_fract_beams = False
        if denominators:
            self.__has_fract_beams = True
            if len(denominators) != len(domains):
                raise ValueError("LEEDEquivalentBeams: Inconsistent number "
                                 f"of denominators. Expected {len(domains)}, "
                                 f"found {len(denominators)}")
            if any(abs(den) < 1e-3 for den in denominators):
                raise ZeroDivisionError("LEEDEquivalentBeams: all "
                                        "denominators must be nonzero.")
            dict_to_fract = beams_dict_numerators_to_fractional
            list_to_fract = beams_list_numerators_to_fractional
            domains = [dict_to_fract(dom, den) if den else dom
                       for dom, den in zip(domains, denominators)]
            extinct_lists = [list_to_fract(ext, den) if den else ext
                             for ext, den in zip(extinct_lists, denominators)]

        try:
            self.add_to_cache(kwargs['which_cache'])
        except RuntimeError:
            # Can't be cached as some hash keys are missing
            pass

        if len(domains) == 1 and kwargs['which_cache'] == 'self':
            # There is nothing to superpose. Skip calculations,
            # And mark it as loaded from cache, as a single
            # LEEDEquivalentBeams has already seen the correct
            # processing all the way.
            self.__loaded_from_cache = True

            # We may have to reprocess beams to fractional, though
            indexed = kwargs['indexed_beams']
            eq_beams_dict = domains[0]
            extinct_lists = extinct_lists[0]
            if denominators and not kwargs['fractional']:
                # Need to turn beams into fractional
                den = denominators[0]
                indexed['beams'] = list_to_fract(indexed['beams'], den)
                eq_beams_dict = dict_to_fract(eq_beams_dict, den)
                extinct_lists = list_to_fract(extinct_lists, den)

            self.__indexed_beams = indexed
            self.__eq_beams_dict = eq_beams_dict
            self.__extinct = extinct_lists

            # NB: this way, the domains_ids are NOT the ones
            # that were originally passed. We assume that the
            # caller will take care of this.
            caller = kwargs.get('caller', None)
            try:  # Do not use isinstance to avoid cyclic imports
                _ = caller.eq_beams_last_config
            except AttributeError:
                _is_struct_domains = False
            else:
                _is_struct_domains = True
            if not caller or not _is_struct_domains:
                warning("LEEDEquivalentBeams: Single domain was already "
                        "calculated, but domain_ids passed were unused. "
                        "They will need to be fixed at the source!",
                        RuntimeWarning)
            return

        beam_groups = self.__build_beam_equivalence_dict(domains)
        beam_groups = self.__sort_beams(beam_groups)
        self.__indexed_beams = self.__index_beams(beam_groups, domains,
                                                  extinct_lists, domain_ids)

    @staticmethod
    def __process_input(domains, kwargs):
        """Process input as per __init__ doc."""
        # Do some checking of the parameters passed:
        # (0) check if the domain_ids passed have the right length
        domain_ids = kwargs.get('domain_ids', range(len(domains)))
        if len(domain_ids) != len(domains):
            raise ValueError("LEEDEquivalentBeams: incompatible number "
                             + f"of domain_ids ({len(domain_ids)}) "
                             + f"found. Expected {len(domains)}")
        kwargs['domain_ids'] = domain_ids

        # (1) see if the domains passed are all instances of
        #     LEEDEquivalentBeams. In this case, treat the
        #     current instance as a superposition of the inputs.
        #     This allows to treat different STRUCTURES, each
        #     on the same bulk, each possibly containing several
        #     SYMMETRY-related domains.
        if all(isinstance(dom, LEEDEquivalentBeams) for dom in domains):
            domains, kwargs = LEEDEquivalentBeams.__preprocess_eqbeams(domains,
                                                                       kwargs)

        # (2) The only other acceptable option is passing a
        #     list of beam-equivalence dictionaries, as well
        #     as the (now mandatory) parameter basis.
        elif all(isinstance(dom, dict) for dom in domains):
            kwargs['which_cache'] = 'dict'
        else:
            raise TypeError("LEEDEquivalentBeams: only list of dictionaries "
                            "or list of LEEDEquivalentBeams are acceptable. "
                            "Not possible to use mixed lists")

        basis = kwargs.get('basis', None)
        if basis is None:
            raise TypeError("LEEDEquivalentBeams: missing "
                            "mandatory basis parameter")
        if np.shape(basis) != (2, 2):
            raise ValueError("LEEDEquivalentBeams: invalid bulk_basis "
                             + f"{basis}. ".replace('\n', '')
                             + "Expected a 2x2 array-like")

        extinct_lists = kwargs.get('extinct_lists', [[]]*len(domains))
        if len(extinct_lists) != len(domains):
            raise ValueError("LEEDEquivalentBeams: incompatible number of "
                             + f"domains ({len(extinct_lists)}) found in "
                             + f"extinct_lists. Expected {len(domains)}")
        kwargs['extinct_lists'] = extinct_lists

        return domains, kwargs

    # TODO: pylint too-complex mccabe=14
    @staticmethod
    def __preprocess_eqbeams(domains, kwargs):
        """Trasform LEEDEquivalentBeams input as if dicts were passed.

        Parameters
        ----------
        domains : list of LEEDEquivalentBeams
            Domains to be combined in superpositions.
        kwargs : dict
            Optional arguments passed during instantiation.

        Raises
        ------
        ValueError
            If the domains have inconsistent bases.

        Returns
        -------
        domains : list of dict
            Beam equivalence dictionaries, each one is
            {beam: set(beams equivalent to beam)} where
            each beam is a 2-tuple-like.
        kwargs : dict
            Modified optional arguments to mimic instance
            creation from a list of dictionaries.
        """
        # Check that the elements passed are consistent with one
        # another, i.e. that the lattice bases are compatible.
        if not same_bases([dom.basis for dom in domains]):
            raise ValueError("LEEDEquivalentBeams: domains have "
                             "different bases, and cannot be combined")

        # Collect the information that is needed to build a
        # hash (if present for all)
        try:  # pylint: disable=too-many-try-statements
            kwargs['angle_key'] = ''.join(dom.hash_dict['angle_key']
                                          for dom in domains)
            kwargs['superlattice'] = tuple(dom.hash_dict['superlattice']
                                           for dom in domains)
        except KeyError:
            pass

        # Take a look in the cache, to see if, by chance, there
        # are already some versions of the same domains, computed
        # as single domains, i.e., with 'domains_id' = (same_id,)
        # that have fractional beams.
        # TODO: is 1000 cached too little/much?
        if len(LEEDEquivalentBeams.__self_cache) < 1000:
            cached_domains = tuple(LEEDEquivalentBeams.__self_cache.values())
            cached_ids = [dom.domain_ids for dom in cached_domains]
            for i, (dom, dom_id) in enumerate(zip(domains,
                                                  kwargs['domain_ids'])):
                if dom.has_fractional_beams:
                    continue

                idx = 0
                while 1:  # There may be multiple ones
                    try:
                        idx = cached_ids.index((dom_id,), idx)
                    except ValueError:
                        # no cached copy with the same ids
                        break
                    cached = cached_domains[idx]
                    cached_dict = {k: (v if k != 'superlattice' else v[0])
                                   for k, v in cached.hash_dict.items()}
                    if (equal_dicts(cached_dict, dom.hash_dict,
                                    ignore_keys=['domain_ids', 'caller'])
                            and cached.has_fractional_beams):
                        domains[i] = cached
                        break
                    idx += 1

        # Then collect the rest of the information
        # needed for creating the superposition
        kwargs['basis'] = domains[0].basis
        kwargs['extinct_lists'] = [dom.extinct for dom in domains]

        # Pre-process the denominators, if passed: if the
        # beams of some domains are already fractional,
        # we should not divide again by the denominator
        denominators = kwargs.get('denominators', None)
        if denominators:
            if len(denominators) != len(domains):
                raise ValueError("Inconsistent number of denominators. "
                                 f"Expected {len(domains)}, found "
                                 f"{len(denominators)}")
            denominators = [None if dom.has_fractional_beams else den
                            for dom, den in zip(domains, denominators)]
        kwargs['denominators'] = denominators

        # If there is a single domain, store its indexed
        # beams and whether the instance has fractional
        # beams, as this will allow us to skip calculations
        # later in the __init__.
        if len(domains) == 1:
            kwargs['indexed_beams'] = domains[0].indexed_beams
            kwargs['fractional'] = domains[0].has_fractional_beams

        # Finally, collect the beams
        domains = [dom.equivalent_beams_dict for dom in domains]
        kwargs['which_cache'] = 'self'

        return domains, kwargs

    @staticmethod
    def clear_cache(which_cache='all'):
        """Clear the selected cache of LEEDEquivalentBeams instances.

        This may be useful in case a completely new bulk structure is
        loaded, in case memory usage starts being an issue.

        Parameters
        ----------
        which_cache : {'all', 'self', 'dict'}
            Which of the two caches is to be cleared

        Raises
        ------
        ValueError
            If which_cache is not one of 'all', 'self', or 'dict'
        """
        if which_cache not in ('all', 'self', 'dict'):
            raise ValueError("LEEDEquivalentBeams: invalid cache "
                             f"{which_cache} cannot be cleared.")
        if which_cache in ('self', 'all'):
            LEEDEquivalentBeams.__self_cache = {}
        if which_cache in ('dict', 'all'):
            LEEDEquivalentBeams.__dict_cache = {}

    @property
    def basis(self):
        """Basis vectors as a 2x2 iterable. a = basis[0], b = basis[1]."""
        return self.__basis

    @property
    def domain_ids(self):
        """Return a tuple of the domain_ids used for construction."""
        if 'domain_ids' not in self.hash_dict:
            raise RuntimeError("Just for me: move call to __update_hash"
                               "earlier on during __init__")
        return self.hash_dict['domain_ids']

    @property
    def extinct(self):
        """List of truly extinct beams, each is a tuple.

        Beams are truly extinct if they do not overlap, or if they
        originate from the superposition of extinct beams
        """
        return self.__extinct

    @property
    def equivalent_beams_dict(self):
        """Return dictionary of equivalent beams.

        Returns
        -------
        dict
            keys : tuple
                beam i
            values : set of tuples
                beams equivalent to beam i, taking the superposition
                of domains into account
        """
        return self.__eq_beams_dict

    @property
    def has_fractional_beams(self):
        """Return whether the beams in self are fractional."""
        return self.__has_fract_beams

    @property
    def hash_dict(self):
        """Return dictionary used for computing the hash."""
        return self.__hash_dict

    @property
    def hash_(self):
        """Compute the hash of the geometry defined during instantiation."""
        if not self.hash_dict:
            return None
        return hash(tuple(self.__hash_dict.values()))

    @property
    def indexed_beams(self):
        """Return a list of beams with indexing.

        When a single instance of LEEDEquivalentBeams was
        used for superposition, the keys 'overlap_domains'
        and 'extinct_domains' are the same as those of
        domains[0], i.e., they do not contain the domain_ids
        passed. The caller is expected to fix this after
        instantiation.

        Returns
        -------
        dict
            keys : {'beams', 'group_indices',
                    'overlap_domains', 'extinct_domains'}
            values : list
                One element per each one of the beams passed
                at instantiation.

        Each element in the values is:
            'beams' : tuple or BeamIndex
                Tuple with only numerators, if no 'denominators'
                keyword was passed at instaniation. Full fractional
                index of type BeamIndex otherwise.
            'group_indices' : int
                Progressive index that uniquely identifies the beam
                equivalence class to which beam belongs.  If the beam
                is fully extinct, its group index is negative.
            'overlap_domains' : list
                List of domains that contribute to this beam.  The
                idenfiers used are those passed at instantiation.
            'extinct_domains' : list
                List of domains, among those that overlap, that
                only contribute to the beam in an extinct fashion.
        """
        return self.__indexed_beams

    @property
    def is_cached(self):
        """Return whether self is stored in the cache."""
        return self.hash_ and (self.hash_ in self.__dict_cache
                               or self.hash_ in self.__self_cache)

    @property
    def was_cached(self):
        """Return whether self was loaded from cache.

        This property is useful to check whether all the properties
        of self are up to date (under the assuption that, if it
        was cached they are), or if further processing may be needed.
        """
        return self.__loaded_from_cache

    def add_to_cache(self, which_cache):
        """Insert self in the correct cache."""
        # If it can be hashed, cache self
        if not self.hash_:
            raise RuntimeError("LEEDEquivalentBeams: not possible to "
                               "cache when hash parameters are missing.")

        if which_cache == 'dict':
            cache = self.__dict_cache
        elif which_cache == 'self':
            cache = self.__self_cache
        else:
            raise ValueError("LEEDEquivalentBeams: invalid "
                             f"cache {which_cache}. Should "
                             "be either 'self' or 'dict'")
        cache[self.hash_] = self

    def load_from_cache(self, which_cache):
        """Load self from one of the caches."""
        if not self.is_cached:
            raise RuntimeError("LEEDEquivalentBeams: self is not in "
                               "the cache. Impossible to load.")
        self.__loaded_from_cache = True

        if which_cache == 'dict':
            cache = self.__dict_cache
        elif which_cache == 'self':
            cache = self.__self_cache
        else:
            raise ValueError("LEEDEquivalentBeams: invalid "
                             f"cache {which_cache}. Should "
                             "be either 'self' or 'dict'")
        cached = cache[self.hash_]
        self.__basis = cached.basis
        self.__eq_beams_dict = cached.equivalent_beams_dict
        self.__extinct = cached.extinct
        self.__indexed_beams = cached.indexed_beams
        self.__has_fract_beams = cached.has_fractional_beams

    @dev_.exec_time
    def __build_beam_equivalence_dict(self, beams_by_domain):
        """Split beams into beam equivalence classes.

        This is the core method of this class. It is run once at
        instantiation, and should not be run ever again.

        It takes a list of dictionaries of beam equivalences, one per
        domain, and combines them into a single dictionary containing
        the equivalence classes corresponding to the superposition of
        the domains.  This is stored into the __eq_beams_dict instance
        attribute. self.__eq_beams_dict has then the form
        {beam: set(beams equivalent to beam)} where each beam has the
        same type as passed, and the equivalence set is the one after
        accounting for superposition.  The information in this class
        attribute is a bit redundant, as the 'values' in the
        dictionary are repeated several times. However, (i) using a
        dictionary gives O(1) lookup, and (2) the redundancy allows one
        to combine an infinte number of LEEDEquivalentBeams to account
        for superposition.  This is the way one can treat symmetry and
        structural domains with the same logics.

        It also returns a list of the equivalence classes found.

        Parameters
        ----------
        beams_by_domain : iterable
            The iterable should support copying, i.e., a generator
            is not approppriate.  The i-th element of the iterable
            is a dictionary of the form:

            key : tuple or BeamIndex
                beam j of domain i
            value : set
                beams (tuple or BeamIndex) equivalent to beam j,
                possibly only considering beams of the domain, and
                certainly without accounting for superposition
                between domains.

        Returns
        -------
        list
            List of all beam equivalence classes.  Each element is a
            set of tuples, containing beams that are equivalent to one
            another as a result of the superposition of the domains.
        """
        # Set up the beams to process taking all keys from all
        # the dictionaries passed.
        beams_to_process = set(itertools.chain.from_iterable(beams_by_domain))
        processed_beams = set()

        # Now figure out which of the beams are symmetry equivalent
        # taking superposition of spots into account:
        #
        # Spots originating from beams that are equivalent within each
        # domain, are equivalent only if all domains contribute with
        # equivalent beams in all spots.
        #
        # For example (symmetry domains on hexagonal bulk):
        #   dom1: (1 0), (-1 0), (1 -1), (-1 1) equivalent
        #   dom2: (1 0), (-1 0), (0 -1), (0  1) equivalent
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
        append_beams = beam_equivalence.append
        while beams_to_process:
            beam = beams_to_process.pop()
            # From each structure, if there is a beam <beam>, take the
            # list of all those equivalent to it.
            # Example: beam = (1 1)
            #   -> take [[(1 1), (1 -1), (-1 -1), (1 -1)],  from p(2x2)
            #            [(1 1), (1 -1)]]                   from c(2x2)
            #
            # NB: beams_dict[beam] is a much smaller set than
            #     processed_beams. The set difference s - t
            #     has time complexity O(len(s)).
            beam_sets = [set(beams_dict[beam]) - processed_beams
                         for beams_dict in beams_by_domain
                         if beam in beams_dict]

            # Take the set intersection, i.e., pick all those in
            # common to all domains that contribute
            common_beams = set.intersection(*beam_sets)

            # Mark these beams as processed.
            #
            # NB: common_beams is a small set. Both s.update(t)
            #     and s.difference_update(t) have O(len(t)) time
            #     complexity .
            processed_beams.update(common_beams)
            beams_to_process.difference_update(common_beams)

            # And set up the corresponding dictionary keys,
            # so that the dictionary looks like
            #     {beam: set(beams equivalent to beam)}
            for common_beam in common_beams:
                self.__eq_beams_dict[common_beam] = common_beams
            # beam_equivalence.append(common_beams)
            append_beams(common_beams)
        return beam_equivalence

    @dev_.exec_time
    # @dev_.profile_lines
    # @dev_.profile_calls()
    def __index_beams(self, beam_groups, beams_by_domain, extinct, domain_ids):
        """Index the input beams.

        Given the (sorted) list of beam groups, assign to each of the
        groups one index. Also, determine which of the original domains
        contribute to the overlapped spots, and which of these
        contributions come from extinct beams. The domains (both in the
        overlapping and in the extinct lists) are labeled according to
        the domain_ids given.

        This method also sets up the __extinct instance attribute,
        which contains all the beams that are truly extinct, as a
        result of the superposition of extinct beams.

        Parameters
        ----------
        beam_groups : iterable
            List of beam groups to be indexed. The sorting of the
            input (as given) is used to assign group indices.
            Used only once, so also a generator works.
        beams_by_domain : iterable
            One entry per domain, containing a list of all the beams
            of the domain. Iterated over multiple times, so a generator
            should not be used. No check is done.
        extinct : iterable
            One entry per domain, containing the beams that are extinct
            in that domain. Iterated over multiple times, so a generator
            should not be used. No check is done.
        domain_ids : iterable
            Contains names to be used for the domains. Typically would
            be integers or strings. Iterated over multiple times, so a
            generator should not be used. No check is done.

        Returns
        -------
        dict
            keys : {'beams', 'group_indices',
                    'overlap_domains', 'extinct_domains'}
            values : list
                One element per each one of the beams originally
                passed.

        Each element in the values is:
            'beams' : same as input.
                Beam identifier
            'group_indices' : int
                Index of the symmetry group to which beam belongs to.
            'overlap_domains' : list
                Domains that contribute to beam.  The names given in
                domains_ids are used.
            'extinct_domains' : list
                Those 'overlap_domains' that contribute to the spot
                with only extinct terms.  The names given in
                domains_ids are used.
        """
        # indexed_beams = []
        indexed_beams = {'beams': [],
                         'group_indices': [],
                         'overlap_domains': [],
                         'extinct_domains': []}
        for i, beams in enumerate(beam_groups):
            for beam in beams:
                # Find which domains overlap
                overlapping = [
                    d for d, domain_beams in enumerate(beams_by_domain)
                    if beam in domain_beams
                    ]
                group_idx = i

                # Figure out whether the beam is extinct in all
                # the structures that overlap, in which case the
                # index goes negative.
                extinct_domains = [d
                                   for d in overlapping
                                   if beam in extinct[d]]
                if len(overlapping) == len(extinct_domains) > 0:
                    group_idx *= -1
                    self.__extinct.append(beam)
                indexed_beams['beams'].append(beam)
                indexed_beams['group_indices'].append(group_idx)
                indexed_beams['overlap_domains'].append(
                    [domain_ids[d] for d in overlapping]
                    )
                indexed_beams['extinct_domains'].append(
                    [domain_ids[d] for d in extinct_domains]
                    )
                # indexed_beams.append(
                #     (beam, group_idx,
                #       [domain_ids[d] for d in overlapping],
                #       [domain_ids[d] for d in extinct_domains])
                #     )
        return indexed_beams

    # @dev_.profile_calls()
    @dev_.exec_time
    def __sort_beams(self, beam_groups):
        """Sort beams given as a list of iterables.

        Parameters
        ----------
        beam_groups : iterable
            Each of the elements in the list is treated as a group
            of symmetry-equivalent beams.

        Returns
        -------
        list of list
            Sorted version of the input, using the criteria defined
            by __sort_hk and self.__sort_energy, i.e., the beam
            groups will appear: (1) in increasing radial position
            (i.e., energy), and (ii) for beams of equal radial
            position, beams with larger h+k index within a group will
            appear first.  There is no criterion for sorting
            beam-equivalence groups one with respect to the other,
            except for the radial position.
        """
        # sort within each equivalence group and by energy
        return sorted((sorted(beams, key=sort_hk) for beams in beam_groups),
                      key=self.__sort_energy)

    def __sort_energy(self, list_of_equivalent_indices):
        """Key-sorting criterion for energies.

        This method can be used as the 'key' optional argument of
        list.sort() or sorted() to sort beams given as (h, k) indices
        in the iterable to be sorted with increasing radial distance
        from the (0, 0) beam.

        Parameters
        ----------
        list_of_equivalent_indices : iterable
            List of indices of beams in a beam-equivalence group.
            The method assumes that the indices passed represent
            beams that are equivalent to one another, i.e., they
            index reciprocal-space vectors which all have the same
            length.

        Returns
        -------
        float
            Length of reciprocal-space vector.
        """
        # This version uses np.dot, and it is ~8% faster than
        # using pure python and calculating explicitly the
        # norm squared

        # Use the first index only, under the assumption that all
        # beams passed are equivalent, i.e., same energy
        g_vector = np.dot(list_of_equivalent_indices[0], self.basis)
        # Sort by energy, i.e., by ||g||^2:
        return g_vector[0]**2 + g_vector[1]**2

    def __update_hash(self, **kwargs):
        """Update hash dictionary with the parameters passed.

        Parameters
        ----------
        basis : iterable
            2x2 iterable of the bulk basis. Converted to a 2x2 tuple.
        superlattice : iterable
            Single superlattice matrix, or list of superlattice
            matrices. When a list, make sure that each of the elements
            is a tuple, so that a hash can be built. Only one
            superlattice matrix per domain is needed.
        domain_ids : iterable
            List of hashable entries that identify the domains.
        angle_key : str
            Identifier that defines the geometry of the primary beam.
        **kwargs : dict
            Other parameters are ignored

        Returns
        -------
        None.

        """
        missing = [key for key in self.hash_keys if key not in kwargs]
        if missing:
            warning(f"LEEDEquivalentBeams: missing {missing}. "
                    "No hash possible")
            return
        # Set up __hash_dict to always have the same order
        self.__hash_dict = {k: kwargs[k] for k in self.hash_keys}
        for k, hash_component in self.__hash_dict.items():
            if k in ('basis', 'superlattice'):
                hash_component = two_by_two_array_to_tuple(hash_component)
            if k in ('basis', 'superlattice', 'domain_ids'):
                self.__hash_dict[k] = tuple(hash_component)
