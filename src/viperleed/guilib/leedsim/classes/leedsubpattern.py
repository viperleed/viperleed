"""Module leedsubpattern of viperleed.guilib.leedsim.classes.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the CachedLEEDSubpatternConstructor and LEEDSubpattern classes,
used for the plotting of a LEED pattern

Author: Michele Riva
Created: 2021-03-13
"""

from collections import defaultdict
import copy

import numpy as np
from matplotlib import (cm as color_maps, colors as mpl_colors)

from viperleed.guilib.helpers import two_d_iterable_to_array

from viperleed.guilib import decorators as dev_


class CachedLEEDSubpatternConstructor:
    """Construct LEEDSubpatterns caching instances."""

    hash_keys = ('leed', 'domains_ids',)

    __cache = {}

    def __new__(cls, leed, **kwargs):
        """Class instance constructor.

        The only purpose of this is to check whether the instance
        that is about to be created is already cached, in which
        case, the cached copy is returned.

        Parameters
        ----------
        leed : LEEDPattern
            LEEDPattern instance that will be used to build
            subpatterns. Used for computing hash.
        domains_ids : iterable, keyword only
            List of identifying numbers for the domains of leed
            for which we construct LEEDSubpatterns. Used
            for computing hash. No caching is possible if this
            is not given.
        **kwargs
            Other kewyword arguments are ignored.

        Returns
        -------
        CachedLEEDSubpatternConstructor
            New instance

        Raises
        ------
        TypeError
            If leed is not a LEEDPattern instance
        ValueError
            If domains_ids is not given or it is None
        """
        try:  # Don't use isinstance to avoid cyclic imports
            _ = leed.bulk, leed.domains, leed.parameters
        except AttributeError:  # Probably not a LEEDPattern
            raise TypeError(f"Invalid 'leed' type {type(leed).__name__}. "
                            "Should be a 'LEEDPattern'")
        domains_ids = kwargs.get('domains_ids', None)
        if domains_ids is None:
            raise ValueError("Mandatory keyword argument domains_ids"
                             "missing. Cannot construct subpatterns.")

        instance = super().__new__(cls)
        hash_dict = {'leed': leed, 'domains_ids': tuple(domains_ids.items())}
        object.__setattr__(instance, f'_{cls.__name__}__hash_dict', hash_dict)
        if cls.__is_cached(instance):
            return cls.__cache[hash(instance)]
        return instance

    @dev_.exec_time
    def __init__(self, leed, **kwargs):
        """Initialize class instance.

        Parameters
        ----------
        leed : LEEDPattern
            LEEDPattern instance that will be used to build
            subpatterns. Used for computing hash.
        domains_ids : iterable, keyword only
            List of identifying numbers for the domains of leed
            for which we construct LEEDSubpatterns. Used
            for computing hash. No caching is possible if this
            is not given.
        **kwargs
            Other kewyword arguments are ignored.

        Returns
        -------
        None.
        """
        # __new__ already checks that the argument is there
        domains_ids = kwargs['domains_ids']

        # The following two lines are completely useless, but prevent
        # static code checkers from going nuts about the 'missing'
        # .__hash_dict instance attribute (that is actually already
        # created in __new__ with the same value)
        hash_dict = {'leed': leed, 'domains_ids': tuple(domains_ids.items())}
        self.__hash_dict = hash_dict

        if self.__is_cached():
            return

        beam_dicts = self.__split_beams_for_subpatterns()

        self.__subpatterns = []
        for domains_id, beams in beam_dicts.items():
            self.__subpatterns.append(LEEDSubpattern(leed, beams,
                                                     domains_id=domains_id))
        self.__cache[hash(self)] = self

    @property
    def patterns(self):
        """Return the subpatterns.

        Returns
        -------
        list
            List of LEEDSubpattern
        """
        return self.__subpatterns

    @staticmethod
    def clear_chache():
        """Completely empty the cache.

        Call this method when the bulk lattice is changed, before
        creating a new LEEDPattern instance.
        """
        CachedLEEDSubpatternConstructor.__cache = {}

    # @dev_.profile_lines
    @dev_.exec_time
    def __split_beams_for_subpatterns(self):
        """Sort out beams into groups appropriate for subpatterns.

        Also converts indices into reciprocal-space vectors

        Returns
        -------
        dict
            One dictionary per subpattern with

            keys : {'__s', '__s__e', struct_id+sym_id, struct_id+sym_id+'__e'}
                identifier for superposed ('__s') or strctural/symmetry
                domains, possibly with their extinct portions (+'__e')
            values : numpy.ndarray
                List of reciprocal-space vectors
        """
        leed = self.__hash_dict['leed']

        # Get the beams from the LEEDPattern. Use the dictionary
        # indexed_beams, with keys {'beams', 'group_indices',
        # 'overlap_domains', 'extinct_domains'}.
        # group_idx < 0 for fully extinct, and
        # overlap_ids = [(struct_id1, (..., sym_ids, ...)), ...]
        # indexed_beams = leed.domains.indexed_beams  # OLD
        indexed_beams = zip(*leed.domains.indexed_beams.values())

        # Go through each beam, and split with the following
        # logics:
        #     * superposed, non-extinct beams belong to their
        #       own subpattern  --> special key: '__s'
        #     * superposed, extinct beams belong to their own
        #       subpattern      --> special key: '__s' + '__e'
        #     * single-domain, non-extinct beams belong to
        #       one subpattern per each domain
        #                       --> key:  struct_id + sym_id
        #     * single-domain, extinct beams belong to
        #       one subpattern per each domain
        #                       --> key:  struct_id + sym_id + '__e'

        subpattern_beams = defaultdict(list)
        for beam_index, group_index, overlap_indices, _ in indexed_beams:
            # Check if multiple domains overlap
            if len(overlap_indices) > 1:
                # Multiple structural domains
                key = '__s'
            elif len(overlap_indices[0][1]) > 1:
                # Multiple symmetry domains of one structural domain
                key = '__s'
            else:
                # Only one symmetry+structural domain.
                # overlap_indices = [(struct_id, (sym_id,)), ]
                idx = overlap_indices[0]
                key = idx[0] + str(idx[1][0])
            if group_index < 0:
                # fully extinct
                key += '__e'
            subpattern_beams[key].append(beam_index)

        # Now get the actual reciprocal-space vectors from
        # the beam indices by taking the dot product with
        # the bulk basis and dividing by the denominator,
        # if needed.  Notice that it is better to do this
        # now rather than earlier, although we are looping
        # through the subpatterns, because this way the
        # beams in each subpattern used for instantiation
        # are already contigous in memory.
        # This way, LEEDSubpattern.transform_beams, which
        # is likely going to be called for plotting, does
        # not need to do any pre-copying of the arrays
        scale = 1
        if not leed.domains.fractional:
            scale = 1/leed.domains[0].denominator_for_bulk_beams
        for k, beams in subpattern_beams.items():
            beams = scale*two_d_iterable_to_array(beams,
                                                  dtype=float,
                                                  shape=(-1, 2))
            subpattern_beams[k] = np.einsum('ij,jk->ik',
                                            beams, leed.bulk.basis)

        return subpattern_beams

    def __hash__(self):
        """Return hash of self."""
        missing = (k for k in self.hash_keys if k not in self.__hash_dict)
        if any(missing):
            return 0
        return hash(tuple(v for v in self.__hash_dict.values()))

    def __is_cached(self):
        """Check if self is already in cache."""
        return hash(self) in self.__cache


class LEEDSubpattern:
    """Sub-pattern for a LEED. Used for plotting with matplotlib.

    A subpattern is a uniform collection of LEED beams with
    specific characteristics for plotting, i.e., a marker,
    a color, and a marker scaling factor (determining how the
    size of a marker scales with energy).  There exist subpatterns
    for the bulk reciprocal lattice, for beams that originate
    from the superposition of multiple domains (structural and/or
    symmetry-induced), and for beams originatig exclusively from
    a single structural+symmetry domain.
    """

    # Eventually one could even pull in the defaults from a config
    # file, probably as a @staticmethod to call the first time
    # the first instance is created?
    defaults = {'__b': {'facecolors': 'none',  # bulk
                        'edgecolors': 'k',
                        'marker': 'o',
                        'alpha': None,
                        'scale': 25.},  # spot area in pts**2
                '__s': {'facecolors': 'gray',  # superposed
                        'edgecolors': 'face',
                        'marker': 'o',
                        'alpha': None,
                        'scale': 1.},
                'single': {'facecolors': 'k',  # 1 domain overall
                           'edgecolors': 'face',
                           'marker': 'o',
                           'alpha': None,
                           'scale': 1.},
                'other': {'facecolors': None,  # 1 domain of many
                          'edgecolors': 'face',
                          'marker': 'o',
                          'alpha': None,
                          'scale': 1.},
                '__e': {'marker': 'x',         # extinct, except bulk
                        'alpha': 0.2,
                        'scale': 2.},
                'cmap': 'gnuplot'}             # multi-domain colors

    __init_keys = ('cmap', 'color', 'alpha')

    def __init__(self, leed, beams, domains_id, **scatter_format):
        """Initialize LEEDSubpattern.

        Parameters
        ----------
        leed : LEEDPattern
            LEED pattern to which this subpattern belongs to.
        beams : numpy.ndarray
            List of unrotated reciprocal-space vectors belonging
            to this subpattern.  Should have shape (2, N) or (N, 2).
            The correct axis is selected automatically.
        domains_id : str
            Domains identifier specifying which domain(s),
            compose this subpattern. It can be:
                * '__b'
                    bulk subpattern
                * '__s'
                    subpattern of superposed beams
                * struct_id + symm_id
                    representing a symmetry domain of a structural
                    domain
            Additionally, each key may contain a '__e' as the last
            3 characters indicating that the subpattern contains
            only extinct beams.
        cmap : str or matplotlib.colors.Colormap, optional
            color map used to pick the colors for each structural
            and/or symmetry domain (default=matplotlib.cm.gnuplot).
            Notice that the cmap will be used to set the colors for
            ALL domains in the LEED pattern, and the colors for the
            domains of this subpattern will be picked from the full
            list. Use facecolors/edgecolors to specify which colors
            should be used only for the domains currently drawn.
            If cmap is given, it overrides facecolors/edgecolors.
        color : iterable, optional
            (R, G, B[, A]) color to be used for the subpattern. Use
            'alpha' for setting the transparency. Default is: (i)
            for single domains in a multi-domain LEED pattern,
            take color from cmap, (ii) gray for superposed beams,
            (iii) black for bulk (used for edge only), and (iv)
            black for a single-domain pattern.
        alpha : float, optional
            Sets the transparency of the colors of markers (1 for
            opaque, 0 for invisible, None is the same as 1).  Takes
            precedence over the value given in color.  Default:
            1 for non-extinct beams, 0.2 for extinct ones.
        **kwargs
            Other keyword arguments are ignored

        Returns
        -------
        None.

        """
        self.__leed = leed
        self.__beams = beams
        self.__transformed_beams = beams
        self.__domains_id = domains_id
        self.__scatter_format = {k: v
                                 for k, v in scatter_format.items()
                                 if k in self.__init_keys}
        self.__all_leed_colors = None
        self.__process_arguments()
        self.__process_color_arguments()
        self.__select_scatter_format()

    @property
    def beams(self):
        """Return the beams in the latest coordinate system.

        Returns
        -------
        numpy.ndarray
            List of reciprocal-space vectors
        """
        return self.__transformed_beams

    @property
    def leed_colors(self):
        """Return the colors for all the domains in the LEED pattern.

        This is independent of whether facecolors and/or edgecolors
        have been given in the constructor.

        Returns
        -------
        colors : dict
            A mapping between domain ids and their colors when multiple
            domains are visualized in the LEED pattern.

            keys : str
                struct_id + symm_id
            values : tuple
                (r, g, b, a) values, each [0, 255]
        """
        if self.__all_leed_colors is not None:
            return self.__all_leed_colors

        # Get the correct domain ids
        struct_ids = self.__leed.domains.domain_ids
        all_symm_ids = [dom.domain_ids for dom in self.__leed.domains]
        keys = []
        for struct_id, symm_ids in zip(struct_ids, all_symm_ids):
            for symm_id in symm_ids:
                keys.append(struct_id + str(symm_id))

        # Get the right colors
        n_domains = len(keys)
        if n_domains == 1:
            colors = (self.defaults['single']['facecolors'],)
        else:
            cmap = self.__scatter_format.pop('cmap', self.defaults['cmap'])
            cmap = color_maps.get_cmap(cmap)
            colors = cmap(np.linspace(0.1, 0.9, n_domains))

        self.__all_leed_colors = dict(zip(keys, colors))
        return self.__all_leed_colors

    @property
    def marker_scaling(self):
        """Return scaling factor for marker.

        Returns
        -------
        float
        """
        return self.__scatter_format['scale']

    @property
    def scatter_plot_settings(self):
        """Return the color/marker settings for scatter-plotting.

        Returns
        -------
        dict
        """
        return {k: v for k, v in self.__scatter_format.items() if k != 'scale'}

    def transform_beams(self, transform_matrix):
        """Apply coordinate transform to the beams.

        Applies the given matrix transformation to the ORIGINAL
        beams passed during construction, NOT to the latest
        transformed beams.  After calling this function, self.beams
        is updated to the latest transformed coordinate system.

        Parameters
        ----------
        transform_matrix : array-like
            2x2 transformation matrix to be applied to the beams.
            The transform is applied on the right, i.e., each of the
            transformed beams will be (bx, by) @ transform_matrix.

        Raises
        ------
        ValueError
            if the transformation is not a 2x2 matrix.
        """
        if np.shape(transform_matrix) != (2, 2):
            raise ValueError("Transformation matrix should be a 2x2"
                             "array-like. Found shape "
                             f"{np.shape(transform_matrix)} instead.")
        self.__transformed_beams = np.dot(self.__beams, transform_matrix)

    def __process_arguments(self):
        """Check and process the initializer input parameters.

        Parameters
        ----------
        leed : LEEDPattern
            LEED pattern to which this subpattern belongs to.
        beams : numpy.ndarray
            List of unrotated reciprocal-space vectors belonging
            to this subpattern.  Should have shape (2, N) or (N, 2).
        domains_ids : iterable of str
            List of domain identifiers specifying which domain(s),
            compose this subpattern. Each one can be:
                * '__b'
                    bulk subpattern
                * '__s'
                    subpattern of superposed beams
                * struct_id+symm_id
                    representing a symmetry domain of a structural
                    domain
            Additionally, each key may contain a '__e' as the last
            three characters indicating that the subpattern contains
            only extinct beams.

        Returns
        -------
        None.

        Raises
        ------
        TypeError
            If any of the input types does not match the specification.
        ValueError
            If the shape of beams is not (2, N) or (N, 2).
        """
        leed = self.__leed
        try:  # Don't use isinstance to avoid cyclic imports
            _ = leed.bulk, leed.domains, leed.parameters
        except AttributeError:  # Probably not a LEEDPattern:
            raise TypeError("Invalid argument 'leed'. Expected 'LEEDPattern' "
                            f"found {type(self.__leed).__name__}")
        if not isinstance(self.__domains_id, str):
            raise TypeError("Invalid argument 'domains_id'. "
                            "Expected 'str', not "
                            f"{type(self.__domains_id).__name__}")
        if not hasattr(self.__beams, '__len__'):
            raise TypeError("Invalid argument 'beams'. Expected "
                            "'numpy.ndarray', not "
                            f"{type(self.__beams).__name__}")
        self.__beams = np.asarray(self.__beams)
        if (len(self.__beams.shape) != 2
                or all(s != 2 for s in self.__beams.shape)):
            raise ValueError("Invalid shape of 'beams'. Expected a 2xN or "
                             f"Nx2 array. Found {self.__beams.shape}")

    def __process_color_arguments(self):
        """Figure out the right colors from the input.

        At the end of this, self.__scatter_format has:
            * only 'cmap' if 'cmap' is appropriate, no 'cmap'
              if inappropriate;
            * 'colors' if no 'cmap', and 'colors' was appropriate;
            * 'alpha' if appropriate

        Returns
        -------
        None.

        Raises
        ------
        TypeError
            If alpha cannot be converted to a number,
            or if color is not a Sequence
        ValueError
            If float(alpha) is not between 0 and 1,
            if color does not have exactly three elements,
            or if any of the three elements is not and
            integer between 0 and 255
        """
        alpha = self.__scatter_format.get('alpha', None)
        if alpha is None:
            self.__scatter_format.pop('alpha', None)
        else:
            try:
                alpha = float(alpha)
            except TypeError as err:
                raise TypeError("Invalid 'alpha'. Expected 'float', found "
                                f"{type(alpha).__name__}") from err
            if alpha < 0 or alpha > 1:
                raise ValueError(f"Invalid 'alpha'={alpha}. Should be "
                                 "between 0 and 1")

        cmap = self.__scatter_format.get('cmap', None)
        color = self.__scatter_format.get('color', None)

        # cmap given:
        if cmap is not None:
            color_maps.get_cmap(cmap)  # raises errors if invalid
            self.__scatter_format.pop('facecolors', None)
            self.__scatter_format.pop('edgecolors', None)
            return

        # cmap not given or None:
        self.__scatter_format.pop('cmap', None)

        #     color not given:
        if color is None:
            self.__scatter_format.pop('color', None)
            return

        #     color given
        if not hasattr(color, '__len__'):
            raise TypeError("Invalid 'color'. Expected list-like, "
                            f"found {type(color).__name__}")
        if len(color) < 3:
            raise ValueError("Invalid 'color'. Expected 3 entries,"
                             f"not {len(color)}")
        color = color[:4]
        if any(c > 255 or c < 0 for c in color):
            raise ValueError("Color components should be integers "
                             "between 0 and 255")
        self.__scatter_format['color'] = tuple(round(c) for c in color)

    def __select_scatter_format(self):
        """Get the scatter-plot format appropriate for the domains.

        Processes the information in self.__scatter_format and
        self.__domains_id to select the correct colors for the
        current subpattern. The correct format is stored in
        self.__scatter_format

        Returns
        -------
        None.
        """
        dom_key = self.__domains_id[:3]
        ext_key = self.__domains_id[-3:]
        if dom_key not in ('__b', '__s'):
            # Need to pick either 'single' or 'other'
            if len(self.leed_colors) == 1:
                dom_key = 'single'
            else:
                dom_key = 'other'
        # Pull in the correct default, using a deepcopy to prevent
        # modification of the class attribute further down
        scatter_format = copy.deepcopy(self.defaults[dom_key])

        # Pull in the default for extinct, if appropriate
        if ext_key == '__e':
            extinct_format = self.defaults['__e']
            if dom_key == '__b':
                scatter_format['alpha'] = extinct_format['alpha']
            else:
                scatter_format.update(extinct_format)

        if dom_key == '__b':
            color_key = 'edgecolors'
        else:
            color_key = 'facecolors'

        # Now pull in the color input
        if 'color' in self.__scatter_format:
            color = self.__scatter_format['color']
        elif dom_key == 'other':
            color = self.leed_colors[self.__domains_id.replace('__e', '')]
        else:
            color = scatter_format[color_key]

        alpha = self.__scatter_format.get('alpha', None)
        color = (mpl_colors.to_rgba(color, alpha),)
        scatter_format[color_key] = color

        self.__scatter_format = scatter_format
