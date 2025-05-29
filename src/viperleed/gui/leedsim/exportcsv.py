"""Module exportcsv of viperleed.gui.leedsim.

This module contains functions needed for exporting the pattern files
needed as input for the ViPErLEED ImageJ plug-ins.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-01-21'
__license__ = 'GPLv3+'

import numpy as np

from viperleed import __version__
from viperleed.gui.classes.beamindex import BeamIndex
from viperleed.gui.helpers import format_floats
from viperleed.gui.helpers import integer_part_length
from viperleed.gui.leedsim.classes.oldleedpatterns import LEEDPattern
from viperleed.gui.leedsim.classes.leedparameters import LEEDParameters


def export_pattern_csv(fnames, leeds, **kwargs):
    """
    Exports the pattern file needed for the ViPErLEED ImageJ plug-in for
    extracting experimental I(V) curves.
    
    Parameters
    ----------
    - fnames: str, or list of str
        Paths to files to be saved. When passing multiple LEED patterns, if only
        one file name is given, the pattern files of each single structure are
        combined there. If multiple paths are given, the first one is used for
        the combined pattern file, all the others for the single-structure
        files. Thus in this case (n+1) paths need to be given, where
        n = len(leeds).
    - leeds: list of dict(s) or viperleed.LEEDPattern(s)
        List determining what has to be exported. All LEED patterns must have
        the same bulk basis. No attempt to make them consistent by a unitary
        transformation (rotate/mirror/swap unit vectors) is made.
    - keyword arguments:
        * source: str
            Path to file (or any other useful reference to it) that was used
            as an input to generate the pattern file. In GUI, this is the .tlm
            file. When called from calc it may be related to the names given
            automatically by the bookkeeper.
        * name: str
            Descriptive text giving a name to the structure, independent of
            file name
        * domains: iterable of (iterable of int or None)
            the 'outer' list should contain as many entries as there are LEED
            parameters. Each entry is a list of indices (0-based) of the
            symmetry-equivalent domains to be exported. This is tricky to
            infer without looking at the GUI. If any of the lists is None, all
            symmetry-equivalent domains are passed. If None or not given, all
            symmetry-equivalent domains from all structures are exported.
    """
    # Check and process parameters:
    source = kwargs.get('source', None)
    name = kwargs.get('name', None)
    all_domains = kwargs.get('domains', [None]*len(leeds))
    
    if isinstance(fnames, (list, tuple, np.ndarray)):
        # will export one combined file, plus as many files as
        # there are patterns
        fname_combined = fnames[0]
        fnames = fnames[1:] if fnames[1:] else [None]*len(leeds)
        if len(fnames) != len(leeds):
            raise ValueError("exportcsv: not enough file paths given for "
                             "exporting multiple LEED patterns")
    elif isinstance(fnames, str):
        fname_combined = fnames
        fnames = [None]*len(leeds)

    if len(all_domains) != len(leeds):
        raise ValueError("exportcsv: The number of domains definition indices "
                         f"({len(all_domains)}) is inconsistent with that of "
                         f"LEED patterns ({len(leeds)}).")
    # TODO: the whole block that follows will be done with a single
    #       call to the new LEEDPattern
    if any(not isinstance(leed, (dict, LEEDPattern, LEEDParameters))
           for leed in leeds):
        raise TypeError("exportcsv: each of the LEED parameters passed must "
                        "be either a leed_parameters dictionary or a "
                        "viperleed.LEEDPattern")
    for i, leed in enumerate(leeds):
        if isinstance(leed, (dict, LEEDParameters)):
            # replace the dictionary entry with a LEEDPattern
            leeds[i] = LEEDPattern(leed)
        # and check that the bulk bases are consistent
        if not np.allclose(leed.bulk_basis, leeds[0].bulk_basis):
            raise ValueError("exportcsv: Incompatible bulk bases found among "
                             "the LEED patterns")

    # TODO: Also this part needs to change. I will do the loop below
    #       only if there are fnames to save, iterating through the
    #       leed.domains
    #       In all cases, I'll then call _format_beams_ on the full LEEDPattern
    #       to get the combination of all structures at once

    # Now go for the actual code
    all_lines = []
    for leed, domains, fname in zip(leeds, all_domains, fnames):               # with the new LEEDPattern, leeds -> leed.domains, so that each iteration brings out a LEEDSymmetryDomains
        lines = _format_beams_(leed,
                               domains=domains,
                               name=name,
                               source=source)
        if fname:  # write to file the single patterns
            with open(fname, 'w+') as f:
                f.write('\n'.join(lines))
        all_lines.append(lines)
    
    # and export the combined one
    combined = _combine_pattern_files_(all_lines)
    with open(fname_combined, 'w+') as f:
        f.write('\n'.join(combined))


def _combine_pattern_files_(all_lines):
    """
    TO BE IMPLEMENTED
    """
    return all_lines[0]


def _format_beams_(leed, **kwargs):
    """
    Format beams BLAH BLAH
    
    Parameters
    ----------
    - leed: LEEDPattern
    - keyword arguments:
        * domains: iterable of int
            Indices (0-based) of the domains to be exported. This is tricky to
            infer without looking at the GUI
        * others are passed to _format_header_
    """

    # Check parameters
    if not isinstance(leed, LEEDPattern):                                       # This will not be a LEEDPattern, but rather a LEEDSymmetryDomains
        raise TypeError("exportcsv: leed must be a "
                        "LEEDPattern instance")

    beams = leed.get_equivalentSpots(kwargs.get('domains', None))
    fractions, groups, overlaps, extinct_domains = zip(*beams)

    # Process the data to get them in a better format for packing each
    # line; also find the correct number of characters for each column
    hk = []  # beam indices
    gg = []  # reciprocal lattice vectors
    dd = []  # overlapping domains

    # lengths is a dictionary that will contain the maximum
    # number of characters appearing in each part of the columns
    lengths = {'numerator': [], 'denominator': [], 'hk_integer': [],
               'g_integer': [], 'domains': []}

    for fract in fractions:
        # construct beam (= hk indices) and g (= reciprocal-space vector
        #                                      in Cartesian coordinates)
        beam = BeamIndex(fract)
        hk.append(beam)
        g = np.dot(beam, leed.bulk_basis)
        gg.append(g)

        # get relevant lengths:
        # 1) fractional representation of beam indices
        num, den = beam.get_format_widths()
        lengths['numerator'].append(num)
        lengths['denominator'].append(den)

        # 2) float representation of beam indices
        lengths['hk_integer'].append(integer_part_length(*beam))

        # 3) and reciprocal-space vector
        lengths['g_integer'].append(integer_part_length(*g))

    # construct the list of overlapping domains
    for (overlap, extinct) in zip(overlaps, extinct_domains):
        overlap_txt = [str(dom) for dom in overlap]
        for e, (s, ext) in enumerate(zip(overlap_txt, extinct)):
            if ext:
                overlap_txt[e] = f"({s})"
        overlap_txt = '+'.join(dom for dom in overlap_txt)
        dd.append(overlap_txt)
        lengths['domains'].append(len(overlap_txt))

    # get the maximum of all lengths, which will determine the formatting
    for key in lengths.keys():
        lengths[key] = max(lengths[key])

    # and account also for the text in the column headers
    lengths['group'] = max(len(str(np.max(groups))), len('group'))
    lengths['domains'] = max(lengths['domains'], len('domain(s)'))

    to_export = [*_format_header_(lengths, leed, **kwargs)]

    for (beam, g, group, doms) in zip(hk, gg, groups, dd):
        n, d = lengths['numerator'], lengths['denominator']
        
        line = (    # list of entries for each column:
            f"{beam:({n},{d})s},"                              # fractional hk
            + f"{beam:{lengths['hk_integer']}f}"[1:-1] + ','   # floating hk
            + format_floats(f"{lengths['g_integer']}f", *g) + ','  # gx, gy
            + f"{group:>{lengths['group']}},"                  # group index
            + f"{doms:>{lengths['domains']}},"                 # domains overlap
            )

        to_export.append(line)

    return to_export


def _format_header_(lengths, leed, **kwargs):                                   # TODO: This needs to handle both leed=LEEDStructuralDomains and leed=LEEDSymmetryDomains
    """
    This function is used internally by export_pattern_file.

    Parameters
    ----------
    - lengths: dict
        Dictionary, generated by _format_beams_ with the following mandatory
        keys, used for formatting the column headers:
        * 'numerator'
        * 'denominator'
        * 'hk_integer'
        * 'g_integer'
        * 'group'
        * 'domains'
    - leed: LEEDPattern
        A full LEEDPattern that is used for exporting.
    - keyword arguments:
        * source: str
            Path to file (or any other useful reference to it) that was used
            as an input to generate the pattern file. In GUI, this is the .tlm
            file. When called from calc it may be related to the names given
            automatically by the bookkeeper.
        * name: str
            Descriptive text giving a name to the structure, independent of
            file name
        * domains: iterable of int
            Indices (0-based) of the domains to be exported. This is tricky to
            infer without looking at the GUI
    Returns
    -------
    - list of str: list of lines to be exported
    
    """
    # 1st Header
    # * ViPErLEED header
    # * reference to source used for exporting (path or other)
    #
    # 2nd Header:
    # * optional structure name
    # * max energy
    # * bulk shape (group); lattice parameters and basis
    # * surf shape (group); lattice parameters
    # * total number of domains, and number of domains exported
    # * for each exported domain: basis and superlattice matrix
    # * finally an uncommented header line for the columns

    # check mandatory parameters
    if not isinstance(leed, LEEDPattern):
        raise TypeError("exportcsv: leed must be a "
                        "LEEDPattern instance")
    needed_keys = ('numerator', 'denominator', 'hk_integer', 'g_integer',
                   'group', 'domains')
    if any((key not in lengths.keys()) for key in needed_keys):
        raise ValueError("exportcsv: some keys are missing from the required "
                         "parameter 'lengths'.")

    # Parse kwargs:
    source = kwargs.get('source', None)
    name = kwargs.get('name', None)
    export_domains = kwargs.get('domains', range(leed.n_domains))

    # Fill up the header
    header = ["# This file was automatically generated by ViPErLEED"
              f" v{__version__}"]
    if source:
        header.append("# It contains the pattern file generated from:\n# "
                      + source)
        header.append('#')
    if name:
        header.append(f"# Structure: {name}\n#")

    header.append(f"# Max. LEED Energy: {leed.max_energy:.1f} eV")
    
    # TODO: add incidence angles

    # TODO: the next block inevitably needs to change, as we have to be
    #       more general to handle multiple structures, so bulk and surface(s)
    #       must be treated separately.
    #       probably easiest: (1) when a LEEDSymmetryDomains is passed, just
    #       wrap it into a 1-element tuple. (2) Define an extra indent if there
    #       is more than one pattern to process. (3) The line "Surface
    #       lattice..." will need to be different if there is a single
    #       structure or multiple ones. Probably something like "Structural
    #       domain {name}: ..." in the latter case
    
    #       LEEDSymmetryDomains has a .bulk property (a Lattice2D), and self[0]
    #       is the first-domain Lattice2D

    for txt in ('Bulk', 'Surface'):
        lattice = leed.reciprocal_lattices[txt[:4].lower()]
        basis = lattice.reciprocal_basis  # TODO: probably better to use real_basis
        shape = lattice.cell_shape
        group = lattice.group.group
        basis_text = _format_lattice_basis_(basis, shape)
        header.append(f"# {txt} lattice: {shape} ({group}); {basis_text}")
        if txt == 'Bulk':
            header.append('#' + ' '*8 + _format_basis_vectors_(basis))

    basis = leed.reciprocal_lattices['bulk'].real_basis

    # TODO: the next part will just be a loop over the structures, using
    #       the extra indent from above
    matrices = leed.superlattices[export_domains]

    header.append(f"#\n# Exporting beams from {len(matrices)} domain(s) "
                  f"out of {leed.n_domains} symmetry-equivalent ones")

    for (dom, m) in zip(export_domains, matrices):
        super_basis = np.dot(m, basis)
        dom_txt = f" Domain {dom + 1} - "
        header.append(f"#{dom_txt}{_format_basis_vectors_(super_basis)}")
        header.append('#' + ' '*len(dom_txt) + 'Superlattice: '
                      + np.array2string(m, separator=' ').replace('\n',''))

    # Prepare the headers of the columns
    # 0) Some description
    header.append(
        '#\n'
        '# * h and k are the surface Miller indices of each LEED spot;'
        ' both fractional\n'
        '#   and floating-point versions are provided.\n'
        '# * gx and gy are the horizontal and vertical components of'
        ' reciprocal lattice\n'
        '#   vectors in AA^(-1), and include a factor of 2*pi.\n'
        '# * Beams with the same absolute value of \'group\' are'
        ' symmetry equivalent (at\n'
        '#   normal incidence); extinct spots have a negative \'group\''
        ' index.\n'
        '# * The domains contributing to each spot are listed in'
        ' \'domain(s)\'. Domains\n'
        '#   contributing with glide-extinct spots are reported in'
        ' parentheses.\n'
        )

    # then all column headers centered with the contents
    n, d = lengths['numerator'], lengths['denominator']
    fract_width = n + d

    # if the indices are fractions, add 1 to account for the "/"
    fract_width += 1 if lengths['denominator'] else 0

    column_headers = (
        f"({'h':^{fract_width}}|{'k':^{fract_width}}),"  # fractional (h|k)
        f"{'h':^{lengths['hk_integer'] + 6}},"           # float h
        f"{'k':^{lengths['hk_integer'] + 6}},"           # float k
        f"{'gx':^{lengths['g_integer'] + 6}},"           # float gx
        f"{'gy':^{lengths['g_integer'] + 6}},"           # float gy
        f"{'group':^{lengths['group']}},"                # group
        f"{'domain(s)':^{lengths['domains']}},"          # domains
        )

    header.append(column_headers)

    return header


def _format_lattice_basis_(basis, shape):
    a, b = np.linalg.norm(basis, axis=1)
    alpha = np.degrees(np.arccos(np.dot(basis[0], basis[1])/(a*b)))

    txt = [f"a = {a:.4f} AA"]

    if shape not in ['Square', 'Hexagonal', 'Rhombic']:
        txt.append(f"b = {b:.4f} AA")
    if shape not in ['Square', 'Rectangular']:
        txt.append(f"alpha = {alpha:.2f} deg")

    return '; '.join(txt)


def _format_basis_vectors_(basis):
    txt = (
        'Basis: '
        + 'a = {}'.format(np.array2string(basis[0], precision=4,
                                          floatmode='fixed', separator=' ',
                                          suppress_small=True))
        + '; '
        + 'b = {}'.format(np.array2string(basis[1], precision=4,
                                          floatmode='fixed', separator=' ',
                                          suppress_small=True))
        )
    return txt

