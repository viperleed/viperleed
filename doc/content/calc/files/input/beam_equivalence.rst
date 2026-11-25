.. include:: /substitutions.rst

.. _beam_equivalence:

================
BEAM_EQUIVALENCE
================

.. versionadded:: 0.15.0

The BEAM_EQUIVALENCE file is an optional input file that defines
explicitly which calculated beams are compared to each experimental
beam in the |R-factor| calculation.

Without a BEAM_EQUIVALENCE file, |calc| finds this correspondence
automatically: each experimental beam (from the
:ref:`EXPBEAMS.csv<EXPBEAMS>` file) is compared to the calculated beam
(from the :ref:`IVBEAMS` file) with the same indices or, failing that,
with a symmetry-equivalent one (see also :ref:`AVERAGE_BEAMS<averagebeams>`
and :ref:`BEAMINCIDENCE`). If a BEAM_EQUIVALENCE file is present, it
**replaces** the automatic correspondence entirely: only the beam
assignments given in the file are used, and none are added from
symmetry.

A BEAM_EQUIVALENCE file is useful, for example,

-   to compare one experimental beam to the average of several
    calculated beams, also when |calc| would not consider them
    symmetry equivalent, e.g., for rotational domains of an overlayer
    that are not equally populated (see :ref:`the example
    below<beam_equivalence_example>`);
-   to check or override the automatic correspondence, for instance if
    it is not what you expect. The correspondence that |calc| used is
    always written to the |SUPP| directory (see
    :ref:`below<beam_equivalence_supp>`), so it can serve as a starting
    point for your own BEAM_EQUIVALENCE file.


Syntax
======

Each line of the file is one **beam assignment**: an experimental beam
on the left-hand side of an equal sign, and the calculated beams
compared to it on the right-hand side, separated by ``+``:

..  code-block:: none

    ! experimental beam = calculated beams
    (1 0) = (1 0) + (-1 0) + (0 1) + (0 -1)
    (1/2 0) = (1/2 0)
    (1, 1) = 0.75*(1, 1) + 0.25*(1, -1)   % with assignment weights

The calculated beams on the right-hand side of a line are averaged
before they are compared with the experimental beam on the left-hand
side.

Beams
    Each beam is written as two indices in round brackets, separated by
    whitespace, a comma, or a vertical bar: ``(1 0)``, ``(1, 0)``, and
    ``(1 | 0)`` are the same beam. Indices may be integers, fractions
    such as ``1/2``, or decimals such as ``0.5`` or ``.5``, and may
    carry a sign (``(+1 -1/2)``). Decimals that differ from a fraction
    with a denominator of at most 99 by less than the rounding to six
    decimals are read as that fraction, as in the :ref:`IVBEAMS` file;
    for example, ``0.333333`` is read as ``1/3``. Numbers in exponential
    notation (e.g., ``5e-1``) are not accepted.

Assignment weights
    A calculated beam may be preceded by a weight and ``*``, e.g.,
    ``0.75*(1 1)`` or ``1/3*(1 0)``. A beam without a weight has weight
    1. Weights are non-negative integers, decimals or fractions; signs
    and exponential notation are not accepted. The weights in each line
    are normalised to sum to one, so ``2*(1 0) + 2*(0 1)`` is the same
    as ``(1 0) + (0 1)``, and ``3*(1 1) + (1 -1)`` is the same as
    ``0.75*(1 1) + 0.25*(1 -1)``. A line whose weights are all equal is
    an unweighted average. A calculated beam with weight zero is left
    out of the line, as if it was not written. See
    :ref:`below<beam_equivalence_weights>` for the |R-factor|
    calculations that currently accept unequal weights.

Comments and empty lines
    Everything after ``!``, ``#``, or ``%`` is a comment. Empty lines
    are ignored.

The file cannot be read, and |calc| stops with an error naming the
file and, where applicable, the line, if

-   a line is not of the form described above (e.g., it has no or more
    than one equal sign, a beam label has not exactly two indices, or a
    weight is negative);
-   an experimental beam is on the left-hand side of more than one
    line;
-   a calculated beam appears more than once on the same right-hand
    side;
-   all the weights of a line are zero;
-   the file contains no beam assignment at all.

The same calculated beam may appear on the right-hand side of several
lines (but see :ref:`below<beam_equivalence_weights>`).


.. _beam_equivalence_example:

Example: rotational domains
===========================

Consider a (2×1) overlayer on a square substrate that forms two
rotational domains, rotated by 90° with respect to each other, with
70 % and 30 % of the surface area. If the beam incidence is
perpendicular to the surface, the |IV| curves of the second domain are
those of the first one with the beam indices rotated by 90°. Instead of
running a :ref:`domain calculation<domain_calculation>`, it is then
sufficient to calculate the first domain only, and to assign its beams
to the experimental beams as follows:

..  code-block:: none

    # Integer-order beams: weighted average of the two domains
    (1 0) = 0.7*(1 0) + 0.3*(0 -1)
    (0 1) = 0.7*(0 1) + 0.3*(1 0)
    (1 1) = 0.7*(1 1) + 0.3*(1 -1)
    # Fractional-order beams: each comes from one domain only
    (1/2 0) = (1/2 0)
    (0 1/2) = (1/2 0)       ! (1/2 0) of the rotated domain
    (1/2 1) = (1/2 1)
    (1 1/2) = (1/2 -1)      ! (1/2 1) of the rotated domain

The intensities of the fractional-order beams scale with the area of
the domain they come from, but the |R factor| does not depend on the
overall scale of each curve (see :ref:`r-factor_calculation`).

Note that this file uses both unequal assignment weights and calculated
beams that appear in more than one line. See
:ref:`below<beam_equivalence_weights>` for the calculations that
currently accept such a file.


.. _beam_equivalence_weights:

Assignment weights and shared calculated beams
==============================================

The |R-factor| calculations currently available in |calc| average the
calculated beams of each experimental beam with equal weights, and
compare each calculated beam to at most one experimental beam.
Therefore, a BEAM_EQUIVALENCE file that contains

-   a line with unequal assignment weights, or
-   a calculated beam on the right-hand side of more than one line

is currently not accepted for any run that computes an |R factor|.
Unless no experimental data are given (no :ref:`EXPBEAMS.csv<EXPBEAMS>`
file), the :ref:`initialization` then stops with an error that names
the cause and lists the offending lines (with normalised weights) if
the :ref:`RUN` parameter contains any of the sections reference
calculation (``1``), structure search (``3``), superpos (``31``),
domain search (``4``), error calculation (``5``), full-dynamic
optimization (``6``), or |R-factor| calculation (``11``, ``12``).
This is the case for both values of
:ref:`R_FACTOR_LEGACY<rfactorlegacy>`. Such a file can still be used,
e.g., with ``RUN = 0`` (initialization only) or for a
:ref:`delta-amplitudes calculation<sec_deltas>` alone (``RUN = 2``).
Weights and shared calculated beams are planned to become available
with the Python |R-factor| implementation.


How the file is used
====================

The BEAM_EQUIVALENCE file is read during the :ref:`initialization`. Its
experimental beams are matched to those of the
:ref:`EXPBEAMS.csv<EXPBEAMS>` file, and its calculated beams to those of
the :ref:`BEAMLIST` and :ref:`IVBEAMS` files:

-   A line whose experimental beam is not in the
    :ref:`EXPBEAMS.csv<EXPBEAMS>` file is ignored, with a warning. If
    this leaves no line at all, |calc| stops with an error (also for
    ``RUN = 0``).
-   An experimental beam in the :ref:`EXPBEAMS.csv<EXPBEAMS>` file that
    is on no left-hand side is not used in the |R-factor| calculation,
    with a warning.
-   If a calculated beam on a right-hand side is not in the
    :ref:`BEAMLIST`, |calc| stops with an error listing these beams.
-   If no :ref:`IVBEAMS` file is given, it is generated from all the
    calculated beams on the right-hand sides (rather than from the
    :ref:`EXPBEAMS.csv<EXPBEAMS>` file). If an :ref:`IVBEAMS` file is
    given but lacks some of these beams, |calc| stops with an error
    listing them. Beams of the :ref:`IVBEAMS` file that are on no
    right-hand side are calculated, but not used in the |R factor|;
    they are listed in the log at the ``INFO`` level.
-   The :ref:`AVERAGE_BEAMS<averagebeams>` parameter has no effect, and
    a warning is issued if it is given in the :ref:`PARAMETERS` file.
-   Experimental or calculated beams are never merged because they are
    symmetry equivalent. Beams with the same indices are still treated
    as duplicates: two experimental beams with the same indices cause a
    warning that stops the run at the default :ref:`HALTING` level.

In a :ref:`domain calculation<domain_calculation>`, the
BEAM_EQUIVALENCE file in the main input folder is used. It assigns the
calculated beams, after averaging over the domains (by their area
fractions), to the experimental beams, exactly as for a single domain.
BEAM_EQUIVALENCE files in the input folders of the individual domains
are ignored.


.. _beam_equivalence_supp:

Beam correspondence in SUPP
===========================

Whenever experimental data are given, |calc| writes the beam
correspondence it uses to the file :file:`SUPP/BEAM_EQUIVALENCE`,
whether it was found automatically or read from a BEAM_EQUIVALENCE file
(in the latter case, without the lines that were ignored). The file is
not written if no experimental beam has a calculated beam assigned.
Indices are written as integers or fractions (e.g., ``(1/2 -1)``), and
weights are left out of lines whose weights are all equal.

Copying this file to the input folder and running again reproduces the
same |R-factor| input. It does not necessarily reproduce every warning
or the same :ref:`HALTING` level. For example, if the
:ref:`EXPBEAMS.csv<EXPBEAMS>` file contains symmetry-equivalent beams,
the automatic correspondence discards them with a warning that stops
the run at the default :ref:`HALTING` level, while with the
BEAM_EQUIVALENCE file they are on no left-hand side and are left out
with a weaker warning.
