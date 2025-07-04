.. _delta-input:

Delta-amplitudes calculation input
==================================

Since the input files for the delta calculations are highly redundant
(one "delta calculation" is performed for each separate delta file),
input is collected in the delta-input file.
The input files listed are of two different categories:

-   PARAM files: Define static fortran array sizes.
    Every different instance of PARAM requires re-compiling the delta
    calculation fortran code, so delta input is grouped by files that
    can run based on the same PARAM files.
    In the delta-input code, PARAM is printed whenever it changes, and
    then proceeded by all input files that can use this PARAM.
-   Input for specific delta calculations:
    Contains the energy range to calculate (from :ref:`THEO_ENERGIES`),
    the surface and bulk unit cell vectors (from :ref:`POSCAR` and
    :ref:`SUPERLATTICE`), and the :ref:`BEAMINCIDENCE`.
    Then contains the entire content of the :ref:`AUXBEAMS` and
    :ref:`PHASESHIFTS` files, which are abbreviated simply as
    [AUXBEAMS] and [PHASESHIFTS] in the delta-input file.
    Then assigns an index of blocks in the :ref:`PHASESHIFTS` file to use,
    and lists all geometric and vibration displacements, as given in the
    :ref:`DISPLACEMENTS` file.
