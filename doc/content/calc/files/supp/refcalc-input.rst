.. include:: /substitutions.rst

.. _refcalc-input:

Reference-calculation input files
=================================

Below is a short list of files used for the :ref:`ref-calc`.

.. _auxnonstruct:

AUXNONSTRUCT
------------
Non-structural information; contains the numbers of the beams to
calculate as listed in the :ref:`BEAMLIST`, and the parameters
:ref:`BEAMINCIDENCE`, :ref:`BULKDOUBLEEPS`, :ref:`BULKDOUBLEITER`,
and :ref:`LMAX`.

.. _auxlatgeo:

AUXLATGEO
---------
Lateral geometry information; contains the energy range to calculate
(from :ref:`THEO_ENERGIES`), the surface and bulk unit-cell vectors
(from :ref:`POSCAR` and :ref:`SUPERLATTICE`), and :ref:`v0_z_onset`.

.. _auxgeo:

AUXGEO
------
Geometry information; this essentially translates the input from :ref:`POSCAR`
and :ref:`VIBROCC`, together with all parameters relevant for interpretation of
these files, into the format expected by TensErLEED. Also contains information
about which layer should produce Tensor output (:ref:`TENSOR_OUTPUT`).

.. _refcalc-fin:

refcalc-FIN
-----------
Contains the input for the TensErLEED reference calculation. Combines (and is
completely redundant with) the input from :ref:`AUXLATGEO`, :ref:`BEAMLIST`,
:ref:`AUXNONSTRUCT`, :ref:`PHASESHIFTS`,  :ref:`AUXBEAMS` and :ref:`AUXGEO`.

The refcalc-FIN file does not have to be in the folder for the reference
calculation to run, but rather is piped directly into TensErLEED from |calc|.

.. todo::
    This was the list before. Check which one is the right order!

    :ref:`AUXNONSTRUCT`, :ref:`AUXLATGEO`, and :ref:`AUXGEO`, as
    well as :ref:`AUXBEAMS` and :ref:`PHASESHIFTS`.

muftin.f
--------
Only used up to TensErLEED v1.61. Fortran source code file. Contains the
subroutine defining muffin-tin potentials, based on the :ref:`MUFTIN`
and :ref:`V0_IMAG` parameters. Gets compiled at runtime, together
with the other refcalc fortran files.

refcalc-PARAM
-------------
Fortran source code file. Initializes fortran array sizes, which are
automatically determined from the rest of the input. Part of fortran
compilation at runtime.
