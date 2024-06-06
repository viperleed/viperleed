.. _refcalc-input:

Reference-calculation input files
=================================

Below is a short list of files used for the
:ref:`Reference Calculation<ref-calc>`.

.. _auxnonstruct:

AUXNONSTRUCT
------------
Non-structural information; contains the numbers of the beams to
calculate as listed in the :ref:`BEAMLIST<BEAMLIST>`, and the
parameters :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`,
:ref:`BULKDOUBLING_EPS<BULKDOUBLEEPS>`,
:ref:`BULKDOUBLING_MAX<BULKDOUBLEITER>`, and :ref:`LMAX<LMAX>`.

.. _auxlatgeo:

AUXLATGEO
---------
Lateral geometry information; contains the energy range to calculate
(from :ref:`THEO_ENERGIES<theo_energies>`), the surface and bulk unit-cell
vectors (from :ref:`POSCAR<POSCAR>`  and :ref:`SUPERLATTICE<SUPERLATTICE>`),
and :ref:`V0_Z_ONSET<INPOTZ>`.

.. _auxgeo:

AUXGEO
------
Geometry information; this essentially translates the input from
:ref:`POSCAR<POSCAR>`  and :ref:`VIBROCC<vibrocc>`, together with all
parameters relevant for interpretation of these files, into the format
expected by TensErLEED. Also contains information about which layer
should produce Tensor output (:ref:`TENSOR_OUTPUT<TOUTPUT>`).

refcalc-FIN
-----------
Compiled reference calculation input file that is actually fed into the
refcalc script. Combines (and is completely redundant with) the input
from AUXNONSTRUCT, AUXLATGEO and AUXGEO (see above), as well as
:ref:`AUXBEAMS<AUXBEAMS>` and :ref:`PHASESHIFTS<PHASESHIFTS>`.

muftin.f
--------
Only used up to TensErLEED v1.61. Fortran source code file. Contains the
subroutine defining muffin-tin potentials, based on the :ref:`V0_REAL<MUFTIN>`
and :ref:`V0_IMAG<v0_imag>`  parameters. Gets compiled at runtime, together
with the other refcalc fortran files.

refcalc-PARAM
-------------
Fortran source code file. Initializes fortran array sizes, which are
automatically determined from the rest of the input. Part of fortran
compilation at runtime.
