.. _lmax:

LMAX
====

LMAX defines the maximum angular momentum number to be used in the spherical
harmonics expansion, used for calculating the scattering matrices within each
layer. In general, this is determined by the
:ref:`PHASESHIFT_EPS<PHASESHIFTMIN>`  parameter, based on the maximum
phase-shift values at a given energy. The LMAX parameter can be used
to instead fix it to a constant value for all energies, or to define
upper and lower bounds.

**Default**: Determined by :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>`, minimum 6

**Syntax**:

::

   LMAX = 9     ! Use a constant LMAX=9 at all energies
   LMAX = 6-11  ! Base LMAX on PHASESHIFT_EPS, but keep between bounds 6 and 11. Equivalent: '6:11' or '6 11'

**Acceptable values**: one or two integer values between 1 and 18

Energy-dependent LMAX values are used in the reference calculation
(TensErLEED version 1.61 and higher). Delta calculations will use
the highest LMAX, filling missing matrix elements (from Tensors with
lower LMAX) with zeroes.

Because the :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>`  uses a relatively low cutoff
("fine") by default, most calculations can be accelerated by setting an upper
bound for LMAX.
