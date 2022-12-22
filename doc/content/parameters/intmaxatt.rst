.. _intmaxatt:

ATTENUATION_EPS
===============

ATTENUATION_EPS (corresponds to FORTRAN parameter TST in beamgen, 
refcalc and search routines) controls which beams will be considered 
or not when plane-wave propagation occurs between neighboring layers. 
Specifically, the program will (i) propagate all beams between two 
neighboring layers, and (ii) keep all those beams that got attenuated by
a factor of less than ATTENUATION_EPS to interact with the next layer.

**TODO**: This may be wrong. I (MR) think that all beams whose intensity
had a relative change I_propagated/I_initial > ATTENUATION_EPS will be 
kept! To be checked in the code.
AI: TensErLEED has conflicting comments about it.
Tt also mentioned in subroutine BEAMS, where it says beams with imaginary out-of-plane wave vector > TST are considered evanescent.

**Default**: ATTENUATION_EPS = 0.001

.. admonition:: Syntax

   ::

      ATTENUATION_EPS = 0.005

**Acceptable values**: 1e-6 < ``attenuation`` < 1 (exception: 
the lower limit is 1e-4 for TensErLEED versions 1.61 and lower).
Typical values < 0.005.


.. note::
  *  This parameter can usually be left to its default value. Tweaking it can help convergence only in those cases in which the minimum interlayer distance goes below 1.0 Å.
  *  The minimum value is 1e-4 up to TensErLEED 1.61 because FORTRAN reads it as a F7.4
