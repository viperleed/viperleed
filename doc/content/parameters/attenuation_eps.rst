.. _attenuation_eps:

ATTENUATION_EPS
===============

ATTENUATION_EPS (corresponds to FORTRAN parameter `TST` in beamgen, refcalc and
search routines) controls which beams will be considered or not when plane-wave
propagation occurs between neighboring layers.
Specifically, the program will (i) propagate all beams between two neighboring
layers, and (ii) keep all those beams that got attenuated by a factor of less
than ATTENUATION_EPS to interact with the next layer.

.. todo:: Add physics background info here.

.. todo:: Would be nice to have advice by Lutz.

**Default**: ATTENUATION_EPS = 0.001

.. admonition:: Syntax

   ::

      ATTENUATION_EPS = 0.005

**Acceptable values**: 1e-6 < ``ATTENUATION_EPS`` < 1


.. note::
  *  This parameter can usually be left to its default value.
     Tweaking (typically increasing) it can help convergence only in those cases
     in which the minimum interlayer distance goes below 1.0 Å.
  *  Smaller values will increase the calculation time as more beams will be
     considered.
  *  Small values of ATTENUATION_EPS can lead to numerical instabilities and
     the calculation of the beam intensities.


Changelog
---------

.. versionchanged:: TensErLEED 1.61
   The lower limit was changed from 1e-4 to 1e-6 (read by FORTRAN as a F7.4).
