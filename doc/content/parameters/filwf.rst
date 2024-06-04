.. include:: /substitutions.rst

.. _filwf:

===========
FILAMENT_WF
===========

FILAMENT_WF defines the workfunction of the emitting cathode (in
electronvolts). This number will be added to the electron energy
when calculating the muffin-tin potentials.

**Default**: FILAMENT_WF = 2.65 (LaB6)

**Syntax**:

::

   FILAMENT_WF = LaB6
   FILAMENT_WF = W
   FILAMENT_WF = <wf>

**Acceptable**: a single floating-point number (``<wf>`` above). Also the
special keys LaB6 (2.65 eV) and W (4.5 eV) are allowed.

**Note**: An erroneous setting of this parameter is not critical since
the search/|R-factor| algorithm will anyway allow a relative shift of
the calculated |IV| curves with respect to the experimental ones to find
the best match. However, this can also be used as a cross-check at the end
of the search: If the the filament workfunction is set correctly, any
significant energy shift in the fit will indicate that something is wrong.
