.. _screen_aperture:

SCREEN_APERTURE
===============

SCREEN_APERTURE defines the aperture (in degrees) of the acceptance cone of the LEED screen. This, in turn, determines the highest-order beams that will be visible in the LEED experiment.

**Default**: SCREEN_APERTURE = 110 (i.e., the one of the ErLEED optics)

**Syntax**:

::

   SCREEN_APERTURE = 95

**Acceptable values**: integer/floating point number between 0 and 180

**Notes**: a beam with indices (h, k) and corresponding in-plane reciprocal vector **g** = h **a**\ \* + k **b**\ \* is visible at energy *E* (in eV) if *θ* < SCREEN_APERTURE / 2, where sin\ *θ* = *ħ g* / √(2 *m*\ :sub:`e` *e* *E*)
