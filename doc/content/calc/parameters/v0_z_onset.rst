.. include:: /substitutions.rst

.. _v0_z_onset:

V0_Z_ONSET
==========

V0_Z_ONSET (original Fortran parameter ASE) defines the *z* position (in
Ångström) of the onset of the inner potential of the solid, i.e., the plane
that separates "vacuum" from the "material". Notice that, in this case,
POSITIVE values will move the plane AWAY from the solid. ``V0_Z_ONSET = 0``
corresponds to the *z* coordinate of the topmost atom.

**Default**: V0_Z_ONSET = 1.0

**Syntax**:

::

   V0_Z_ONSET = 0.8

**Acceptable values**: any non-negative floating-point number. Typically
:math:`0 \leq`  ``V0_Z_ONSET``  :math:`\leq 2`.

**Notes:** In most cases the default value should be accurate enough. The
only exceptions are those cases in which the topmost layer is severely rough
(e.g., sparse adsorbates with large free areas in between). In this case, it
might be a good idea to try out a few different values: moving the plane closer
to the surface (i.e., decreasing V0_Z_ONSET) might lead to an improvement.
