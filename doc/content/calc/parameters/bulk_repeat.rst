.. include:: /substitutions.rst

.. _bulk_repeat:

BULK_REPEAT
===========

BULK_REPEAT defines how the bulk layers (as defined by :ref:`LAYER_CUTS`
and :ref:`N_BULK_LAYERS`) are repeated after reaching the bottom of the
:ref:`POSCAR` file.
If the |c| vector of the POSCAR is oriented parallel to a bulk repeat
vector, BULK_REPEAT can be defined simply as a bulk "layer thickness"
(which may be equivalent to step height).
Note that the :ref:`BULK_LIKE_BELOW` parameter offers an
easy way to detect the bulk repeat unit automatically.
See also :numref:`fig_Fe2O3_layers` for an illustration of the bulk repeat
vector and the layers.

.. note::

   If :ref:`BULK_LIKE_BELOW` is defined, BULK_REPEAT and :ref:`N_BULK_LAYERS`
   can be automatically detected from the :ref:`POSCAR` file.
   The :ref:`PARAMETERS` file will automatically be updated accordingly.

**Default:** Detected automatically from :ref:`POSCAR` and
:ref:`N_BULK_LAYERS`: If as many non-bulk layers as bulk
layers are "bulk-like" (i.e. unrelaxed), automatically detects the minimal
repeat vector. Otherwise, assumes that the bulk repeat vector is parallel
to |c| and defines BULK_REPEAT as a bulk layer thickness. This is done by
detecting the cartesian distance in |z| direction between bottom atom in
bottom bulk layer and bottom atom in bottom non-bulk layer.

**Allowed values:** positive float, or any three float values if vector

**Syntax:**

::

   ! options for defining bulk repeat vectors along c:
   BULK_REPEAT = 1.428         ! bulk layer thickness along z in cartesian coordinates
   BULK_REPEAT = z(1.428)      ! same as above
   BULK_REPEAT = c(0.1173)     ! bulk layer thickness as a fraction of the unit cell vector perpendicular to the surface

   ! explicitly defined bulk repeat vectors:
   BULK_REPEAT = [-2.5175   2.3358  3.6820]     ! repeat vector in cartesian coordinates
   BULK_REPEAT = xyz[-2.5175   2.3358  3.6820]  ! same as above
   BULK_REPEAT = abc[0.5 0.5 0.1173]            ! fractional coordinates instead of cartesian

Depending on your initial POSCAR file, three cases can be distinguished:

-  If enough layers above the bulk are sufficiently bulk-like, BULK_REPEAT can
   safely be left undefined. It will be detected automatically by finding the
   shortest repeat vector which translates the bottom-most non-bulk layers such
   that they match the bulk layers.
-  If there are too few unrelaxed layers above the bulk, but the |c|
   unit-cell vector of the :ref:`POSCAR` file is oriented as a
   bulk repeat vector, BULK_REPEAT can also be left undefined, as long as
   at least the bottom-most non-bulk atom is still sufficiently "bulk-like".
   In this case, BULK_REPEAT will be calculated simply by assuming that the
   bulk layers should repeat along |c| with the same spacing as in the
   :ref:`POSCAR`.
-  The only case where BULK_REPEAT *must* be defined manually is when the |c|
   unit-cell vector is not parallel to a bulk repeat vector, and there are too
   few bulk-like layers to detect a repeat vector automatically.

When defining BULK_REPEAT manually as a single float value, a definition in
*cartesian coordinates* will be interpreted as a bulk thickness, i.e., along z.
The bulk layers will be repeated along c such that their total "thickness"
(including spacing) equals the defined value. If the value is given as
``c(value)``, i.e. in fractional coordinates, the bulk layers will simply
be repeated along c, with the value again defining the total "thickness"
of a bulk repeat unit.

When BULK_REPEAT is defined in vector form, the vector does not have to be
parallel to c. The given repeat vector defines how each individual atom should
be displaced to get the repeat unit.
The vector can point either up or down, but note that it is defined in
coordinates as in the :ref:`POSCAR` file, that is with |c|
*pointing out of the surface*.

**Note:** If BULK_REPEAT is not defined, the calculated value will be written
to the :ref:`PARAMETERS` file during the initialization to ensure
that in future runs, the bulk repeat value is conserved even if the bottom-most
non-bulk layers are varied.
