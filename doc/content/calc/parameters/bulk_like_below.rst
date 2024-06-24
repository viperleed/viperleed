.. _bulk_like_below:

BULK_LIKE_BELOW
===============

BULK_LIKE_BELOW is a helper parameter that can be used to automatically
determine :ref:`BULK_REPEAT`, :ref:`N_BULK_LAYERS` and :ref:`LAYER_CUTS`.
It specifies a fraction of **c** in the :ref:`POSCAR` unit cell below
which the structure is bulklike, i.e., unrelaxed (if coming from DFT), or
sufficiently close to unrelaxed to fall within :ref:`sym_eps`.
Generally, this is the most easy way to detect and define the bulk if at least
two unrelaxed bulk repeat units are present in the :ref:`POSCAR` file.
BULK_LIKE_BELOW does not need to be specified if the bulk is manually defined
through :ref:`N_BULK_LAYERS`, :ref:`LAYER_CUTS`, and :ref:`BULK_REPEAT`.

**Default:** No default, no auto-detection of bulk if undefined

**Allowed values:** positive float between 0 and 1

**Syntax:**

::

   BULK_LIKE_BELOW = 0.22

BULK_LIKE_BELOW will only be used once during initialization, then
automatically commented out in PARAMETERS for future runs. Instead,
the automatically detected :ref:`N_BULK_LAYERS`, :ref:`LAYER_CUTS`,
and :ref:`BULK_REPEAT` will be written into the PARAMETERS file
explicitly. This ensures stability of the bulk definition in case
some of the layers that are initially bulk-like are then varied
during optimization.

The bulk unit cell is detected by taking only the slab below BULK_LIKE_BELOW
and finding the smallest translation vector (towards vacuum) with a **c**
component that preserves the structure. Only translation vectors with **z**
components of less than or equal to the thickness of the slab below
BULK_LIKE_BELOW are considered. The minimal such translation vector
is by definition the minimal :ref:`BULK_REPEAT` vector.
By default, :ref:`N_BULK_LAYERS` is 1, but if a further cut with layer
spacing of at least 1.2 Å is possible within the bulk repeat unit,
then :ref:`N_BULK_LAYERS` will be chosen as 2 and :ref:`LAYER_CUTS`
will be set accordingly.

Note that for the automatic detection to work, the bulk structure has
to repeat at least once below BULK_LIKE_BELOW. If there are not enough
unrelaxed layers, you can try to choose BULK_LIKE_BELOW more liberally
and increase :ref:`sym_eps` to still detect the bulk repeat unit
automatically. However, be aware that the numerical precision of the
detected :ref:`BULK_REPEAT` vector is then limited by the relaxations
of your structure, and may need to be corrected manually.

If :ref:`BULK_REPEAT` is already defined, BULK_LIKE_BELOW will
be *ignored* (with a warning).

.. todo::
    Add examples of one correct and one incorrect definition with figures.
