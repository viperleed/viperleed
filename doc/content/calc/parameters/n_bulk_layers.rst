.. _n_bulk_layers:

N_BULK_LAYERS
=============

N_BULK_LAYERS defines how many layers, counting from the bottom, form the
repeat unit of the bulk structure. Layers are defined by the :ref:`LAYER_CUTS`
parameter.

.. note::

   The :ref:`BULK_LIKE_BELOW` parameter offers an easy way to detect the
   bulk repeat unit automatically, which will set both N_BULK_LAYERS and
   :ref:`LAYER_CUTS`. If :ref:`BULK_LIKE_BELOW` is defined and the bulk is
   detected, N_BULK_LAYERS will be updated in the :ref:`PARAMETERS` file
   during initialization.

**Default:** 1

**Allowed values:** 1, 2

**Syntax:**

::

   N_BULK_LAYERS = 2

If N_BULK_LAYERS = 2, the bottom two layers will be repeated, alternating
between the two.
If N_BULK_LAYERS = 1, only the bottom layer will be repeated.
See also :numref:`fig_Fe2O3_layers` and :ref:`this example<example_Fe2O3>` for
an illustration of the layers and the bulk repeat unit.

The layers defined as bulk should form a bulk unit cell, or a larger bulk
repeat unit.
The restriction to only one or two layers is given by TensErLEED.
If your smallest bulk unit cell consists of more than two layers, use
the :ref:`LAYER_CUTS` parameter to manually combine those layers to
only one or two layers, and set N_BULK_LAYERS accordingly. An example for such
a case is the ABAC stacking sequence of the Ni\ :sub:`3`\ Ti (D0\ :sub:`24`)
structure. There, one should define the two bottommost layer pairs (AB and AC)
as single layers by cutting at the appropriate heights (using a list and
restricting cutting by interlayer distance to higher layers in LAYER_CUTS),
and ``N_BULK_LAYERS = 2``.

By default, as many layers *above* the bulk as are present *within* the bulk 
are used to determine the bulk repeat vector, see :ref:`BULK_REPEAT`. If not 
enough bulk-like layers are found, :ref:`BULK_REPEAT` will instead be defined 
as parallel to the POSCAR **c** vector, using only the z position of the 
bottommost non-bulk atom. See the page on the :ref:`BULK_REPEAT` parameter 
for more complicated cases.

Correct identification of the bulk layers and bulk repeat vectors can be
checked by looking at the :ref:`POSCAR_bulk and POSCAR_bulk_appended<POSCAR>`
files, created during initialization.

**See also**: :ref:`BULK_LIKE_BELOW`, another way to specify which layers of 
the input file are used as bulk.
