.. _tensor_output:

TENSOR_OUTPUT
=============

TENSOR_OUTPUT defines whether Tensor output is required for non-bulk 
layers (starting from the top) in the refcalc. Layers are defined by 
the :ref:`LAYER_CUTS`  parameter.

**Default:** 1 for all layers

**Allowed values:** 0, 1, T, F, True, False

**Syntax:** 1 = TRUE, 0 = FALSE

::

   TENSOR_OUTPUT = 1 1 0 0 1

OR (same result):

::

   TENSOR_OUTPUT = 2*1 2*0 1

OR:

::

   TENSOR_OUTPUT = False   ! disables tensor output for all layers

Give a space-separated list of values (0 or 1), or specify to repeat values
several times using ``times*value``. Layers are numbered from top to bottom,
so the first value given in TENSOR_OUTPUT will be assigned to the top layer,
the second value to the first subsurface layer, and so on.

If only one value is specified, this value will be applied to all
layers; so, ``TENSOR_OUTPUT = True`` explicitly turns tensor output
on, and ``TENSOR_OUTPUT = False`` disables output for all layers.

**Notes**:

-  Only those layers for which tensors are output will be accessible for the
   structural optimization!
-  Setting this parameter to 0 for some layers will speed up the calculation
-  This parameter can generally be left undefined, unless: (1) you want to run
   the very final reference calculation, and have no intention of optimizing
   any longer; (2) you want to create a thick 'pseudobulk' whose geometry will
   not be optimized
-  Notice that defining this parameter requires you to know already how many
   layers you are going to split your slab into. Therefore, it's a good idea
   to check layer numbering in the :ref:`POSCAR` file after
   initialization to make sure you are addressing the correct layers.
