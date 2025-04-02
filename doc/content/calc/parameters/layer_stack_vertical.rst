.. include:: /substitutions.rst

.. _layer_stack_vertical:

LAYER_STACK_VERTICAL
====================

LAYER_STACK_VERTICAL defines whether the interlayer vectors of non-bulk layers
in the :ref:`AUXGEO` file (input to the underlying TensErLEED program) are
parallel to z (if True) or parallel to the |c| vector of the input unit
cell (if False). In this file, the atom coordinates per layer will change
accordingly.

This parameter is for debugging only, it does not affect the result of the
calculation.

**Default**: LAYER_STACK_VERTICAL = True (layer stacking along z, all single
atom coordinates given relative to lowest non-bulk layer)

**Syntax**:

::

   LAYER_STACK_VERTICAL = False

**Acceptable values** (not case sensitive): ``True``, ``False``, ``T``, ``F``,
``Z`` *(=True)*, ``C`` *(=False)*

``LAYER_STACK_VERTICAL = True`` makes it easier to recognize
the overall arrangement of atoms from the atom coordinates.
However, ``LAYER_STACK_VERTICAL = False`` can make mostly
equivalent layers more easily comparable in the AUXGEO file.
