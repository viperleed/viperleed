.. include:: /substitutions.rst

.. _utilities:

Utilities
---------

ViPErLEED includes several additional utilities for file processing and
organization. These are installed as part of the Python package, and can
be called from the command line via ``viperleed poscar`` or ``viperleed util``
followed by the utility name.

-  :ref:`POSCAR utilities<poscar_utils>`: A set of utilities for quick POSCAR
   manipulation. Examples include deleting parts of cell, structure symmetry
   detection or setting flags for use with :term:`DFT` software. These are
   called via ``viperleed poscar`` followed by the utility name.
-  :ref:`rearrange_phaseshifts<rearrange_phaseshifts>`: A simple utility for
   taking an existing :ref:`PHASESHIFTS` file and duplicating
   or re-arranging blocks.
-  :ref:`aux_to_exp`: Transforms :ref:`AUXEXPBEAMS` format files (input for
   TensErLEED) to :ref:`EXPBEAMS.csv<EXPBEAMS>`  format.

.. tip::
    If you find yourself using a utility frequently, you can add an alias to
    your ``.bashrc`` file.


.. seealso::
   :ref:`bookkeeper`: A helper utility built into |calc| that stores files
   from previous runs into a :file:`history` folder, and keeps track of
   executions in the :file:`history.info` file.

.. toctree::

   utilities/poscar_utils
   utilities/rearrange_phaseshifts
   utilities/aux_to_exp
