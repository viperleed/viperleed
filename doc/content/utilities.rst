.. _utilities:

Utilities
---------

ViPErLEED includes several small utilities for file processing and organization. Proper embedding of these into a shared interface is pending; the following list is a loose summary and may not be complete.

-  :ref:`Bookkeeper<bookkeeper>`: Helper utility for the TensErLEED Manager tleedm. Sorts files from previous runs into a "history" folder, and keeps track in the history.info file.
-  Combine-POSCAR, POSCAR_get_bulk_repeat: Very simple POSCAR editing scripts. **TODO** documentation
-  ModifyPhaseshifts: Simple utility for taking an existing PHASESHIFTS file and duplicating and/or re-arranging blocks. **TODO** usage (how to call it)
-  AUXEXPBEAMS_TO_EXPBEAMS: Transforms :ref:`AUXEXPBEAMS<AUXEXPBEAMS>`  format files (input for TensErLEED) to :ref:`EXPBEAMS.csv<EXPBEAMS>`  format. **TODO** usage (how to call it)

.. toctree:: 
   :maxdepth: 1
   :hidden:

   Atomic Simulation Environment API<utilities/aseapi>
   Bookkeeper<utilities/bookkeeper>