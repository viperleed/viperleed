.. _utilities:

Utilities
---------

ViPErLEED includes several small utilities for file processing and organization. 
Proper embedding of these into a shared interface is pending; the following list is a loose summary and may not be complete.

-  :ref:`Bookkeeper<bookkeeper>`: Helper utility for the TensErLEED Manager tleedm. Sorts files from previous runs into a "history" folder, and keeps track in the history.info file.
-  :ref:`Combine-POSCAR, POSCAR_get_bulk_repeat<poscar_combine_repeat>`: Very simple POSCAR editing scripts.
-  :ref:`ModifyPhaseshifts<modify_phaseshifts>`: Simple utility for taking an existing PHASESHIFTS file and duplicating and/or re-arranging blocks.
-  :ref:`AUXEXPBEAMS_TO_EXPBEAMS<aux_to_exp>`: Transforms :ref:`AUXEXPBEAMS<AUXEXPBEAMS>`  format files (input for TensErLEED) to :ref:`EXPBEAMS.csv<EXPBEAMS>`  format.
   **TODO Alex, Florian** usage (how to call it)

.. toctree:: 
    utilities/bookkeeper
    utilities/poscar_combine_repeat
    utilities/modify_phaseshifts
    utilities/aux_to_exp
