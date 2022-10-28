.. _list_input_files:

List of input files
-------------------

ViPErLEED has a number of input and control files, some of which are optional, 
depending on the work-segment (:ref:`see above<work-segments>`).

+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| File                                 | Function                                          | Required                                                          |
+======================================+===================================================+===================================================================+
| :ref:`PARAMETERS<PARAMETERS>`        | What to execute, and how to interpret other input | Always required                                                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`POSCAR<POSCAR>`                | Structural information                            | Always required                                                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`VIBROCC<VIBOCCIN>`             | Vibrational amplitudes and site occupations       | Always required, can be generated automatically                   |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`IVBEAMS<IVBEAMS>`              | Which beams to calculate                          | Always required, can be generated automatically from EXPBEAMS.csv |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`DISPLACEMENTS<DISPLACEMENTS>`  | What to vary during the search                    | Required for delta calculations and search                        |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`EXPBEAMS.csv<EXPBEAMS>`        | Experimental I(V) curves                          | Required for R-factor calculations and search                     |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+
| :ref:`PHASESHIFTS<PHASESHIFTS>`      | Phase shifts per energy \* site \* L              | Generated automatically if needed                                 |
+--------------------------------------+---------------------------------------------------+-------------------------------------------------------------------+

.. toctree::
   :glob:

   files/input/*
