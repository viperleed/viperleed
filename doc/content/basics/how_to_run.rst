.. _how_to_run:

==================================
Directory Structure and How To Run
==================================

In ViPErLEED, each calculation must have its own directory. 
The input and output files have fixed names, see :ref:`the list of files<list_input_files>`. 
Below, we give an example directory tree with the files needed to start a LEED :math:`I(V)` calculation in ViPErLEED.
See also the ViPErLEED examples section for some of working calculations to help you get started.

Minimum input
-------------

To set up a ViPErLEED calculation first create a source directory (in this example ``my_surface``) and place all input files inside.

.. code-block:: console
    :caption: Minimum input directory tree

    my_surface
    ├── EXPBEAMS.csv
    ├── POSCAR
    ├── PARAMETERS
    └── job.py

The minimum information required to start a :ref:`refercence calculation<ref-calc>` is contained in three files:
-   :ref:`EXPBEAMS.csv file<expbeams>` contains the experimentally measured LEED :math:`I(V)` curves.
    Unless specified otherwise in :ref:`PARAMETERS<parameters>`, ViPErLEED will also use the information in the :ref:`EXPBEAMS.csv file<expbeams>` to set energy rangesa and choose which beams to calculate.
-   :ref:`POSCAR file<poscar>` contains the reference surface structure.
    ViPErLEED will determine the applicable symmetry from the :ref:`POSCAR file<poscar>`. See also the :ref:`symmetry settings<symmetry_settings>` and :ref:`input structure settings<input_structure_settings>`.
-   :ref:`PARAMETERS<parameters>` contains the settings for the calculation (see :ref:`the list of parameters<paramname>`).
    If no :ref:`VIBROCC file<vibrocc>` is given, :ref:`PARAMETERS<parameters>` needs to contain values for :ref:`T_EXPERIMENT<t_experiment>` and :ref:`T_DEBYE<t_debye>`.
-   The ``job.py`` script. **TODO**

To run also a :ref:`delta amplitudes calculation<sec_deltas>` and a :ref:`structure search<sec_search>`, you additionally need to provide a :ref:`DISPLACEMENTS<displacements>` file that contains the requested pertubations of the structure.

.. note:: 
    Setting up a :ref:`domain calculation<domains>` with multiple surface structures, requires a slightly different directory tree.
    See the :ref:`domain calculation page<domains>`.

work directory and output files
-------------------------------

A large number of files are created in the directory that tleedm is executed in.
Exemplary job scripts are provided.
These generally create a "work" directory, copy input files there, execute tleedm, and then copy the relevant output files back to the data directory.
For this purpose, tleedm also creates a :ref:`manifest` file that lists the relevant output files which should be copied back.




After the first run, an ``OUT`` directory is created that contains the output files.