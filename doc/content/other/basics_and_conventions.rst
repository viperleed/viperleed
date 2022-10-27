.. _basics:

Basics and conventions
----------------------

Installation
^^^^^^^^^^^^

**TODO** Needs an own docu page. Python 3.8. Info on Fortran compiler currently on FORTRAN_COMP page should be moved there. Does the Windows version require Windows-Subsystem für Linux? !!

Directory structure and file names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each calculation must have its own directory. 
The input and output files have fixed names, see :ref:`the list of files<_list_input_files>`. 
In case of automated multiple search runs (which can be specified in the :ref:`DISPLACEMENTS<DISPLACEMENTS>`  file), 
tleedm creates a “workhistory” directory and moves a snapshot of all input and output files that may be relevant and may get overwritten into a subfolder there. **TODO** Give an example tree for a LEED project!!

Running the program
^^^^^^^^^^^^^^^^^^^

A large number of files are created in the directory that tleedm is executed in. Exemplary job scripts are provided. These generally create a "work" directory, copy input files there, execute tleedm, and then copy the relevant output files back to the data directory. For this purpose, tleedm also creates a "manifest" file that lists the relevant output files which should be copied back.

Coordinate system
^^^^^^^^^^^^^^^^^

The lattice vectors describing the unit cell are named **a**, **b**, and **c**. **a** and **b** must be in the xy plane, and **c** have a positive z component. The +z direction is outwards from the surface.

Units
^^^^^

All distances and vibration amplitudes are in Ångström, energies in eV.!!
