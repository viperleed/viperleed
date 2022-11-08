.. _how_to_run:

==================================
Directory Structure and How To Run
==================================

A large number of files are created in the directory that tleedm is executed in.
Exemplary job scripts are provided.
These generally create a "work" directory, copy input files there, execute tleedm, and then copy the relevant output files back to the data directory.
For this purpose, tleedm also creates a :ref:`manifest` file that lists the relevant output files which should be copied back.

.. code-block:: console
    :caption: Minimum input directory tree

    my_surface
    ├── EXPBEAMS.csv
    ├── PARAMETERS
    ├── job.py
    └── POSCAR


After the first run, an ``OUT`` directory is created that contains the output files.