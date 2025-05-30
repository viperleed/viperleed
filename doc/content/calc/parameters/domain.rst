.. include:: /substitutions.rst

.. _domain:

DOMAIN
======

The DOMAIN parameter defines a file source for optimization using
domains, i.e. when multiple structures are optimized at the same time.
See also the page about :ref:`Domain calculations<domain_calculation>`.

**Default:** None

Syntax
------

::

   DOMAIN 1x1 = ~/LEED/Fe2O3_1x1_supercell2x1       ! absolute path to files for domain 1, named 1x1
   DOMAIN 2x1 = ~/LEED/Fe2O3_2x1                    ! absolute path to files for domain 2, named 2x1
   ! alternative:
   DOMAIN 1x1 = Tensors/Tensors_003.zip    ! relative path with Tensors_003.zip as a file
   ! or:
   DOMAIN 1x1 = Tensors/Tensors_003        ! relative path; will first check if the directory exists, but otherwise also accept Tensors_003.zip as a file

A unique name should be defined for each domain on the left-hand side. If no
names are defined, the domains will be numbered instead. The right-hand side
accepts paths to either directories or to Tensors.zip files.
Paths can be absolute or relative. When given as a relative path, the directory
in which |calc| was started is considered first. Then, the temporary
directory in which calculations are executed (see also :ref:`how_to_run` and
the ``--work`` :ref:`command-line argument<cli_calc>`) is searched. 
Input files from the different domains will then be compared to determine if a
new reference calculation is needed. If the source is a Tensors.zip file, the
input files saved in the .zip archive will be compared.
If they are compatible (i.e. same number and order of beams, same unit cell,
same energy range etc.), then no further reference calculation is executed,
and the Tensor files are used as-is. If the right-hand side is a directory
instead, the program will first check for existing Tensors files in that
directory, and check them as described above. If no Tensors files are found,
or if the Tensors files are incompatible, the target path will be searched
for the required input files, and a new reference calculation will be executed.

The target directories or Tensor archives
must contain the following input files:

-   :ref:`PARAMETERS` (most parameters from the subfolders will be ignored,
    but those relevant to interpreting the respective POSCAR files are
    required)
-   :ref:`POSCAR`
-   :ref:`VIBROCC`

In Tensor archives, an :ref:`IVBEAMS` file is also required in order to check
compatibility with other sources. If a :ref:`PHASESHIFTS` file is present and
compatible with the input, it will be used; otherwise, a new PHASESHIFTS file
will be generated automatically.

.. versionchanged:: 0.13.0
    Relative ``DOMAIN`` paths also consider the folder in which |calc| was
    started. In earlier versions, only the work directory was searched.
