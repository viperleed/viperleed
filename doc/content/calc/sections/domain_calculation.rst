.. include:: /substitutions.rst

.. _domain_calculation:

===================
Domain calculations
===================

If multiple structures coexist on the sample, then the experimental beams are
an incoherent superposition of the beam sets diffracted from each of the
different structural domains. This means that the different structures must be
optimized at the same time to optimize the |R factor| between calculated and
experimental |IV| curves. This page covers how the structure of the input files
for ViPErLEED changes when domains are present.

.. tip::
    While domain calculations can be run from scratch (i.e., including the
    :ref:`reference calculations<ref-calc>`) from two sets of input files,
    a way to get more user control is to first execute reference calculations
    separately for the different structures, then start the domain calculation
    from the :file:`Tensors_<index>.zip` :ref:`files<tensorszip>` produced
    by those reference calculations. This way, one can verify the correct
    execution of the :ref:`initialization` for each domain, that the detected
    symmetries are sensible, that the theoretical beam sets for the different
    domains are as expected, etc., before launching into the simultaneous
    optimization of the structures.

.. note::
    :term:`TensErLEED` requires the different structures to be calculated
    on the same surface supercell, i.e. with the same unit cell in the
    :ref:`POSCAR` file, and the same :ref:`SUPERLATTICE` parameter for
    each domain. Since changing the unit-cell size may change the detected
    symmetry, the :ref:`SYMMETRY_CELL_TRANSFORM` parameter allows defining
    a base unit cell and linking redundant atoms by translational symmetry.

In addition to the surface unit cell and the :ref:`SUPERLATTICE` parameter, the
following settings must also be the same in each of the reference calculations
(it is recommended to fix the :ref:`PARAMETERS` explicitly):

-  :ref:`THEO_ENERGIES`
-  :ref:`LMAX`
-  :ref:`BEAMINCIDENCE`
-  :ref:`IVBEAMS`

Running the domain calculations
-------------------------------

To execute a domain calculation, set up a base folder (:file:`my_domain_calc`
in :numref:`list_domains_directories`) containing all the normal input files,
except for the ones specifically concerned with the structures,
i.e., **without** :ref:`POSCAR` and :ref:`VIBROCC`. Likewise, the
:ref:`PARAMETERS` file should contain **no** parameters concerned with
interpretation of :file:`POSCAR` or :file:`VIBROCC`, such as, e.g.,
:ref:`BULK_REPEAT`, :ref:`ELEMENT_MIX`, :ref:`SITEDEF`, etc.; if any such
parameter is present, it will be ignored. Finally, in the :file:`PARAMETERS`
file, **do** define the :ref:`DOMAIN` parameter once for each of the domains
that should be included. The :ref:`DOMAIN` parameter can point to an absolute
or relative path from which the input data for a given domain should be
fetched. The path may point to a complete :file:`Tensors_<index>.zip`
:ref:`file<tensorszip>` (e.g., :file:`my_domain_2/Tensors/Tensors_005.zip` in
:numref:`list_domains_directories`) or to a folder containing the
:ref:`usual input files<list_input_files>` for a reference calculation
(:ref:`EXPBEAMS` and :ref:`DISPLACEMENTS` files in the subfolder are ignored),
as those in :file:`my_domain_1` of :numref:`list_domains_directories`. If the
target path is a directory in which previous ViPErLEED calculations
have been executed, |calc| will check whether there is a :file:`Tensors`
folder, and fetch the highest-number :file:`Tensors_<index>.zip` file.
For more information, see the :ref:`DOMAIN` page. Use the :ref:`DOMAIN_STEP`
parameter to define the step width for domain-area variations.

.. _list_domains_directories:
.. code-block:: console
    :caption: Example directory tree for a domain calculation.

    my_domain_calc/
    ├── Domain_1/         <-- Created by calc at initialization
    │   ├── OUT/          <-- Created by calc at end
    │   │   ├── POSCAR_OUT
    │   │   ├── VIBROCC_OUT
    │   │   └── ...
    │   ├── SUPP/          <-- Created by calc at end
    │   │   ├── POSCAR_bulk
    │   │   └── ...
    │   ├── Tensors/       <-- Copied by calc from /my_domain_1/
    │   │   └── Tensors_027.zip
    │   ├── PARAMETERS     <-- Copied by calc from /my_domain_1/
    │   ├── POSCAR         <-- Copied by calc from /my_domain_1/
    │   └── VIBROCC        <-- Copied by calc from /my_domain_1/
    │
    ├── Domain_another/    <-- Created by calc at initialization
    │   ├── OUT/           <-- Created by calc at end
    │   │   ├── POSCAR_OUT
    │   │   ├── VIBROCC_OUT
    │   │   └── ...
    │   ├── SUPP/          <-- Created by calc at end
    │   │   ├── POSCAR_bulk
    │   │   └── ...
    │   └── Tensors/       <-- Copied by calc from /my_domain_2/
    │       └── Tensors_005.zip
    │
    ├── OUT/
    │   ├── Rfactor_analysis_superpos.pdf
    │   ├── Search_progress.pdf
    │   └── ...
    ├── calc-<timestamp>.log
    ├── EXPBEAMS.csv
    ├── DISPLACEMENTS
    └── PARAMETERS


    my_domain_1/
    ├── Tensors/
    │   ├── Tensors_001.zip
    │   ├── ...
    │   └── Tensors_027.zip
    ├── PARAMETERS
    ├── POSCAR
    └── VIBROCC


    my_domain_2/
    └── Tensors/
        └── Tensors_005.zip


You can then run |calc| :ref:`as usual<cli_calc>`. In the base directory,
a subfolder is created for each domain (folders :file:`Domain_1` and
:file:`Domain_another` in :numref:`list_domains_directories`). Input
files for each domain are copied there. All structure-specific output
files (e.g., :ref:`POSCAR_OUT<poscar>`, :ref:`VIBROCC_OUT<vibrocc>`)
will go to these subfolders, *not* to the original paths from which
the inputs were fetched. Output concerning all domains taken together
will go to the main folder (e.g., the :ref:`searchprogresspdf` file
and the |R-factor| :ref:`plots<Rfactorplots>` after the :ref:`super_pos`).

To specify which segments should be run, either use the :ref:`RUN` parameter
as usual, or set ``RUN = 4`` as a shorthand for a domain calculation. This
will be interpreted as ``RUN = 1-3`` or ``RUN = 2-3``, depending on whether
the input files are compatible :file:`Tensors.zip` files or whether a reference
calculation is needed, respectively. For ``RUN = 4``, reference calculations
will only be executed for the domains that need them. Specify ``RUN = 1-3``
explicitly to rerun reference calculations for all domains. However, as
discussed above, it is recommended to run the reference calculations separately
beforehand for better control, and specify ``RUN = 2-3`` explicitly here.

.. warning::

  The :ref:`bookkeeper<bookkeeper>` functionality is only partially implemented
  for domain calculations.
  The bookkeeper will archive and clean up the top level directory as usual, but
  the domain-specific directories will not be cleaned up.
  To preserve the domain-specific output files, you must manually run the 
  bookkeeper in each of the domain directories using the command
  ``viperleed bookkeeper --archive``.
  To clean the directories and remove old `_ori` and `.log` files, run the 
  bookkeeper with the ``-clear`` flag.

The DISPLACEMENTS file for domains
----------------------------------

Instead of specifying :ref:`DISPLACEMENTS` in each of the input subfolders,
DISPLACEMENTS are defined in the main folder for all domains.
The syntax is similar to the way consecutive searches are specified, with an
extra header line specifying which domain is being addressed, e.g., for a
calculation with two domains called ``1x1`` and ``2x1``:

::

   == SEARCH z

     == DOMAIN 1x1

       = GEO_DELTA
       * L(1) z = -0.1 0.1 0.02
       * L(2) z = -0.05 0.05 0.025

     == DOMAIN 2x1

       = GEO_DELTA
       * L(1-2) z = -0.1 0.1 0.02
       * L(3) z = -0.05 0.05 0.025

.. note:: Indentation is allowed, but does not affect the functionality.