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

.. _run_domains:

Running the domain calculations
-------------------------------

To execute a domain calculation, set up a base folder (\ :file:`my_domain_calc`
in :numref:`list_domains_inputs`) containing all the normal input files,
except for the ones specifically concerned with the structures,
i.e., **without** :ref:`POSCAR` and :ref:`VIBROCC`. Likewise, the
:ref:`PARAMETERS` file should contain **no** parameters concerned with
interpretation of |POSCAR| or |VIBROCC|, such as, e.g., :ref:`BULK_REPEAT`,
:ref:`ELEMENT_MIX`, :ref:`SITEDEF`, etc.; if any such parameter is present,
|calc| will raise an error.
Finally, in the |PARAMETERS| file, **do define** the :ref:`DOMAIN`
parameter once for each of the domains that should be included. The
:ref:`DOMAIN` parameter can point to an absolute or relative path from
which the input data for a given domain should be fetched. The path may
point to a complete :file:`Tensors_<index>.zip` :ref:`file<tensorszip>`
(e.g., :file:`my_domain_2/Tensors/Tensors_005.zip` in
:numref:`list_domains_inputs`) or to a folder containing the
:ref:`usual input files<list_input_files>` for a reference calculation
(:ref:`EXPBEAMS` and :ref:`DISPLACEMENTS` files in the subfolder are ignored),
as those in :file:`my_domain_1` of :numref:`list_domains_inputs`. If
the target path is a directory in which previous ViPErLEED calculations
have been executed, |calc| will check whether there is a :file:`Tensors`
folder, and fetch the highest-number :file:`Tensors_<index>.zip` file.
For more information, see the :ref:`DOMAIN` page. Use the :ref:`DOMAIN_STEP`
parameter to define the step width for domain-area variations.

.. _list_domains_inputs:
.. code-block:: console
    :caption: Example directory tree with inputs for a domain calculation.

    my_domain_calc/
    ├── my_domain_1/       <-- Use input files, requires reference calculation
    │   ├── PARAMETERS
    │   ├── POSCAR
    │   ├── VIBROCC
    │   └── DISPLACEMENTS  <-- Not used
    ├── my_domain_2/       <-- Use most recent Tensors
    │   ├── Tensors/
    │   │   ├── ...
    │   │   └── Tensors_005.zip
    │   ├── PARAMETERS     <-- Use Tensors_005 instead of this
    │   ├── POSCAR         <-- Use Tensors_005 instead of this
    │   └── VIBROCC        <-- Use Tensors_005 instead of this
    ├── EXPBEAMS.csv
    ├── DISPLACEMENTS
    └── PARAMETERS


You can then run |calc| :ref:`as usual<cli_calc>`.
:numref:`list_domains_outputs` shows the tree structure at the end of the
domain calculation. All structure-specific output files (e.g., :ref:`poscar`,
:ref:`vibrocc`) are collected in the respective domain subfolders (i.e.,
:file:`my_domain_1` and :file:`my_domain_2`). Instead, output concerning all
domains taken together will go to the main :file:`my_domain_calc` folder (e.g.,
the :ref:`searchprogresspdf` file and the |R-factor| :ref:`plots<Rfactorplots>`
after the :ref:`super_pos`).

.. _list_domains_outputs:
.. code-block:: console
    :caption:
        Same directory tree as in :numref:`list_domains_inputs` after execution
        of a domain calculation.

    my_domain_calc/
    ├── my_domain_1/       <-- Use input files, requires reference calculation
    │   ├── OUT/           <-- Created by calc at end
    │   │   ├── POSCAR
    │   │   ├── VIBROCC
    │   │   └── ...
    │   ├── SUPP/          <-- Created by calc at end
    │   │   ├── POSCAR_bulk
    │   │   └── ...
    │   ├── Tensors/       <-- Created by calc at end
    │   │   └── Tensors_001.zip
    │   ├── PARAMETERS
    │   ├── POSCAR
    │   ├── VIBROCC
    │   └── DISPLACEMENTS  <-- Not used
    │
    ├── my_domain_2/       <-- Use most recent Tensors
    │   ├── OUT/           <-- Created by calc at end
    │   │   ├── POSCAR
    │   │   ├── VIBROCC
    │   │   └── ...
    │   ├── SUPP/          <-- Created by calc at end
    │   │   ├── POSCAR_bulk
    │   │   └── ...
    │   ├── Tensors/
    │   │   ├── ...
    │   │   └── Tensors_005.zip
    │   ├── PARAMETERS
    │   ├── POSCAR
    │   └── VIBROCC
    │
    ├── OUT/
    │   ├── Rfactor_analysis_superpos.pdf
    │   ├── Search_progress.pdf
    │   └── ...
    ├── viperleed-calc-<timestamp>.log
    ├── EXPBEAMS.csv
    ├── DISPLACEMENTS
    └── PARAMETERS

To specify which segments should be run, either use the :ref:`RUN` parameter
as usual, or set ``RUN = 4`` as a shorthand for a domain calculation. This
will be interpreted as ``RUN = 2-3`` or ``RUN = 1-3``, depending on whether
the input files are compatible :file:`Tensors.zip` files or whether a reference
calculation is needed, respectively. For ``RUN = 4``, reference calculations
will only be executed for the domains that need them. Specify ``RUN = 1-3``
explicitly to rerun reference calculations for all domains. However, as
discussed above, it is recommended to run the reference calculations separately
beforehand for better control.

.. warning::
    The :ref:`bookkeeper<bookkeeper>` functionality is only partially
    implemented for domain calculations. Only the :file:`history` folder and
    :file:`history.info` file for the root directory (\ :file:`my_domain_calc`
    in :numref:`list_domains_outputs`) are handled automatically.
    The domain-specific subfolders (i.e., :file:`my_domain_1` and
    :file:`my_domain_2` in :numref:`list_domains_outputs`) will not be
    processed. To preserve the domain-specific output files, you must manually
    run the |bookkeeper| in each of the domain subfolders using the command
    ``viperleed bookkeeper --archive``. Then, to clean the directories
    and remove old ``*_ori`` and ``*.log`` files, run the |bookkeeper|
    with the ``--clear`` flag in each of the domain subfolders. If you
    do not run the |bookkeeper| (in ``--archive`` mode) in the subfolders,
    the results of a structure optimization (especially, files |POSCAR| and
    |VIBROCC|) **will be lost**: the next calculation will **start from the**
    **same inputs as the previous one**.

.. versionchanged:: 0.13.0
    In earlier versions of |calc|, the results of the calculations from each
    domain in a domain calculation would **not be copied back from the work**
    **folder**.


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


.. _domain_inputs_outside:

Fetching domain inputs from other filesystem locations
------------------------------------------------------

The :ref:`DOMAIN` parameter also supports the specification of other locations
from which inputs for each domain should be collected.
:numref:`list_domains_inputs_abs` shows an example of such a filesystem: inputs
for the domains defined in :file:`my_domain_calc/PARAMETERS` are to be fetched
from the :file:`domain_1_somewhere_else` and :file:`domain_2_somewhere_else`
folders. This section describes the differences that such a layout entails, as
compared to the simpler one in :numref:`list_domains_inputs`.

.. _list_domains_inputs_abs:
.. code-block:: console
    :caption:
        Example directory tree for a domain calculation with inputs
        for each domain stored at arbitrary locations on the filesystem.
        A better structure for the input tree is displayed in
        :numref:`list_domains_inputs`.

    my_domain_calc/
    ├── EXPBEAMS.csv
    ├── DISPLACEMENTS
    └── PARAMETERS

    domain_1_somewhere_else/
    ├── PARAMETERS
    ├── POSCAR
    ├── VIBROCC
    └── DISPLACEMENTS  <-- Not used

    domain_2_somewhere_else/
    └── Tensors/
        ├── ...
        └── Tensors_005.zip

In this case, |calc| will create new subfolders of :file:`my_domain_calc`,
uniquely named from either the user-given name or a progressive index, as
:file:`Domain_<name-or-index>`. Only the necessary input files for each domain
are copied there from the source directories.
:numref:`list_domains_outputs_abs` shows the structure of the
:file:`my_domain_calc` folder after execution. Notice that file
:file:`DISPLACEMENTS` was not copied from :file:`domain_1_somewhere_else`,
and only the most recent tensor file (\ :file:`Tensors_005.zip`) was
collected from :file:`domain_2_somewhere_else`.

.. _list_domains_outputs_abs:
.. code-block:: console
    :caption:
        Example directory tree at the end of a domain calculation with inputs
        fetched from arbitrary locations on the filesystem.

    my_domain_calc/
    ├── Domain_1/          <-- Created by calc at initialization
    │   ├── OUT/           <-- Created by calc at end
    │   │   ├── POSCAR
    │   │   ├── VIBROCC
    │   │   └── ...
    │   ├── SUPP/          <-- Created by calc at end
    │   │   ├── POSCAR_bulk
    │   │   └── ...
    │   ├── Tensors/       <-- Created by calc at end
    │   │   └── Tensors_001.zip
    │   ├── PARAMETERS     <-- Copied by calc from /domain_1_somewhere_else/
    │   ├── POSCAR         <-- Copied by calc from /domain_1_somewhere_else/
    │   └── VIBROCC        <-- Copied by calc from /domain_1_somewhere_else/
    │
    ├── Domain_another/    <-- Created by calc at initialization
    │   ├── OUT/           <-- Created by calc at end
    │   │   ├── POSCAR
    │   │   ├── VIBROCC
    │   │   └── ...
    │   ├── SUPP/          <-- Created by calc at end
    │   │   ├── POSCAR_bulk
    │   │   └── ...
    │   └── Tensors/       <-- Copied by calc from /domain_2_somewhere_else/
    │       └── Tensors_005.zip
    │
    ├── OUT/
    │   ├── Rfactor_analysis_superpos.pdf
    │   ├── Search_progress.pdf
    │   └── ...
    ├── viperleed-calc-<timestamp>.log
    ├── EXPBEAMS.csv
    ├── DISPLACEMENTS
    └── PARAMETERS


Execution then proceeds as described in :ref:`run_domains`, treating
:file:`my_domain_calc` as an independent folder. No further reference
is made to the original paths where the inputs were collected from.
In particular, all structure-specific output files will go to the
newly created subfolders, *not* to the original paths.

Additionally, the ``DOMAIN`` parameters in the |PARAMETERS| file of
the main directory (i.e., :file:`my_domain_calc/PARAMETERS`, as well
as :file:`my_domain_calc/OUT/PARAMETERS`) are automatically modified
to point to the newly created subfolders of :file:`my_domain_calc`
(i.e., :file:`Domain_1` and :file:`Domain_another` in
:numref:`list_domains_outputs_abs`).

.. versionchanged:: 0.13.0
    The main |PARAMETERS| file is updated such that each ``DOMAIN`` points
    to the subfolders of the directory in which |calc| was executed.
