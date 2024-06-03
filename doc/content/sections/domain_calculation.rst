.. _domain_calculation:

.. include:: /substitutions.rst

===================
Domain calculations
===================

If multiple structures coexist on the sample, then the experimental beams 
are a superposition from the beam sets diffracted from each of the different 
structural domains. This means that the different structures must be optimized 
at the same time to optimize the R-factor between theoretical and experimental 
I(V)-curves. This page covers how the input file structure for ViPErLEED 
changes when domains are present.

While domain calculations can be run from scratch (i.e., including the 
reference calculations) from two sets of input files, a way to get more 
user control is to first execute reference calculations separately for 
the different structures, then start the domain calculation from the 
Tensors.zip files from those pre-calculations. This way, you can convince 
yourself that initializations have run correctly, that the detected 
symmetries are sensible, that the theoretical beam sets for the different 
domains are as expected, etc., before launching into the simultaneous 
optimization of the structures.

Note that TensErLEED requires the different structures to be calculated 
on the same surface supercell, i.e. with the same unit cell in the 
:ref:`POSCAR file<POSCAR>`, and the same :ref:`SUPERLATTICE<SUPERLATTICE>` 
parameter for each domain. Since changing the unit cell size may change the 
detected symmetry, the :ref:`SYMMETRY_CELL_TRANSFORM<SYMMETRY_CELL_TRANSFORM>`  
parameter allows defining a base unit cell and linking redundant atoms by 
translational symmetry.

In addition to the surface unit cell and 
:ref:`the SUPERLATTICE parameter<SUPERLATTICE>`, the following settings 
must also be the same in each of the reference calculations (it is recommended 
to fix the :ref:`PARAMETERS<PARAMETERS>` explicitly):

-  :ref:`THEO_ENERGIES<theo_energies>` 
-  :ref:`LMAX<LMAX>` 
-  :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>` 
-  :ref:`IVBEAMS<IVBEAMS>` 

Running the domain calculations
-------------------------------

To execute a domain calculation, set up a base folder containing all the normal 
input files, except for the ones specifically concerned with the structures, 
i.e., **without** POSCAR and VIBROCC. Likewise, the PARAMETERS file should 
contain **no** parameters concerned with interpretation of POSCAR or VIBROCC, 
such as, e.g., :ref:`BULK_REPEAT<BULK_REPEAT>`, :ref:`ELEMENT_MIX<ELSPLIT>`, 
:ref:`SITE_DEF<SITEDEF>`, etc.; if any such parameters are present, they will 
be ignored here. Finally, **do** define the :ref:`DOMAIN<DOMAIN>` parameter in 
the PARAMETERS file once for each domain that should be included. The DOMAIN 
parameter can point to an absolute or relative path from which the input data 
for a given domain should be fetched; input data can be a complete Tensors.zip 
file, or the path can point to a folder containing the usual input files for a 
reference calculation. If the target path is a directory in which previous 
ViPErLEED calculations have been executed, the program will check whether there 
is a Tensors folder and fetch the highest-number Tensors.zip file.
For more information, see the :ref:`DOMAIN<DOMAIN>`  page. Use the 
:ref:`DOMAIN_STEP<DOMAIN_STEP>`  parameter to define the step width 
for domain area variations.

Then, run |calc| :ref:`as usual<cli_calc>`. A subfolder will be created for 
each domain, in which the input files will be placed. All structure-specific 
output files (e.g. POSCAR_OUT, VIBROCC_OUT) will go to these subdirectories, 
*not* the original paths from which the inputs were fetched; output concerning 
all domains taken together will go to the main folder (e.g., the 
:ref:`Search-progress.pdf<searchprogresspdf>` file and 
:ref:`Rfactor plots<Rfactorplots>` after the Superpos calculation).

To specify which segments should be run, either use the :ref:`RUN<RUN>`  
parameter as usual, or set ``RUN = 4`` as a shorthand for a domain 
calculation. This will be interpreted as ``RUN = 1-3`` or ``RUN = 2-3``, 
depending on whether the input files are compatible Tensors.zip files or 
whether a reference calculation is needed. For ``RUN = 4``, reference 
calculations will only be executed for the domains that need them; 
specify ``RUN = 1-3`` explicitly to re-run reference calculations 
for all domains. However, as discussed above, it is recommended you 
run the reference calculations separately beforehand for better control, 
and specify ``RUN = 2-3`` explicitly here.

.. warning:: 
  In the current version, automatic clean-up after domain calculations is 
  implemented only in a rudimentary manner. Domain-specific output files 
  are not copied out from the work folder, and the bookkeeper ``--cont`` 
  functionality will not work.

The DISPLACEMENTS file for domains
----------------------------------

Instead of specifying :ref:`DISPLACEMENTS<DISPLACEMENTS>` in each of the input 
subdirectories, DISPLACEMENTS are defined in the main folder for all domains. 
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

.. note:: Indentation is allowed, but does not affect function.