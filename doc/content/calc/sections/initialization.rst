.. include:: /substitutions.rst

.. _initialization:

==============
Initialization
==============


The initialization is the first section of |calc| that is run. It is
**always** inserted at the beginning of any |calc| calculation, even if
not explicitly specified via the :ref:`RUN` parameter. To test out
symmetry recognition, generation of phase shifts, detection of bulk, etc.,
the initialization can be invoked without a subsequent calculation by
specifying ``RUN = 0``.

.. note::
    The settings from the :ref:`PARAMETERS` file and, in a single-domain 
    calculation, the structure from the :ref:`POSCAR` file are read and 
    interpreted **before** the initialization. This is important for 
    the :ref:`ASE API<aseapi>` where structure and settings can be 
    passed programmatically to |calc| in the form of ``ase.Atoms`` objects.

A large number of tasks and checks are performed during initialization.
The major steps are listed below in order of execution.

1.  Create the :file:`SUPP/original_inputs` directory in which all files used
    to start the calculation are stored. This way it is possible to look up
    the used settings, even if, for example, the PARAMETERS file was altered
    by the user during the run. (Notice that all potential input files are
    copied to :file:`SUPP/original_inputs`, whether they will be used for
    the calculation or not.)
2.  If found, read an :ref:`EXPBEAMS.csv or EXBEAMS<expbeams>` file containing
    experimental data. If :ref:`theo_energies` is not defined explicitly, pick
    a slightly expanded version of the energy range of the experimental data
    for use in all calculations.
3.  Check whether the calculation involves
    :ref:`multiple structural domains<domain_calculation>`. If so, perform
    steps 4–9 separately for each domain.
4.  For a domain calculation, read the individual :ref:`PARAMETERS` and
    :ref:`POSCAR` files for each domain.
5.  If necessary, determine the symmetry and minimal unit cell of the
    input structure.
6.  Generate a new :ref:`POSCAR` file with the symmetry applied.
    Write the original structure to :ref:`poscar_oricell`.
7.  If not specified, try to determine the
    :ref:`bulk repeat vector<BULK_REPEAT>` and the
    :ref:`bulk plane group<symmetrybulk>`. Write :ref:`poscar_bulk`
    and :ref:`POSCAR_bulk_appended<poscar_bulk>`.
8.  Check whether a :ref:`PHASESHIFTS` file is present and
    consistent with the structure and settings. If not, and a phase-shifts
    calculation utility is available, generate new phase shifts automatically.
9.  Generate the :ref:`BEAMLIST`, :ref:`patterninfo`, and, if not provided by 
    the user, :ref:`IVBEAMS` files.
