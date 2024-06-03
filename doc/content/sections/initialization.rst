.. include:: /substitutions.rst

.. _initialization:

==============
Initialization
==============


The Initialization is the first section of |calc| that is run. It is
**always** inserted at the beginning of any |calc| calculation, even if
not explicitly specified via the :ref:`RUN parameter<run>`. To test out
symmetry recognition, phaseshift generation, etc. the initialization can
be invoked without a subsequent calculation by specifying ``RUN = 0``.

A large number of tasks and checks will be performed during initialization.
The major steps are listed below in order of execution.

.. note::
    The structure from the :ref:`POSCAR file<poscar>` and the settings from the
    :ref:`PARAMETERS file<parameters>` are read and interpreted **before** the
    initialization. This is important for the :ref:`ASE API<aseapi>` where
    structure and settings can be passed programmatically to |calc| in the form
    of ``ase.Atoms`` objects.

1.  Check whether a :ref:`domain calculation<domain_calculation>` is being
    performed.
    If so, perform all following steps and calculations separately for
    each domain.
#.  Look for an :ref:`EXPBEAMS.csv or EXBEAMS file<expbeams>`
    containing experimental data and read it if found.
    If experimental data is found and :ref:`THEO_ENERGIES<theo_energies>` is
    not defined explicitly, the experimental energy range will be used for
    all calculations.
#.  If necessary, determine the symmetry and minimal unit cell of the
    input structure.
#.  Generate a new :ref:`POSCAR file<poscar>` with the symmetry applied.
    The original structure will be written to
    :ref:`POSCAR_oricell<poscar_oricell>`.
#.  If not specified, try to determine the
    :ref:`bulk repeat vector<BULK_REPEAT>` and the bulk plane
    group. Following this, :ref:`poscar_bulk` and
    :ref:`POSCAR_bulk_appended<poscar_bulk>` will be written.
#.  Check whether a :ref:`PHASESHIFTS file<phaseshifts>` is present and
    consistent with the structure and settings. If not (and a
    phaseshifts utility is available), new phaseshifts will be
    generated automatically.
#.  Generate the :ref:`BEAMLIST<beamlist>`,
    :ref:`PatternInfo.tld<patterninfo>`, and
    :ref:`IVBEAMS<ivbeams>` files (:ref:`IVBEAMS<ivbeams>` is only generated if
    not provided by the user).
#.  Create the directory ``original_inputs`` in which all files used to
    start the calculation are stored.
    This way it is possible to look up the used settings, even if, for example,
    the PARAMETERS file was altered by the user during the run.
