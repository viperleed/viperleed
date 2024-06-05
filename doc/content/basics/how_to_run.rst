.. include:: /substitutions.rst

.. _how_to_run:

==================================
Directory structure and how to run
==================================

In ViPErLEED, each calculation must have its own directory.
The input and output files have fixed names (case sensitive!),
see :ref:`the list of files<list_input_files>`. Below, we give an
example directory tree with the files needed to start a |LEED-IV|
calculation in ViPErLEED.
See the :ref:`installation instructions<installation>` for a guide on
how to install ViPErLEED. See also the ViPErLEED examples section for
some working calculations to help you get started.

Minimum input
=============

To set up a ViPErLEED calculation, first create a source directory
(in this example ``my_surface``) and place all input files inside.

.. code-block:: console
    :caption: Minimum input directory tree

    my_surface
    ├── EXPBEAMS.csv [and/or IVBEAMS]
    ├── POSCAR
    └── PARAMETERS

The minimum information required to start
a  :ref:`reference calculation<ref-calc>`
is contained in three files:

-   :ref:`EXPBEAMS.csv<expbeams>`: This file contains the experimentally
    measured |LEED-IV| curves in :term:`CSV` format.
    Unless specified otherwise in :ref:`PARAMETERS<parameters>`, ViPErLEED
    will use the information in the :ref:`EXPBEAMS.csv file<expbeams>` to
    set energy ranges (:ref:`THEO_ENERGIES<theo_energies>`) and choose which
    beams should be written to the output file (:ref:`file IVBEAMS<ivbeams>`).
    Alternatively, you can also provide the :ref:`IVBEAMS file<ivbeams>`
    directly.
-   :ref:`POSCAR<poscar>`: This file contains the structure for the
    :ref:`reference calculation<ref-calc>`. ViPErLEED will determine the
    applicable symmetry from the :ref:`POSCAR file<poscar>`. See also the
    :ref:`symmetry settings<symmetry_settings>` and
    :ref:`input-structure settings<input_structure_settings>`.
-   :ref:`PARAMETERS<parameters>`: This file contains the settings for the
    calculation (see :ref:`the list of parameters<paramname>`). If no
    :ref:`VIBROCC file<vibrocc>` is given, :ref:`PARAMETERS<parameters>`
    needs to contain values for :ref:`T_EXPERIMENT<t_experiment>` and
    :ref:`T_DEBYE<t_debye>`.


To run also a :ref:`delta amplitudes calculation<sec_deltas>` and/or a
:ref:`structure search<sec_search>`, you additionally need to provide a
:ref:`DISPLACEMENTS file<displacements>` that contains the requested
perturbations of the structure.

.. note::
    Setting up a :ref:`domain calculation<domain_calculation>` with multiple
    coexisting surface structures, requires a slightly different directory
    tree. See the :ref:`domain-calculation page<domain_calculation>`.


Starting the calculation
========================

Once you have set up the input files you are ready to start the calculation.

You can then start a ViPErLEED calculation by invoking ``viperleed calc``
(or |calc|) via the command line. A list of all available command-line
options for |calc| can be found :ref:`here<cli_calc>`.

A typical call may look like this:

.. code-block:: console

    $ tensorleed_path="path/to/tensorleed_dir"
    $ work_path="path/to/work_dir"
    $
    $ viperleed calc -w $work_path -t $tensorleed_path

Here, the work directory ``work_dir`` is the directory where the calculation
will be executed and all temporary files will be stored.
``work_dir`` will be created if it does not yet exist.
The ``tensorleed_path`` is the path to the TensErLEED source code.
If the ``-t`` option is not given, ViPErLEED will look for the TensErLEED
source code under the ``VIPERLEED_TENSORLEED`` :term:`environment variable`.


.. tip::
    As a calculation can take a long time, it is recommended to start the
    calculation using `nohup <https://en.wikipedia.org/wiki/Nohup>`__ or
    in a `tmux <https://github.com/tmux/tmux/wiki>`__ session. This way,
    the calculation will not be aborted if the user is logged out (or the
    connection of an ``ssh`` session breaks).

HPC systems
-----------

If you are running ViPErLEED on an :term:`HPC` system with a workload scheduler
such as `slurm <https://slurm.schedmd.com/documentation.html>`__, make sure to
load the required compilers, :term:`MPI` implementations, and Python
packages/environment in the submission-script (system-dependent, e.g.,
``module load mpiifort``).

Such a submission script usually contains details on the requested hardware
(e.g., declared via ``#SBATCH`` in slurm) and instructions on which precompiled
packages to make available.
Below, you find an example for a submission script for the
`Vienna Scientific Cluster (VSC-4) <https://vsc.ac.at//home/>`__, which
uses the slurm workload manager.
The script first loads the required Intel compilers and
:term:`conda` distribution, before executing ViPErLEED
using the :ref:`viperleed command<cli_calc>`.


.. literalinclude :: /_static/example_job_script.txt
   :language: bash
   :caption: Example submission script for the job script.

.. _dir_organization_output:

Output organization
===================

A large number of files are created in the work directory is executed in
defined by the :ref:`-w option<cli_calc>` of |calc|. ViPErLEED will copy
input files there, run all requested calculations, and then copy the
relevant output files back to the input directory. For this purpose,
ViPErLEED also creates a :ref:`manifest` file that lists the relevant
output files which will be copied back.

The directory tree after a run may look something like this:

.. code-block:: console
    :caption: Normal output directory tree

    my_surface
    ├── EXPBEAMS.csv
    ├── POSCAR
    ├── POSCAR_user
    ├── PARAMETERS
    ├── IVBEAMS
    ├── VIBROCC
    ├── PHASESHIFTS
    ├── DISPLACEMENTS
    ├── work
    │   ├── manifest
    │   └── ...
    ├── OUT
    │   ├── THEOBEAMS.csv
    │   ├── Rfactor_analysis_refcalc.pdf
    │   └── ...
    ├── SUPP
    │   ├── POSCAR_bulk
    │   └── ...
    ├── Tensors
    │   └── Tensors_001.zip
    ├── Deltas
    │   └── Deltas_001.zip
    └── calc-$timestamp.log

ViPErLEED will create the additional input files :ref:`IVBEAMS<ivbeams>`,
:ref:`BEAMLIST<beamlist>`, :ref:`PHASESHIFTS<phaseshifts>`,
and :ref:`VIBROCC<vibrocc>` under certain conditions; see
the respective pages for details.

The original :ref:`POSCAR file<poscar>` is renamed to ``POSCAR_user`` while
the new ``POSCAR`` contains the structure as interpreted by ViPErLEED. For
details see the page on the :ref:`POSCAR file<poscar>`.

After the first run, an ``OUT`` directory is created that contains the output
files, see the :ref:`list of output files<output_files>` for details.
ViPErLEED further produces additional :ref:`supplementary files<supp_files>`
that are required during execution, that contain intermediate results or that
may be of interest for debugging purposes.
These files are stored in the ``SUPP`` subfolder.

If a :ref:`refercence calculation<ref-calc>` is run with
:ref:`Tensor output<toutput>`, a ``Tensors`` directory will
be created that stores the :ref:`tensor files<tensorszip>`.
Similarly, if a :ref:`delta-amplitudes<sec_deltas>` calculation is run,
a ``Deltas`` directory will be created that contains the resulting
:ref:`delta files<deltaszip>`.

In case of automated multiple search runs (which can be specified in the
:ref:`DISPLACEMENTS<DISPLACEMENTS>` file), |calc| creates a ``workhistory``
directory.
A snapshot of all input and output files that may be relevant and may get
overwritten will be moved into a subfolder of ``workhistory``.
