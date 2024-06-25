.. include:: /substitutions.rst

.. |slurm|  replace:: :program:`slurm`


.. _how_to_run:

=====================================
How to run and directory organization
=====================================

In ViPErLEED, each set of calculations for one system must have its own
directory. This is because the input and output files have **case sensitive**,
fixed names, see :ref:`the list of files<list_input_files>`. A typical approach
is to create one such directory for each experimental dataset and each time
a new structural model (i.e., :ref:`POSCAR` file) is considered.

:numref:`list_minimum_dir_inputs` shows an example directory tree
with the files needed to start a |LEED-IV| calculation in ViPErLEED.
See the :ref:`installation instructions<installation>` for a guide on
how to install ViPErLEED. See also the :ref:`examples` section for
some working calculations to help you get started.

Minimum input
=============

To set up a ViPErLEED calculation, first create a source directory
(in this example ``my_surface``) and place all input files inside.

.. _list_minimum_dir_inputs:
.. code-block:: console
    :caption:
        Minimum input directory tree for
        a :ref:`reference calculation<ref-calc>`.

    my_surface/
    ├── IVBEAMS or EXPBEAMS.csv
    ├── PARAMETERS
    └── POSCAR

The minimum information required to start
a :ref:`reference calculation<ref-calc>`
is contained in three files:

-   :ref:`IVBEAMS` file: This file contains information on
    which beams should be written to the output. It can be generated
    automatically starting from an :ref:`EXPBEAMS.csv<expbeams>` file
    that contains the |LEED-IV| curves measured experimentally. If an
    :ref:`EXPBEAMS.csv<expbeams>` file is provided, |calc| will use it
    also to decide the energy range to calculate (\ :ref:`THEO_ENERGIES`),
    unless specified otherwise in :ref:`PARAMETERS`.
-   :ref:`PARAMETERS`: This file contains the settings for the calculation
    (see :ref:`the list of parameters<paramname>`). If no :ref:`VIBROCC`
    file is given, the :ref:`PARAMETERS` file must contain values for
    :ref:`T_EXPERIMENT` and :ref:`T_DEBYE`.
-   :ref:`POSCAR`: This file contains the structure for the
    :ref:`reference calculation<ref-calc>`. ViPErLEED will determine the
    applicable symmetry from the :ref:`POSCAR` file. See also the
    :ref:`symmetry settings<symmetry_settings>` and
    :ref:`input-structure settings<input_structure_settings>`.


To run also a :ref:`delta amplitudes calculation<sec_deltas>` or a
:ref:`structure search<sec_search>`, you additionally need to provide a
:ref:`DISPLACEMENTS` file that contains the requested
perturbations of the structure.

.. note::
    A :ref:`domain calculation<domain_calculation>` , that is, with multiple
    coexisting surface structures, requires a slightly different directory
    tree. See the :ref:`domain-calculation section<domain_calculation>`.


Starting the calculation
========================

Once you have set up the input files you are ready to start the calculation.
You can start it by invoking ``viperleed calc`` via the command line. A list
of all available command-line options can be found :ref:`here<cli_calc>`.

A typical call may look like this:

.. tab-set::

    .. tab-item:: Linux, macOS, WSL

        .. code-block:: bash

            tensorleed_path="path/to/tensorleed_dir"
            work_path="path/to/work_dir"

            viperleed calc -w $work_path -t $tensorleed_path

        .. tip::
            As |LEED-IV| calculations can take a long time, it is a good
            idea to start them in a "detached" manner, for example using
            `nohup <https://en.wikipedia.org/wiki/Nohup>`__ or in a
            `tmux <https://github.com/tmux/tmux/wiki>`__ session. This
            way, the calculation will not be aborted if the user logs
            out (or the connection of an ``ssh`` session breaks).

    .. tab-item:: Windows

        .. code-block:: bat

            SET tensorleed_path="path\to\tensorleed_dir"
            SET work_path="path\to\work_dir"

            viperleed calc -w %work_path% -t %tensorleed_path%

The directory at ``work_path`` is the one where the calculation will be
executed and all temporary files will be stored. It is created if it does not
exist. The ``tensorleed_path`` is the path to the :ref:`install_tensorleed`.
If the ``-t`` option is not given, ViPErLEED looks for the TensErLEED source
code under the path specified by the :envvar:`VIPERLEED_TENSORLEED`
:term:`environment variable` (see also :ref:`set_envvar`).


HPC systems
-----------

If you are running ViPErLEED on a high-performance computing
(\ :term:`HPC`) system with a workload scheduler such as
`slurm <https://slurm.schedmd.com/documentation.html>`__,
make sure to load the required compilers, :term:`MPI`
implementations, and Python packages/environment in
the submission script (e.g., via ``module load mpiifort``).

Such a submission script usually contains details on the requested hardware
(e.g., declared via ``#SBATCH`` in |slurm|\ ) and instructions on which
precompiled packages to make available.
:numref:`list_slurm_example` shows an example for a submission script for
the `Vienna Scientific Cluster (VSC-4) <https://vsc.ac.at//home/>`__, which
uses the |slurm| workload manager. The script first loads the required Intel
compilers and :term:`conda` distribution, before executing ViPErLEED using
the :ref:`viperleed command<cli_calc>`.

.. _list_slurm_example:
.. literalinclude :: /_static/example_job_script.txt
   :language: bash
   :caption:
        Example submission script for running |calc| using
        |slurm|, as on the Vienna Scientific Cluster.

.. _dir_organization_output:

Directory organization
======================

|calc| executes calculations in a work directory, distinct from the directory
containing the input files. The path to such work directory can be specified
via the ``-w`` :ref:`command-line option<cli_calc>`. There are two main reasons
for running the calculations in a dedicated directory: (i) |calc| creates
a large number of files to handle :term:`TensErLEED` and its input/output,
and (ii) some :term:`HPC` systems do not allow execution of jobs in the
"user space" that typically contains the input files.

Upon starting execution, |calc| copies all the input files into the work
directory, runs all requested calculations, and then copies the relevant
output files back to the input directory. For this purpose, it also creates
a :ref:`manifest` file that lists the files and directories to be copied back.

:numref:`list_organization_output` shows an example of the directory tree
after a run. The :file:`my_work` directory may be on a different file-system
path (including a different drive, a network path, or a remote server) or be
a subfolder of :file:`my_surface`.

.. _list_organization_output:
.. code-block:: console
    :caption: Typical directory tree after a |calc| run.

    my_surface/
    ├── Deltas/
    │   ├── Deltas_001.zip
    │   └── ...
    ├── OUT/
    │   ├── THEOBEAMS.csv
    │   ├── Rfactor_analysis_refcalc.pdf
    │   └── ...
    ├── SUPP/
    │   ├── POSCAR_bulk
    │   └── ...
    ├── Tensors/
    │   ├── Tensors_001.zip
    │   └── ...
    ├── workhistory/
    │   └── ...
    ├── calc-<timestamp>.log
    ├── DISPLACEMENTS        <-- input
    ├── EXPBEAMS.csv         <-- input
    ├── IVBEAMS              <-- input, or created by calc
    ├── PARAMETERS           <-- input, may be edited by calc
    ├── PHASESHIFTS          <-- input, or created by calc
    ├── POSCAR               <-- input, may be edited by calc
    ├── POSCAR_user          <-- input
    └── VIBROCC              <-- input, or created by calc

    my_work/
    ├── manifest
    └── ...

|calc| will create the additional input files :ref:`IVBEAMS`,
:ref:`PHASESHIFTS`, and :ref:`VIBROCC` if not provided by the
user; see the respective sections for details.

During the very first run, the original :ref:`POSCAR` file is renamed to
:file:`POSCAR_user`, while the new :file:`POSCAR` contains the structure
as interpreted by |calc|. See the section on the :ref:`POSCAR` file for
more details.

Information about the progress of the calculation and about errors
in the user input are printed to the terminal and recorded in the
:ref:`log_files_calc` file. It is always a good idea to check the
log file thoroughly for ``WARNING`` and ``ERROR`` messages.

The ``OUT`` directory (created automatically) contains the results of the
calculation, see the :ref:`list of output files<output_files>` for details.
|calc| also produces some :ref:`supplementary files<supp_files>`, stored
in the ``SUPP`` directory. These files contain intermediate results or may
be of interest for debugging purposes.

If a :ref:`refercence calculation<ref-calc>` is run with
:ref:`Tensor output<TENSOR_OUTPUT>`, a ``Tensors`` directory will
be created that stores the :ref:`tensor files<tensorszip>`.
Similarly, if a :ref:`delta-amplitudes<sec_deltas>` calculation is run,
a ``Deltas`` directory will be created that contains the resulting
:ref:`delta files<deltaszip>`.

In case of automated multiple-search runs (which can be specified in the
:ref:`DISPLACEMENTS` file), |calc| creates a ``workhistory`` directory.
It contains one subfolder for each intermediate run, storing input and
output files that will be overwritten during subsequent runs.

.. note:
    For :ref:`multi-domain calculations<domain_calculation>` the input structure
    will be different, as separate directories are used for the inputs of each
    domain. See the :ref:`domain-calculation section<domain_calculation>` for 
    more details.

.. todo:: Add also a few words about bookkeper here.
