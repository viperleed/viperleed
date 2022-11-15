.. _how_to_run:

==================================
Directory Structure and How To Run
==================================

In ViPErLEED, each calculation must have its own directory. 
The input and output files have fixed names, see :ref:`the list of files<list_input_files>`. 
Below, we give an example directory tree with the files needed to start a LEED :math:`I(V)` calculation in ViPErLEED.
See also the ViPErLEED examples section for some working calculations to help you get started.

Minimum input
=============

To set up a ViPErLEED calculation, first create a source directory (in this example ``my_surface``) and place all input files inside.

.. code-block:: console
    :caption: Minimum input directory tree

    my_surface
    ├── IVEBAMS [and/or EXPBEAMS.csv]
    ├── POSCAR
    ├── PARAMETERS
    └── job.py [or job.sh]

The minimum information required to start a :ref:`refercence calculation<ref-calc>` is contained in three files:

-   :ref:`EXPBEAMS.csv file<expbeams>` contains the experimentally measured LEED-:math:`I(V)` curves.
    Unless specified otherwise in :ref:`PARAMETERS<parameters>`, ViPErLEED will also use the information in the :ref:`EXPBEAMS.csv file<expbeams>` to set energy ranges (:ref:`THEO_ENERGIES<theo_energies>`) and choose which beams should be written to the output file (:ref:`file IVBEAMS<ivbeams>`).
    Alternatively, you can also provide the :ref:`IVBEAMS file<ivbeams>` directly.
-   :ref:`POSCAR file<poscar>` contains the reference surface structure.
    ViPErLEED will determine the applicable symmetry from the :ref:`POSCAR file<poscar>`. See also the :ref:`symmetry settings<symmetry_settings>` and :ref:`input-structure settings<input_structure_settings>`.
-   :ref:`PARAMETERS<parameters>` contains the settings for the calculation (see :ref:`the list of parameters<paramname>`).
    If no :ref:`VIBROCC file<viboccin>` is given, :ref:`PARAMETERS<parameters>` needs to contain values for :ref:`T_EXPERIMENT<t_experiment>` and :ref:`T_DEBYE<t_debye>`.
-   :ref:`job.py / job.sh<job_script>` is the entry point for the ViPErLEED calculation.
    Defines the paths to the ViPErLEED source code and the desired ``work`` directory.

To run also a :ref:`delta amplitudes calculation<sec_deltas>` and a :ref:`structure search<sec_search>`, you additionally need to provide a :ref:`DISPLACEMENTS file<displacements>` that contains the requested perturbations of the structure.

.. note:: 
    Setting up a :ref:`domain calculation<domain_calculation>` with multiple coexisting surface structures, requires a slightly different directory tree.
    See the :ref:`domain-calculation page<domain_calculation>`.

Starting the calculation
========================

Once you have set up the input files you are ready to start the calculation.
If you are running ViPErLEED from Python (recommend), make sure that all tleedm dependencies are in the Python PATH, i.e. make sure the correct conda or virtual environment is loaded, otherwise this will raise an Error.
You can then start a ViPErLEED calculation by invoking the :ref:`job script<job_script>` via the command line.

.. code-block:: console
    
    $ src_path="path/to/source_dir"
    $ wrk_path="path/to/work_dir"
    $
    $ python3 job.py -s $src_path -w $wrk_path

You can also set the source and work directory path directly in the job script, rather than giving them as command line arguments.

If you are running using a pre-packaged version of tleedm, you can start the calculation by running the job shell script ``job.sh``. Make sure to edit the source and work path in the script beforehand.

.. code-block:: console
    
    $ ./job.sh

.. tip:: 
    As a ViPErLEED calculation can take a long time, it is recommended to start the calculation using `nohup <https://en.wikipedia.org/wiki/Nohup>`__ or in a `tmux <https://github.com/tmux/tmux/wiki>`__ session. This way, the calculation will not be aborted if the user is logged out (or the connection of an ``ssh`` session breaks).

If you are running ViPErLEED on an HPC system with a workload scheduler such as `slurm <https://slurm.schedmd.com/documentation.html>`__, make sure to load the required compilers, :term:`MPI` implementations and Python packages in the submission-script.

**TODO** job submission script example?

.. _dir_organization_output:

Output organization
===================

A large number of files are created in the directory that tleedm is executed in.
The :ref:`job script<job_script>` defines the path to a ``work`` directory (typically just a subdirectory of the source directory ``my_surface``) that will be used during the calculation.
ViPErLEED will  copy input files there, execute tleedm, and then copy the relevant output files back to the data directory.
For this purpose, tleedm also creates a :ref:`manifest` file that lists the relevant output files which should be copied back.

The directory tree after a run may look something like this:

.. code-block:: console
    :caption: Normal output directory tree

    my_surface
    ├── EXPBEAMS.csv
    ├── POSCAR
    ├── PARAMETERS
    ├── job.py
    ├── IVBEAMS
    ├── VIBROCC
    ├── PHASESHIFTS
    ├── DISPLACEMENTS
    ├── work
    │   ├── manifest
    │   └── ...
    ├── OUT
    │   ├── THEOBEAMS.csv
    │   └── ...
    ├── SUPP
    │   ├── POSCAR_bulk
    │   └── ...
    └── tleedm-$timestamp.log

ViPErLEED will create the additional input files :ref:`IVBEAMS<ivbeams>`, :ref:`BEAMLSIT<beamlist>`, :ref:`PHASESHIFTS<phaseshifts>`, and :ref:`VIBROCC<viboccin>` under certain conditions; see the respective pages for details.

After the first run, an ``OUT`` directory is created that contains the output files, see the :ref:`list of output files<output_files>` for details.
ViPErLEED further produces additional :ref:`supplementary files<supp_files>` that are required during execution, that contain intermediate results or that may be of interest for debugging purposes.
These files are stored in the ``SUPP`` subfolder.

If a :ref:`refercence calculation<ref-calc>` is run with :ref:`Tensor output<toutput>`, a ``Tensors`` directory will be created that stores the :ref:`tensor files<tensorszip>`.
Similarly, if a :ref:`delta-amplitudes<sec_deltas>` calculation is run, a ``Deltas`` directory will be created that contains the resulting :ref:`delta files<deltaszip>`.

In case of automated multiple search runs (which can be specified in the :ref:`DISPLACEMENTS<DISPLACEMENTS>` file), tleedm creates a ``workhistory`` directory and moves a snapshot of all input and output files that may be relevant and may get overwritten into a subfolder there.