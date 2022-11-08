.. _installation:

============
Installation
============


**TODO** Needs an own docu page. 
Python 3.8. Info on Fortran compiler currently on FORTRAN_COMP page should be moved there. 
Does the Windows version require Windows-Subsystem für Linux? !!

Below, we provide a short guide on how to install all required components of tleedm.

**TODO**: Tleedm refactor so we can install it!!!

Python
======

**TODO**

A conda environment and a pip environment are available.


Fortran compilers
=================

The tleedm (TensErLEED manager) package acts as a wrapper and feature extension to the :term:`TensErLEED` package.
It requires TensErLEED source files to be present and will compile them (with task-specific adjustments) *at run-time*.
This requires the presence of a suitable :term:`Fortran` 77 & 90 compiler on the system.
Unlike the original version of TensErLEED by Blum and Heinz :cite:p:`blumFastLEEDIntensity2001a`, all TensErLEED versions supported by ViPErLEED (TensErLEED >= 1.6) also require :term:`BLAS` and :term:`LAPACK` dependencies to be available.
ViPErLEED supports :term:`gfortran` from the GNU Compiler Collection (:term:`gcc`) and the Intel Fortran compiler :term:`ifort`.
ViPErLEED will default to using :term:`ifort` if available.
Use the parameter :ref:`FORTRAN_COMP<fortran_comp>` to adjust this behavior.

The :ref:`structure search section<sec_search>`, which is the computationally most expensive part of ViPErLEED and TensErLEED, supports compillation and execution with :term:`MPI`.
To use the :term:`MPI` version of TensErLEED, you need to also install an :term:`MPI` implementation and the :term:`MPI` compiler corresponding to your Fortan compiler.
We recommend using Open MPI on Linux and MacOS.
The MPI compiler for :term:`gfortran` this is :term:`mpifort`, for :term:`ifort` it is :term:`mpiifort` (sic!).

If you are running ViPErLEED on a :term:`HPC` system, appropriate Fortran compilers and a prefered :term:`MPI` implementation are likely already installed.
Please consult the documentation for your system and the administrators of details regarding their usage.


.. note:: 

    -  If you are running on an Intel-processor based system, we recommend using ``ifort``. It is known from experience to give better performance for TensErLEED.
    -  Using the :term:`MPI` version of TensErLEED is not strictly required, but **highly** recommended.
       Execution times for the :ref:`structure search<sec_search>` may be significantly higher without :term:`MPI`.
       A working MPI implementation is necessary to make use of multi-processing in the :ref:`structure search section<sec_search>`, even if you are working on a single processor.


tleedm can run on Linux, MacOS and Microsoft Windows, but the installation of the compilers in particular differs significantly for each system.


``gfortran`` and ``mpifort``
----------------------------

Linux
#####

First, using your distributions package-manager, update the package list and install the newest version of :term:`gfortran`.
In this manual, we use ``apt-get``, the standard package-manager for Debian based distributions.\ [#]_


.. code-block:: console

    $ sudo apt-get update
    $ sudo apt-get install gfortran -y

The compiler can be invoced with the ``gfortran`` command.
You can show the version and check if :term:`gfortran` was installed properly using:

.. code-block:: console
    
    $ gfortran --version

In addition to :term:`gfortran`, we also need to install the :term:`BLAS` and :term:`LAPACK` libraries, as they are required for the :ref:`reference calculation section<ref-calc>`:

.. code-block:: console
    
    $ sudo apt-get install libblas-dev liblapack-dev

Next install Open MPI (or alternatively another MPI implementation of your choosing) to make ``mpirun`` available:

.. code-block:: console
    
    $ sudo apt-get install openmpi-bin

Finally, install the :term:`gfortran` MPI wrapper ``mpifort``:

.. code-block:: console

    $ sudo apt-get install libopenmpi-dev


MacOS
#####

.. note:: 
    Newer Macs using "Apple Silicon" ARM-based chips are incompatible with the Intel compilers.
    Use :term:`gfortran` and :term:`mpifort` instead.

For running under MacOS, it is recommened to first install a package manager such as `brew <https://brew.sh>`__.
This will also install the XCode Command Line Tools which are required for installing most other components.

Using the ``brew`` command, you can then easily install gfortran and the Open MPI implementation (automatically including ``mpifort``).

.. code-block:: console

    $ brew install gfortran
    $ brew install open-mpi

There is no need to install :term:`BLAS` and :term:`LAPACK`, as MacOS already ships with these libraries pre-installed.


``ifort`` and ``mpiifort``
----------------------------

**TODO**

As a first step, update the package index:

.. code-block:: console

    $ sudo apt-get update

Then follow the instructions 

For ViPErLEED you need the Intel Base Toolkit (``intel-basekit``) and the Intel HPC Toolkit (``intel-hpckit``).

Windows
#######

To run tleedm and TensErLEED under Windows, we recommend using the :term:`Windows Subsystem for Linux<WSL>`.
Follow the `instructions by Microsoft to install the WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`__.
With the :term:`WSL` installed, you can follow the same instructions as provided below for Linux.
Running natively on Windows is possible (:ref:`see below<native_windows>`), but experimental and *not recommended*.


.. _native_windows:

Natively running on (64-bit) Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here are some notes on which steps are needed to run (tested up to refcalc) natively on Windows (test only from python source), i.e., get a working Fortran compiler with LAPACK/BLAS
The notes below are for gfortran (gcc), and for the very basic, unoptimized LAPACK/BLAS versions.
Hence, execution of the code will be rather slow.

-  Install `MSys2 <https://www.msys2.org/>`__, which then installs MinGW, then open the MSys2 shell.
-  Update MSys2 running ``pacman -Syu``
-  Install gfortran and other useful stuff via ``pacman -S mingw-w64-x86_64-toolchain``
-  Add the ``<path_to_mingw_installation>/mingw64/bin`` path to your ``Path`` environment variable (this way, calling gfortran from shell will find the one just installed with no need to explicitly passing the whole path)
-  Install dev tools with ``pacman -S base-devel``
-  Install cmake with ``pacman -S mingw-w64-x86_64-cmake``
-  Install git with ``pacman -S git``
-  Clone the LAPACK git repository with ``git clone https://github.com/msys2/MINGW-packages.git`` This is the 'basic', unoptimized version. There are ways to also build better versions (see `here <https://icl.cs.utk.edu/lapack-for-windows/lapack/>`__).
-  Move to LAPACK directory with ``cd MINGW-packages/mingw-w64-lapack``
-  Build LAPACK and BLAS pacakges with ``makepkg-mingw`` Should ``curl`` complain about some certificates, you can also `download <http://www.netlib.org/lapack/>`__ the LAPACK/BLAS source code as a ``.tar.gz`` archive. Take the version that ``curl`` complains about, and place the archive in the package folder (which you can find in ``<path_to_mingw_installation>/home/<user_name>/MINGW-packages/mingw-w64-lapack``). This build will take quite a while.
-  Install LAPACK/BLAS packages with ``pacman -U mingw-w64-x86_64-lapack-<REPLACE_WITH_VERSION>.pkg.tar.zst`` Note: archive may have a different suffix. Run ``ls`` in the same folder to check the correct name.

You can then test the LAPACK installation with:

.. code-block:: console

   cd ~
   wget http://www.math.ucla.edu/~wotaoyin/software/lapack_test.cpp  # download
   g++ lapack_test.cpp -llapack -o lapack_test     # build
   ./lapack_test                                   # run

For actually running, set :ref:`FORTRAN_COMP<fortran_comp>` as follows:

::

   FORTRAN_COMP = 'gfortran -O2 -std=legacy'        # -std=legacy makes it work for Fortran77
   FORTRAN_COMP post = '-llapack -lblas -lpthread'  # NOTE: order of LAPACK and BLAS is important!



Directory structure and file names
==================================


In case of automated multiple search runs (which can be specified in the :ref:`DISPLACEMENTS<DISPLACEMENTS>`  file), tleedm creates a “workhistory” directory and moves a snapshot of all input and output files that may be relevant and may get overwritten into a subfolder there.









.. code-block:: console

    normal

    ├── EXPBEAMS.csv
    ├── PARAMETERS
    ├── work          
    │   ├── area.py
    │   └── bboxinout.py
    ├── pywps.cfg          
    ├── requirements.txt
    └──  server.py          



::

    after

    ├── EXPBEAMS.csv
    ├── PARAMETERS
    ├── job.py
    ├── work
    ├── OUT
    └── POSCAR



.. [#] For other distributions have a look at e.g. this tutorial `<https://fortran-lang.org/en/learn/os_setup/install_gfortran/>`__