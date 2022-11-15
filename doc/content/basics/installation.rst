.. _installation:

============
Installation
============



Below, we provide a short guide on how to install all required components of tleedm.

**TODO**: Tleedm refactor so we can install it!!!

Python
======

**TODO !!! - needs package refactor first**

A conda environment and a pip environment are available.


Fortran compilers
=================

The tleedm (TensErLEED manager) package acts as a wrapper and feature extension to the :term:`TensErLEED` package.
It requires TensErLEED source files to be present and will compile them (with task-specific adjustments) *at run-time*.
This requires the presence of a suitable :term:`Fortran` 77 & 90 compiler on the system.
Unlike the original version of TensErLEED by Blum and Heinz :cite:p:`blumFastLEEDIntensity2001a`, all TensErLEED versions supported by ViPErLEED (TensErLEED >= 1.6) also require :term:`BLAS` and :term:`LAPACK` dependencies to be available.
By default ViPErLEED supports :term:`gfortran` from the GNU Compiler Collection (:term:`gcc`) and the Intel Fortran compiler :term:`ifort`.
ViPErLEED will default to using :term:`ifort` if available.
Use the parameter :ref:`FORTRAN_COMP<fortran_comp>` to adjust this behavior.
You can also use the :ref:`FORTRAN_COMP<fortran_comp>` parameter to use any other compiler installed on your system.

The :ref:`structure-search section<sec_search>`, which is the computationally most expensive part of ViPErLEED and TensErLEED, supports compilation and execution with :term:`MPI`.
To use the :term:`MPI` version of TensErLEED, you need to also install an :term:`MPI` implementation and the :term:`MPI` compiler corresponding to your Fortran compiler.
We recommend using Open MPI on Linux and MacOS.
The MPI compiler for :term:`gfortran` is :term:`mpifort`, for :term:`ifort` it is :term:`mpiifort` (sic!).

If you are running ViPErLEED on a :term:`HPC` system, appropriate Fortran compilers and a preferred :term:`MPI` implementation are likely already installed.
Please consult the documentation for your system and the administrators of details regarding their usage.

:term:`tleedm` can run on Linux, MacOS and Microsoft Windows, but the installation of the compilers in particular differs significantly for each system.

.. note:: 

    -  If you are running on an Intel-processor-based system, we recommend using ``ifort``. It is known from experience to give better performance for TensErLEED.
    -  Using the :term:`MPI` version of TensErLEED is not strictly required, but **highly** recommended.
       Execution times for the :ref:`structure search<sec_search>` may be significantly higher without :term:`MPI`.
       A working MPI implementation is necessary to make use of multi-processing in the :ref:`structure-search section<sec_search>`, even if you are working on a single processor.



``ifort`` and ``mpiifort``
----------------------------

.. _ifort_linux:

Linux
#####

Installation of the Intel compilers and :term:`MPI` implementation for Linux can be performed using a few shell commands.
In this manual, we use ``apt``, the standard package-manager for Debian based distributions.
For installation instructions with other package-managers see the `guides by Intel <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top.html>`__.

As a first step, update the package index:

.. code-block:: console

    $ sudo apt update && sudo apt upgrade

Then follow the `instructions by Intel to add the Intel oneAPI repository <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html#apt>`__.
Following this, you can install the required packages with the package-manager.
For ViPErLEED you need the Intel Base Toolkit (``intel-basekit``) and the Intel HPC Toolkit (``intel-hpckit``):

.. code-block:: console

    $ sudo apt install intel-basekit -y
    $ sudo apt install intel-hpckit -y

.. note:: The toolkits are multiple GB in size and will take a while to download and install.

After installation, we still need to configure the system and add the compilers to our path (see also `here <https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-hpc-linux/top/before-you-begin.html#before-you-begin>`__).
First, we need to make sure required build tools (such as Cmake) are present:

.. code-block:: console

    $ sudo apt install cmake pkg-config build-essential -y

Then, we finally need to configure the Intel one API installation such that it is discovered by by our environment.
For this, we need to source the file `/opt/intel/oneapi/setvars.sh` which sets the required :term:`CLI` arguments.
We recommend you do this by adding the following line to the end of your shell startup script (usually `~/.bashrc`):

.. code-block:: console

    . /opt/intel/oneapi/setvars.sh

Afterwards, the required compilers should be available for use.
You can check if :term:`ifort` is present using:

.. code-block:: console

    $ which ifort

If the result is a path, it means that the shell knows the compiler exists.
You can do the same check with `mpirun` and `mpiifort` to check that they are properly configured as well.

macOS
#####

.. warning::
    Newer Macs using "Apple Silicon" ARM-based chips are incompatible with the Intel compilers (since they don't use Intel chips).
    Use :term:`gfortran` and :term:`mpifort` instead.

To install the Intel oneAPI Toolkits under macOS please follow `the guide provided by Intel <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-macos/top.html>`__.
As for Linux, you will need to install the Intel Base Toolkit and the Intel HPC Toolkit.

Windows
#######

.. warning::
    To run tleedm and TensErLEED under Windows, we recommend using the :term:`Windows Subsystem for Linux<WSL>` (available starting from Windows 10).
    Follow the `instructions by Microsoft to install the WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`__.
    With the :term:`WSL` installed, you can follow the same instructions as provided in `the Linux section<ifort_linux>`.
    Running natively on Windows is possible (:ref:`see below<native_windows>`), but experimental and *not recommended*.

To install the Intel oneAPI Toolkits under Windows please follow `the guide provided by Intel <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-windows/top.html>`__.
As for Linux, you will need to instal the Intel Base Toolkit and the Intel HPC Toolkit.


``gfortran`` and ``mpifort``
----------------------------

Below, we provide a simple guide on how to install the GNU Fortran compiler :term:`gfortran`\ [#.], the Open MPI implementation and the :term:`gfortran` MPI wraper :term:`mpifort`.


Linux
#####

First, using your distributions package-manager, update the package list and install the newest version of :term:`gfortran`.
In this manual, we use ``apt``, the standard package-manager for Debian based distributions.\ [#]_


.. code-block:: console

    $ sudo apt update
    $ sudo apt install gfortran -y

The compiler can be invoked with the ``gfortran`` command.
You can show the version and check if :term:`gfortran` was installed properly using

.. code-block:: console
    
    $ gfortran --version

In addition to :term:`gfortran`, you also need to install the :term:`BLAS` and :term:`LAPACK` libraries.

.. code-block:: console
    
    $ sudo apt install libblas-dev liblapack-dev

Next install Open MPI (or alternatively another MPI implementation of your choosing) to make ``mpirun`` available:

.. code-block:: console
    
    $ sudo apt install openmpi-bin

Finally, install the :term:`gfortran` MPI wrapper ``mpifort``:

.. code-block:: console

    $ sudo apt install libopenmpi-dev


macOS
#####


For running under MacOS, it is recommended to first install a package manager such as `brew <https://brew.sh>`__.
This will also install the XCode Command Line Tools which are required for installing most other components.

Using the ``brew`` command, you can then easily install gfortran and the Open MPI implementation (automatically including ``mpifort``).

.. code-block:: console

    $ brew install gfortran
    $ brew install open-mpi

There is no need to install :term:`BLAS` and :term:`LAPACK`, as MacOS already ships with these libraries pre-installed.

.. warning:: 
    If the XCode Command Line Tools are not installed before you install :term:`gfortran`, you will get an error stating that the ``-lSystem`` library is not available.
    If this happens, make sure to first install the XCode Command Line Tools via
    
    .. code-block:: console

        $ xcode-select --install

    and then reinstall :term:`gfortran`:

    .. code-block:: console

        $ brew reinstall gfortran

Windows
#######

.. warning::
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
-  Update MSys2 running

   .. code-block:: console

        $ pacman -Syu

-  Install gfortran and other useful stuff via

   .. code-block:: console

        $ pacman -S mingw-w64-x86_64-toolchain

-  Add the ``<path_to_mingw_installation>/mingw64/bin`` path to your ``%PATH%`` environment variable (this way, calling gfortran from shell will find the one just installed with no need to explicitly passing the whole path).

-  Install dev tools, cmake and git  with

   .. code-block:: console

       $ pacman -S base-devel
       $ pacman -S mingw-w64-x86_64-cmake
       $ pacman -S git

-  Clone the LAPACK git repository with

   .. code-block:: console

        $ git clone https://github.com/msys2/MINGW-packages.git

   This is the 'basic', unoptimized version.
   There are ways to also build better versions (see `here <https://icl.cs.utk.edu/lapack-for-windows/lapack/>`__).
-  Move to LAPACK directory with 
   
   .. code-block:: console

        $ cd MINGW-packages/mingw-w64-lapack
-  Build LAPACK and BLAS pacakges with 

   .. code-block:: console

        $ makepkg-mingw

   Should ``curl`` complain about some certificates, you can also `download <http://www.netlib.org/lapack/>`__ the LAPACK/BLAS source code as a ``.tar.gz`` archive.
   Take the version that ``curl`` complains about, and place the archive in the package folder (which you can find in ``<path_to_mingw_installation>/home/<user_name>/MINGW-packages/mingw-w64-lapack``).
   This build will take quite a while.

-  Install LAPACK/BLAS packages with
   
   .. code-block:: console

        $ pacman -U mingw-w64-x86_64-lapack-<REPLACE_WITH_VERSION>.pkg.tar.zst

    Note, the archive may have a different suffix.
    Run ``ls`` in the same folder to check the correct name.

You can then test the LAPACK installation with:

.. code-block:: console

   $ cd ~
   $ wget http://www.math.ucla.edu/~wotaoyin/software/lapack_test.cpp  # download
   $ g++ lapack_test.cpp -llapack -o lapack_test     # build
   $ ./lapack_test                                   # run

For actually running, set :ref:`FORTRAN_COMP<fortran_comp>` as follows:

**TODO** Michele: is -std=legacy required on native Windows?

::

   FORTRAN_COMP = 'gfortran -O2 -std=legacy'        # -std=legacy makes it work for Fortran77
   FORTRAN_COMP post = '-llapack -lblas -lpthread'  # NOTE: order of LAPACK and BLAS is important!


To compile the static files described :ref:`below<static_compile>`, go into ``viperleed/tensorleed`` and call:

.. code-block:: console

   gfortran beamgen_source/beamgen.v1.7.f -o beamgen.v1.7 -Ofast -fno-finite-math-only
   gfortran eeasisss_code/modified/imported_routines.f90 eeasisss_code/modified/eeasisss.f90 -o EEASiSSS.x -Ofast -fno-finite-math-only
   del "*.mod"



.. _static_compile:

Compiling static files
======================

In addition to the TensErLEED source code, which is compiled *at run-time*, ViPErLEED needs a few auxilary scripts that need compiling before a calculation can be started.
These can be compiled automatically using a provided Makefile.
To do this, go into the ``tensorleed`` folder in the ``viperleed`` directory and call:

.. code-block:: console

   $ make all

The Makefile uses the ``gfortran`` compiler by default, if you only have ``ifort`` installed, change the variable ``gcomp`` in the first line of the Makefile accordingly.


.. [#] See also `here <https://fortran-lang.org/en/learn/os_setup/install_gfortran/>`__ for a guide on how to install gfortran on various operating systems.

.. [#] For other distributions have a look at, for example, this tutorial `<https://fortran-lang.org/en/learn/os_setup/install_gfortran/>`__.