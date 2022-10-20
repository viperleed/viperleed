.. _fortran_comp:

FORTRAN_COMP
============

FORTRAN_COMP defines the fortran compile statement that will be used to compile the non-precompiled parts of the TensErLEED program. Statements that go *after* the file name can be added with the ``post`` flag.

The **search** uses mpirun for parallelization, and therefore requires an MPI compiler. If mpiifort or mpifort are present, they will be used automatically. Optimization can be customized with the ``mpi`` flag. If mpiifort/mpifort or mpirun are not found on the system, the search will be compiled with the standard FORTRAN_COMP compiler and executed *without* parallelization.

**Default:** if ifort is present:

::

   FORTRAN_COMP = 'ifort -O2 -I/opt/intel/mkl/include'
   FORTRAN_COMP post = '-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl'

if ifort is not present, but gfortran is:

::

   FORTRAN_COMP = 'gfortran -O2'
   FORTRAN_COMP post = '-llapack -lpthread'

Additionally, if mpirun is present:

::

   FORTRAN_COMP mpi = 'mpiifort -Ofast'

Or, if mpiifort is not present, but mpifort is:

::

   FORTRAN_COMP mpi = 'mpifort -Ofast -no-pie'

**Syntax:**

::

   FORTRAN_COMP = 'ifort -O3 -I/opt/intel/mkl/include'
   FORTRAN_COMP post = '-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl'
   FORTRAN_COMP mpi = 'mpiifort -O2'

OR

::

   FORTRAN_COMP = ifort
   FORTRAN_COMP = gfortran
   FORTRAN_COMP mpi = mpifort
   FORTRAN_COMP mpi = mpiifort

**Acceptable values**: Either ``ifort`` or ``gfortran`` (or ``mpifort`` for ``mpiifort`` for ``FORTRAN_COMP mpi``) without quotation marks (this will add the default flags listed above for both ``FORTRAN_COMP`` and ``FORTRAN_COMP post``), or any string, with quotation marks on the left and right.

This compile statement will be used for all fortran code compiled while the program runs. Any optimization statements should be included in the string. Needless to say, the fortran compiler needs to be installed on the computer that executes the script.

Whether the default flags work depends not only on the ifort or gfortran packages being present, but also on additional packages and the local library structure. Therefore, explicitly declaring ``FORTRAN_COMP`` and ``FORTRAN_COMP post`` as strings is generally preferable.

.. note:: 
   -  **mpifort**: mpifort is only a "wrapper" compiler, which effectively calls an underlying back-end fortran compiler (usually gfortran) with a number of flags suitable for mpi scripts. If there are problems or you want to see exactly what mpifort is doing, use ``mpifort -showme`` to see what the ``mpifort`` command will be interpreted as (see also the `mpifort man page <https://www.open-mpi.org/doc/v4.0/man1/mpifort.1.php>`__). If this does not work well, you may want to try replacing the ``FORTRAN_COMP mpi = mpifort`` with the output of the ``-showme`` command as an explicit string, fixing paths if necessary.
   -  **mpifort**: if compiling with mpifort causes an error like :literal:`Symbol `time' causes overflow in R_X86_64_PC32 relocation` in the search log (use :ref:`LOG_SEARCH<LOG_SEARCH>`  to produce such a log), this can be resolved by using mpifort with the additional flag ``-no-pie`` (set by default by ViPErLEED for the ``mpifort`` option, see above).

.. warning::
   -  **gfortran/mpifort**: When using aggressive (``-Ofast``) optimization flags, checks for NaNs and +/-Inf values are disabled by the compiler. This poses no known problems for TensErLEED up to at least v.1.73, but it could lead to unexpected behavior in the future. Use the flag ``-fno-finite-math-only`` to re-enable these checks.

**TODO**: We have to add instructions on how to install the compilers!

Natively running on (64-bit) Windows
------------------------------------

Here are some notes on which steps are needed to run (tested up to refcalc) natively on Windows (test only from python source), i.e., get a working Fortran compiler with LAPACK/BLAS. The notes below are for gfortran (gcc), and for the very basic, unoptimized LAPACK/BLAS versions. Hence, execution of the code will be rather slow.

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

.. code-block:: bash

   cd ~
   wget http://www.math.ucla.edu/~wotaoyin/software/lapack_test.cpp  # download
   g++ lapack_test.cpp -llapack -o lapack_test     # build
   ./lapack_test                                   # run

For actually running, set ``FORTRAN_COMP`` as follows:

::

   FORTRAN_COMP = 'gfortran -O2 -std=legacy'        # -std=legacy makes it work for Fortran77
   FORTRAN_COMP post = '-llapack -lblas -lpthread'  # NOTE: order of LAPACK and BLAS is important!
