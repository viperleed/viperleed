.. include:: /substitutions.rst

.. _fortran_comp:

FORTRAN_COMP
============

FORTRAN_COMP defines the fortran compile statement that will be used to compile
the non-precompiled parts of the TensErLEED program. Statements that go *after*
the file name can be added with the ``post`` flag.

The :ref:`structure search<sec_search>` can use mpirun for parallelization.
This is optional, but highly recommended. If you use :term:`MPI` you require
an MPI compiler. If ``mpiifort`` or ``mpifort`` are present, they will be used
automatically. Optimization can be customized with the ``mpi`` flag.
If ``mpiifort``/``mpifort`` or mpirun are not found on the system, the search
will be compiled with the standard FORTRAN_COMP compiler and executed *without*
parallelization.

See the :ref:`installation section<installation>` for details on how to install
the the Fortran compilers, the :term:`MPI` wrappers and ``mpirun``.

**Default:** if :term:`ifort` is present:

.. tab-set::

  .. tab-item:: Linux, macOS, Windows Subsystem for Linux
    :sync: unix

    ::

       FORTRAN_COMP = 'ifort -O2 -I/opt/intel/mkl/include'
       FORTRAN_COMP post = '-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -traceback'

  .. tab-item:: Windows
    :sync: win

    ::

       FORTRAN_COMP = 'ifort /O2 /I"%MKLROOT%/include"'
       FORTRAN_COMP post = '/Qmkl:parallel /traceback'

    Where ``%MKLROOT%`` is the path that is automatically set as
    an environment variable by Intel oneAPI® when it is initialized.


if :term:`ifort` is not present, but gfortran is:

::

   FORTRAN_COMP = 'gfortran -O2'
   FORTRAN_COMP post = '-llapack -lpthread -lblas -fbacktrace'

Additionally, if :program:`mpirun` or, on Windows, :program:`mpiexec` are
present:

.. tab-set::

  .. tab-item:: Linux, macOS, Windows Subsystem for Linux
    :sync: unix

    ::

       FORTRAN_COMP mpi = 'mpiifort -Ofast'

  .. tab-item:: Windows
    :sync: win

    ::

       FORTRAN_COMP mpi = 'cmd mpiifort /Ofast'

    On Windows, :program:`mpiifort` is a :file:`.bat` file, not an executable.
    Python's subprocess cannot directly run it. This is why it is prepended by
    ``cmd`` (i.e., the Windows Command Prompt).


Or, if :program:`mpiifort` is not present, but :program:`mpifort` is:

::

   FORTRAN_COMP mpi = 'mpifort -Ofast -fallow-argument-mismatch'

.. warning::
   If you are using an older version of gfortran packaged with GCC 9 or
   earlier, you need to remove the compiler flag ``-fallow-argument-mismatch``.
   |calc| should do this automatically if it can detect the compiler version.
   Should this fail, you can explicitly specify in :ref:`PARAMETERS`:

   ::

      FORTRAN_COMP mpi = 'mpifort -Ofast'

   This is necessary, because type-checks were made stricter in GCC 10,
   making ``-fallow-argument-mismatch`` mandatory to compile unaltered
   TensErLEED. However, earlier versions of GCC and gfortran may not
   recognize the flag.

.. versionchanged:: 0.14.1
    Added correct defaults for ``ifort``/``mpiifort`` on Windows. On
    earlier versions, Windows users needed to manually set ``FORTRAN_COMP``.


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

**Acceptable values**: Either ``ifort`` or ``gfortran`` (or
``mpifort``/``mpiifort`` for ``FORTRAN_COMP mpi``) (this will add the default
flags listed above for both ``FORTRAN_COMP`` and ``FORTRAN_COMP post``),
or a string containing the full compile statement, including optimization
flags, libraries, etc. Quotation marks are optional and will be stripped
from the string.

This compile statement will be used for all Fortran code compiled while
the program runs. Any optimization statements should be included in the
string. Needless to say, the Fortran compiler needs to be installed on
the computer that executes the script.

Whether the default flags work depends not only on the ifort or gfortran
packages being present, but also on additional packages and the local
library structure. Therefore, explicitly declaring ``FORTRAN_COMP`` and
``FORTRAN_COMP post`` as strings is generally preferable.

.. note::
   -  **mpifort**: :term:`mpifort` is only a "wrapper" compiler, which
      effectively calls an underlying back-end fortran compiler (usually
      gfortran) with a number of flags suitable for mpi scripts. If there
      are problems or you want to see exactly what mpifort is doing, use
      ``mpifort -showme`` to see what the ``mpifort`` command will be
      interpreted as (see also the
      `mpifort man page <https://www.open-mpi.org/doc/v4.0/man1/mpifort.1.php>`__).
      If this does not work well, you may want to try replacing the
      ``FORTRAN_COMP mpi = mpifort`` with the output of the ``-showme``
      command as an explicit string, fixing paths if necessary.
   -  **mpifort**: if compiling with mpifort causes an error like
      :literal:`Symbol `time' causes overflow in R_X86_64_PC32 relocation`
      in the search log (use :ref:`LOG_SEARCH` to produce such
      a log), this can be resolved by using mpifort with the additional flag
      ``-no-pie``.

.. warning::
   -  **gfortran/mpifort**: When using aggressive (``-Ofast``) optimization
      flags, checks for NaNs and +/-Inf values are disabled by the compiler.
      This poses no known problems for TensErLEED up to at least v.1.7.3, but
      it could lead to unexpected behavior in the future. Use the flag
      ``-fno-finite-math-only`` to re-enable these checks.
