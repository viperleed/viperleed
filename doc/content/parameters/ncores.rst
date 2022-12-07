.. _ncores:

N_CORES
=======

N_CORES defines the number of cores available for TensErLEED execution.

**Default:** Determine available cores on system automatically (might fail)

**Allowed values:** positive integer

**Syntax:**

::

   N_CORES = 16

The TensErLEED search program is executed using ``mpirun``, with the number of cores specified by N_CORES. If mpirun / mpiifort are not present on the computer running tleedm, no parallelization will be performed in the search (see also :ref:`FORTRAN_COMP<FORTRAN_COMP>`). Reference and Delta calculations instead use python multiprocessing to execute multiple single-threaded TensErLEED processes, so those can run without a special fortran compiler.
