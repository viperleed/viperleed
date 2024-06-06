.. include:: /substitutions.rst

.. _ncores:

N_CORES
=======

N_CORES defines the number of cores available for TensErLEED execution.

**Default:** Determine available cores on system automatically (might fail)

**Allowed values:** positive integer

**Syntax:**

::

   N_CORES = 16

.. warning::

      Some :term:`BLAS` and :term:`LAPACK` libraries (e.g. Intel MKL, which is
      used by default with ifort) support use of multiple threads.
      This can conflict with the parallelization used by ViPErLEED and
      **massively** reduce performance (we have observed a factor of **up to
      50**).

      To avoid this, restrict the number of multithreading threads by setting
      the corresponding environment variable (e.g. ``$MKL_NUM_THREADS`` for
      Intel MKL):

      .. code-block:: bash

         export MKL_NUM_THREADS=2

      Ideally, :math:`{\mathrm{N\_CORES} \times \mathrm{THREADS}}` should be
      set to the number of available hyperthreading cores on the system with
      the N_CORES as high as possible without running into memory limitations.


The TensErLEED search program is executed using ``mpirun``, with the number
of cores specified by N_CORES. If mpirun / mpiifort are not present on the
computer running |calc|, no parallelization will be performed in the search
(see also :ref:`FORTRAN_COMP<FORTRAN_COMP>`).
Reference and Delta calculations instead use python multiprocessing to execute
multiple TensErLEED processes, so those can run without an MPI compiler.
