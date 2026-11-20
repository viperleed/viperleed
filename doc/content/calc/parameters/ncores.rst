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

.. note::

      Some :term:`BLAS` and :term:`LAPACK` libraries (e.g., OpenBLAS, or
      Intel MKL, which is used by default with ifort) support use of
      multiple threads. This can conflict with the parallelization used
      by ViPErLEED and **massively** reduce performance (we have observed
      a factor of **up to 50**).

      For this reason, |calc| runs each reference- and delta-amplitude
      calculation with the environment variables ``OPENBLAS_NUM_THREADS``,
      ``OMP_NUM_THREADS``, and ``MKL_NUM_THREADS`` set to 1.

      You can control threading manually (e.g., on a computing cluster)
      by setting ``OMP_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``,
      ``GOTO_NUM_THREADS``, or ``MKL_NUM_THREADS`` in your environment
      to a non-empty value. If you set ``OMP_NUM_THREADS``, which all
      libraries read, |calc| leaves all of these variables as they are.
      A variable of a specific library (``OPENBLAS_NUM_THREADS`` or
      ``GOTO_NUM_THREADS`` for OpenBLAS, ``MKL_NUM_THREADS`` for Intel
      MKL) applies only to that library: |calc| keeps your value, and
      sets the variable of the other library and ``OMP_NUM_THREADS``
      to 1. In both cases, the log file says which variables you set.
      For example, to use two OpenBLAS threads:

      .. code-block:: bash

         export OPENBLAS_NUM_THREADS=2

      If you do so, :math:`{\mathrm{N\_CORES} \times \mathrm{THREADS}}`
      should not exceed the number of available hyperthreading cores on
      the system.

.. versionchanged:: 0.15.0
   |calc| sets the number of BLAS threads for reference- and
   delta-amplitude calculations to 1, unless you set it yourself.
   Previously, the number of threads had to be restricted by hand.


The TensErLEED search program is executed using ``mpirun``, with the number
of cores specified by N_CORES. If mpirun / mpiifort are not present on the
computer running |calc|, no parallelization will be performed in the search
(see also :ref:`FORTRAN_COMP`).
Reference and Delta calculations instead use python multiprocessing to execute
multiple TensErLEED processes, so those can run without an MPI compiler.
