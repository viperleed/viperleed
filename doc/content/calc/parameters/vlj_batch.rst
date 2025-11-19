.. _vlj_batch:

VLJ_BATCH
=========

VLJ_BATCH sets the batch sizes for the parallelized
:ref:`tensor-LEED<tensor_leed>` optimization using the viperleed-jax plugin.
The calculation can be parallelized over multiple energies and/or atoms.
The choice of batch sizes can have a significant impact on the performance, with
larger batch sizes generally leading to better performance, but also higher
memory usage.

**Default:** ``energies -1, atoms -1``


**Syntax:**

::

   ! setting the algorithm to use
   VLJ_BATCH = energies 80, atoms 20 ! parallelize over 80 energies and 20 atoms
   VLJ_BATCH = energies 20           ! parallelize over 20 energies and all atoms

The batch sizes can be set for the energies and atoms separately. The
default value of ``-1`` means that viperleed-jax will parallelize over all
energies or atoms, respectively.

Generally, larger batch sizes will lead to better performance, as calculations
are more efficient when they can be vectorized over larger arrays. However,
this also means that the memory usage will increase, as more data needs to be
kept in memory at once. If you run into memory issues, you can try reducing the
batch sizes.

.. tip::

   In most calculations, atoms constitute the "inner" loop, i.e., it is better
   to keep a larger batch size for the atoms than for the energies. If you need
   to decrease memory usage, try to reduce the batch size for the energies
   first.

   When choosing batch sizes, it is preferable to use equal fractions of the
   total number, e.g. for 208 energy steps, you could use energy batch sizes
   of 26, 52, 104, or 208.