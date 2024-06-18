.. _deltaszip:

Deltas
======

Delta directories Deltas\_\ *XYZ* and archives Deltas\_\ *XYZ*.zip use the
same Tensor index *XYZ* as the :ref:`Tensor files<Tensorszip>`. In fact,
**the index XYZ of a Deltas\_\ XYZ.zip archive indicates that it was **
**created based on the Tensors with the same index, Tensors\_\ XYZ.zip**.
Delta files are named using the pattern 'DEL\_\ *N*\ \_\ *El*\ \_\ *k*',
where *N* is the number of the atom (same numbering as in
:ref:`POSCAR`), *El* is the chemical element being addressed by
this Delta file (may differ from POSCAR element through
:ref:`ELSPLIT`), and *k* is a running index
to distinguish multiple Delta files generated for the same
atom and element (see below).

When a delta-amplitudes calculation starts and a Deltas\_\ *XYZ* folder or
Deltas\_\ *XYZ*.zip archive already exist for the given Tensors\_\ *XYZ*,
ViPErLEED will read the previous Delta files to check whether they match
the requested parameter variations. New Delta files will only be calculated
if no file with the requested set of variations exists yet. If this is the
case, the old Delta file will be kept nevertheless, and the index *k* will
simply increment by 1 to distinguish the files. Therefore, several Delta
files can accumulate for each atom if consecutive searches are performed,
especially when searches are looped.

As with the Tensors, at the end of the calculation, the directory
Deltas\_\ *XYZ* is packed into an archive Deltas\_\ *XYZ*.zip, which
is moved into the 'Deltas' directory.
