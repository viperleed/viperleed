.. _tensor_index:

TENSOR_INDEX
============

TENSOR_INDEX defines which Tensors_*.zip (and corresponding Deltas_*.zip) file
from the Tensors (and Deltas) folder should be used. This can be used to start
a delta calculation or a search from an earlier state of the project. Note that
to yield sensible results, several parameters have to be the same as were used
to produce the tensor files. The input files used to produce a given set of
Tensors can be found in the Tensors_*.zip file.

**Default:** uses the highest tensor index found in the Tensors folder.

**Allowed values:** positive integer

**Syntax:**

::

   TENSOR_INDEX = 3

.. note::
  * Setting TENSOR_INDEX is pointless if the run starts with a reference
    calculation, as this would produce an automatically determined new index.
    Similarly, if a reference calculation occurs later during the job, a new
    tensor index would be produced at that point.
  * The given index is taken at face value, and if the required tensor files
    (and for immediate search, delta files) are not found, the job will stop.
  * Note that if your job script deletes the old work folder, and doesn't copy
    all tensor and delta files (e.g., copying only the newest ones), you need
    to either adapt your job script accordingly, or move the required .zip
    files to the work folder yourself, for TENSOR_INDEX to work.
