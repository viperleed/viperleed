.. _zip_compression_level:

ZIP_COMPRESSION_LEVEL
=====================

ZIP_COMPRESSION_LEVEL sets the compression level applied to
ZIP archives storing the :ref:`tensor<tensorszip>` and
:ref:`delta-amplitude files<deltaszip>`.


**Default:** 2

**Allowed values:** Integer values 0-9

**Syntax:**

::

   ZIP_COMPRESSION_LEVEL = 5 ! use compression level 5
   COMPRESSION_LEVEL = 5     ! allowed synonym
   COMPRESSION = 5           ! allowed synonym

   ZIP_COMPRESSION_LEVEL = 0 ! do not compress files


The compression level is a trade-off between file size and compression time.
Higher compression levels reduce the size of stored :ref:`tensor<tensorszip>`
and :ref:`delta-amplitude files<deltaszip>` but ViPErLEED will take longer to
create the archives. You may want to tweak this parameter if you have limited
disk space available or if you are working with large unit cells where
compression takes a significant amount of time.

If ZIP_COMPRESSION_LEVEL is set to ``0``, files will be stored in
*uncompressed* archives. Compression time will be almost negligible,
but the files will be much larger.


.. figure:: /_static/zip_compression_levels_plot.png
   :width: 70%
   :align: center

   Benchmark results for a ~250MB :ref:`tensor file<tensorszip>`: Compressed
   file size and compression time vs. compression level. Average over 10
   repetitions.


ViPErLEED uses the Python
`zipfile <https://docs.python.org/3/library/zipfile.html>`__ module for
creating ZIP archives and passes the value set for ZIP_COMPRESSION_LEVEL
via the ``compresslevel`` parameter.
