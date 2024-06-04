.. include:: /substitutions.rst

.. _ivbeams:

IVBEAMS
=======

The IVBEAMS file lists the beams for which output should be generated.
These are the beams for which the R-factors will be calculated, as well
as those on which the structural optimization will run.

If no IVBEAMS file is found, but the file with experimental |IV| curves
(:ref:`EXPBEAMS file<EXPBEAMS>`) is present, an IVBEAMS file will be
generated from the EXPBEAMS file.
Note that this will only happen once, so later iterations of the
program will use the IVBEAMS file created in the first run, if it is
not deleted manually.

The first line is treated as a header and can contain a comment, all
other lines contain two numbers indicating a beam. Input format can be
integer, float, or fractions as ``a/b``:

..  code-block:: none


   list of beams to output (input form: integer, float or fractions as 'a/b')
   1 0
   0 -1
   1 1
   2 0
   0 2
   0.5 1
   1 1/2
   1/3 1
   1/3 -1

The order of the beams is not important, and if one beam is given twice,
it will be ignored the second time.

Since everything to the right of the first two columns is ignored,
copying beams from the :ref:`BEAMLIST<beamlist>`  file will also work.
For example, the following is a valid IVBEAMS file:

..  code-block:: none


   list of beams to output (input form: integer, float or fractions as 'a/b')
      1.00000   0.00000  1  1          E =     5.9289  NR.   4
      1.00000   1.00000  1  1          E =    11.0329  NR.   6
      1.00000  -1.00000  1  1          E =    11.0329  NR.   9
      0.00000   2.00000  1  1          E =    20.4162  NR.  10
      0.00000  -2.00000  1  1          E =    20.4162  NR.  11
      2.00000   0.00000  1  1          E =    23.7154  NR.  12
      1.00000   2.00000  1  1          E =    26.3450  NR.  14
      1.00000  -2.00000  1  1          E =    26.3450  NR.  17
      2.00000   1.00000  1  1          E =    28.8195  NR.  18
      2.00000  -1.00000  1  1          E =    28.8195  NR.  21
      2.00000   2.00000  1  1          E =    44.1316  NR.  22
      2.00000  -2.00000  1  1          E =    44.1316  NR.  25
      ...

Likewise, round brackets and vbars are ignored even if they appear on
the left, so input of format ``(h k)`` or ``(h | k)`` is also accepted.
