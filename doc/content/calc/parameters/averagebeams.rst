.. include:: /substitutions.rst

.. _averagebeams:

=============
AVERAGE_BEAMS
=============

AVERAGE_BEAMS can be used to turn averaging of symmetry-equivalent beams
off, or to make it follow a different averaging scheme than would be 
dictated by :ref:`BEAMINCIDENCE`. This can be useful to 
e.g. enforce "normal-incidence averaging" for slightly off-normal 
incidence, or to force calculation of multiple equivalent theoretical 
beams at normal incidence.
However, use of this parameter is **not recommended** for the general 
case.

**Default**: Use averaging as given by :ref:`BEAMINCIDENCE` 

**Syntax**:

::

   AVERAGE_BEAMS = off     ! ignore symmetry equivalence entirely. Alias: 'none', 'false'
   AVERAGE_BEAMS = all     ! pretend that beam incidence is perpendicular (equivalent to setting '0 0'). Alias: 'perpendicular'
   AVERAGE_BEAMS = 0 0     ! same as above
   AVERAGE_BEAMS = THETA 0.3, PHI 10.1
   AVERAGE_BEAMS = 0.3 10.1

**Acceptable values**: ``off``, ``none``, ``false`` / ``all``, 
``perpendicular`` / beam incidence specified by ``theta`` and 
``phi``, following the same syntax as specified in 
:ref:`BEAMINCIDENCE` 

See page on :ref:`BEAMINCIDENCE`  for details on the 
definition of ``theta`` and ``phi``.

If AVERAGE_BEAMS is disabled entirely (i.e. ``off``), then 
:ref:`IVBEAMS`  and :ref:`EXPBEAMS.csv<EXPBEAMS>` 
are taken at face value. Only beams with the same label will 
be recognized as equivalent and compared in the |R-factor| calculation. 
If :ref:`IVBEAMS`  is not present, the generated 
:ref:`IVBEAMS`  file will contain the same beams as the 
:ref:`EXPBEAMS.csv<EXPBEAMS>`  file. Note that this affects not only 
spot equivalences due to mirror or rotation symmetry, but also disables 
equivalence based on multiple rotational domains that may be present on 
the surface. Averaging over rotational domains, but not over 
mirror/rotation symmetries can be achieved by defining an arbitrary, 
non-zero ``theta``, and a ``phi`` that is not parallel to any mirror 
plane.
