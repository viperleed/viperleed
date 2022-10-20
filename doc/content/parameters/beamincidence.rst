.. _beamincidence:

BEAM_INCIDENCE
==============

BEAM_INCIDENCE defines the incidence angles (in degrees) of the electron beam on the surface. The polar angle ``theta`` is measured from the surface normal, the azimuthal angle ``phi`` is positive counterclockwise with ``phi=0`` **corresponding to the positive x axis. TODO: Perhaps more intuitive if we take as relative to a1??**, as defined in the :ref:`POSCAR<POSCAR>`  file.

**Default**: BEAM_INCIDENCE = THETA 0, PHI 0

**Syntax**:

::

   BEAM_INCIDENCE = THETA 0.3, PHI 10.1
   BEAM_INCIDENCE = 0.3 10.1

**Acceptable values**: -90 ≤ ``theta`` ≤ 90, 0 ≤ ``phi`` < 360. All numbers are considered floats. Negative values for ``theta`` will internally be corrected to positive by adding or subtracting 180° **should be 360°** from ``phi``.

``theta`` and ``phi`` represent tilt and azimuthal angles, respectively, in degrees. Notice that if the flags THETA and PHI are not specified, only the first two floats are considered: the first is taken as the tilt angle theta, the second as the azimuth phi. If you want to use ranges, the flags THETA and PHI are mandatory.

**Notes**:

-  In general, unless the experiment was performed at large off-normal incidence (>2°), keeping the default value should lead to the correct optimized geometry, but the R factor will be worse than for an experiment at normal incidence. In case the incidence angle is significantly off, one needs to measure a simple system (clean, unreconstructed metal) with the same experimental settings, to determine first the correct angle of incidence to be used for more complicated situations.
-  Beam incidence optimization is commonly one of the last refinement steps of the LEED fit.
-  Even for normal incidence, during the last polishing one can even use a purposely off-normal BEAM_INCIDENCE (``theta`` different from zero) and average the resulting almost-equivalent beams in order to account for the fact that the electron beam has a finite aperture angle. (For averaging, see the :ref:`AVERAGE_BEAMS<AVERAGEBEAMS>`  parameter). This commonly leads to ``theta`` values in the order of 0.3 -- 0.5°. This option makes sense only if the R factor is very low and the surface has sufficiently high symmetry (at least threefold rotation symmetry or mirror/glide planes in two directions), so that averaging simultates incoming beams from several azimuthal directions. The azimuth ``phi`` of the incoming beam should be chosen such that it is midway between the azimuth values of two beams that are symmetry-equivalent at normal incidence.

**TODO**: have a decent default for ``phi``?. One would like to run a refcalc such that the beam impinges midway between mirror/glide planes.

**TODO**: it might be worth implementing ranges at a later stage. However, if there are ranges, this requires to (i) run one refcalc at each (theta,phi), (ii) run delta+search -> get best R, and compare with other (theta,phi).

Syntax would then also allow:

::

   BEAM_INCIDENCE = THETA thetamin thetamax thetastep PHI phimin phimax phistep

0≤\ ``theta``\ ≤90, 0≤\ ``phi``\ <360; 0≤\ ``thetamin``\ <``thetamax``\ ≤90, 0≤\ ``phimin``\ <``phimax``\ <360, 0<``thetastep``\ ≤\ ``thetamax`` – ``thetamin``, 0<``phistep``\ ≤\ ``phimax`` – ``phimin``. **AMI**: Isn't this redundant now that we have FD calculations? In my limited experience, this works reasonably well, even without Delta+Search for each angle combination.
