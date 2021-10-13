C  Tensor LEED subroutines for reference & delta amplitude calculation
C  v1.72, LH,MR,AMI 13.10.2021
C  for use with ref-calc.f v1.72, delta.f v1.7+
C
C  as described in
C
C  V. Blum, K. Heinz, submitted to Comp. Phys. Comm. (2000).
C
C  Please cite this reference when using this program.
C
C  Welcome to the Erlangen Tensor LEED program package. Permission is granted
C  to use these programs freely, and to modify them as you may find appropriate.
C  However, we cannot take any responsibility for possible bugs or other
C  limitations you may encounter. If you find a bug or worthwhile improvement,
C  please notify us at
C
C                       tleed@fkp.physik.uni-erlangen.de
C
C  so we can make the updated code available to everyone. Likewise, in re-
C  distributing the code please refer back to us first to ensure the latest version
C  of the package is passed on.
C
C**********************************************************************************
C
C  Please read the comment in ref-calc.f, v1.72 .
C
C  Most subroutines are derived from those described in
C
C  Van Hove / Tong, "Surface Crystallography by LEED" (Springer, Berlin 1979),
C
C  with modifications for Tensor LEED where appropriate.
C  Initial modifications were made by P.J. Rous for Tensor LEED,
C  further ones by N. Bickel, W. Oed, R. Doell, U. Loeffler, M. Kottcke and others.
C
C  For this version, subroutines were modified resulting in some
C  loss of generality compared to original VHT subroutines. E.g., the MPERTI
C  subroutine is no longer compatible with the current version.
C
C  One word of advice: If you find a bug in these subroutines, do not ignore
C  it or just patch it somehow. Look into the code, find the bug and fix it
C  once and for all (and let others using the code know). It'll help a lot
C  of people avoid hitting the same trap over and over again.
C  Another one: Please write your comments in English!!
C
C  Subroutines are included in alphabetical order

!-------------------------------------------------------------------
!  SUBROUTINE GET_TOPLAY_PROPAGATORS computes propagators for
!  the complex amplitudes of each beam:
!    WK(i, 1)
!      beam i to/from 'vacuum' from/to topmost surface layer (i.e.,
!      where the solid starts).
!    WK(i, 2)/WK(i, 3)
!      beam i going towards (2) and coming from (3) the underlying stack
!
!    Parameters
!    ----------
!    WK : complex, OUTPUT
!        The propagators are stored here.
!    PQ : real
!        Lateral momentum for all the beams in the beam list in
!        atomic units (i.e., the reciprocal-space vectors). Only
!        the first N are to be considered.
!    N : integer
!        Total number of beams included in the calculation at the
!        current energy
!    INTL_VEC : real
!        Interlayer vector (z, x, y) between top layer and the
!        underlying stack
!
!    From common
!    -----------
!    E : real
!        Current energy in the solid (including real part of inner
!        potential)
!    AK2, AK3 : real
!        Lateral momentum of incident beam in vacuum
!    VPI : real, unused
!        Imaginary part of the inner potential
!    VPIO, VPIS : real
!        Imaginary part of the inner potential in the top layer
!        and in the rest of the stack
!    VO : real
!        Offset of real part of inner potential between top layer
!        and rest of the stack (VO == E_rest - E_top)
!    FR : real
!        Fraction of the phase shift due to perpendicular propagation
!        that should use the inner potential of the 'rest of the stack'.
!        1 - FR will be the contribution of the top layer. When using
!        VPIS == VPIO and VO == 0 this has no effect.
!    ASE : real
!        Onset 'height' of inner potential, i.e., where the solid starts
!        with respect to the geometrical position of the topmost layer.
!
!    Author: Michele Riva
!    Added: 2021-08-26

      SUBROUTINE  GET_TOPLAY_PROPAGATORS (WK, PQ, N, INTL_VEC)

      INTEGER I, N
      COMPLEX WK, XX, YY, IU
      DIMENSION WK(N, 3)
      COMPLEX SQRT, EXP
      REAL AK, BK2, BK3, PQ, INTL_VEC, X
      DIMENSION INTL_VEC(3), PQ(2, N)

!    ! Ideally we do not want to use VPIS/VPIO/VO as they give weird
!    ! results according to Lutz, but rather VPI. This would also make
!    ! XX and YY entirely equal, and would simplify the code a bit.
!    ! For now keep it this way for consistency with the code that was
!    ! in ADREF2T. This may perhaps turn out to become handy for 'rough'
!    ! surfaces, although the fact that a difference is assumed only
!    ! between the very top layer and what's underneath may be a bit
!    ! of a poor approximation.
      COMMON  E, AK2, AK3    !, VPI
      COMMON  /ADS/ FR, ASE, VPIS, VPIO, VO, VV  ! VV necessary (although unused) as this is a named COMMON block from main

      IU = CMPLX(0.0,1.0)
      AK = 2.0 * E

      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)

!       YY is k_perp in topmost layer
        YY = CMPLX(AK - 2.0 * VO - BK2 * BK2 - BK3 * BK3,
     +             -2.0 * VPIO + 0.000001)
        YY = SQRT(YY)
!       XX is k_perp in the rest of the stack
        XX = CMPLX(AK - BK2 * BK2 - BK3 * BK3,
     +             -2.0 * VPIS + 0.000001)
        XX = SQRT(XX)

!       X is the phase shift due to in-plane propagation between top layer
!       and stack underneath
        X = BK2 * INTL_VEC(2) + BK3 * INTL_VEC(3)

        WK(I, 1) = EXP(IU * (YY * ASE))
        WK(I, 2) = EXP(IU * (YY * (1.0 - FR) * INTL_VEC(1)
     +                       + XX * FR * INTL_VEC(1) + X))
        WK(I, 3) = EXP(IU * (YY * (1.0 - FR) * INTL_VEC(1)
     +                       + XX * FR * INTL_VEC(1) - X))
      ENDDO
      RETURN
      END

!-------------------------------------------------------------------
!  SUBROUTINE GET_TOPLAY_MULTSCATT computes the LU factorization of
!  the complex multiple-scattering matrix between the top layer and
!  the rest of the stack.
!  The LU part of the matrix is returned into MULT_SCATT, while PIVOT
!  stores the permutation pivots of the LU decomposition.
!  The return values should be used in a subsequent call to CGETRS.
!  The call to CGETRS is then equivalent to computing
!  AMP1 = MULT_SCATT*AMP2,
!  with MULT_SCATT(i, j) corresponding to the multiple scattering of
!  beam j (exiting the top layer) into beam i, which originates from
!  an infinite number of reflections between the rest of the stack
!  (RREST_ABOVE) and the top layer itself (RTOP_BELOW) through all
!  beams. Beam i is then right below the top layer.
!
!  The matrix is computed without approximations relying on the fact
!  that (see Neumann Series)
!
!          (Id - A)^(-1) = SUM(A**k) .
!
!  Thus, SUM(A**k) can be calculated by inverting Id - A.
!
!  Author: Michele Riva
!  Added: 2021-08-27
!
!  Parameters
!  ----------
!  RTOP_BELOW : complex
!     Matrix of reflection for beams reaching the top layer from below
!  RREST_ABOVE : complex
!     Matrix of reflection for beams reaching the rest of the stack
!     from above
!  WK : complex
!     Propagators between top layer and rest, as calculated via
!     GET_TOPLAY_PROPAGATORS
!  MULT_SCATT : complex
!     The multiple scattering matrix output
!  N : integer
!     Total number of beams at the current energy

      SUBROUTINE  GET_TOPLAY_MULTSCATT(RTOP_BELOW, RREST_ABOVE, WK,
     +                                 MULT_SCATT, INV_PIVOT, N)

      INTEGER  N
      COMPLEX  RTOP_BELOW, RREST_ABOVE, WK, MULT_SCATT, RU, CZ
      DIMENSION RTOP_BELOW(N, N), RREST_ABOVE(N, N)
      DIMENSION WK(N, 3)
      DIMENSION MULT_SCATT(N, N)
      INTEGER   INV_PIVOT, INV_INFO   ! For matrix inversion
      DIMENSION INV_PIVOT(N)

      RU = CMPLX(1.0, 0.0)            ! Complex one
      CZ = CMPLX(0.0, 0.0)            ! Complex zero

!     Prepare the matrix (Id - SINGLE_SCATTERING) that will later be inverted

      DO J = 1, N        ! outer loop on second index
        DO I = 1, N
          MULT_SCATT(I, J) = CZ
          DO K = 1, N
            MULT_SCATT(I, J) =
     +          MULT_SCATT(I, J)
     +          - (WK(J, 2)                ! Propagate j toward rest
     +             * RREST_ABOVE(K, J)     ! Reflect j at rest into k
     +             * WK(K, 3)              ! Propagate k back toward top
     +             * RTOP_BELOW(I, K))     ! And reflect k at top into i
          ENDDO
        ENDDO
        MULT_SCATT(J, J) = MULT_SCATT(J, J) + RU
      ENDDO

!     Gaussian-elimination step of matrix inversion (LU decomposition)
      CALL CGETRF(N, N, MULT_SCATT, N, INV_PIVOT, INV_INFO)

      RETURN
      END

!-------------------------------------------------------------------
!  SUBROUTINE ADREF2T computes the complex amplitudes for beams
!  as a result of stacking the top layer above the rest of the slab.
!  The rest of the slab already includes the substrate bulk and all
!  surface layers except the topmost one.
!  Amplitudes are calculated for ONE INCIDENT BEAM from (1) reflection
!  and transmission matrices of the top layer for beams coming from above
!  and below, (2) reflection matrix of the rest, (3) multiple-scattering
!  matrix between top and rest, and (4) propagators.
!  Only the amplitudes necessary to output Tensors (AMP0, AMP1, AMP2) are
!  calculated for a generic beam. In addition, the outgoing wavefield XI
!  is calculated when beam (0, 0) is the current beam.
!
!  2021-08-27 Michele Riva: Reimplement subroutine and CHANGE ITS
!  SIGNATURE. The previous signature was:
!     (RA, TA, RAM, TAM, RB, ST, WK, U, INT, XI, PQ, N, EMACH, S,
!      AMP0, AMP1, AMP2, PSQ1, PSQ2)
! and was called in ref-calc with
!     (RAMP, TAPP, RAPM, TAMM, REABOV, S1, WK, XS, JNT, XI, PQ, NT,
!      EMACH, ASV, APLUS, AMINUS, ADOWN, PSQ1, PSQ2)
!
!  Parameters
!  -----------
!  RA, TA : complex
!      Reflection and transmission matrices for top layer and beams
!      coming from above. RA(i, j) is reflection of incident beam j
!      (from above) into beam i (moving up). TA(i, j) is transmission
!      of incident beam j (from above) into beam i (going down).
!  RAM, TAM : complex
!      Reflection and transmission matrices for top layer and beams
!      coming from below. RAM(i, j) is reflection of incident beam j
!      (from below) into beam i (going down). TAM(i, j) is transmission
!      of incident beam j (from below) into beam i (going up). RAM is
!      not used anymore (its effect is already included in the multiple
!      scattering matrix ST), but is kept in the signature for symmetry.
!  RB : complex
!      Reflection matrix for the rest of the slab for beams impinging
!      from above (i.e., from the region between top layer and rest)
!  ST : complex
!      Multiple-scattering matrix between top layer and the rest.
!      ST(i, j) is scattering of beam j -- exiting the top layer toward
!      the rest of the slab, right below the top layer -- into beam i,
!      which is also exiting the top layer propagating toward the rest
!      of the slab. ST should e calculated (once per each 'unknown' beam)
!      with GET_TOPLAY_MULTSCATT
!  WK : complex
!      Propagators for beams. WK(i, 1) propagates beam i from/to the
!      vacuum to/from the geometrical 'height' of the topmost layer.
!      WK(i, 2)/WK(i, 3) propagate beam i down/up from/to the top
!      layer to/from the rest of the slab. WK should be calculated
!      (once per each 'unknown' beam) with GET_TOPLAY_PROPAGATORS
!  XI : complex, OUTPUT
!      The wavefield in vacuum when the (0, 0) beam impinges on the
!      whole stack from vacuum will be stored in here. While the amplitudes
!      AMP0, AMP1, and AMP2 are always computed for all the requested
!      outgoing beams (needed for Tensors), XI is computed only once, i.e.,
!      when calculating incidence of the (0, 0) beam (IIN == 1).
!  N : integer
!      Total number of beams in the current beam set (i.e., at the
!      current energy)
!  AMP0, AMP1, AMP2 : complex, OUTPUT
!      Wavefield amplitudes as per following scheme:
!      ---- VACUUM ---------------------------------------------------
!                         |                           ^ XI
!                         v AMP0                      |
!      ---- TOP LAYER ------------------------------------------------
!                                  |           ^ AMP1
!                                  v AMP2      |
!      ---------------------------------------------------------------
!      //// REST OF SLAB /////////////////////////////////////////////
!  IIN : integer
!      Index of current beam in the beam set

      SUBROUTINE  ADREF2T (RA, TA, RAM, TAM, RB, ST, ST_PIVOT, WK, XI,
     +                     N, AMP0, AMP1, AMP2, IIN)

!   Alex: moved ST_PIVOT from complex to int -> should be?
      INTEGER    I, J, N, IIN, ST_PIVOT
      COMPLEX    RA, TA, RAM, TAM, RB, ST, WK, XI, CZ
      DIMENSION  RA(N, N), TA(N, N), RAM(N, N), TAM(N,N), RB(N, N)
      DIMENSION  ST(N, N), ST_PIVOT(N), XI(N)
      DIMENSION  WK(N, 3)
      COMPLEX    AMP0(N), AMP1(N), AMP2(N)
      COMPLEX    temp

      CZ = CMPLX(0.0, 0.0)

!     2021-06-27 Michele Riva: clean up and reimplement.
!     Now uses the propagators (WK) and the multiple-scattering matrix
!     computed (only once for each 'unknown' beam) via GET_TOPLAY_PROPAGATORS
!     and GET_TOPLAY_MULTSCATT.
!     AMP2(I) is the complex amplitude of beam I, propagating towards the
!     'rest of the slab' and right above it, when beam IIN comes in with unit
!     amplitude from the vacuum.

!     Construct the amplitudes of the incoming beams, right below the top
!     layer: beam IIN produces contributions in all other beams. In the
!     meanwhile, also zero the other amplitudes.
      temp = WK(IIN, 1)             ! Current beam from vacuum
      DO I = 1, N
        AMP2(I) = temp              ! Current beam from vacuum
     +            * TA(I, IIN)      ! Transmitted into beam I at top layer
        AMP1(I) = CZ
        AMP0(N) = CZ
      ENDDO

!     Now each beam is multiple-scattered between the top and rest. AMP2(I)
!     contains then the amplitude of beam I right below the top layer.
      call CGETRS('N', N, 1, ST, N, ST_PIVOT, AMP2, N, INFO)

!     The quantity we're interested in is the amplitude right above the rest,
!     though --> Propagate it down.
      DO I = 1, N
        AMP2(I) = AMP2(I) * WK(I, 2)
      ENDDO

!     AMP1(I) is the complex amplitude of beam I, propagating towards the
!     top layer from below and right below it, after reflection with the
!     rest of the slab of AMP2.
      DO J = 1, N                     ! Outer loop on second index
        temp = AMP2(J)                ! Beam J right above the 'rest'
        DO I = 1, N
          AMP1(I) = AMP1(I)
     +              + (temp           ! Beam J right above the 'rest'
     +                 * RB(I, J)     ! Reflected at the rest into beam I
     +                 * WK(I, 3))    ! Propagated right below the top layer
        ENDDO
      ENDDO

!     AMP0 contains the actual complex amplitudes of all beams right above
!     the top layer and coming from the vacuum side, when beam IIN comes
!     from the vacuum. All amplitudes are zero (zeroed above), except for
!     the one of IIN that is just propagated from vacuum to solid.

      AMP0(IIN) = WK(IIN, 1)

!     Now, only for IIN == 1, i.e., when calculating the (0, 0) beam, we also
!     compute the total wave-field exiting the solid for all beams, as this
!     is what is used to eventually calculate the IV curves (with a call to
!     RINT). This is stored into XI. We do the calculation only once, as we
!     only need the actual intensities when the incident beam is the (0, 0)
!     one. However, the calculation of the amplitudes done above is necessary
!     for computing the Tensors, and must be done for each calculated beam.
      IF (IIN.EQ.1) THEN
        DO I = 1, N
          ! First get the straight reflection of the (0, 0) beam at
          ! the top layer, right above it
          XI(I) = WK(IIN, 1)     ! (0, 0) beam from vacuum to top layer
     +            * RA(I, IIN)   ! Reflected at top layer into beam I

          ! Then prepare the contributions of the scattering from below
          temp = CZ
          DO J = 1, N
            temp = temp
     +             + (AMP1(J)      ! Beam J right below the top layer
     +                * TAM(I, J)) ! Transmitted into beam I while going
!                                    through the top layer from below
          ENDDO
          XI(I) = (temp + XI(I)) * WK(I, 1)  ! Finally sum up, and propagate to vacuum
        ENDDO
      ENDIF
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE BEAMS selects those beams from the input list that are
!  needed at the current energy, based on the parameter TST which limits
!  the decay of plane waves from one layer to the next (the interlayer
!  spacing must already be incorporated in TST)
!
!   KNBS = No. of input beam sets
!   KNB = No. of beams in each input beam set
!   SPQ, SPQF = List of input beams (reciprocal lattice vectors G)
!   KNT = Total no. of input beams
!   NPU  Indicates for which beams I(V) curves should be given at the end
!        of the refcalc ('punched')
!   NPUN = No. of beams for which intensities are to be output
!   AK2, AK3 = Components of incident wave vector parallel to surface
!              must be set to AK21, AK31 in main
!   E = Current energy
!   TST = Criterion for beam selection: beams whose out-of-plane
!         wave vector (i.e., perpendicular to the surface) is
!         imaginary with imaginary part larger than sqrt(TST) are
!         considered evanescent. The interlayer spacing must already
!         be incorporated in TS
!   NB, PQ, PQF, NPUC, MPU, NT
!         same as KNB, SPQ, SPQF, NPU, NPUN, KNT after selection of
!         the non-evanescent beams
!   NP = Largest no. of beams selected in any of the beam sets
!
!   Edit 2021-09-10: remove symmetry codes
      SUBROUTINE BEAMS(KNBS, KNB, SPQ, SPQF, KNT, NPU, NPUN, AK2, AK3,
     &                 E, TST, NB, PQ, PQF, NPUC, MPU, NT, NP)
      INTEGER KNBS, KNBJ, NT, NP, MPU, J, N, K, KK, MPU1, NPUN, I
      INTEGER KNB, NB, NPU, NPUC
      REAL E, AK2, AK3, SPQ, TST, PQ, PQF, SPQF
      DIMENSION KNB(KNBS), SPQ(2, KNT), SPQF(2, KNT), NPU(NPUN),
     &          NB(KNBS), PQ(2, KNT), PQF(2, KNT), NPUC(NPUN)

      KNBJ = 0
      NT = 0
      NP = 0
      MPU = 0
      DO J = 1, KNBS
        N = KNB(J)
        NB(J) = 0
        DO K = 1, N
          KK = K + KNBJ
          MPU1 = MPU + 1
          IF ((MPU1.LE.NPUN).AND.(KK.EQ.NPU(MPU1))) THEN
            MPU = MPU + 1
            NPUC(MPU) = 0
          ENDIF
          IF ((2.0 * E - (AK2 + SPQ(1, KK))**2 - (AK3 + SPQ(2, KK))**2)
     &        .GE.(-TST)) THEN  ! beam propagates
            NB(J) = NB(J) + 1
            NT = NT + 1
            IF ((MPU1.LE.NPUN).AND.(KK.EQ.NPU(MPU1))) THEN
              NPUC(MPU) = NT
            ENDIF
            DO I = 1, 2
              PQ(I, NT) = SPQ(I, KK)
              PQF(I, NT) = SPQF(I, KK)
            ENDDO
          ENDIF
        ENDDO
        KNBJ = KNBJ + KNB(J)
        NP = MAX0(NP, NB(J))
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
C  Subroutine Bessel generates the spherical Bessel functions by direct
C  evaluation of Pendry's expression (A.27). Since the argument is complex,
C  downward recursion is said to be more stable. For the range of interest
C  here, I have tested downward recursion vs. the current subroutine, and
C  find no deviations whatsoever. This one is more straightforward though
C  less elegant. VB.

      subroutine BESSEL(BJ, Z, N1)

C  note N1 is highest angular momentum for which BJ are calculated - here
C  2 * LMAX + 1. For low LMAX, the sums required when calculating temperature-
C  dependent phase shifts do not converge fast enough. For high LMAX, 2*LMAX+1
C  is a bit of an overkill, I guess ...

      complex Z
      complex*16 BJ(N1)

      complex*16 PRE,TERM,SUM,ZSQ
      integer K,L

C  begin computation

      ZSQ = Z*Z/2.

C  L is indeed the angular momentum in question. PRE is the L-dependent
C  pre-factor.

      PRE = 1.

      DO L = 0,N1-1,1

C  K is index of series expansion.

        TERM = 1.
        SUM = 1.
        K = 1

C  repeat the loop defined by 100 ...

 100    CONTINUE

          TERM = - TERM * (ZSQ/FLOAT(K*(2*L+2*K+1)))
          SUM = SUM + TERM

          K = K+1

C  ... until TERM really small!
        IF (CABS(CMPLX(TERM)).gt.1.E-16) GO TO 100

C  evaluate j_L(Z), i.e. BJ(L+1)

        BJ(L+1) = PRE * SUM

        PRE = PRE * Z / (2.*FLOAT(L+1)+1.0)
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
COMMENT FUNCTION BLM PROVIDES THE INTEGRAL OF THE PRODUCT
C       OF THREE SPHERICAL HARMONICS, EACH OF WHICH CAN BE
C       EXPRESSED AS A PREFACTOR TIMES A LEGENDRE FUNCTION.
C       THE THREE PREFACTORS ARE LUMPED TOGETHER AS FACTOR
C       'C'; AND THE INTEGRAL OF THE THREE LEGENDRE FUNCTIONS
C       FOLLOWS GAUNT'S SUMMATION SCHEME SET OUT BY SLATER
C       ATOMIC STRUCTURE, VOL1, 309,310
C  AUTHOR  PENDRY
      FUNCTION  BLM (L1, M1, L2, M2, L3, M3, LMAX, NFAC, FAC)

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL BLM
      DOUBLE PRECISION FAC
      DIMENSION FAC(NFAC)
      FLOAT(I)=DFLOAT(I)
   40 FORMAT(28H INVALID ARGUMENTS FOR BLM  ,6(I3,1H,))
      PI = 3.14159265
      IF (M1+M2+M3)  420, 350, 420
  350 IF (L1-LMAX-LMAX)  360, 360, 530
  360 IF (L2-LMAX)  370, 370, 530
  370 IF (L3-LMAX)  380, 380, 530
  380 IF (L1-IABS(M1))  530, 390, 390
  390 IF (L2-IABS(M2))  530, 400, 400
  400 IF (L3-IABS(M3))  530, 410, 410
  410 IF (MOD(L1+L2+L3,2))  420, 430, 420
  420 BLM = 0.0
      RETURN
  430 NL1 = L1
      NL2 = L2
      NL3 = L3
      NM1 = IABS(M1)
      NM2 = IABS(M2)
      NM3 = IABS(M3)
      IC = (NM1 + NM2 + NM3)/2
      IF (MAX0(NM1,NM2,NM3)-NM1)  470, 470, 440
  440 IF (MAX0(NM2,NM3)-NM2)  450, 450, 460
  450 NL1 = L2
      NL2 = L1
      NM1 = NM2
      NM2 = IABS(M1)
      GO TO 470
  460 NL1 = L3
      NL3 = L1
      NM1 = NM3
      NM3 = IABS(M1)
  470 IF (NL2-NL3)  480, 490, 490
  480 NTEMP = NL2
      NL2 = NL3
      NL3 = NTEMP
      NTEMP = NM2
      NM2 = NM3
      NM3 = NTEMP
  490 IF (NL3-IABS(NL2-NL1))  500, 510, 510
  500 BLM = 0.0
      RETURN
C
C      CALCULATION OF FACTOR 'A'
C
  510 IS = (NL1 + NL2 + NL3)/2
      IA1 = IS - NL2 - NM3
      IA2 = NL2 + NM2
      IA3 = NL2 - NM2
      IA4 = NL3 + NM3
      IA5 = NL1 + NL2 - NL3
      IA6 = IS - NL1
      IA7 = IS - NL2
      IA8 = IS - NL3
      IA9 = NL1 + NL2 + NL3 + 1
      A=((-1.0)**IA1)/FAC(IA3+1)*FAC(IA2+1)/FAC(IA6+1)*FAC(IA4+1)
      A=A/FAC(IA7+1)*FAC(IA5+1)/FAC(IA8+1)*FAC(IS+1)/FAC(IA9+1)
C
C      CALCULATION OF SUM 'B'
C
      IB1 = NL1 + NM1
      IB2 = NL2 + NL3 - NM1
      IB3 = NL1 - NM1
      IB4 = NL2 - NL3 + NM1
      IB5 = NL3 - NM3
      IT1 = MAX0(0, - IB4) + 1
      IT2 = MIN0(IB2,IB3,IB5) + 1
      B = 0.
      SIGN = ( - 1.0)**(IT1)
      IB1 = IB1 + IT1 - 2
      IB2 = IB2 - IT1 + 2
      IB3 = IB3 - IT1 + 2
      IB4 = IB4 + IT1 - 2
      IB5 = IB5 - IT1 + 2
      DO 520 IT = IT1, IT2
      SIGN =  - SIGN
      IB1 = IB1 + 1
      IB2 = IB2 - 1
      IB3 = IB3 - 1
      IB4 = IB4 + 1
      IB5 = IB5 - 1
      BN=SIGN/FAC(IT)*FAC(IB1+1)/FAC(IB3+1)*FAC(IB2+1)
      BN=BN/FAC(IB4+1)/FAC(IB5+1)
  520 B=B+BN
C
C      CALCULATION OF FACTOR 'C'
C
      IC1 = NL1 - NM1
      IC2 = NL1 + NM1
      IC3 = NL2 - NM2
      IC4 = NL2 + NM2
      IC5 = NL3 - NM3
      IC6 = NL3 + NM3
      CN = FLOAT((2 * NL1 + 1) * (2 * NL2 + 1) * (2 * NL3 + 1))/PI
      C=CN/FAC(IC2+1)*FAC(IC1+1)/FAC(IC4+1)*FAC(IC3+1)
      C=C/FAC(IC6+1)*FAC(IC5+1)
      C = (SQRT(C))/2.
C
C
      BLM = (( - 1.0)**IC) * A * B * C
      RETURN

530   WRITE(6,40)L1,M1,L2,M2,L3,M3
      STOP

      END

C-----------------------------------------------------------------------
C FUNCTION CA PERFORMS SAME PURPOSE FOR CAAA AS BLM DOES FOR CELMG
      FUNCTION  CA (L1, MA1, L2, M2, L3, MA3)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      REAL CA,FACT
      DOUBLE PRECISION FACT, PREF, PREFA, PREFB, PREFT
      M1 =  - MA1
      M3 =  - MA3
      IF ((IABS(M1) .GT. L1) .OR. (IABS(M2) .GT. L2) .OR. (IABS(M3) .GT.
     1 L3))  GO TO 1570
      IF ((M1+M2) .NE. (-M3))  GO TO 1570
      IF ( .NOT. ((L3 .GE. IABS(L1-L2)) .AND. (L3 .LE. (L1+L2))))  GO
     1TO 1570
      IF (MOD(L1+L2+L3,2) .NE. 0)  GO TO 1570
      IF ((L1 .EQ. 0) .OR. (L2 .EQ. 0) .OR. (L3 .EQ. 0))  GO TO 1580
      DELT = FACT(L1 + L2 - L3) * FACT(L2 + L3 - L1) * FACT(L3 + L1 -
     1L2)/FACT(L1 + L2 + L3 + 1) * (( - 1)**(L1 - L2 - M3 + ((L1 + L2 +
     2L3)/2)))
      PREF = FACT(L1 + M1) * FACT(L1 - M1) * FACT(L2 + M2) * FACT(L2 -
     1M2) * FACT(L3 - M3) * FACT(L3 + M3)
      PREFA = SQRT(PREF) * (10.D0**(L1 + L2 + L3))
      PREFB = FACT((L1 + L2 + L3)/2)/(FACT((L1 + L2 - L3)/2) * FACT((L2
     1+ L3 - L1)/2) * FACT((L3 + L1 - L2)/2)) * (0.1)
      PREFT = DELT * PREFA * PREFB
      IGT = MIN0(L2 + M2,L1 - M1,L1 + L2 - L3)
      ILST = MAX0(L1 + M2 - L3, - M1 + L2 - L3,0)
      NUM = IGT - ILST + 1
      SUM = 0.0
      DO 1560 IT = 1, NUM
      ITA = ILST + IT - 1
      SUM = SUM + (( - 1.0)**ITA)/(FACT(ITA) * FACT(L3 - L2 + ITA + M1)
     1* FACT(L3 - L1 + ITA - M2) * FACT(L1 + L2 - L3 - ITA) * FACT(L1 -
     2ITA - M1) * FACT(L2 - ITA + M2)) * (10.D0**( - L1 - L2 - L3))
 1560 CONTINUE
      YINT = ((2.0 * FLOAT(L1) + 1.0) * (2.0 * FLOAT(L2) + 1.0) * (2.0 *
     1 FLOAT(L3) + 1.0)/(12.56637))**(0.5) * PREFT * SUM
      GO TO 1590
 1570 YINT = 0.0
      GO TO 1590
 1580 YINT = ((2.0 * FLOAT(L1) + 1.0) * (2.0 * FLOAT(L2) + 1.0) * (2.0 *
     1 FLOAT(L3) + 1.0)/(12.56637))**(0.5) * ( - 1.0)**((IABS(M1) +
     2IABS(M2) + IABS(M3))/2)/(FLOAT(IABS(L1) + IABS(L2) + IABS(L3)) +
     31.0)
 1590 CA = (( - 1.0)**(MA1 + MA3)) * YINT
      RETURN
      END
C---------------------------------------------------------------------
C  SUBROUTINE CAAA COMPUTES CLEBSCH-GORDON COEFFICIENTS FOR USE BY
C  SUBROUTINE GHD IN THE SAME ORDER AS GHD USES THEM
      SUBROUTINE  CAAA (CAA,NCAA,LMMAX)
      DIMENSION  CAA(NCAA)
      II = 1
      DO 60 I = 1, LMMAX
      L1 = INT(SQRT(FLOAT(I - 1)+0.00001))                              261179
      M1 = I - L1 - L1 * L1 - 1
      DO 50 J = 1, LMMAX
      L2 = INT(SQRT(FLOAT(J - 1)+0.00001))                              261179
      M2 = J - L2 - L2 * L2 - 1
      M3 = M2 - M1
      IL = IABS(L1 - L2)
      IM = IABS(M3)
      LMIN = MAX0(IL,IM + MOD(IL + IM,2))
      LMAX = L1 + L2
      LMIN = LMIN + 1
      LMAX = LMAX + 1
      IF (I-J) 20,10,20
   10 IF (M1) 30,30,50
   20 IF (L1.LT.L2.OR.(L1.EQ.L2.AND.(IABS(M1).LT.IABS(M2))))
     1GO TO 50
   30 DO 40 ILA = LMIN, LMAX, 2
      LA = ILA - 1
      CAA(II) = CA(L1,M1,L2,M2,LA,M3)
   40 II = II + 1
   50 CONTINUE
   60 CONTINUE
      II=II-1
C     WRITE(6,70)II
   70 FORMAT(27H NO. OF ELEMENTS IN CAAA   ,I5)
      RETURN
      END

C-----------------------------------------------------------------------
COMMENT ROUTINE TO TABULATE THE CLEBSCH-GORDON TYPE COEFFICIENTS
C       CLM AND ELM, FOR USE WITH SUBROUTINES XM AND XMT.
C       THE NON-ZERO VALUES ARE TABULATED FIRST FOR (L2+M2) AND
C       (L3+M3) ODD; AND THEN FOR (L2+M2) AND (L3+M3) EVEN -
C       THE SAME SCHEME AS THAT BY WHICH THEY ARE RETRIEVED IN
C       SUBROUTINES XM AND XMT.
C  AUTHOR  PENDRY
      SUBROUTINE  CELMG (CLM,NLM,YLM,FAC2,NN,FAC1,N,LMAX,NFAC,FAC)
C     DIMENSION  CLM(NLM),           YLM(NN), FAC2(NN), FAC1(N)

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL BLM,  CLM(NLM),           YLM(NN), FAC2(NN), FAC1(N)
      DOUBLE PRECISION FAC
      DIMENSION FAC(NFAC)

      PI = 3.14159265
      L2MAX = LMAX + LMAX
      NF=4*LMAX+1

      FAC(1) = 1.0
      DO 340 I=1,NF
  340 FAC(I + 1) = DBLE(I) * FAC(I)

C
COMMENT THE ARRAY YLM IS FIRST LOADED WITH SPHERICAL
C       HARMONICS, ARGUMENTS THETA=PI/2.0, FI=0.0
      LM = 0
      CL = 0.0
      A = 1.0
      B = 1.0
      ASG = 1.0
      LL = L2MAX + 1
COMMENT MULTIPLICATIVE FACTORS REQUIRED
      DO 240 L = 1, LL
      FAC1(L) = ASG * SQRT((2.0 * CL + 1.0) * A/(4.0 * PI * B * B))
      CM =  - CL
      LN = L + L - 1
      DO 230 M = 1, LN
      LO = LM + M
      FAC2(LO) = SQRT((CL + 1.0 + CM) * (CL + 1.0 - CM)/((2.0 * CL + 3.
     10) * (2.0 * CL + 1.0)))
  230 CM = CM + 1.0
      CL = CL + 1.0
      A = A * 2.0 * CL * (2.0 * CL - 1.0)/4.0
      B = B * CL
      ASG =  - ASG
  240 LM = LM + LN
COMMENT FIRST ALL THE YLM FOR M=+-L AND M=+-(L-1) ARE
C       ARE CALCULATED BY EXPLICIT FORMULAE
      LM = 1
      CL = 1.0
      ASG =  - 1.0
      YLM(1) = FAC1(1)
      DO 250 L = 1, L2MAX
      LN = LM + L + L + 1
      YLM(LN) = FAC1(L + 1)
      YLM(LM + 1) = ASG * FAC1(L + 1)
      YLM(LN - 1) = 0.0
      YLM(LM + 2) = 0.0
      CL = CL + 1.0
      ASG =  - ASG
  250 LM = LN
COMMENT USING YLM AND YL(M-1) IN A RECURRENCE RELATION
C       YL(M+1) IS CALCULATED
      LM = 1
      LL = L2MAX - 1
      DO 270 L = 1, LL
      LN = L + L - 1
      LM2 = LM + LN + 4
      LM3 = LM - LN
      DO 260 M = 1, LN
      LO = LM2 + M
      LP = LM3 + M
      LQ = LM + M + 1
      YLM(LO) =  - (FAC2(LP) * YLM(LP))/FAC2(LQ)
  260 CONTINUE
  270 LM = LM + L + L + 1
C
      K = 1
      II = 0
  280 LL = LMAX + II
      DO 310 IL2 = 1, LL
        L2 = IL2 - II
        M2 =  - L2 + 1 - II
        DO 310 I2 = 1, IL2
          DO 300 IL3 = 1, LL
            L3 = IL3 - II
            M3 =  - L3 + 1 - II
            DO 300 I3 = 1, IL3
              LA1 = MAX0(IABS(L2 - L3),IABS(M2 - M3))
              LB1 = L2 + L3
              LA11 = LA1 + 1
              LB11 = LB1 + 1
              M1 = M3 - M2
              DO 290 L11 = LA11, LB11, 2
                L1 = LA11 + LB11 - L11 - 1
                L = (L3 - L2 - L1)/2 + M3
                M = L1 * (L1 + 1) - M1 + 1
                ALM = (( - 1.0)**L) * 4.0 * PI *
     +                BLM(L1,M1,L2,M2,L3, - M3,LMAX,NFAC,FAC)
                CLM(K) = YLM(M) * ALM
  290         K = K + 1
  300     M3 = M3 + 2
  310 M2 = M2 + 2
      IF (II)  320, 320, 330
  320 II = 1
      GO TO 280
330   CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
C  Copypos copies global array SUBPOS (containing atomic positions in
C  composite layer) to local array POS (used in RTINV), ensuring that
C  the RTINV array dimensions are correct.

      SUBROUTINE  COPYPOS(POS,SUBPOS,NLTYPE,NNSUB,NLAY,ILTYPE)

      INTEGER NLTYPE,MNSUB,NLAY
      REAL POS,SUBPOS
      DIMENSION POS(NLAY,3),SUBPOS(NLTYPE,NNSUB,3)

      DO IPOS = 1,3
        DO ILAY = 1,NLAY
          POS(ILAY,IPOS) = SUBPOS(ILTYPE,ILAY,IPOS)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
COMMENT CPPP TABULATES THE FUNCTION PPP(I1,I2,I3), EACH ELEMENT
C       CONTAINING THE INTEGRAL OF THE PRODUCT OF THREE LEGENDRE
C       FUNCTIONS P(I1),P(I2),P(I3). THE INTEGRALS ARE CALCULATED
C       FOLLOWING GAUNT'S SUMMATION SCHEME SET OUT BY SLATER
C       ATOMIC STRUCTURE, VOL1, 309,310
C  PPP IS USED BY SUBROUTINE PSTEMP IN COMPUTING TEMPERATURE-DEPENDENT
C  PHASE SHIFTS
C  AUTHOR  PENDRY
C  DIMENSION 42 REQUIRES N1+N2+N3-1 = 4*LMAX+2 .LE. 42
      SUBROUTINE  CPPP (PPP, N1, N2, N3, NFAC,F)
      DOUBLE PRECISION F
      DIMENSION  PPP(N1,N2,N3), F(NFAC)
      F(1) = 1.0
      DO 370 I = 1, NFAC-1
  370 F(I + 1) = F(I) * DBLE(I)
      DO 460 I1 = 1, N1
      DO 460 I2 = 1, N2
      DO 460 I3 = 1, N3
      IM1 = I1
      IM2 = I2
      IM3 = I3
      IF (I1-I2)  380, 390, 390
  380 IM1 = I2
      IM2 = I1
  390 IF (IM2-I3)  400, 420, 420
  400 IM3 = IM2
      IM2 = I3
      IF (IM1-IM2)  410, 420, 420
  410 J = IM1
      IM1 = IM2
      IM2 = J
  420 A = 0.0
      IS = I1 + I2 + I3 - 3
      IF (MOD(IS,2)-1)  430, 460, 430
  430 IF (IABS(IM2-IM1)-IM3+1)  440, 440, 460
  440 SUM = 0.0
      IS = IS/2
      SIGN = 1.0
      DO 450 IT = 1, IM3
      SIGN =  - SIGN
      IA1 = IM1 + IT - 1
      IA2 = IM1 - IT + 1
      IA3 = IM3 - IT + 1
      IA4 = IM2 + IM3 - IT
      IA5 = IM2 - IM3 + IT
  450 SUM = SUM - SIGN * F(IA1) * F(IA4)/(F(IA2) * F(IA3) * F(IA5) * F(
     1IT))
      IA1 = 2 + IS - IM1
      IA2 = 2 + IS - IM2
      IA3 = 2 + IS - IM3
      IA4 = 3 + 2 * (IS - IM3)
      A =  - ( - 1.0)**(IS - IM2) * F(IA4) * F(IS + 1) * F(IM3) * SUM/(
     1F(IA1) * F(IA2) * F(IA3) * F(2 * IS + 2))
  460 PPP(I1,I2,I3) = A
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE CXMTXT SOLVES A SET OF LINEAR EQUATIONS.                 ! << TODO Called inside TAUMAT with INOPT == -1. Could it be replaced with CGERF/CGERS?
C     A = INPUT MATRIX.  A  IS DESTROYED DURING THE COMPUTATION AND
C         IS REPLACED BY THE UPPER TRIANGULAR MATRIX RESULTING FROM THE
C         GAUSSIAN ELIMINATION PROCESS (WITH PARTIAL PIVOTING).
C     NC = FIRST DIMENSION OF A IN CALLING ROUTINE. NC.GE.NR.
C     NR = ORDER OF A
C     NSYS = NO. OF SYSTEMS TO BE SOLVED.  IF INVERSE OPTION IS CHOSEN
C            NSYS MUST BE AT LEAST AS LARGE AS NR .
C            COEFFICIENT MATRIX MUST BE STORED IN A(I,J), I=1,NR, J=1,NR
C            CONSTANT VECTORS MUST BE STORED IN A(I,NR+1), I=1,NR:
C            A(I,NR+2), I=1,NR  .... A(I,NR+NSYS), I=1,NR.
C            RESULT OVERWRITES CONSTANT VECTORS
C     NTOT = NR + NSYS
C     MARK = SINGULARITY INDICATOR (MARK=1 FOR SINGULAR A)
C     DET = DET(A)
C     INOPT = -1 FOR SYSTEM SOLN. AND DET
C              0 FOR DET ONLY
C             +1 FOR INVERSE AND DET
C     X is auxiliary array carried through only to ensure proper dimensions
C     DIM. OF X ARRAY MUST BE AT LEAST AS LARGE AS FIRST DIM. OF A ARRAY
      SUBROUTINE  CXMTXT (A, NC,NR, NSYS, NTOT, MARK, DET, INOPT, X)
      REAL  AMAX, QZ
      COMPLEX  A, X
      COMPLEX  AGG, DET, CONST, SIGN, TEMP
      DIMENSION A(NC, NTOT), X(NC)
  300 FORMAT(5X,6HDETC= ,2E20.5)
  310 FORMAT(//, 1X, 16HSINGULAR MATRIX , //)
C
C     PRESET PARAMETERS
      SIGN = (1.0E + 00,0.0E + 00)
      MARK = 0
      IFLAG = INOPT
      N = NR
      NPL = N + 1
      NMI = N - 1
      NN = N + N
      NPLSY = N + NSYS
      IF (IFLAG)  950, 950, 920
C
C     INVERSE OPTION - PRESET AUGMENTED PART TO I
  920 DO 930 I = 1, N
      DO 930 J = NPL, NN
  930 A(I,J) = (0.0E + 00,0.0E + 00)
      DO 940 I = 1, N
      J = I + N
  940 A(I,J) = (1.0E + 00,0.0E + 00)
      NPLSY = NN
C
C     TRIANGULARIZE A
  950 DO 1020 I = 1, NMI
      IPL = I + 1
C     DETERMINE PIVOT ELEMENT
      MAX = I
      AMAX = CABS(A(I,I))
      DO 970 K = IPL, N
      QZ = CABS(A(K,I))
      IF (AMAX-QZ)  960, 970, 970
  960 MAX = K
      AMAX = CABS(A(K,I))
  970 CONTINUE
      IF (MAX-I)  980, 1000, 980
C     PIVOTING NECESSARY - INTERCHANGE ROWS
  980 DO 990 L = I, NPLSY
      TEMP = A(I,L)
      A(I,L) = A(MAX,L)
  990 A(MAX,L) = TEMP
      SIGN =  - SIGN
C     ELIMINATE A(I+1,I)---A(N,I)
 1000 DO 1020 J = IPL, N
      TEMP = A(J,I)
      QZ = CABS(TEMP)
      IF(QZ .LT. 1.0E-10) GO TO 1020
      CONST =  - TEMP/A(I,I)
      DO 1010 L = I, NPLSY
 1010 A(J,L) = A(J,L) + A(I,L) * CONST
 1020 CONTINUE
C
C     COMPUTE VALUE OF DETERMINANT
      TEMP = (1.0E + 00,0.0E + 00)
      DO 1030 I = 1, N
      AGG = A(I,I)
      QZ = CABS(AGG)
      IF(QZ .GT. 1.0E-10) GO TO 1030
C     MATRIX SINGULAR
C     WRITE(6,310)
      MARK = 1
      GO TO 1040
 1030 TEMP = TEMP * AGG
      DET = SIGN * TEMP
C     WRITE(6,300) DET
C
C     EXIT IF DET ONLY OPTION
 1040 IF (IFLAG)  1050, 1160, 1050
C     CHECK FOR INVERSE OPTION OR SYSTEMS OPTION
 1050 IF (IFLAG-1)  1070, 1060, 1070
C     INVERSE OPTION - ABORT IF A IS SINGULAR
 1060 IF (MARK-1)  1070, 1160, 1070
C
C     BACK SUBSTITUTE TO OBTAIN INVERSE OR SYSTEM SOLUTION(S)
 1070 DO 1150 I = NPL, NPLSY
      K = N
 1080 TEMP = A(K,I)                                                     100589
      IF (K-N)  1090, 1110, 1090
 1090 DO 1100 J = KPL, N
 1100 TEMP = TEMP - A(K,J) * X(J)                                       100589
 1110 X(K) = TEMP/A(K,K)                                                100589
      IF (K-1)  1120, 1130, 1120
 1120 KPL = K
      K = K - 1
      GO TO 1080
C     PUT SOLN. VECT. INTO APPROPRIATE COLUMN OF A
 1130 DO 1140 L = 1, N
 1140 A(L,I) = X(L)
 1150 CONTINUE
C
 1160 RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE DBGT_MOD COMPUTES THE REFLECTION MATRIX FOR A PAIR OF LAYERS   ! TODO: optimize?
!  FOR INCIDENCE FROM ONE SIDE ONLY. LAYER DOUBLING IS USED. NO
!  OVERWRITING OF INPUT TAKES PLACE.
!  NOTE  SUBROUTINE DBLG PRODUCES ADDITIONNALLY THE REFLECTION FROM
!  THE OTHER SIDE OF THE LAYERS AND THE TRANSMISSIONS IN BOTH DIRECTIONS
!   RA1,TA1,RA2,TA2= REFLECTION AND TRANSMISSION MATRICES FOR FIRST
!    LAYER (ON THE SIDE OF THE INCIDENT BEAMS). 1 AND 2 REFER TO
!    INCIDENCE ON FIRST LAYER FROM SIDE OPPOSITE SECOND LAYER AND SIDE
!    TOWARDS SECOND LAYER, RESP.
!   RB1= REFLECTION MATRIX OF SECOND LAYER (E.G. SUBSTRATE).
!   RAB= OUTPUT REFLECTION MATRIX OF COMBINED LAYER.
!   ASD= INTERLAYER VECTOR FROM FIRST LAYER TO SECOND LAYER.
!   N=NO. OF BEAMS USED.
!   NA= OFFSET IN LIST PQ OF PARTICULAR SUBSET OF BEAMS USED.
!   PQ= LIST OF BEAMS G.
!   NT= TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
!   S1,S2,XS,PP,IPL= WORKING SPACES.
!   NP=DIMENSION FOR WORKING SPACES, NP.GE.N.
!   EMACH= MACHINE ACCURACY.
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY.
!   VPI= IMAGINARY PART OF ENERGY.
!   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!
!  ** SPECIAL VERSION FOR TENSOR-LEED
!  ** THIS SUBROUTINE PRODUCES ALL MATRICES NEEDED FOR CALCULATION OF
!  ** AMPLITUDES AT THE CURRENT LAYER IN -X DIRECTION  AND THOSE AT
!  ** THE LAYER BELOW IN +X DIRECTION BY MULTIPLYING THEM WITH THE
!  ** AMPLITUDES OF CURRENT LAYER IN +X DIRECTION (E.G. PRODUCED
!  ** BY SUBROUTINE ADREF1T). USE MULTAMP FOR THIS PURPOSE.
!  ** DIAGRAM OF ADDITIONAL INFORMATION:
!  **
!  ** LAYER   INTERL VEC   DIFFR MATRICES     AMPLITUDES  MATRICES
!  ** --------------------------------------------------------------
!  **   >> N-1 TOP LAYERS
!  **          ASE         + RA1,TA1            + AMP0
!  **   N ----------------------------------------------------------
!  **          ASD         - RA2,TA2            - AMP0  *  S1
!  **                      + RAB                + AMP0  *  S2
!  **   N+1 --------------------------------------------------------
!  **   >> PREVIOUS CALCULATION BY LAYER-DOUBLING
!  ** --------------------------------------------------------------
!  **
!  ** TO CALCULATE AMPLITUDES OF DEEPER LAYERS USE DBGT INSTEAD OF DBG
!  ** WHICH WILL PRODUCE THE MATRICES TO BE MULTIPLIED WITH AMP'S
!  ** CALCULATED BY ADREF1T
!  **
      SUBROUTINE DBGT_MOD(RA1, TA1, RA2, TA2, RB1, RAB, ASD, N,                     **T
     &                    PQ, S1, S2, XS, PP, IPL)
      COMPLEX RAB, RA1, TA1, RB1
      COMPLEX  S1, S2, PP, XS, XX, CZ, IU
      COMPLEX TEMP
      COMPLEX RA2(N, N), TA2(N, N)
      DIMENSION RA1(N, N), TA1(N, N), RB1(N, N), RAB(N, N),
     &          S1(N, N), S2(N, N)
      DIMENSION PP(N, 2), XS(N), PQ(2, N), IPL(N)
      DIMENSION  ASD(3)
      COMMON  E, AK2, AK3, VPI

      CZ = CMPLX(0.0,0.0)
      IU = CMPLX(0.0,1.0)
      AK = 2.0 * E

!     COMPUTE INTERLAYER PROPAGATORS PP
      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)
        XX = CMPLX(AK - BK2 * BK2 - BK3 * BK3, - 2.0 * VPI + 0.000001)
        XX = SQRT(XX) * ASD(1)           ! Phase shift due to perpendicular propagation
        X = BK2 * ASD(2) + BK3 * ASD(3)  ! Phase shift due to parallel propagation
        PP(I, 1) = EXP(IU * (XX + X))
        PP(I, 2) = EXP(IU * (XX - X))
      ENDDO

!     PRODUCE multiple-scattering MATRIX TO BE INVERTED
      DO K = 1, N
        TEMP = PP(K, 1)
        DO J = 1, N
          XX = CZ
          DO L = 1, N
            XX = XX - (RA2(J, L)
     +                 * PP(L, 2)
     +                 * RB1(L, K)
     +                 * TEMP)
          ENDDO
          S1(J, K) = XX
        ENDDO
        S1(K, K) = S1(K, K) + 1.0
      ENDDO

!     GAUSSIAN ELIMINATION STEP OF MATRIX INVERSION
      CALL CGETRF(N, N, S1, N, IPL, INFO)

!     SUBSTITUTION STEP OF MATRIX INVERSION
      CALL CGETRS('N', N, N, S1, N, IPL, TA1, N, INFO)

      DO K = 1, N
        DO L = 1, N
          XS(L) = PP(L, 1) * TA1(L, K)
        ENDDO
        DO J = 1, N  ! This loop can be replaced with a call to CGEMV??
          S1(J, K) = CZ
          DO L = 1, N
            S1(J, K) = S1(J, K) + RB1(J, L) * XS(L)
          ENDDO
        ENDDO

        DO L = 1, N
          S2(L, K) = XS(L)
          S1(L, K) = PP(L, 2) * S1(L, K)
        ENDDO
        DO J = 1, N  ! TODO: This loop, with the outer K, can be replaced with a call to CGEMM('N', 'N', N, N, N, (1.0, 0.0), TA2, N, S1, N, (1.0, 0.0), RAB, N), after assigning RAB = RA1
          RAB(J, K) = RA1(J, K)
          DO L = 1, N
            RAB(J, K) = RAB(J, K) + TA2(J, L) * S1(L, K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

!---------------------------------------------------------------------------
!  SUBROUTINE DBLGAS_MOD computes the combined reflection/transmission
!  matrices of two, possibly different, bulk layers from their
!  respective reflection/transmission matrices. Both reflection and
!  transmission may be different for incidence toward the solid and
!  away from it.
!
!  Edits: MRiva 2021-09-15 -- remove GOTOs, remove useless arguments
!         NT, NP (both == N) and EAMCH (unused). Rework code to use
!         effectively LAPACK (BLAS) routine CGEMM for matrix-matrix
!         multiplication. Get rid of S3, S4, and XS, unused as a
!         result of the reworking of the code.
!
!   RAMP,TAPP= REFLECTION AND TRANSMISSION MATRICES OF FIRST LAYER FOR
!    INCIDENCE FROM SIDE OF FIRST LAYER (OVERWRITTEN BY REFLECTION AND
!    TRANSMISSION MATRICES FOR INCIDENCE FROM SIDE OF FIRST LAYER),
!    i.e., beams propagating toward the solid
!   RAPM,TAMM= REFLECTION AND TRANSMISSION MATRICES OF FIRST LAYER FOR
!    INCIDENCE FROM SIDE OF SECOND LAYER, i.e., beams propagating away
!    from the solid.
!   RBPM,TBMM= REFLECTION AND TRANSMISSION MATRICES OF SECOND LAYER FOR
!    INCIDENCE FROM SIDE OF SECOND LAYER (OVERWRITTEN BY REFLECTION AND
!    TRANSMISSION MATRICES FOR INCIDENCE FROM SIDE OF SECOND LAYER),
!    i.e., beams propagating toward the solid
!   RBMP,TBPP= REFLECTION AND TRANSMISSION MATRICES OF SECOND LAYER FOR
!    INCIDENCE FROM SIDE OF FIRST LAYER, i.e., beams propagating away
!    from the solid.
!   ASD= INTERLAYER VECTOR FROM FIRST LAYER TO SECOND LAYER.
!   N=TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
!   PQ= LIST OF BEAMS G.
!   S1,S2,IPL= WORKING SPACES.
!   PP= interlayer propagators. PP(...,1) is propagation toward the
!       solid, PP(...,2) away from it. Calculated and used internally.
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY.
!   VPI= IMAGINARY PART OF ENERGY.
!   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
      SUBROUTINE DBLGAS_MOD(RAMP, TAPP, RAPM, TAMM,
     &                      RBMP, TBPP, RBPM, TBMM,
     &                      ASD, N, PQ, S1, S2, PP,
     &                      IPL)
      COMPLEX RAMP, TAPP, RAPM, TAMM, RBPM, TBMM, RBMP, TBPP
      COMPLEX  S1, S2, PP, XX, CZ, IU, C_ONE,
     &         C_M_ONE, P_UP, P_DOWN
      DIMENSION  RAMP(N, N), TAPP(N, N), RAPM(N, N), TAMM(N, N)
      DIMENSION  RBPM(N, N), TBMM(N, N), RBMP(N, N), TBPP(N, N)
      DIMENSION  S1(N, N), S2(N, N)
      DIMENSION  PP(N, 2), PQ(2, N), IPL(N)
      DIMENSION  ASD(3)
      COMMON  E, AK2, AK3, VPI

      CZ = CMPLX(0.0, 0.0)
      IU = CMPLX(0.0, 1.0)
      C_ONE = CMPLX(1.0, 0.0)
      C_M_ONE = CMPLX(-1.0, 0.0)
      AK = 2.0 * E

!     Compute interlayer propagators PP
      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)
        XX = CMPLX(AK - BK2 * BK2 - BK3 * BK3, - 2.0 * VPI + 0.000001)
        XX = SQRT(XX) * ASD(1)            ! out-of-plane, propagation phase shift
        X = BK2 * ASD(2) + BK3 * ASD(3)   ! in-plane, propagation phase shift
        PP(I, 1) = EXP(IU * (XX + X))     ! Beam i propagating toward solid
        PP(I, 2) = EXP(IU * (XX - X))     ! Beam i propagating away from solid
      ENDDO

!     Modify some of the reflection/transmission matrices to also
!     include propagation: RAPM(i, j), RBMP(i, j), TAMM(i, j)
!     and TBPP(i, j) are now:
!
!                                                       ^ i
!                                                  TAMM |
!     --<<TOP>>------------------------------------------------------
!                  j ^ RAPM |                         j ^
!                    |      v i                         |
!                                     |      ^ i               |
!                                   j v RBMP |               j v
!     --<<BOTTOM>>---------------------------------------------------
!                                                              | TBPP
!                                                              v i
!
!     Include for all the propagation of the incoming beam (j) from
!     the other layer. This means that the matrices become then:
!
!                                                       ^ i
!                                                  TAMM |
!     --<<TOP>>------------------------------------------------------
!                    ^ RAPM |       j |                 ^    j |
!                    |      v i       |                 |      |
!                    |                |      ^ i        |      |
!                  j |                v RBMP |        j |      v
!     --<<BOTTOM>>---------------------------------------------------
!                                                              | TBPP
!                                                              v i
      DO J = 1, N
        P_UP = PP(J, 2)
        P_DOWN = PP(J, 1)
        DO I = 1, N
          RAPM(I, J) = RAPM(I, J) * P_UP
          RBMP(I, J) = RBMP(I, J) * P_DOWN
          TAMM(I, J) = TAMM(I, J) * P_UP
          TBPP(I, J) = TBPP(I, J) * P_DOWN
        ENDDO
      ENDDO

!     Prepare first multiple-scattering matrix S1, used for calculating
!     reflection/transmission for incidence from above. Multiple
!     scattering is computed exactly from S1 = Id - X, since
!     (Id - X)**(-1) == sum_n X**n. Each single-pass scattering X(i, j)
!     is beam j coming from right below top layer, scattered at bottom,
!     propagated up, and scattered at top again. This means
!     X(i, j) = sum_k [RBMP(k, j)     j from above, incl. propagation, scattered into k at bottom
!                     * RAPM(i, k)]   k then propagated to top and scattered there into i, moving down
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S1 = Id - X
     &           RAPM, N, RBMP, N,
     &           CZ, S1, N)
      DO I = 1, N
        S1(I, I) = C_ONE + S1(I, I)  ! Id - X, diagonal
      ENDDO

!     LU decomposition for matrix inversion
      CALL CGETRF(N, N, S1, N, IPL, INFO)

!     Make TAPP = inv(S1)*TAPP, i.e., TAPP(i, j) is then beam j from
!     above, transmitted through top layer, and multiple-scattered
!     into beam i, which is right below the top layer, propagating down
      CALL CGETRS('N', N, N, S1, N, IPL, TAPP, N, INFO)

!     Prepare to compute the total reflection from above (stored in RAMP).
!     First, propagate-reflect the transmitted beam at the bottom layer
!     by taking the product S1(i, j) = sum_k RBMP(i, k) TAPP(k, j),
!     corresponding to incoming beam j from above top layer, and outgoing
!     beam i right above bottom layer (propagating up).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, RBMP, N, TAPP, N, CZ, S1, N)

!     Then propagate-transmit this at the top layer by multiplying
!     TAMM(i, k) with S1(k, j). Also, add to this 'reflection due to
!     multiple scattering in between' the 'specular' reflection RAMP(i, j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TAMM, N, S1, N, C_ONE, RAMP, N)

!     Finally, also compute the total transmission from above top layer
!     to below bottom one: Propagate outgoing beam k of TAPP(k, j) to
!     bottom layer, and transmit it there. Both done by multiplying
!     with TBPP(i, k), which already includes the propagation.
!     Notice that the calculation is done in S1, and not directly in
!     TAPP, as the CZ==0 term would zero it as a buffer before it can
!     be used in the calculation. TAPP = S1 is done later at the end.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TBPP, N, TAPP, N, CZ, S1, N)


!     Prepare second multiple-scattering matrix, used for calculating
!     reflection/transmission for incidence from below. Multiple
!     scattering is computed exactly from S2 = Id - X, since
!     (Id - X)**(-1) == sum_n X**n. Each single-pass scattering is
!     beam j coming from right below top layer, scattered at bottom,
!     propagated up, and scattered at top again. This means
!     X(i,j) = sum_k [RAPM(k, j)      j from below, incl. propagation, scattered into k at top
!                     * RBMP(i, k)]   k then propagated to bottom and scattered there into i, propagating up
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S2 = Id - X
     &           RBMP, N, RAPM, N,
     &           CZ, S2, N)
      DO I = 1, N
        S2(I, I) = C_ONE + S2(I, I)  ! Id - X, diagonal
      ENDDO

!     LU decomposition for matrix inversion
      CALL CGETRF(N, N, S2, N, IPL, INFO)

!     Make TBMM = inv(S2)*TBMM, i.e., TBMM(i,j) is then beam j from
!     below, transmitted through bottom layer, and multiple-scattered
!     into beam i, which is above the bottom layer, propagating up
      CALL CGETRS('N', N, N, S2, N, IPL, TBMM, N, INFO)

!     Prepare to compute the total reflection from below (stored in RBPM).
!     First, propagate-reflect the transmitted beam at the top layer
!     by taking the product S2(i, j) = sum_k RAPM(i, k) TBMM(k, j),
!     corresponding to incoming beam j from below bottom layer, and outgoing
!     beam i right below top layer (propagating down).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, RAPM, N, TBMM, N, CZ, S2, N)

!     Then propagate-transmit it at the bottom layer by multiplying
!     TBPP(i, k) with S2(k, j). Also, add to this 'reflection due to
!     multiple scattering in between' the 'specular' reflection RBPM(i, j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TBPP, N, S2, N, C_ONE, RBPM, N)

!     Finally, also compute the total transmission from below bottom
!     layer to above top one: Propagate outgoing beam k of TBMM(k, j) to
!     top layer, and transmit it there. Both done by multiplying
!     with TAMM(i, k), which already includes the propagation.
!     Notice that the calculation is done in S2, and not directly in
!     TBMM, as the CZ==0 term would zero it as a buffer before it can
!     be used in the calculation. TBMM = S2 is done later at the end.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TAMM, N, TBMM, N, CZ, S2, N)

      DO K = 1, N
        DO J = 1, N
          TAPP(J, K) = S1(J, K)
          TBMM(J, K) = S2(J, K)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE DEBWAL COMPUTES DEBYE-WALLER FACTORS FOR
C  DIFFRACTION FROM BEAM G(PRIME) INTO BEAM G (AND ITS SYMMETRICAL
C  COUNTERPARTS) FOR, IN GENERAL, ANISOTROPIC ATOMIC VIBRATION,
C  INCLUDING ZERO-TEMPERATURE VIBRATION.
C   NG= NO. OF BEAMS CORRESPONDING TO G.
C   G= SET OF SYMMETRICAL BEAMS G.
C   GP= INCIDENT BEAM G(PRIME).
C   E= ENERGY.
C   VPI= OPTICAL POTENTIAL
C   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C   T= ACTUAL TEMPERATURE.
C   T0= REFERENCE TEMPERATURE FOR VIBRATION AMPLITUDES.
C   DRX= RMS VIBRATION AMPLITUDE PERPENDICULAR TO SURFACE.
C   DRY= RMS VIBRATION AMPLITUDE PARALLEL TO SURFACE.
C   D04= FOURTH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE
C    (ISOTROPIC).
C   EDW= OUTPUT DEBYE-WALLER FACTOR (EDW(1,I) FOR REFLECTION, EDW(2,I)
C    FOR TRANSMISSION).
      SUBROUTINE  DEBWAL (NG, G, GP, E, VPI, AK2, AK3, T, T0, DRX, DRY,
     1D04, EDW)
      COMPLEX SQRT
      COMPLEX EDW
      DIMENSION  GP(2), G(2,12), EDW(2,NG)

      write(6,*) "Welcome to subroutine Debwal. If you really wish to"
      write(6,*) "compute anisotropic vibrations, you should consider"
      write(6,*) "using Tensor LEED. VB, 12.06.98."

      A1 = GP(1) + AK2
      A2 = GP(2) + AK3
      CC = REAL(SQRT(CMPLX(2.0 * E - A1 * A1 - A2 * A2, - 2.0 * VPI)))
      A1 = G(1,1) + AK2
      A2 = G(2,1) + AK3
      DD = REAL(SQRT(CMPLX(2.0 * E - A1 * A1 - A2 * A2, - 2.0 * VPI)))
C  C IS PERPENDICULAR COMPONENT OF SCATTERING VECTOR FOR REFLECTION,
C  D IS SAME FOR TRANSMISSION
      C = CC + DD
      D = CC - DD
      DO 330 I = 1, NG
      A1 = GP(1) - G(1,I)
      A2 = GP(2) - G(2,I)
C  D2 IS MEAN-SQUARE VIBRATION AMPLITUDE PARALLEL TO SURFACE AT ACTUAL
C  TEMPERATURE
      D2 = DRY * DRY * T/T0
C  ZERO-TEMPERATURE VIBRATION IS NOW MIXED IN
      D2 = SQRT(D2 * D2 + D04)
      A1 = (A1 * A1 + A2 * A2) * D2
C  D2 IS SAME AS ABOVE, BUT FOR PERPENDICULAR COMPONENTS
      D2 = DRX * DRX * T/T0
      D2 = SQRT(D2 * D2 + D04)
      EDW(1,I) = CMPLX(EXP( - 0.166667 * (A1 + C * C * D2)),0.0)
  330 EDW(2,I) = CMPLX(EXP( - 0.166667 * (A1 + D * D * D2)),0.0)
      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE DECIDE determines which parts of the computation are
C  necessary for the current beam NEXIT and sets calculational
C  quantities appropriately (replaces old cryptographic masterpiece
C  wavenew). Note ICRSYM is not retained.

C  criteria: NEXIT = 0 -> do normal calculation
C            NEXIT > 0
C              evanescent beam -> skip complete calculation
C              non-evanescent beam
C                does (incident!) beam belong to those (exit!) beams for
C                which layer diffraction matrices were last calculated?
C                yes -> skip fd calculation up to ADREF2T
C                no (e.g. for first TLEED beam in case of off-normal incidence)
C                    -> perform entire calculation using this beam as inc.

      SUBROUTINE DECIDE(NEXIT,KNOWN,EMERGE,PSQ1,PSQ2,AK2,AK3,AK21,AK31,
     +                  KSQ,AK2M,AK3M,NT0,RBR1,RBR2)

C  NEXIT is current beam number
C  EMERGE is 1 if current beam NEXIT propagates in vacuum, 0 otherwise
C  KNOWN is 1 if scattering matrices for current beam have already been
C        calculated, 0 otherwise
C  PSQ1,PSQ2 is lateral momentum of current beam relative to last beam for
C        which layer diffraction matrices were computed (used in ADREF2T)
C  AK2,AK3 is lateral momentum of (true) incident beam
C  AK21,AK31 is lateral momentum of incident beam for which scattering
C        matrices were last computed (set to current value if computation required)
C  KSQ   is square of total momentum of wave
C  AK2M,AK3M is array containing negative lateral momentum of all NT0 Tensor LEED beams
C  NT0   is number of TLEED beams
C  RBR1,RBR2 are superlattice reciprocal vectors

      INTEGER NEXIT
      INTEGER EMERGE,KNOWN
      REAL PSQ1,PSQ2,AK2,AK3,AK21,AK31,KSQ
      REAL AK2M,AK3M
      DIMENSION AK2M(NT0),AK3M(NT0)
      INTEGER NT0
      REAL RBR1,RBR2
      DIMENSION RBR1(2),RBR2(2)

C  local variables

      REAL DET1,DET2

C  end declarations

      IF (NEXIT.eq.0) THEN

c  incident beam emerges, and scattering is certainly not yet known.

        EMERGE = 1
        KNOWN  = 0

        PSQ1=0.
        PSQ2=0.

        AK21 = AK2
        AK31 = AK3

        write(6,*) "Computation for incident beam started."

      ELSE

        IF ( (KSQ-AK2M(NEXIT)**2.-AK3M(NEXIT)**2.) .gt. 0. ) THEN

          EMERGE = 1

          write(6,*) "Beam NEXIT = ", NEXIT," can propagate in vacuum",
     +               " - computation required."

C  AK21,AK31 are lateral components of last incident beam for which entire
C  scattering was computed. Normally, that would be incident beam. In case
C  of off-normal incidence, recomputation is necessary for the first outgoing
C  beam as well. Only for beams not present in the reference structure at all,
C  another computation is necessary - this feature is currently not implemented.

          DET1 = (AK21-AK2M(NEXIT))*RBR2(2)
     +         - (AK31-AK3M(NEXIT))*RBR2(1)
          DET2 = (AK31-AK3M(NEXIT))*RBR1(1)
     +         - (AK21-AK2M(NEXIT))*RBR1(2)

          DET1 = DET1/(RBR1(1)*RBR2(2)-RBR1(2)*RBR2(1))
          DET2 = DET2/(RBR1(1)*RBR2(2)-RBR1(2)*RBR2(1))

          DET1 = ABS(ABS(DET1)-FLOAT(IFIX(ABS(DET1)+1.E-3)))
          DET2 = ABS(ABS(DET2)-FLOAT(IFIX(ABS(DET2)+1.E-3)))

C  now, if DET1 or DET2 .ne. 0, the current beam can not be related to the
C  one for which the last computation was performed by a superlattice vector,
C  i.e. a new computation is required
C  s.w. non-perpendicular incidence can be handled now

          IF ((DET1.gt.1.E-3).or.(DET2.gt.1.E-3).or.
     +        ((NEXIT.eq.1).and.((AK2.ne.0.).or.(AK3.ne.0.)))) THEN
            KNOWN = 0
            AK21 = AK2M(NEXIT)
            AK31 = AK3M(NEXIT)
            PSQ1 = 0.
            PSQ2 = 0.
            write(6,*) "Scattering matrices not yet known",
     +                 " - recompute entire scattering."
          ELSE
            KNOWN = 1
            PSQ1 = AK21 - AK2M(NEXIT)
            PSQ2 = AK31 - AK3M(NEXIT)
C            write(6,*) "AK21=",AK21,"AK2M=", AK2M(NEXIT)
C            write(6,*) "AK31=",AK31,"AK3M=", AK3M(NEXIT)
            write(6,*) "Scattering matrices known from previous beams",
     +                 " - only compute incoming amplitudes."
          END IF

        ELSE

          EMERGE = 0
          KNOWN = 0

          write(6,*) "Beam NEXIT = ", NEXIT," can't propagate at ",
     +               "current energy - skip computation."

        END IF

      END IF

      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE DELIMIT WRITES -1 TO OUTPUTFILE FOR AMPLITUDES
C
      SUBROUTINE DELIMIT(IFILE,IFORM)
C
      IF (IFORM.EQ.0) THEN
           WRITE (IFILE) -1
      ELSE
           WRITE (IFILE,10) -1
   10 FORMAT (I5) ! changed to I5 in v1.71
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
C  FUNCTION FACT COMPUTES FACTORIAL(L)/10**L, USING AN ASYMPTOTIC
C  EXPANSION FOR L.GT.4. USED IN CA.
C  expansion works extremely well, much higher than L=15
      DOUBLE PRECISION FUNCTION  FACT (L)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      IF (L .GT. 4)  GO TO 1600
      IF (L .EQ. 0)  FACT = 1.0
      IF (L .EQ. 1)  FACT = 0.1
      IF (L .EQ. 2)  FACT = 0.02
      IF (L .EQ. 3)  FACT = 6.0 * 0.001
      IF (L .EQ. 4)  FACT = 24.0 * 0.0001
      RETURN
 1600 X = L + 1
      FACT = DEXP( - X) * (10.0D0**((X - 0.5D0) * DLOG10(X) - (X - 1.
     10D0))) * (DSQRT(6.283185307179586D0)) * (1.0 + (1.0/(12.0 * X)) +
     2(1.0/(288.0D0 * (X**2))) - (139.0D0/(51840.0D0 * (X**3))) - (571.
     30D0/(2488320.0D0 * (X**4))))
      RETURN
      END

C-----------------------------------------------------------------------
C      ROUTINE FINAL1 EVALUATES THE AMPLITUDE OF SPHERICAL
C      WAVES INCIDENT UPON THE ATOM AT THE ORIGIN OF THE
C      TOP LAYER GIVEN THE AMPLITUDE OF PLANE WAVES INCIDENT
C      ON EITHER SIDE OF THE TOP LAYER.
C
C  A0LM(LMMAX): INITIALLY INCIDENT SPH WAVE AMPLITUDES.
C  ALM(LMMAX) : FINALLY   INCIDENT SPH WAVE AMPLITUDES.
C
C  APLUS(NT)  : PLANE WAVES INCIDENT FROM -X (I.E. WITH
C               POSITIVE KZ- KZPLUS)
C  AMINUS(NT) : PLANE WAVES INCIDENT FROM +X (I.E. WITH
C               NEGATIVE KZ- KZMINUS)
C
C  CYLM(NT,LMMAX): SPHERICAL HARMONICS FOR EACH BEAM WITH
C                  POSITIVE KZ.FOR NEGATIVE KZ MULTIPLY BY
C                  (-1.)**(L+M)
C
C  XODST(LOD,LOD): TRANSPOSE OF ODD BLOCK OF MATRIX (1-X)
C  XEVST(LEV,LEV): TRANSPOSE OF EVEN BLOCK OF MATRIX (1-X)
C
C  IPLO(LOD),IPLE(LEV) : ELIMINATION INFORMATION FROM ZGE.
C
      SUBROUTINE FINAL1(ALM,A0LM,APLUS,AMINUS,
     1CYLM,XODST,XEVST,IPLO,IPLE,LMAX,LMMAX,
     1LOD,LEV,NT,LX,AEV,AOD,EMACH)
C
C
      COMPLEX AOD(LOD),AEV(LEV)
      COMPLEX ALM(LMMAX),A0LM(LMMAX)
      COMPLEX APLUS(NT),AMINUS(NT),CYLM(NT,LMMAX)
      COMPLEX XODST(LOD,LOD),XEVST(LEV,LEV)
      COMPLEX PRE,CZ,CI,ATEMP
      INTEGER IPLO(LOD),IPLE(LEV),LX(LMMAX)
C
C
      PI=4.0*ATAN(1.0)
      CI=CMPLX(0.0,1.0)
      CZ=CMPLX(0.0,0.0)
C
C
C   EVALUATE THE INTIALLY INCIDENT SPHERICAL WAVE AMPLITUDES
C   INCIDENT UPON THE ATOM AT THE ORIGIN OF THE TOP LAYER
C   PRODUCED BY THE INCIDENT PLANE WAVES APLUS & AMINUS
C
      DO 1 L=0,LMAX
      DO 1 M=-L,L
C
      LM=(L+1)**2-L+M
      LMM=LM-2*M
C
      PRE=4.0*PI*(CI**L)
      AL=((-1.)**L)
      AM=((-1.)**M)
      A0LM(LM)=CZ
C
      DO 1 IG=1,NT
C
      ATEMP=PRE*(AM*APLUS(IG)+AL*AMINUS(IG))*CYLM(IG,LMM)
C
   1  A0LM(LM)=A0LM(LM)+ATEMP
C
  32  FORMAT(I4,2E14.4)
C
C
C
C  EVALUATE FINAL SPHERICAL WAVE AMPLITUDES ALM FROM A0LM
C  BY MULTIPLYING ON THE LEFT BY (1-X) TRANSPOSED. I.E.
C  ALM=A0LM*(1-X) = (1-X)**T * A0LM
C
C  FIRST SORT A0LM INTO ODD & EVEN PARTS.
C
C
      DO 2 N=1,LEV
      AEV(N)=A0LM(LX(N))
    2  CONTINUE
C  2  WRITE(6,93)N,AEV(N)
C
      DO 3 N=(LEV+1),LMMAX
      AOD(N-LEV)=A0LM(LX(N))
   3  CONTINUE
C  3  WRITE(6,93)N,AOD(N)
  93  FORMAT(I4,2E14.4)
C
C
C  GAUSSIAN ELLIMINATE TRANSPOSED (1-X) MATRIX       ! Michele: calls can be replaced with CGETRF /TODO
C
      CALL ZGE(XODST,IPLO,LOD,LOD,EMACH)
      CALL ZGE(XEVST,IPLE,LEV,LEV,EMACH)
C
C  INVERT (1-X) & MULTIPLY BY A0LM ON THE RIGHT      ! Michele: calls can be replaced with CGETRS? Check what ZTU does differently than ZSU. Looks like it solves the system A**T X = B for X. /TODO
C
      CALL ZTU(XODST,IPLO,AOD,LOD,LOD,EMACH)
      CALL ZTU(XEVST,IPLE,AEV,LEV,LEV,EMACH)
C
C  STORE AOD & AEV
C
C
      DO N=1,LEV
          ALM(LX(N))=AEV(N)
      ENDDO
      DO N=(LEV+1),LMMAX
          ALM(LX(N))=AOD(N-LEV)
      ENDDO
C
C
      RETURN
      END

C---------------------------------------------------------------------
C  AUTHOR P. ROUS, MODIFIED W. OED 121190
      SUBROUTINE FINALOV(TSTORE,ALM,APLUS,AMINUS,NT,NLAYER,
     1LMN,LMMAX,CAF,LMAX,LXM,E,VPI)
C
      INTEGER LXM(LMMAX)
      COMPLEX ALM(LMMAX),CSUM,RSUM,AK,CI
      COMPLEX TSTORE(2,LMN,NT),CAF(LMAX+1)
      COMPLEX APLUS(NT), AMINUS(NT)
C
      CI=CMPLX(0.0,1.0)
      AK=-0.5/SQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
C
      IN=(NLAYER-1)*LMMAX
      K=0
      DO 1 L=0,LMAX
      DO 1 M=-L,L
      K=K+1
      KLM=LXM(K)
      CSUM=CMPLX(0.0,0.0)
      DO 2 JGP=1,NT
2     CSUM = CSUM + TSTORE(1,IN+KLM,JGP)*APLUS(JGP)+
     &TSTORE(2,IN+KLM,JGP)*AMINUS(JGP)
C
      IF( CABS(CAF(L+1)).GE.1.0E-06) THEN
      CSUM=(CI**L)*CSUM/(CAF(L+1)*AK)
      ELSE
      CSUM=CMPLX(0.0,0.0)
      ENDIF
C
      ALM(K)=CSUM

 1    CONTINUE
C
      RETURN
      END

!  TODO: Alex: remove FMAT_OLD since FMAT now exists?
C-----------------------------------------------------------------------
COMMENT FMAT_OLD CALCULATES THE VALUES OF THE SUM FLMS(JS,LM),
C       OVER LATTICE POINTS OF EACH SUBLATTICE JS, WHERE
C       LM=(0,0),(1,-1),(1,1),.....
C       NOTE  FOR ODD (L+M), FLMS IS ZERO
C  AUTHOR  PENDRY
C  DIMENSIONS 80 AND 41 REQUIRE LMAX.LE.10. AND IDEG.LE.4
C   FLMS= OUTPUT LATTICE SUMS.
C   V,JJS= INPUT FROM SUBROUTINE SLIND.
C   NL= NO. OF SUBLATTICES
C   NLS= ACTUAL NO. OF SUBLATTICE SUMS DESIRED (E.G. 1 OR NL).
C   DCUT= CUTOFF DISTANCE FOR LATTICE SUM.
C   IDEG= DEGREE OF SYMMETRY OF LATTICE (IDEG-FOLD ROTATION AXIS).
C    NOTE  DO NOT USE IDEG=1. IDEG=3 PREFERABLE OVER IDEG=6.
C   LMAX= LARGEST VALUE OF L.
C   KLM=(2*LMAX+1)*(2*LMAX+2)/2.
C  IN COMMON BLOCKS
C   E= CURRENT ENERGY.
C   AK  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C   VPI= IMAGINARY PART OF ENERGY.
C   BR1,BR2,RAR1,RAR2,NL1,NL2  NOT USED.
C   AR1,AR2= BASIS VECTORS OF SUPERLATTICE.
COMMENT FMAT CALCULATES THE VALUES OF THE SUM FLMS(JS,LM),
C       OVER LATTICE POINTS OF EACH SUBLATTICE JS, WHERE
C       LM=(0,0),(1,-1),(1,1),.....
C       NOTE  FOR ODD (L+M), FLMS IS ZERO
C  AUTHOR  PENDRY
C  DIMENSIONS 80 AND 41 REQUIRE LMAX.LE.10. AND IDEG.LE.4
CVB These dimensions are no longer set excplicitly; only restriction
C   is IDEG <= 6 now, from ANC and ANS
C   FLMS= OUTPUT LATTICE SUMS.
C   V,JJS= INPUT FROM SUBROUTINE SLIND.
C   NL= NO. OF SUBLATTICES
C   NLS= ACTUAL NO. OF SUBLATTICE SUMS DESIRED (E.G. 1 OR NL).
C   DCUT= CUTOFF DISTANCE FOR LATTICE SUM.
C   IDEG= DEGREE OF SYMMETRY OF LATTICE (IDEG-FOLD ROTATION AXIS).
C   NOTE  DO NOT USE IDEG=1. IDEG=3 PREFERABLE OVER IDEG=6.
C   LMAX= LARGEST VALUE OF L.
C   KLM=(2*LMAX+1)*(2*LMAX+2)/2.
C  IN COMMON BLOCKS
C   E= CURRENT ENERGY.
C   AK  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C   VPI= IMAGINARY PART OF ENERGY.
C   BR1,BR2,RAR1,RAR2,NL1,NL2  NOT USED.
C   AR1,AR2= BASIS VECTORS OF SUPERLATTICE.
      SUBROUTINE  FMAT_OLD (FLMS, V, JJS, NL, NLS, DCUT, IDEG, LMAX,KLM,
     +                  SCC,SA)
      COMPLEX  FLMS, SCC, SA, RTAB, CZERO, CI, KAPPA, SC, SD, SE, Z,
     1ACS, ACC, RF
      COMPLEX SQRT,EXP,COS,SIN
      DIMENSION  FLMS(NL,KLM), V(NL,2), JJS(NL,IDEG), BR1(2)
      DIMENSION  BR2(2), SCC(IDEG,4*LMAX+1),
     +           SA(2*LMAX*IDEG), ANC(6), ANS(6), RTAB(4)
      DIMENSION  AK(2),AR1(2),AR2(2),RAR1(2),RAR2(2),R(2)
      COMMON  E, AK, VPI
      COMMON  /SL/BR1, BR2, AR1, AR2, RAR1, RAR2, NL1, NL2
      PI = 3.14159265
      CZERO = CMPLX(0.0,0.0)
      CI = CMPLX(0.0,1.0)
      KAPPA = CMPLX(2.0 * E, - 2.0 * VPI + 0.000001)
      KAPPA = SQRT(KAPPA)
      AG=SQRT(AK(1)*AK(1)+AK(2)*AK(2))
C
COMMENT ANC,ANS AND SA ARE PREPARED TO BE USED IN THE SUM
C       OVER SYMMETRICALLY RELATED SECTORS OF THE LATTICE
      L2MAX = LMAX + LMAX
      LIM = L2MAX * IDEG
      LIML = L2MAX + 1
      ANG = 2.0 * PI/FLOAT(IDEG)
      D = 1.0
CDIR$ NOVECTOR                                                            170389
      DO 10 J = 1, IDEG
      ANC(J) = COS(D * ANG)
      ANS(J) = SIN(D * ANG)
   10 D = D + 1.0
CDIR$ VECTOR                                                              170389
      D = 1.0
      DO 20 J = 1, LIM
      SA(J) = EXP( - CI * D * ANG)
   20 D = D + 1.0
      DO 30 J = 1, NL
      DO 30 K = 1, KLM
   30 FLMS(J,K) = CZERO
COMMENT THE LATTICE SUM STARTS.THE SUM IS DIVIDED INTO ONE
C       OVER A SINGLE SECTOR,THE OTHER (IDEG-1) SECTORS
C       ARE RELATED BY SYMMETRY EXCEPT FOR FACTORS
C       INVOLVING THE DIRECTION OF R
C  THE RANGE OF SUMMATION IS LIMITED BY DCUT
      D = SQRT(AR1(1) * AR1(1) + AR1(2) * AR1(2))
      LI1=INT(DCUT/D)+2                                                 020780
      D = SQRT(AR2(1) * AR2(1) + AR2(2) * AR2(2))
      LI2=INT(DCUT/D)+2                                                 020780
      DCUT2=DCUT*DCUT                                                   020780
C  ONE SUBLATTICE AT A TIME IS TREATED IN THE FIRST SECTOR
      DO 160 JS = 1, NLS
      LI11 = LI1
      LI22 = LI2
      ASST = 0.0
      ADD = 1.0
      ANT =  - 1.0
      ADR1=V(JS,1)*COS(V(JS,2))
      ADR2=V(JS,1)*SIN(V(JS,2))
C  SHIFT POINT ADR1,2 BY MULTIPLES OF AR1 AND AR2                       020780
C  INTO SUPERLATTICE UNIT CELL NEAR ORIGIN, DEFINED BY THE LIMITS       020780
C  (-0.001 TO 0.999)*AR1 AND (0.001 TO 1.001)*AR2                       020780
      DET=AR1(1)*AR2(2)-AR1(2)*AR2(1)                                      .
      B1=(ADR1*AR2(2)-ADR2*AR2(1))/DET
      B2=(AR1(1)*ADR2-AR1(2)*ADR1)/DET
      BP1=AMOD(B1,1.)
      IF (BP1.LT.-.001) BP1=BP1+1.
      IF (BP1.GT.0.999) BP1=BP1-1.
      BP2=AMOD(B2,1.)
      IF (BP2.LT.0.001) BP2=BP2+1.
      IF (BP2.LT.0.001) BP2=BP2+1.
      ADR1=BP1*AR1(1)+BP2*AR2(1)                                           .
      ADR2=BP1*AR1(2)+BP2*AR2(2)                                        020780
   40 AST =  - 1.0
   60 AN1 = ANT
      DO 130 I1 = 1, LI11
      AN1 = AN1 + ADD
      AN2 = AST
      DO 130 I2 = 1, LI22
      AN2 = AN2 + 1.0
COMMENT R=THE CURRENT LATTICE VECTOR IN THE SUM
C       AR=MOD(R)
C       RTAB(1)=-EXP(I*FI(R))
      R(1)=AN1*AR1(1)+AN2*AR2(1)+ADR1
      R(2)=AN1*AR1(2)+AN2*AR2(2)+ADR2
      AR=R(1)*R(1)+R(2)*R(2)                                            020780
      IF (AR.GT.DCUT2) GO TO 130                                        020780
      AR=SQRT(AR)                                                       020780
      RTAB(1) =  - CMPLX(R(1)/AR,R(2)/AR)
      ABC = 1.0
      ABB = 0.0
      IF (AG-1.0E-4)  80, 80, 70
   70 ABC = (AK(1) * R(1) + AK(2) * R(2))/(AG * AR)
      ABB = ( - AK(2) * R(1) + AK(1) * R(2))/(AG * AR)
   80 SC = CI * AG * AR
COMMENT SCC CONTAINS FACTORS IN THE SUMMATION DEPENDENT ON
C       THE DIRECTION OF R. CONTRIBUTIONS FROM SYMMETRICALLY
C       RELATED SECTORS CAN BE GENERATED SIMPLY AND ARE
C       ACCUMULATED FOR EACH SECTOR, INDEXED BY THE SUBSCRIPT J.
C       THE SUBSCRIPT M IS ORDERED M=(-L2MAX),(-L2MAX+1)....
C       (+L2MAX)
      DO 90 J = 1, IDEG
      AD = ABC * ANC(J) - ABB * ANS(J)
      SD = EXP(SC * AD)
      SCC(J,LIML) = SD
      MJ = 0
      SE = RTAB(1)
CDIR$ SHORTLOOP                                                           170389
      DO 90 M = 1, L2MAX
      MJ = MJ + J
      MP = LIML + M
      MM = LIML - M
      SCC(J,MP) = SD * SA(MJ)/SE
      SCC(J,MM) = SD * SE/SA(MJ)
   90 SE = SE * RTAB(1)
      Z = AR * KAPPA
      ACS = SIN(Z)
      ACC = COS(Z)
COMMENT RTAB(3)=SPHERICAL HANKEL FUNCTION OF THE FIRST KIND,L=0
C       RTAB(4)=SPHERICAL HANKEL FUNCTION OF THE FIRST KIND,L=1
      RTAB(3) = (ACS - CI * ACC)/Z
      RTAB(4) = ((ACS/Z - ACC) - CI * (ACC/Z + ACS))/Z
      AL = 0.0
COMMENT THE SUMMATION OVER FACTORS INDEPENDENT OF THE
C       DIRECTION OF R IS ACCUMULATED IN FLM, FOR EACH
C       SUBLATTICE INDEXED BY SUBSCRIPT JSP.  THE SECOND
C       SUBSCRIPT ORDERS L AND M AS  (0,0),(1,-1),(1,1),(2,-2),
C       (2,0),(2,2)...
      JF = 1
      DO 120 JL = 1, LIML
      RF = RTAB(3) * CI
      JM = L2MAX + 2 - JL
      DO 110 KM = 1, JL
C  CONSIDER THE CORRESPONDING LATTICE POINTS IN THE OTHER SECTORS AND
C  GIVE THEIR CONTRIBUTION TO THE APPROPRIATE SUBLATTICE
CDIR$ NOVECTOR                                                            170389
      DO 100 J = 1, IDEG
      JSP = JJS(JS,J)
  100 FLMS(JSP,JF) = FLMS(JSP,JF) + SCC(J,JM) * RF
CDIR$ VECTOR                                                              170389
      JF = JF + 1
      JM = JM + 2
  110 CONTINUE
COMMENT SPHERICAL HANKEL FUNCTIONS FOR HIGHER L ARE
C       GENERATED BY RECURRENCE RELATIONS
      ACS = (2.0 * AL + 3.0) * RTAB(4)/Z - RTAB(3)
      RTAB(3) = RTAB(4)
      RTAB(4) = ACS
      AL = AL + 1.0
  120 CONTINUE
  130 CONTINUE
COMMENT SPECIAL TREATMENT IS REQUIRED IF IDEG=2
C  TWO SECTORS REMAIN TO BE SUMMED OVER
      IF (IDEG-2)  140, 140, 160
  140 IF (ASST)  150, 150, 160
  150 ASST = 1.0
      ADD =  - 1.0
      AST=-1.0                                                          020780
      IF (BP2.GT.0.999) AST=-2.                                         020780
      ANT = 0.0
      GO TO 60
  160 CONTINUE
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE FMAT CALCULATES THE VALUES OF THE SUM FLMS(JS,LM),
!       OVER LATTICE POINTS OF EACH SUBLATTICE JS, WHERE
!       LM=(0,0),(1,-1),(1,1),.....
!       NOTE  FOR ODD (L+M), FLMS IS ZERO
!  AUTHOR  PENDRY
!
!  Edited: 2021-09-23 MRiva. Cleanup: get rid of GOTOs, add ENDDOs
!
!  DIMENSIONS 80 AND 41 REQUIRE LMAX.LE.10. AND IDEG.LE.4
!   FLMS= OUTPUT LATTICE SUMS.
!   V,JJS= INPUT FROM SUBROUTINE SLIND.
!   NL= NO. OF SUBLATTICES
!   NLS= ACTUAL NO. OF SUBLATTICE SUMS DESIRED (E.G. 1 OR NL).
!   DCUT= CUTOFF DISTANCE FOR LATTICE SUM.
!   IDEG= DEGREE OF SYMMETRY OF LATTICE (IDEG-FOLD ROTATION AXIS).
!    NOTE  DO NOT USE IDEG=1. IDEG=3 PREFERABLE OVER IDEG=6.
!   LMAX= LARGEST VALUE OF L.
!   KLM=(2*LMAX+1)*(2*LMAX+2)/2.
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY.
!   AK  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!   VPI= IMAGINARY PART OF ENERGY.
!   BR1,BR2,RAR1,RAR2,NL1,NL2  NOT USED.
!   AR1,AR2= BASIS VECTORS OF SUPERLATTICE.
!
! COMMENT FMAT CALCULATES THE VALUES OF THE SUM FLMS(JS,LM),
!       OVER LATTICE POINTS OF EACH SUBLATTICE JS, WHERE
!       LM=(0,0),(1,-1),(1,1),.....
!       NOTE  FOR ODD (L+M), FLMS IS ZERO
!  AUTHOR  PENDRY
!  DIMENSIONS 80 AND 41 REQUIRE LMAX.LE.10. AND IDEG.LE.4
!VB These dimensions are no longer set explicitly; only restriction
!   is IDEG <= 6 now, from ANC and ANS
!   FLMS= OUTPUT LATTICE SUMS.
!   V,JJS= INPUT FROM SUBROUTINE SLIND.
!   NL= NO. OF SUBLATTICES
!   NLS= ACTUAL NO. OF SUBLATTICE SUMS DESIRED (E.G. 1 OR NL).
!   DCUT= CUTOFF DISTANCE FOR LATTICE SUM.
!   IDEG= DEGREE OF SYMMETRY OF LATTICE (IDEG-FOLD ROTATION AXIS).
!   NOTE  DO NOT USE IDEG=1. IDEG=3 PREFERABLE OVER IDEG=6.
!   LMAX= LARGEST VALUE OF L.
!   KLM=(2*LMAX+1)*(2*LMAX+2)/2.
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY.
!   AK  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!   VPI= IMAGINARY PART OF ENERGY.
!   BR1,BR2,RAR1,RAR2,NL1,NL2  NOT USED.
!   AR1,AR2= BASIS VECTORS OF SUPERLATTICE.
      SUBROUTINE  FMAT (FLMS, V, JJS, NL, NLS, DCUT, IDEG, LMAX, KLM,
     &                  SCC, SA)
      INTEGER ASST ! Used as flag when IDEG=2 to run the sum twice
      COMPLEX  FLMS, SCC, SA, RTAB, CZERO, CI, KAPPA, SC, SD, SE, Z,
     &         ACS, ACC, RF
      COMPLEX SQRT, EXP, COS, SIN
      DIMENSION  FLMS(NL,KLM), V(NL,2), JJS(NL, IDEG), BR1(2)
      DIMENSION  BR2(2), SCC(IDEG, 4 * LMAX + 1),
     &           SA(2*LMAX*IDEG), ANC(6), ANS(6), RTAB(4)
      DIMENSION  AK(2),AR1(2),AR2(2),RAR1(2),RAR2(2),R(2)
      COMMON  E, AK, VPI
      COMMON  /SL/BR1, BR2, AR1, AR2, RAR1, RAR2, NL1, NL2
      PI = 3.14159265
      CZERO = CMPLX(0.0,0.0)
      CI = CMPLX(0.0,1.0)
      KAPPA = CMPLX(2.0 * E, - 2.0 * VPI + 0.000001)
      KAPPA = SQRT(KAPPA)
      AG = SQRT(AK(1) * AK(1) + AK(2) * AK(2))
!
!     ANC, ANS and SA are prepared to be used in the sum
!     over symmetrically related sectors of the lattice
      L2MAX = LMAX + LMAX
      LIM = L2MAX * IDEG
      LIML = L2MAX + 1
      ANG = 2.0 * PI/FLOAT(IDEG)
      D = 1.0
CDIR$ NOVECTOR                                                            170389
      DO J = 1, IDEG
        ANC(J) = COS(D * ANG)
        ANS(J) = SIN(D * ANG)
        D = D + 1.0
      ENDDO
CDIR$ VECTOR                                                              170389
      D = 1.0
      DO J = 1, LIM
        SA(J) = EXP( - CI * D * ANG)
        D = D + 1.0
      ENDDO

      DO K = 1, KLM
        DO J = 1, NL
          FLMS(J, K) = CZERO
        ENDDO
      ENDDO
!       THE LATTICE SUM STARTS. THE SUM IS DIVIDED INTO ONE
!       OVER A SINGLE SECTOR, THE OTHER (IDEG-1) SECTORS
!       ARE RELATED BY SYMMETRY EXCEPT FOR FACTORS
!       INVOLVING THE DIRECTION OF R
!     THE RANGE OF SUMMATION IS LIMITED BY DCUT
      D = SQRT(AR1(1) * AR1(1) + AR1(2) * AR1(2))
      LI1 = INT(DCUT/D) + 2                                             020780
      LI11 = LI1  ! This is just a copy, never edited.
      D = SQRT(AR2(1) * AR2(1) + AR2(2) * AR2(2))
      LI2 = INT(DCUT/D) + 2                                             020780
      DCUT2 = DCUT * DCUT                                               020780
      DET = AR1(1) * AR2(2) - AR1(2) * AR2(1)                                      .
!     ONE SUBLATTICE AT A TIME IS TREATED IN THE FIRST SECTOR
      DO JS = 1, NLS
        LI22 = LI2
        ASST = 0 ! could be made an integer, used as flag only
        ADD = 1.0
        ANT =  - 1.0
        ADR1 = V(JS, 1) * COS(V(JS, 2))
        ADR2 = V(JS, 1) * SIN(V(JS, 2))
!       SHIFT POINT ADR1,2 BY MULTIPLES OF AR1 AND AR2                       020780
!       INTO SUPERLATTICE UNIT CELL NEAR ORIGIN, DEFINED BY THE LIMITS       020780
!       (-0.001 TO 0.999)*AR1 AND (0.001 TO 1.001)*AR2                       020780
        B1 = (ADR1 * AR2(2) - ADR2 * AR2(1))/DET
        B2 = (AR1(1) * ADR2 - AR1(2) * ADR1)/DET
        BP1 = AMOD(B1, 1.)
        IF (BP1.LT.-.001) THEN
          BP1 = BP1+1.
        ELSEIF (BP1.GT.0.999) THEN
          BP1 = BP1-1.
        ENDIF
        BP2 = AMOD(B2, 1.)
        IF (BP2.LT.0.001) THEN
          BP2=BP2+1.
        ENDIF
        ADR1 = BP1 * AR1(1) + BP2 * AR2(1)                                           .
        ADR2 = BP1 * AR1(2) + BP2 * AR2(2)                                020780
        AST =  - 1.0
   60   CONTINUE  ! Continue from here once again later if IDEG <= 2
        AN1 = ANT
        DO I1 = 1, LI11    ! First in-plane direction
          AN1 = AN1 + ADD
          AN2 = AST
          DO I2 = 1, LI22  ! Second in-plane direction
            AN2 = AN2 + 1.0
!           R=THE CURRENT LATTICE VECTOR IN THE SUM
!           AR=MOD(R)
!           RTAB(1)=-EXP(I*FI(R))
            R(1) = AN1 * AR1(1) + AN2 * AR2(1) + ADR1
            R(2) = AN1 * AR1(2) + AN2 * AR2(2) + ADR2
            AR = R(1) * R(1) + R(2) * R(2)                                020780
            IF (AR.LE.DCUT2) THEN                                         020780
!             Lattice sum goes on, not yet reached the limit
              AR = SQRT(AR)                                               020780
              RTAB(1) =  - CMPLX(R(1)/AR, R(2)/AR)
              ABC = 1.0
              ABB = 0.0
              IF (AG.GT.1.0E-4) THEN  ! off-normal beam
                ABC = (AK(1) * R(1) + AK(2) * R(2))/(AG * AR)
                ABB = ( - AK(2) * R(1) + AK(1) * R(2))/(AG * AR)
              ENDIF
              SC = CI * AG * AR
!             SCC CONTAINS FACTORS IN THE SUMMATION DEPENDENT ON
!             THE DIRECTION OF R. CONTRIBUTIONS FROM SYMMETRICALLY
!             RELATED SECTORS CAN BE GENERATED SIMPLY AND ARE
!             ACCUMULATED FOR EACH SECTOR, INDEXED BY THE SUBSCRIPT J.
!             THE SUBSCRIPT M IS ORDERED M=(-L2MAX),(-L2MAX+1)....
!             (+L2MAX)
              DO J = 1, IDEG
                AD = ABC * ANC(J) - ABB * ANS(J)
                SD = EXP(SC * AD)
                SCC(J, LIML) = SD
                MJ = 0
                SE = RTAB(1)
CDIR$ SHORTLOOP                                                           170389
                DO M = 1, L2MAX
                  MJ = MJ + J
                  MP = LIML + M
                  MM = LIML - M
                  SCC(J, MP) = SD * SA(MJ)/SE
                  SCC(J, MM) = SD * SE/SA(MJ)
                  SE = SE * RTAB(1)
                ENDDO   ! M loop, i.e., over angular momentum L**2
              ENDDO     ! J loop, i.e., over rotational symmetry order IDEG
              Z = AR * KAPPA
              ACS = SIN(Z)
              ACC = COS(Z)
!             RTAB(3)=SPHERICAL HANKEL FUNCTION OF THE FIRST KIND,L=0
!             RTAB(4)=SPHERICAL HANKEL FUNCTION OF THE FIRST KIND,L=1
              RTAB(3) = (ACS - CI * ACC)/Z
              RTAB(4) = ((ACS/Z - ACC) - CI * (ACC/Z + ACS))/Z
              AL = 0.0
!             THE SUMMATION OVER FACTORS INDEPENDENT OF THE
!             DIRECTION OF R IS ACCUMULATED IN FLM, FOR EACH
!             SUBLATTICE INDEXED BY SUBSCRIPT JSP.  THE SECOND
!             SUBSCRIPT ORDERS L AND M AS  (0,0),(1,-1),(1,1),(2,-2),
!             (2,0),(2,2)...
              JF = 1
              DO JL = 1, LIML
                RF = RTAB(3) * CI
                JM = L2MAX + 2 - JL
                DO KM = 1, JL
!                 CONSIDER THE CORRESPONDING LATTICE POINTS IN THE OTHER SECTORS AND
!                 GIVE THEIR CONTRIBUTION TO THE APPROPRIATE SUBLATTICE
CDIR$ NOVECTOR                                                            170389
                  DO J = 1, IDEG
                    JSP = JJS(JS, J)
                    FLMS(JSP, JF) = FLMS(JSP, JF) + SCC(J, JM) * RF
                  ENDDO
CDIR$ VECTOR                                                              170389
                  JF = JF + 1
                  JM = JM + 2
                ENDDO  ! loop over angular momentum M
!               SPHERICAL HANKEL FUNCTIONS FOR HIGHER L ARE
!               GENERATED BY RECURRENCE RELATIONS
                ACS = (2.0 * AL + 3.0) * RTAB(4)/Z - RTAB(3)
                RTAB(3) = RTAB(4)
                RTAB(4) = ACS
                AL = AL + 1.0
              ENDDO   ! Again sum over angular momentum L**2
            ENDIF   ! ||R||**2 <= DCUT**2
          ENDDO  ! loop I2, i.e., second in-plane direction
        ENDDO    ! loop I1, i.e., first in-plane direction

!       SPECIAL TREATMENT IS REQUIRED IF IDEG=2
!       TWO SECTORS REMAIN TO BE SUMMED OVER
        IF ((IDEG.LE.2).AND.(ASST.EQ.0)) THEN
          ASST = 1
          ADD =  - 1.0
          AST = -1.0                                                      020780
          IF (BP2.GT.0.999) AST = -2.                                     020780
          ANT = 0.0
          GO TO 60
       ENDIF
      ENDDO ! JS loop, i.e., over sublattices
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE GHD COMPUTES DIRECT LATTICE SUMS FOR (LM)-SPACE
C  PROPAGATORS GH BETWEEN TWO SUBPLANES (HAVING BRAVAIS LATTICES) OF
C  A COMPOSITE LAYER. THE SUBPLANES MAY BE COPLANAR.
C  FOR QUANTITIES NOT EXPLAINED BELOW SEE GHMAT, MPERTI OR MTINV
C   IZ= SERIAL NO. OF CURRENT INTERPLANAR VECTOR DRL.
C   IS=1 FOR PROPAGATION FROM FIRST TO SECOND SUBPLANE.
C   IS=2 FOR PROPAGATION FROM SECOND TO FIRST SUBPLANE.
C   GH= OUTPUT INTERPLANAR PROPAGATOR.
C   S= WORKING SPACE (LATTICE SUM).
C   LMS= (2*LMAX+1)**2.
C   Y= WORKING SPACE (SPHERICAL HARMONICS).
C   L2M= 2*LMAX+1.
C   K0= COMPLEX MAGNITUDE OF WAVEVECTOR.
C   DCUT= CUTOFF RADIUS FOR LATTICE SUM.
C   CAA= CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CAAA.
C   NCAA= NO. OF C.-G. COEFFICIENTS IN CAA.
C   LXM= PERMUTATION OF (LM) SEQUENCE FROM SUBROUTINE LXGENT.
C   LAY,PQ- SEE GHMAT.
C   H is auxiliary array HGHD used only here, carried through for var. dimensions
      SUBROUTINE GHD(IZ,IS,GH,LMG,LMMAX,S,LMS,Y,L2M,DRL,NLAY2,
     1               K0,DCUT,CAA,NCAA,LXM,LAY,PQ,H)                       111181
      COMPLEX H(L2M)
      COMPLEX GH(LMG,LMMAX),Y(L2M,L2M),S(LMS)
      COMPLEX RU,CI,CZ,K0,FF,Z,Z1,Z2,Z3,ST
      DIMENSION DRL(NLAY2,3),ARA1(2),ARA2(2),ARB1(2),ARB2(2)
      DIMENSION RBR1(2),RBR2(2),CAA(NCAA)
      DIMENSION V(3),LXM(LMMAX),PQ(2)                                     111181
      COMMON E,AK2,AK3,VPI
      COMMON /SL/ ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2

      RU=(1.0,0.0)
      CZ=(0.0,0.0)
      CI=(0.0,1.0)
      DCUT2=DCUT*DCUT

      DO I=1,LMS
        S(I)=CZ
      ENDDO

      V(1)=DRL(IZ,1)                                                    100589
      V(2)=DRL(IZ,2)                                                    100589
      V(3)=DRL(IZ,3)                                                    100589
C  TURN INTERPLANAR VECTOR AROUND IF IS=2
      IF (IS.EQ.2) THEN
        V(1) = -V(1)                                                      100589
        V(2) = -V(2)                                                      100589
        V(3) = -V(3)                                                      100589
      ENDIF
C
C  START OF TWO 1-DIMENSIONAL LATTICE LOOPS FOR SUMMATION OVER 1 QUADRAN
40    NUMR=0
      JJ1=0
50    JJ1=JJ1+1
      JJ2=0
60    JJ2=JJ2+1
      NOR=0
      J1=JJ1-1
      J2=JJ2-1
C  START OF LOOP OVER QUADRANTS !TODO clean up: remove statement numbers
      DO 140 KK=1,4
        IF (KK.EQ.1) THEN
70        NR1=J1
          NR2=J2
        ELSEIF (KK.EQ.2) THEN
80        IF (J1.EQ.0.AND.J2.EQ.0) GO TO 150
          IF (J2.EQ.0) GO TO 140
          NR1=J1
          NR2=-J2
        ELSEIF (KK.EQ.3) THEN
90        IF (J1.EQ.0) GO TO 140
          NR1=-J1
          NR2=J2
        ELSE
100       IF (J1.EQ.0.OR.J2.EQ.0) GO TO 140
          NR1=-J1
          NR2=-J2
        ENDIF
C
 110    IF (LAY.EQ.1) THEN                                                111181
          PX=NR1*ARB1(1)+NR2*ARB2(1)                                          .
          PY=NR1*ARB1(2)+NR2*ARB2(2)
        ELSE
          PX=NR1*ARA1(1)+NR2*ARA2(1)                                          .
          PY=NR1*ARA1(2)+NR2*ARA2(2)                                          .
        ENDIF
        X1=(PX+V(2))**2+(PY+V(3))**2                                      111181

C  CUTOFF OF LATTICE SUMMATION AT RADIUS DCUT
        IF (X1.LE.DCUT2) THEN
          NOR=1
          NUMR=NUMR+1
          Z1=EXP(-CI*(PX*(AK2+PQ(1))+PY*(AK3+PQ(2))))                       111181
          X2=SQRT(X1+V(1)**2)
          X1=SQRT(X1)
          Z2=K0*X2
          Z=CMPLX(V(1)/X2,0.0)
          FY=0.0
          IF (ABS(X1).GE.1.E-6) THEN
            CFY=(PX+V(2))/X1
            IF (ABS(ABS(CFY)-1.).GT.1.E-6) THEN
              FY= ACOS(CFY)
            ELSE
              IF (CFY.LT.0.0) FY=3.14159265
            ENDIF
            IF ((PY+V(3)).LT.0.0) FY=-FY
          ENDIF
C  COMPUTE REQUIRED BESSEL FUNCTIONS AND SPHERICAL HARMONICS
 118      CALL SB(Z2,H,L2M)
          CALL SH(L2M,Z,FY,Y)
C
          ST=RU
          DO 130 L=1,L2M
            ST=ST*CI
            Z3=ST*H(L)*Z1
            L1=L*L-L
C  Compute lattice sum into S                                                 100589
            DO 120 M=1,L                                                      100589
              S(L1+M)=S(L1+M)+Z3*Y(L,M)                                       100589
120         CONTINUE                                                          100589
            DO 121 M=1,L-1                                                    100589
              S(L1-M+1)=S(L1-M+1)+Z3*Y(M,L)*(-1.)**MOD(M,2)                   100589
121         CONTINUE                                                          100589
130       CONTINUE
        ENDIF
140   CONTINUE  ! End of KK loop over quadrants

!  TODO: replace GO TO
150   IF (NOR.EQ.1) GO TO 60
      IF (JJ2.EQ.1) GO TO 160
      IF (JJ2.NE.1) GO TO 50
160   CONTINUE

C
C  PRINT NUMBER OF LATTICE POINTS USED IN SUMMATION
C     WRITE(6,170)NUMR
C170   FORMAT(28H NO. OF LATT. PTS. FOR GH   ,I5)
      FF=-8.0*3.14159265*K0

C  USE SUBROUTINE GHSC TO MULTIPLY LATTICE SUM INTO CLEBSCH-GORDON
C  COEFFICIENTS
      CALL GHSC(IZ,IS,GH,LMG,LMMAX,S,LMS,CAA,NCAA,FF,LXM,NLAY2)
      RETURN
      END

C-----------------------------------------------------------------------
C  Subroutine GHMAT computes (LM)-space interplanar propagators GH for
C  the composite layer treated by subroutine MTINV. Depending on the
C  interplanar spacing, either a reciprocal-space summation or a direct
C  space summation is performed (The latter in routine GHD). Using
C  information provided by subroutine SRTLAY, GHMAT reuses (does not
C  recompute) existing values of GH provided previously, computes new values
C  not previously produced and copies these new values into GH's that are
C  identical (instead of computing the latter seperately). Ordering of the
C  elements for each GH is arranged as follows for elemenst referenced by
C  the following (LM) pairings
C  (00),(2-2),(20),(22),(4-4),(4-2),...,(1-1),(11),(3-3),(31),...
C  ,(10),(3-2),(30),(32),(5-4),...,(2-1),(21),(4-3),(4-1),(41),...
C  The GH(I,J)(I,J=Subplane indices) are stacked on top of each other
C  in columnar matrix GH, so that GH(I,J) is found in the K-th position
C  from the top, where K=MGH(I,J).MGH(I,J).
C
C Parameter List;
C ===============
C
C  GH           =   MATRIX CONTAINING ALL INTERPLANAR PROPAGATORS FOR THE
C                   COMPOSITE LAYER CONSIDERED BY MPERTI OR MTINV. THE
C                   INDIVIDUAL MATRICES ARE STACKED VERTICALLY IN GH.
C  LMG          =   2*NLAY2R*LMMAX.
C  LMMAX        =   (LMAX+1)**2.
C  MGH          =   MATRIX CONTAINING KEY TO POSITION OF INDIVIDUAL GH^S IN THE
C                   MATRIX GH  MGH(I,J) IS SEQUENCE NUMBER OF GH(I,J) IN
C                   COLUMNAR MATRIX GH.
C  NLAY          =  NO. OF SUBPLANES CONSIDERED.
C  NUGH          =  LIST OF THOSE GH^S THAT MUST BE COMPUTED.
C  NGEQ          =  LIST OF THOSE GH^S THAT CAN BE COPIED FROM FRESHLY COMPUTED
C                   GH^S.
C  NGOL          =  LIST OF THOSE GH^S THAT CAN BE COPIED FROM PREVIOUSLY
C                   PRODUCED GH^S.
C  NLAY2         =  NLAY*(NLAY-1)/2.
C  TST           =  INPUT QUANTITY FOR DETERMINING NO. OF POINTS REQUIRED IN
C                   RECIPROCAL LATTICE SUMMATION.
C  TEST,Y1,Y,S   =  WORKING SPACE.
C  L2M           =  2*LMAX+1.
C  LM            =  LMAX+1.
C  LMS           =  (2*LMAX+1)**2.
C  DRL           =  SET OF INTERPLANAR VECTORS.
C  TV            =  AREA OF UNIT CELL OF EACH SUBPLANE.
C  LXM           =  PERMUTATION OF (LM) SEQUENCE.
C  LEV           =  (LMAX+1)*(LMAX+2)/2.
C  DCUT          =  CUTOFF RADIUS FOR LATTICE SUMMATION.
C  CAA           =  CLEBSCH-GORDAN COEFFICIENTS.
C  NCAA          =  NO. OF CLEBSCH-GORDAN COEFFICIENTS.
C  LAY         =   1 IF CURRENT CALCULATION REFERS TO OVERLAYER.
C                  0 IF CURRENT CALCULATION REFERS TO SUBSTRATE.
C  PQ          =   ADDITIONAL OFFSET INTRODUCED FOR UNKNOWN PURPOSES
C  HGHD          = working space for subroutine GHD
C
C In Common Blocks;
C =================
C
C  E             =  CURRENT ENERGY.
C  AK2,AK3       =  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C  VPI           =  IMAGINARY PART OF ENERGY.
C  ARA1,ARA2     =  BASIS VECTORS OF SUBSTRATE LAYER LATTICE.
C  ARB1,ARB2     =  BASIS VECTORS OF SUPERLATTICE.
C  RBR1,RBR2     =  RECIPROCAL LATTICE OF SUPERLATTICE.
C  NL1,NL2       =  SUPERLATTICE CHARACTERIZATION CODES.
C
C =========================================================================
C
      SUBROUTINE GHMAT(GH,LMG,LMMAX,MGH,NLAY,NUGH,NGEQ,NGOL,NLAY2,TST,
     & TEST,Y1,L2M,Y,LM,S,LMS,DRL,TV,LXM,LEV,DCUT,CAA,NCAA,LAY,HGHD,PQ)
C
      DIMENSION MGH(NLAY,NLAY),NUGH(NLAY2),NGEQ(NLAY2),NGOL(NLAY2)
      DIMENSION RBR1(2),RBR2(2),DRL(NLAY2,3),TEST(NLAY2),LXM(LMMAX)
      DIMENSION ARA1(2),ARA2(2),ARB1(2),ARB2(2),CAA(NCAA)
      DIMENSION RXR1(2),RXR2(2),PQ(2)
      COMPLEX GH(LMG,LMMAX),Y(LM,LM),Y1(L2M,L2M),S(LMS)
      COMPLEX CI,CZ,KPRG,K0,Z,T1,T1P,T2,T2P,T3,BS,CS,CFAC
      COMPLEX HGHD(L2M)
C
      COMMON E,AK2,AK3,VPI
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
C
      CZ=(0.0,0.0)
      CI=(0.0,1.0)
      IF (LAY.EQ.1) THEN  ! overlayer
         RXR1(1)=RBR1(1)
         RXR1(2)=RBR1(2)
         RXR2(1)=RBR2(1)
         RXR2(2)=RBR2(2)
      ELSE                ! bulk
         PI=3.14159265
         ATV=2.0*PI/TV
         RXR1(1)=ARA2(2)*ATV
         RXR1(2)=-ARA2(1)*ATV
         RXR2(1)=-ARA1(2)*ATV
         RXR2(2)=ARA1(1)*ATV
      ENDIF
      K0=SQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      CFAC=-16.0*(3.14159265)*(3.14159265)*CI/TV
      NLYLM=NLAY2*LMMAX
C
C  COPY OLD VALUES OF GH ONTO GH WHERE APPROPRIATE
C
      DO 6 I=1,NLAY2
         IF (NGOL(I).NE.0) THEN
            IS=(I-1)*LMMAX
            IT=(NGOL(I)-1)*LMMAX
            IU=IS+NLYLM
            IV=IT+NLYLM
            DO 1461 K=1,LMMAX
               DO 2 J=1,LMMAX
                  GH(IS+J,K)=GH(IT+J,K)
                  GH(IU+J,K)=GH(IV+J,K)
2              CONTINUE
1461        CONTINUE
         ENDIF
6     CONTINUE
C
C  INITIALIZE GH IN PREPARATION FOR NEW VALUES
C
      DO 8 I=1,NLAY2
         IF (NUGH(I).NE.0) THEN
            IS=(I-1)*LMMAX
            IT=IS+NLYLM
            DO 1462 K=1,LMMAX
               DO 7 J=1,LMMAX
                  GH(IS+J,K)=CZ
                  GH(IT+J,K)=CZ
7              CONTINUE
1462        CONTINUE
         ENDIF
8     CONTINUE
C
C  TSTS= ESTIMATED NO. OF POINTS IN DIRECT LATTICE SUM
C
      TSTS=DCUT*DCUT*3.14159265/TV
      DO 10 IZ=1,NLAY2
         IF (NUGH(IZ).NE.0) THEN
            IF (ABS(DRL(IZ,1)).GT.0.001) THEN
C
C  AKP2= ESTIMATED NO. OF POINTS IN RECIPROCAL LATTICE SUM
C
               FACT1=ALOG(TST)/DRL(IZ,1)
               AKP2=(2.0*E+FACT1*FACT1)*TV/(4.*3.1415926)
C
C  SKIP DIRECT LATTICE SUM, IF RECIPROCAL LATTICE SUM FASTER (BUT NUMBER
C  OF REC. LATT. POINTS IS TO BE RESTRICTED DUE TO A CONVERGENCE
C  PROBLEM)
C
               IF ((TSTS.GE.2.0*AKP2).AND.(AKP2.LT.80.0)) GOTO 10
            ENDIF
C
C  PRODUCE GH(I,J) AND GH(J,I) WITH DIRECT LATTICE SUM FOR TWO
C  PROPAGATION DIRECTIONS
C
            DO 9 K=1,2
               CALL GHD(IZ,K,GH,LMG,LMMAX,S,LMS,Y1,L2M,DRL,NLAY2,K0,
     &                  DCUT,CAA,NCAA,LXM,LAY,PQ,HGHD)
9           CONTINUE
C
C  GH(I,J) AND GH(J,I) NOW NO LONGER NEED TO BE COMPUTED
C
            NUGH(IZ)=0
         ENDIF
10    CONTINUE

      TSTS=0.0
      DO 1 I=1,NLAY2
!         IF (.NOT.((NUGH(I).EQ.0).OR.(ABS(DRL(I,1)).LE.0.001))) THEN
         IF ((NUGH(I).NE.0).AND.(ABS(DRL(I,1)).GT.0.001)) THEN
            DRL(I,1)=DRL(I,1)
C
C  TEST(I) WILL SERVE AS CUTOFF IN RECIPROCAL LATTICE SUM
C  Notice that the convergence of the planar sum is not
C  related to the convergence of the RFS scheme controlled by TST
C
CC            TEST(I)=(ALOG(TST)/DRL(I,1))*(ALOG(TST)/DRL(I,1))
            TEST(I)=8.*ABS(ALOG(.002)/DRL(I,1))
            TSTS=AMAX1(TEST(I),TSTS)
         ENDIF
1     CONTINUE

      IF (TSTS.GT.0.00001) THEN
C
C  START OF TWO 1-DIMENSIONAL SUMMATION LOOPS IN ONE QUADRANT OF
C  RECIPROCAL SPACE
C
         NCOUNT=0
         NUMG=0
         JJ1=0
1171     JJ1=JJ1+1
         JJ2=0
1172     JJ2=JJ2+1
         NOG=0
         J1=JJ1-1
         J2=JJ2-1
C
C  START OF LOOP OVER QUADRANTS
C
         DO 1370 KK=1,4
            IF (KK.EQ.2) THEN       ! 2nd quadrant
               IF (J1.EQ.0.AND.J2.EQ.0) GOTO 1380
               IF (J2.EQ.0) GOTO 1370
               NG1=J1
               NG2=-J2
            ELSEIF (KK.EQ.3) THEN   ! 3rd quadrant
               IF (J1.EQ.0) GOTO 1370
               NG1=-J1
               NG2=J2
            ELSEIF (KK.NE.4) THEN   ! 4th quadrant
               NG1=J1
               NG2=J2
            ELSEIF (J1.EQ.0.OR.J2.EQ.0) THEN
               GOTO 1370
            ELSE                    ! 1st quadrant
               NG1=-J1
               NG2=-J2
            ENDIF
C
C  CURRENT RECIPROCAL LATTICE POINT
C
            GX=NG1*RXR1(1)+NG2*RXR2(1)
            GY=NG1*RXR1(2)+NG2*RXR2(2)
            GKX=GX+AK2+PQ(1)
            GKY=GY+AK3+PQ(2)
            GK2=GKX*GKX+GKY*GKY
C
C  TEST FOR CUTOFF
C
            KPRG=SQRT(CMPLX(2.0*E-GK2,-2.0*VPI+0.000001))
            AKP2=AIMAG(KPRG)
            IF (AKP2.LE.(TSTS)) THEN
               NUMG=NUMG+1
               NOG=1
               Z=KPRG/K0
               FY=0.0
               IF (GK2.GT.1.0E-8) THEN
                  CFY=GKX/SQRT(GK2)
                  IF (ABS(ABS(CFY)-1.).LE.1.E-6) THEN
                     IF (CFY.LT.0.0) FY=3.14159265
                  ELSE
                     FY=ACOS(CFY)
                  ENDIF
                  IF (GKY.LT.0.0) FY=-FY
               ENDIF
C
C  FIND APPROPRIATE SPHERICAL HARMONICS
C
               CALL SH(LM,Z,FY,Y)
C
C  START OF LOOP OVER INTERPLANAR VECTORS
C
               DO 1360 IZ=1,NLAY2
C
C  SKIP IF NEW GH(I,J) NOT NEEDED
C
                  IF (NUGH(IZ).NE.0) THEN
C
C  SKIP IF THIS CONTRIBUTION OUTSIDE CUTOFF
C
                     IF (AKP2.LE.(TEST(IZ))) THEN
                        T2=GKX*DRL(IZ,2)+GKY*DRL(IZ,3)
                        T2P=-T2
                        T3=KPRG*ABS(DRL(IZ,1))
                        T1=(EXP(CI*(T2+T3))/KPRG)*CFAC
                        T1P=(EXP(CI*(T2P+T3))/KPRG)*CFAC
C
C  PRODUCE A LIMITED SET OF MATRIX ELEMENTS
C
                        DO 1350 I=1,LMMAX
                           IP=LXM(I)
                           NII=IP+(IZ-1)*LMMAX
                           NNI=NII+NLYLM
                           L1=INT(SQRT(FLOAT(I-1)+0.00001))
                           M1=I-L1-L1*L1-1
                           IF (M1.GT.0) THEN
                              BS=Y(M1,L1+1)
                           ELSE
                              BS=(-1.)**(MOD(M1,2))*Y(L1+1,-M1+1)
                           ENDIF
                           DO 1340 J=1,LMMAX
                              JP=LXM(J)
                              L2=INT(SQRT(FLOAT(J-1)+0.00001))
                              M2=J-L2-L2*L2-1
                              IF (I.EQ.J) THEN
                                 IF (M1.GT.0) GOTO 1340
                              ELSEIF (L1.LT.L2.OR.(L1.EQ.L2.AND.(IABS
     &                         (M1).LT.IABS(M2)))) THEN
                                 GOTO 1340
                              ENDIF
                              IF (M2.GE.0) THEN
                                 CS=Y(L2+1,M2+1)
                              ELSE
                                 CS=(-1.)**(MOD(M2,2))*Y(-M2,L2+1)
                              ENDIF
C
C  GH(I,J) AND GH(J,I) ARE PRODUCED TOGETHER
C
                              GH(NNI,JP)=GH(NNI,JP)+T1P*BS*CS
                              GH(NII,JP)=GH(NII,JP)+(-1.)**(MOD(L1+M1+
     &                         L2+M2,2))*T1*BS*CS
1340                       CONTINUE
1350                    CONTINUE
                     ENDIF
                  ENDIF
1360           CONTINUE
            ENDIF
1370     CONTINUE  ! end of quadrant loop

1380     IF (NOG.EQ.1) GOTO 1172
         IF (JJ2.NE.1) GOTO 1171
         IF (JJ2.EQ.1) NCOUNT=NCOUNT+1
         IF (NCOUNT.LE.3) GOTO 1172
C
C  FILL IN MISSING MATRIX ELEMENTS FROM ELEMENTS JUST PRODUCED, USING
C  SYMMETRY RELATIONS
C
         DO 1450 IZ=1,NLAY2
            IF (NUGH(IZ).NE.0) THEN
               DO 1440 I=1,LMMAX
                  IP=LXM(I)
                  NII=IP+(IZ-1)*LMMAX
                  NNI=NII+NLYLM
                  L1=INT(SQRT(FLOAT(I-1)+0.00001))
                  M1=I-L1-L1*L1-1
                  IM1=I-2*M1
                  IMP=LXM(IM1)
                  DO 1430 J=1,LMMAX
                     JP=LXM(J)
                     L2=INT(SQRT(FLOAT(J-1)+0.00001))
                     M2=J-L2-L2*L2-1
                     JM2=J-2*M2
                     JMP=LXM(JM2)
                     NJM=JMP+(IZ-1)*LMMAX
                     NNM=NJM+NLYLM
                     SIGN=(-1.)**MOD(M1+M2,2)
                     IF (I.EQ.J) THEN
                        IF (M1.LE.0) GOTO 1430
                     ELSEIF (L1.GE.L2.AND.(L1.NE.L2.OR.(IABS(M1).GE.
     &                       IABS(M2)))) THEN
                        GOTO 1430
                     ENDIF
                     GH(NNI,JP)=SIGN*GH(NNM,IMP)
                     GH(NII,JP)=SIGN*GH(NJM,IMP)
1430              CONTINUE
1440           CONTINUE
               NUGH(IZ)=0
            ENDIF
1450     CONTINUE
      ENDIF
C
C  PRINT NUMBER OF RECIPROCAL LATTICE POINTS (BEAMS) USED IN SUMMATION
C
C  COPY NEW MATRICES GH(I,J) ONTO OTHERS THAT ARE IDENTICAL
C
      DO 5 IZ=1,NLAY2
         IF (NGEQ(IZ).NE.0) THEN
            II=(IZ-1)*LMMAX
            IIN=II+NLYLM
            IN=(NGEQ(IZ)-1)*LMMAX
            INN=IN+NLYLM
            DO 1463 I=1,LMMAX
               DO 4 J=1,LMMAX
                  GH(II+I,J)=GH(IN+I,J)
                  GH(IIN+I,J)=GH(INN+I,J)
4              CONTINUE
1463        CONTINUE
         ENDIF
5     CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE GHSC COMPUTES THE (LM)-SPACE INTERPLANAR PROPAGATORS GH
C  FROM DIRECT LATTICE SUMS PRODUCED IN SUBROUTINE GHD AND CLEBSCH-
C  GORDON COEFFICIENTS FROM SUBROUTINE CAAA.
C  FOR QUANTITIES NOT EXPLAINED BELOW SEE SUBROUTINES GHD AND GHMAT.
C   IZ= SERIAL NO. OF CURRENT INTERPLANAR VECTOR DRL.
C   IS=1 FOR PROPAGATION FROM FIRST TO SECOND SUBPLANE.
C   IS=2 FOR PROPAGATION FROM SECOND TO FIRST SUBPLANE.
C   FF= PREFACTOR OF GH.
C   LXM= PERMUTATION OF (LM) SEQUENCE
      SUBROUTINE GHSC(IZ,IS,GH,LMG,LMMAX,S,LMS,CAA,NCAA,FF,LXM,NLAY2)
      COMPLEX GH(LMG,LMMAX),S(LMS),FF
      DIMENSION CAA(NCAA),LXM(LMMAX)
      II=1
      NLYLM=LMMAX*NLAY2
C  PRODUCE A LIMITED SET OF MATRIX ELEMENTS
      DO 1350 I=1,LMMAX
      IP=LXM(I)
      NII=IP+(IZ-1)*LMMAX
      IF (IS.EQ.2) NII=NII+NLYLM
      L1=INT(SQRT(FLOAT(I-1)+0.00001))                                  261179
      M1=I-L1-L1*L1-1
      DO 1340 J=1,LMMAX
      JP=LXM(J)
      L2=INT(SQRT(FLOAT(J-1)+0.00001))                                  261179
      M2=J-L2-L2*L2-1
      M3=M2-M1
      IL=IABS(L1-L2)
      IM=IABS(M3)
      LMIN=MAX0(IL,IM+MOD(IL+IM,2))
      LMAX=L1+L2
      LMIN=LMIN+1
      LMAX=LMAX+1
      IF (I-J) 1290,1280,1290
1280  IF (M1) 1300,1300,1340
1290  IF (L1.LT.L2.OR.(L1.EQ.L2.AND.(IABS(M1).LT.IABS(M2))))
     1GO TO 1340
1300  DO 203 ILA=LMIN,LMAX,2
      LA=ILA-1
      CC=CAA(II)
      II=II+1
203   GH(NII,JP)=GH(NII,JP)+CC*S(LA*LA+LA+M3+1)
      GH(NII,JP)=FF*GH(NII,JP)
1340  CONTINUE
1350  CONTINUE
C
C  FILL IN MISSING MATRIX ELEMENTS FROM ELEMENTS JUST PRODUCED
      DO 1440 I=1,LMMAX
      IP=LXM(I)
      NII=IP+(IZ-1)*LMMAX
      IF (IS.EQ.2) NII=NII+NLYLM
      L1=INT(SQRT(FLOAT(I-1)+0.00001))                                  261179
      M1=I-L1-L1*L1-1
      IM1=I-2*M1
      IMP=LXM(IM1)
      DO 1430 J=1,LMMAX
      JP=LXM(J)
      L2=INT(SQRT(FLOAT(J-1)+0.00001))                                  261179
      M2=J-L2-L2*L2-1
      JM2=J-2*M2
      JMP=LXM(JM2)
      NJM=JMP+(IZ-1)*LMMAX
      IF (IS.EQ.2) NJM=NJM+NLYLM
      SIGN=(-1.)**MOD(M1+M2,2)
      IF (I.NE.J) GO TO 1410
      IF (M1.LE.0) GO TO 1430
      GO TO 1420
1410  IF (L1.GE.L2.AND.(L1.NE.L2.OR.(IABS(M1).GE.IABS(M2))))
     1GO TO 1430
1420  GH(NII,JP)=SIGN*GH(NJM,IMP)
1430  CONTINUE
1440  CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE HEAD writes header for fd. output spectrum
C  enjoy the comments :-)

      SUBROUTINE HEAD(IFILE,TITLE,NPUN,NPU,KNT,SPQF,KSYM)

      CHARACTER*80 TITLE
      INTEGER NPUN,KNT
      INTEGER NPU,KSYM
      DIMENSION NPU(NPUN),KSYM(2,KNT)
      REAL SPQF
      DIMENSION SPQF(2,KNT)

  283 FORMAT(I5,2F10.5,I3)
C LH  Adjusted to allow for beams numbers until 99999
  284 FORMAT(I5,2F10.5,I3)

      WRITE(IFILE,'(A80)') TITLE
C  PUNCH A CARD WITH NPUN
      WRITE(IFILE,283)NPUN
      DO K=1,NPUN
      N=NPU(K)
C  PUNCH A CARD FOR EACH BEAM FOR WHICH PUNCH OUTPUT IS DESIRED (THIS
C  PASSES THE BEAM IDENTITY (0,0),(1,0),(0,1),ETC. TO THE DECK OF CARDS)
      WRITE(IFILE,284)N,(SPQF(I,N),I=1,2),KSYM(1,N)
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE LXGENT GENERATES THE NEEDED RELATIONSHIPS BETWEEN
C  DIFFERENT ORDERING SEQUENCES OF THE (L,M) PAIRS (L.LE.LMAX, ABS(M)
C  .LE.L). THREE ORDERINGS ARE CONSIDERED, THE 'NATURAL' ONE (N), THE
C  'COPLANAR' ONE (C) AND THE 'SYMMETRIZED' ONE (S). THEY ARE TABULATED
C  BELOW FOR THE CASE OF LMAX=4.
C
C     I=  0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2
C         1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
C
C  N  L=  0 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4
C     M=  0-1 0 1-2-1 0 1 2-3-2-1 0 1 2 3-4-3-2-1 0 1 2 3 4
C
C  C  L=  0 1 1 2 2 2 3 3 3 3 4 4 4 4 4 1 2 2 3 3 3 4 4 4 4
C     M=  0-1 1-2 0 2-3-1 1 3-4-2 0 2 4 0-1 1-2 0 2-3-1 1 3
C                   L+M=EVEN           *      L+M=ODD
C
C  S  L=  0 2 2 2 4 4 4 4 4 1 1 3 3 3 3 1 3 3 3 2 2 4 4 4 4
C     M=  0-2 0 2-4-2 0 2 4-1 1-3-1 1 3 0-2 0 2-1 1-3-1 1 3
C             L=EVEN       *  L=ODD    * L=ODD *  L=EVEN
C             M=EVEN       *  M=ODD    * M=EVEN*   M=ODD
C
C  TO DESCRIBE THE RELATIONSHIPS, A PARTICULAR PAIR (L,M) IN A PARTICU-
C  LAR SEQUENCE SHALL BE REPRESENTED HERE BY N(I), C(I) OR S(I) (I=1,
C  LMMAX)  E.G. N(10)=(3,-3). THE RELATIONSHIPS LX,LXI,LT,LXM GENERATED
C  BY LXGENT ARE NOW DEFINED BY
C   LX   N(LX(I))= S(I).
C   LXI  C(I)= S(LXI(I)) IF I.LE.LEV,
C              S(LXI(I)+LEV) IF I.GT.LEV.
C   LT   CSM(N(LT(I)))= S(I) WHERE CSM MEANS  CHANGE THE SIGN OF M.
C   LXM  N(I)= S(LXM(I)).
C
C   LX,LXI,LT,LXM= OUTPUT PERMUTATIONS OF (L,M) SEQUENCE.
C   LMAX= LARGEST VALUE OF L.
C   LMMAX= (LMAX+1)**2.
      SUBROUTINE LXGENT(LX,LXI,LT,LXM,LMAX,LMMAX)
      DIMENSION LX(LMMAX),LXM(LMMAX),LXI(LMMAX),LT(LMMAX)
      LEV=(LMAX+1)*(LMAX+2)/2
      LEE=(LMAX/2+1)**2
      LEO=((LMAX+1)/2+1)*((LMAX+1)/2)+LEE
      LOE=((LMAX-1)/2+1)**2+LEO
      LL=0
      L1=0
      LT1=0
      L=-1
1     L=L+1
      M=-L
2     LL=LL+1
      IF (MOD(L+M,2)-1) 3,6,3
3     LT1=LT1+1
      LT(LL)=LT1
      IF (MOD(L,2)-1) 4,5,4
4     L1=L1+1
      LX(L1)=LL
      GO TO 9
5     LEE=LEE+1
      LX(LEE)=LL
      GO TO 9
6     LEV=LEV+1
      LT(LL)=LEV
      IF (MOD(L,2)) 7,8,7
7     LEO=LEO+1
      LX(LEO)=LL
      GO TO 9
8     LOE=LOE+1
      LX(LOE)=LL
9     M=M+1
      IF (L-M) 10,2,2
10    IF (L-LMAX) 1,11,11
11    DO 12 L=1,LMMAX
      L1=LX(L)
      LT1=LT(L1)
12    LXI(LT1)=L
      L1=LEE+1
      DO 13 L=L1,LMMAX
13    LXI(L)=LXI(L)-LEE
      L1=LMAX+1
      JLM=1
      DO 16 L=1,L1
      LL=L+L-1
      DO 16 M=1,LL
      JLP=JLM+LL+1-2*M
      DO 14 LT1=1,LMMAX
      IF (JLM-LX(LT1)) 14,15,14
14    CONTINUE
15    LT(LT1)=JLP
      LXM(JLM)=LT1
16    JLM=JLM+1
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE MFOLD_SIMPLE PRODUCES, FOR SUBROUTINE RSMF, THE INDIVIDUAL
!  REFLECTION AND TRANSMISSION MATRIX ELEMENTS FOR A BRAVAIS-LATTICE
!  LAYER. A DEBYE-WALLER FACTOR IS INCLUDED EXPLICITLY, IF REQUESTED
!  (IT.NE.0). AN ORIGIN SHIFT TO A SYMMETRY AXIS OR PLANE IS INCLUDED
!  FOR UP TO FOUR DIFFERENT SHIFTS (YIELDING FOUR DIFFERENT RESULTS).
!  SYMMETRY-INDUCED FOLDING OF THE MATRIX ELEMENTS IS ALSO CARRIED OUT
!  HERE.
!
!  This is a simplified version that removes support for registry
!  shifts ans symmetry
!
!  FOR QUANTITIES NOT EXPLAINED BELOW SEE SUBROUTINE MSMF.
!   JG= INDEX OF CURRENT SCATTERED BEAM.
!   LM= LMAX+1.
!   YLM= RESULT OF (1-X)**(-1.)*Y*T, FROM MSMF.
!   CYLM= SET OF SPHERICAL HARMONICS, FROM MSMF.
!   RT= OUTPUT MATRIX ELEMENTS. 1st index is reflection, second is transmission

      SUBROUTINE  MFOLD_SIMPLE (JG, LM, YLM, CYLM, LMMAX, NT, RT)

      INTEGER JG, LM
      COMPLEX  CYLM, YLM, CZ, RU, CI, RA, TA, ST, SM, SL, CY, CTR, CTT
      COMPLEX  R, T, RT
      DIMENSION RT(2) ! 1 is reflection, 2 is transmission
      DIMENSION  CYLM(NT, LMMAX), YLM(LMMAX)

!     Unclear why here we use the NOVECTOR that disables vectorization
!     of the L loop. Probably worth testing?
CDIR$ NOVECTOR
      CZ = (0.0, 0.0)
      RU = (1.0, 0.0)
      CI = (0.0, 1.0)

      R = CZ
      T = CZ
!     START SUMMATION OVER (L,M)
      ST = RU
      SL = RU
      JLM = 1
      L1 = 1
      DO L = 1, LM
        SM = SL
        LL = L + L - 1
        L1 = L1 + LL
        DO M = 1, LL
          MM = L1 - M
          CY = CYLM(JG, MM)
          CTT = YLM(JLM) * ST
          CTR = CTT * SL
          CTT = CTT * SM
          R = R + CTR * CY
          T = T + CTT * CY
          JLM = JLM + 1
          SM =  - SM
        ENDDO
        SL =  - SL
        ST = ST * CI
      ENDDO

      RT(1) = R
      RT(2) = T
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE MFOLT_SIMPLE computes the probability of reflecting and
!  transmitting a beam into another for a composite layer.
!  This simplified version removes support for beam symmetry codes
!  and registry shifts. Uses input from RTINV_SIMPLE.
!
!  This stripped down version was created on 2021-09-02 by Michele Riva
!
!  Parameters
!  ----------
!  JG : real
!      Index (within this beam set) of scattered beam for which
!      reflection/transmission are to be computed
!  NA : real
!      Offset of the current beam set within PQ
!  AK2, AK3: real
!      In-plane components of wave vector of primary beam
!  CYLM : complex array, shape (NT, LMMAX)
!      Spherical harmonics for all beams (first index) at all
!      (L, M) combinations.
!  N : INTEGER
!      Total number of beams in the current beam set (same for both
!      incident and scattered)
!  LMMAX : INTEGER
!      Maximum number of (L, M) combinations
!  PQ : real array, shape (2, NT)
!      List of all beams, i.e., reciprocal lattice vectors
!  TS : complex vector, length LMN == NLAY*LMMAX
!      Multiple-scattered amplitudes from each of the (L, M)
!      combinations for each of the sublayers in the composite layer.
!      indices 1...LMMAX contain the (L, M) combinations for layer 1.
!      The other layers follow. This is calculated inside RTINV_SIMPLE.
!  LMN : INTEGER
!      == NLAY*LMMAX. Total number of (L, M) combinations for
!      all sublayers.
!  RG : complex array, shape (2, NLAY, NT)
!      Plane-wave propagators between subplanes. RG(1, ...)
!      corresponds to propagation into the solid; RG(2, ...)
!      towards vacuum;
!      RG(..., ILAY, IBEAM) corresponds to propagation of beam IBEAM
!      of the current beam set up to sublayer ILAY.
!  RG_R, RG_T : complex vectors, length NLAY
!      Contributions to reflection (_R) and transmission (_T) of
!      pure propagation to layer ILAY, relative to incidence/exit
!      at the extremes of the composite layer. These are calculated
!      and used internally in MFOLT_SIMPLE.
!  NLAY : INTEGER
!      Total number of sublayers in the current composite layer.
!  JGP : INTEGER
!      Index of incident beam for which reflection/transmission
!      are computed
!  LXM : integer vector, length LMMAX
!      Permutations of (L, M) sequence
!  RT_OUT : complex, OUTPUT
!      Reflection/transmission for the composite layer. Reflection
!      is RT_OUT(1), transmission is RT_OUT(2)
!  INC : {-1, +1}
!      Direction of the incident beam for which reflection and
!      transmission should be computed. +1 corresponds to
!      an incident beam moving toward the solid, -1 away from it.
!  POSS : real array, shape (NLAY, 3)
!      Atomic positions of each of the Bravais sublayers, sorted
!      from the surface towards the inside of the solid. Order
!      of coordinates is: POSS(..., 1) = z (perp, to surface),
!      POSS(..., 2/3) = x/y (in-plane).
      SUBROUTINE  MFOLT_SIMPLE (JG, NA, AK2, AK3, CYLM, N, LMMAX, PQ,
     +                          NT, TS, LMN, RG, RG_R, RG_T, NLAY, JGP,
     +                          LXM, RT_OUT, INC, POSS)
      INTEGER JG, NA, N, LMMAX, NT, LMN, NLAY, JGP, INC,  ! Input
     +        JGA                                         ! Internal use. TODO: update
      COMPLEX  CYLM, CZ, RU, CI, RT_OUT, ST, SM,
     +         CY, CTR, CTT, R, T, CR, CT, CS
      COMPLEX  EDW_R, EDW_T
      COMPLEX  EXP
      COMPLEX  TS(LMN), RG(2, NLAY, N)
      COMPLEX  RG_R(NLAY), RG_T(NLAY)
      REAL     K_IN, K_OUT          ! K_IN/K_OUT are total wave vectors of incident and scattered beams
      DIMENSION  GP(2), PQ(2, NT)   ! GP is the reciprocal-space part of the wave vector of the incident beam
      DIMENSION  CYLM(NT, LMMAX)
      DIMENSION  LXM(LMMAX), POSS(NLAY, 3), RT_OUT(2)
      DIMENSION  K_IN(2), K_OUT(2)

      REAL BR, BT

      COMMON /MFB/ GP, LM

CDIR$ NOVECTOR
      CZ = (0.0, 0.0)
      RU = (1.0, 0.0)
      CI = (0.0, 1.0)

      JGA = JG + NA
      K_OUT(1) = PQ(1, JGA) + AK2
      K_OUT(2) = PQ(2, JGA) + AK3
      K_IN(1) = GP(1) + AK2
      K_IN(2) = GP(2) + AK3

      IF (INC.EQ.1) THEN
         I1 = 1
         I2 = NLAY
      ELSE
         I1 = NLAY
         I2 = 1
      ENDIF

!     BR is the total in-plane propagation phase shift due to reflection
!     at the first layer encountered during propagation (I1): both K_IN
!     and K_OUT at first.
      BR =   (K_IN(1) - K_OUT(1)) * POSS(I1, 2)
     +     + (K_IN(2) - K_OUT(2)) * POSS(I1, 3)

!     BT is the total in-plane propagation phase shift due to transmission
!     between the first (I1) and last (I2) layers encountered: K_IN
!     at first, K_OUT at last
      BT =    K_IN(1)  * POSS(I1, 2) + K_IN(2)  * POSS(I1, 3)
     +     - (K_OUT(1) * POSS(I2, 2) + K_OUT(2) * POSS(I2, 3))

!     EDW_R/T are the in-plane propagators for reflection and
!     transmission
      EDW_R = EXP(CI * BR)
      EDW_T = EXP(CI * BT)

!     Prepare coefficients taking into account the effect of
!     propagation on the scattered amplitudes. These will then
!     be later combined with the relevant L,M contributions.
!        RG_R is effect on REFLECTION, RG_T effect on TRANSMISSION
!     For RG: RG(1, ...) is propagation towards solid, RG(2, ...) away
      DO J = 1, NLAY
        IF (INC.EQ.1) THEN
          ST = RG(1, J, JGP) / RG(1, 1, JGP)          ! Fraction of incident beam (JGP) reaching plane J
!                                                       (relative to the one reaching first layer) due to
!                                                       propagation. Both propagating down.
          RG_R(J) = ST
     +              * RG(2, J, JG) / RG(2, 1, JG)     ! Fraction of scattered beam (JG) leaving layer J
!                                                       (relative to the one reaching first layer) due to
!                                                       propagation. Both propagating up.
          RG_T(J) = ST
     +              * RG(1, NLAY, JG) / RG(1, J, JG)  ! Fraction of scattered beam (JG) leaving last layer
!                                                       (relative to the one reaching layer J) due to
!                                                       propagation. Both propagating down.
        ELSE
          ST = RG(2, NLAY, JGP) / RG(2, J, JGP)       ! Fraction of incident beam (JGP) reaching last
!                                                       layer relative to the one at plane J, due to
!                                                       propagation. Both propagating up.
          RG_R(J) = ST
     +              * RG(1, NLAY, JG) / RG(1, J, JG)  ! Fraction of scattered beam (JG) leaving last layer
!                                                       relative to the one reaching plane J, due to
!                                                       propagation. Both propagating up.
          RG_T(J) = ST
     +              * RG(2, J, JG) / RG(2, 1, JG)     ! Fraction of scattered beam (JG) leaving plane J
!                                                       relative to the one leaving first layer, due to
!                                                       propagation. Both propagating down.
        ENDIF
      ENDDO

!     Combine the pure effect of propagation with the (L,M)-dependent
!     scattering probabilities.
!       Some optimization from older code: it used to include products
!       of type SM*SL, with
!           SL = -1*(-1)**L,
!       and
!           SM = -SL * (-1)**M = (-1)**(L+M)
!       Thus, SM*SL == -1 * (-1)**(2L+M) == -1*(-1)**M. This quantity
!       is now called SM.
      R = CZ
      T = CZ
      JLM = 1
      L1 = 1
      DO L = 1, LM
        SM = RU
        LL = L + L - 1
        L1 = L1 + LL
        L2 = L1 - 2*L
        DO M = 1, LL
          KLM = LXM(JLM)
!         Sum over all sublayers
          CR = CZ
          CT = CZ
          DO J = 1, NLAY
            CS = TS((J - 1) * LMMAX + KLM)  ! Scattering at layer J with angular momentum numbers (L, M)
            CR = CR + RG_R(J) * CS          ! Including the propagation part on reflection
            CT = CT + RG_T(J) * CS          ! And the propagation part on transmission
          ENDDO

!         Finally include spherical harmonics and some more prefactors
          CY = CYLM(JG, L2 + M)
          IF (INC.EQ.1) THEN
            CTR = SM
            CTT = RU
          ELSE
            CTR = RU
            CTT = SM
          ENDIF
          R = R + CTR * CY * CR
          T = T + CTT * CY * CT
          JLM = JLM + 1
          SM =  - SM
        ENDDO
      ENDDO

      RT_OUT(1) = R * EDW_R
      RT_OUT(2) = T * EDW_T

      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE MULTAMP_OPT multiplies one vector into a matrix.
!  Optimized for Fortran memory access - by MR v1.72
!
      SUBROUTINE MULTAMP_OPT(INAMP,OUTAMP,MATRIX,N)
      COMPLEX INAMP(N),OUTAMP(N),MATRIX(N,N)
      COMPLEX CZ, TEMP

      CZ = (0.,0.)
      DO I=1,N
        OUTAMP(I) = CZ
      ENDDO

      DO J=1,N
        TEMP = INAMP(J)
        DO I=1,N
          OUTAMP(I) = OUTAMP(I) + MATRIX(I,J) * TEMP
        ENDDO
      ENDDO
      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE OPENOUT OPENS THE OUTPUTFILE FOR THE AMPLITUDES
C  AND SUPPLIES ALL NEEDED GLOBAL INFORMATION
C
      SUBROUTINE OPENOUT(IFILE,FILENAM,IFORM)
      CHARACTER*(*) FILENAM
C
      IF (IFORM.EQ.0) THEN
      OPEN (IFILE,FILE=FILENAM,FORM='UNFORMATTED',STATUS='NEW')
      ELSE
      OPEN (IFILE,FILE=FILENAM,FORM='FORMATTED')
      ENDIF
      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE OUTAMP WRITES THE OUTPUT AMPLITUDES TO FILE# IFILE
C
      SUBROUTINE OUTAMP(IFILE,IFORM,NEXIT,PQ1,PQ2,ALM,LMMAX,AK21,AK31)
      COMPLEX ALM(LMMAX)
C
      IF (IFORM.EQ.0) THEN
C  UNFORMATTED OUTPUT
      WRITE (IFILE) NEXIT
      WRITE (IFILE) PQ1,PQ2,AK21,AK31,ALM
      ELSE
C  FORMATTED OUTPUT
      WRITE (IFILE,10) NEXIT
      WRITE (IFILE,20) PQ1,PQ2,AK21,AK31
      WRITE (IFILE,20) (ALM(I),I=1,LMMAX)
   10 FORMAT (I5)
   20 FORMAT(5E12.6)
      ENDIF
      RETURN
      END

C---------------------------------------------------------------------
C  OUTXIST STORES THE AMPLITUDES OF THE INTEGER ORDER
C  BEAMS AFTER THEY HAVE BEEN CALCULATED FOR THE INCIDENT BEAM
C
C  XI(NT): AMPLITUDES OF BEAMS FOR CURRENT ENERGY
C  XIST(NT0): ORIGINAL AMPLITUDES OF EACH EXIT BEAM
C                FOR EACH ENERGY: =0 IF IG=FRACTIONAL
C                ORDER BEAM.
C
      SUBROUTINE OUTXIST(IFILE,IFORM,E,PQF,PQFEX,NT0,NT,XI,XIST,
     1                    L1,CAF)
C
      DIMENSION PQF(2,NT),PQFEX(2,NT0)
      COMPLEX XI(NT), XIST(NT0), CAF(L1)
      COMMON /ADS/ ASE,VPIS,VPIO,VV ! note: V0 removed!

C   START LOOP OVER EXIT BEAMS.

       DO 1 IG=1,NT0

         XIST(IG)=CMPLX(0.0,0.0)

         DO 1 JG=1,NT

           P1=ABS(PQF(1,JG)-PQFEX(1,IG))
           P2=ABS(PQF(2,JG)-PQFEX(2,IG))

   1       IF((P1+P2).LT.0.000001) XIST(IG)=XI(JG)


C  NOW WRITE ALL AMPLITUDES TO FILE # IFILE

      IF (IFORM.EQ.0) THEN
C  UNFORMATTED OUTPUT
        WRITE (IFILE) E,VPIS,VPIO,VV,L1,CAF,XIST
      ELSE
C  FORMATTED OUTPUT
        WRITE(IFILE,20) E,VPIS,VPIO,VV
        WRITE(IFILE,10)L1
        WRITE(IFILE,20)(CAF(I),I=1,L1)
        WRITE(IFILE,20)(XIST(I),I=1,NT0)
   10   FORMAT (I5)
   20   FORMAT (5E16.10)
      ENDIF

      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE PSTEMP INCORPORATES THE THERMAL VIBRATION EFFECTS IN THE
C  PHASE SHIFTS, THROUGH A DEBYE-WALLER FACTOR. ISOTROPIC VIBRATION
C  AMPLITUDES ARE ASSUMED.
C   PPP= CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CPPP.
C   N3= NO. OF INPUT PHASE SHIFTS.
C   N2= DESIRED NO. OF OUTPUT TEMPERATURE-DEPENDENT PHASE SHIFTS.
C   N1= N2+N3-1.
C   DR0= 4TH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDES.
C   DR= ISOTROPIC RMS VIBRATION AMPLITUDE AT REFERENCE TEMPERATURE T0.
C   T0= ARBITRARY REFERENCE TEMPERATURE FOR DR.
C   TEMP= ACTUAL TEMPERATURE.
C   E= CURRENT ENERGY (REAL NUMBER).
C   PHS= INPUT PHASE SHIFTS.
C   DEL= OUTPUT (COMPLEX) PHASE SHIFTS.
C   11.01.95 UL: uses BESSEL from TLEED package to calculate spherical
C   bessel-functions. Yield better convergence for higher vibration amplitudes
C   than original calculation scheme                                     110195
      SUBROUTINE  PSTEMP (PPP, N1, N2, N3, DR0, DR, T0, TEMP, E, PHS,
     1DEL,CTAB,SUM,BJ)
      COMPLEX  DEL, SUM, CTAB
      COMPLEX*16 BJ(N1)
      COMPLEX  Z, CI, CS, CL
      COMPLEX EXP,LOG
      DIMENSION  PPP(N1,N2,N3), PHS(N3), DEL(N2), SUM(N2)
      DIMENSION  CTAB(N3)

      PI = 3.14159265
      CI = CMPLX(0.0,1.0)
      DO 170 J = 1, N2
  170 DEL(J) = (0.0,0.0)
      ALFA = DR * DR * TEMP/T0
      ALFA = 0.166667 * SQRT(ALFA * ALFA + DR0)
      FALFE =  - 4.0 * ALFA * E
      IF (ABS(FALFE)-1.0E-3)  180, 200, 200
  180 DO 190 J = 1, N3
  190 DEL(J) = CMPLX(PHS(J),0.0)
      GO TO 360
C
COMMENT BJ(N1) IS LOADED WITH SPHERICAL BESSEL FUNCTIONS OF
C       THE FIRST KIND; ARGUMENT Z
  200 Z = FALFE * CI
c
c     use subroutine BESSEL
      call bessel(BJ,Z,N1)
c
  270 CS = CMPLX(1.0,0.0)
      FL = 1.0
      DO 280 I = 1, N1

      BJ(I) = EXP(FALFE) * FL * CS * BJ(I)

      FL = FL + 2.0
  280 CS = CS * CI
C
      FL = 1.0
      DO 290 I = 1, N3
      CTAB(I) = (EXP(2.0 * PHS(I) * CI) - (1.0,0.0)) * FL
  290 FL = FL + 2.0
C
      ITEST = 0
      LLLMAX = N2
      FL = 1.0
      DO 350 LLL = 1, N2
      SUM(LLL) = CMPLX(0.0,0.0)
      DO 300 L = 1, N3
      LLMIN = IABS(L - LLL) + 1
      LLMAX = L + LLL - 1
      DO 300 LL = LLMIN, LLMAX
  300 SUM(LLL) = SUM(LLL) + PPP(LL,LLL,L) * CTAB(L) * BJ(LL)

cvb  now, sum(LLL) is already the temperature-dependent t-matrix we were
C    looking for. It is next converted to a temp-dependent phase shift, only
C    to be converted back right after the pstemp call in tscatf. Kept for
C    the sake of compatibility with van Hove / Tong book only.

      DEL(LLL) =  - CI * LOG(SUM(LLL) + (1.0,0.0))/(2.0,0.0)
      ABSDEL = CABS(DEL(LLL))
      IL = LLL - 1
      IF (ABSDEL-1.0E-2)  320, 310, 310
  310 ITEST = 0
      GO TO 350
  320 IF (ITEST-1)  340, 330, 340
  330 LLLMAX = LLL
      GO TO 360
  340 ITEST = 1
  350 FL = FL + 2.0
  360 RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE READGEO reads structural (and some calculational) information
C  for the reference surface in question.
C  VB, 06/98
C

      SUBROUTINE READGEO(NEL,NSITE,NLTYPE,NNSUB,NSTACK,NBRAV,
     +                   CONC,VIB,NSUB,LBRAV,LCOMP,LATT,STYPE,
     +                   SUBPOS,TSLAB,ASA,TOPLAYB,BOTLAYB,ASBULK,LTYPE,
     +                   LDIST,TENS,LAYFILE,ASB,NL)

C  include global quantities

      INCLUDE "GLOBAL"

C  declare variables

C  dimensions:
C  NEL   : number of different elements in phase shifts
C  NSITE : number of different site occupations (vibrational, chemical)
C  NLTYPE: number of different layer types
C  NSTACK: number of layer stacking steps (layers on top of bulk)

C  NSITEH,NLTYPEH,NSTACKH are only read to check consistency of input
C  NBRAV : number of Bravais layer types
C  NCOMP : number of composite layer types (for consistency only)

      integer NEL,NSITE,NLTYPE,NSTACK,NNSUB
      integer NSITEH,NLTYPEH,NSTACKH,NBRAV,NCOMP

C  CONC  : concentration of element IEL on site type ISITE
C  VIB   : vibrational amplitude of element IEL on site type ISITE
C  NSUB  : number of (Bravais) subplanes in layer type ILTYPE
C  LBRAV : layer number ILTYPE is Bravais layer no. LBRAV(ILTYPE) (0 if comp.)
C  LCOMP : layer number ILTYPE is composite layer no. LCOMP(ILTYPE) (0 if Brav.)
C  LATT  : lattice type of layer type ILTYPE ( 1 = overlayer, 2 = substrate)
C  STYPE : site type of layer type ILTYPE, subplane ISUB
C  SUBPOS: 3D position of layer type ILTYPE, subplane ISUB
C  TSLAB : 0: bulk calculation required, 1: no bulk calc. required
C  ASA   : interlayer vector between different bulk units (subras)
C  TOPLAYB: layer type of top "bulk" layer (dblgas)
C  BOTLAYB: layer type of bottom "bulk" layer (dblgas)
C  ASBULK: interlayer vector between top and bottom layer of bulk unit
C          (different from ASA!)
C  LTYPE : layer type of stacked layer no. ISTACK
C          note reverse order: ISTACK = 1 is top layer, ISTACK = NSTACK ist
C          lowest layer
C  LDIST : interlayer vector between layer ISTACK and lower layer ISTACK+1
C  TENS  : tensor output for layer ISTACK desired?

      REAL CONC,VIB
      DIMENSION CONC(NSITE,NEL),VIB(NSITE,NEL)
      INTEGER NSUB,LBRAV,LCOMP,LATT,STYPE
      DIMENSION NSUB(NLTYPE),LBRAV(NLTYPE),LCOMP(NLTYPE)
      DIMENSION LATT(NLTYPE),STYPE(NLTYPE,NNSUB)
      REAL SUBPOS
      DIMENSION SUBPOS(NLTYPE,NNSUB,3)
      INTEGER TSLAB,TOPLAYB,BOTLAYB
      REAL ASA,ASBULK
      DIMENSION ASA(3),ASBULK(3)
      INTEGER LTYPE
      DIMENSION LTYPE(NSTACK)
      REAL LDIST
      DIMENSION LDIST(NSTACK,3)
      INTEGER TENS
      DIMENSION TENS(NSTACK)
      CHARACTER*40 LAYFILE
      DIMENSION LAYFILE(NSTACK,NNSUB)

C  ASB interlayer vector with min. interlayer distance. (to be used in the
C      TST convergence criterion)
C  NL  number of superlattice sublattices in substrate lattice

      REAL ASB
      DIMENSION ASB(3)
      INTEGER NL

C  first, read information for different atom types (distinguished by element,
C  vibrational amplitude, and stoichiometry)

      READ(5,*)
      READ(5,*)
      READ(5,*)

      READ(5,'(I3)') NSITEH
      IF (NSITEH.ne.NSITE) THEN
        write(6,*) "NSITE .ne. MNSITE - please correct!"
        STOP
      END IF

      DO ISITE = 1,NSITE,1

        READ(5,*)

        DO IEL = 1,NEL,1
            ! AMI: note change F7.4 -> F9.4 in v1.71 by LH
          READ(5,'(2F9.4)') CONC(ISITE,IEL),VIB(ISITE,IEL)
          VIB(ISITE,IEL) = VIB(ISITE,IEL)/BOHR

        ENDDO

      ENDDO

C  next, read information on different layer types to be included

      READ(5,*)
      READ(5,*)
      READ(5,*)
      READ(5,'(I3)') NLTYPEH
      IF (NLTYPEH.ne.NLTYPE) THEN
        write(6,*) "NLTYPE .ne. MNLTYPE - please correct!"
        STOP
      END IF

C  KBRAV, KCOMP are only counters

      KBRAV = 0
      KCOMP = 0

      DO ILTYPE = 1,NLTYPE

        READ(5,*)

        READ(5,'(I3)') LATT(ILTYPE)

        READ(5,'(I3)') NSUB(ILTYPE)

        IF (NSUB(ILTYPE).gt.NNSUB) THEN
          write(6,*) "Number of subplanes in layer type ", ILTYPE,
     +               "exceeds dimension MNSUB - please correct!"
          STOP
        ELSE IF (NSUB(ILTYPE).le.0) THEN
          write(6,*) "Number of subplanes in layer type ", ILTYPE,
     +               "less than 1 ? - please correct!"
          STOP
        END IF

        IF (NSUB(ILTYPE).eq.1) THEN
          KBRAV = KBRAV+1
          LBRAV(ILTYPE)=KBRAV
          LCOMP(ILTYPE)=0
        ELSE
          KCOMP = KCOMP+1
          LCOMP(ILTYPE)=KCOMP
          LBRAV(ILTYPE)=0
        END IF

        DO ISUB = 1,NSUB(ILTYPE)

          READ(5,'(I3,3F9.4)') STYPE(ILTYPE,ISUB), ! used to be 3F7.4 pre v1.7
     +                        (SUBPOS(ILTYPE,ISUB,IPOS),IPOS=1,3)

          IF ((STYPE(ILTYPE,ISUB).lt.1)
     +        .or.
     +        (STYPE(ILTYPE,ISUB).gt.NSITE)) THEN
            write(6,*) "Site type",STYPE(ILTYPE,ISUB)," in sublayer ",
     +                 ISUB," of layer ",ILTYPE,":"
            write(6,*) "No such site defined - please correct!"
            STOP
          END IF

          DO IPOS=1,3
            SUBPOS(ILTYPE,ISUB,IPOS)=SUBPOS(ILTYPE,ISUB,IPOS)/BOHR
          ENDDO

        ENDDO

        IF ((ISUB.eq.1).and.
     +      ((SUBPOS(ILTYPE,1,1).ne.0.).or.
     +       (SUBPOS(ILTYPE,1,2).ne.0.).or.
     +       (SUBPOS(ILTYPE,1,3).ne.0.)    ))   THEN

          write(6,*)
     +    "For a Bravais layer, a subplane position makes no sense and"
          write(6,*)
     +    "should be (0,0,0). To make sure, please correct."
          STOP

        END IF

      ENDDO

      IF (KBRAV.ne.NBRAV) THEN
        write(6,*) "Number of Bravais layers in actual input does not"
        write(6,*) "match dimension MNBRAV - please correct!"
        STOP
      END IF

C  now gather info for layer doubling

      READ(5,*)
      READ(5,*)
      READ(5,*)

      READ(5,'(I3)') TSLAB

      IF (TSLAB.eq.0) THEN

C  ASA IS THE SUBSTRATE INTERLAYER VECTOR (LINKING REFERENCE POINTS
C  ON SUCCESSIVE LAYERS)(ANGSTROM)
        READ(5,'(3F9.4)')(ASA(I),I=1,3) ! used to be 3F7.4 pre v1.7
        READ(5,'(I3)') TOPLAYB
        READ(5,'(I3)') BOTLAYB
        READ(5,'(3F9.4)') (ASBULK(I),I=1,3) ! used to be 3F7.4 pre v1.7

        IF ((TOPLAYB.lt.1).or.(TOPLAYB.gt.NLTYPE)) THEN
          write(6,*) "Top bulk layer is illegal layer type.",
     +               " Please correct."
          STOP
        END IF

        IF ((BOTLAYB.lt.1).or.(BOTLAYB.gt.NLTYPE)) THEN
          write(6,*) "Bottom bulk layer is illegal layer type.",
     +               " Please correct."
          STOP
        END IF

        DO  I=1,3
          ASA(I)=ASA(I)/BOHR
          ASBULK(I)=ASBULK(I)/BOHR
        ENDDO

        ASB(1) = ASBULK(1)

      ELSE IF (TSLAB.eq.1) THEN
        write(6,*) "No bulk calculation by subref required!"
        READ(5,*)
        READ(5,*)
        READ(5,*)
        READ(5,*)
      ELSE
        write(6,*) "TSLAB should be 0 or 1 - please correct!"
        STOP
      END IF

C  next, read info for layer stacking

      READ(5,*)
      READ(5,*)
      READ(5,*)

      READ(5,'(I3)') NSTACKH
      IF (NSTACKH.ne.NSTACK) THEN
        write(6,*) "NSTACK .ne. MNSTACK - please correct!"
        STOP
      END IF

      DO ISTACK = NSTACK,1,-1

        READ(5,'(I3,3F9.4)') LTYPE(ISTACK),(LDIST(ISTACK,I),I=1,3,1) ! used to be 3F7.4 pre v1.7

        IF ((LTYPE(ISTACK).lt.1).or.(LTYPE(ISTACK).gt.NLTYPE)) THEN
          write(6,*) "Stacked layer number ",ISTACK,
     +               "has illegal layer type. Please correct."
          STOP
        END IF

        DO IPOS = 1,3,1
          LDIST(ISTACK,IPOS) = LDIST(ISTACK,IPOS)/BOHR
        ENDDO

        IF (ABS(LDIST(ISTACK,1)).le.ASB(1)) THEN
          ASB(1) = LDIST(ISTACK,1)
        END IF

        READ(5,'(I3)') TENS(ISTACK)

        IF (TENS(ISTACK).eq.0) THEN
          write(6,*) "Tensor output for layer no. ",ISTACK,
     +               " not required."
        ELSE IF (TENS(ISTACK).eq.1) THEN

          IF ((NL.gt.1).and.(LATT(LTYPE(ISTACK)).eq.2)) THEN
            write(6,*) "In the presence of a non-trivial superlattice,"
            write(6,*) "Tensor LEED is currently not possible for a"
            write(6,*) "substrate type layer."
            STOP
          END IF

          write(6,*) "Tensor output for layer no. ",ISTACK,
     +               " required. Filename(s): "

          DO ISUB=1,NSUB(LTYPE(ISTACK)),1

            READ(5,'(A20)') LAYFILE(ISTACK,ISUB)
            WRITE(6,*) "Layer ",ISTACK,", subplane ",ISUB,": ",
     +                 LAYFILE(ISTACK,ISUB)

          ENDDO

        ELSE
          write(6,*) "Tensor LEED flag for layer ",ISTACK,
     \               " not 0 or 1 - please correct!"
          STOP
        END IF

      ENDDO

      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE READIN READS IN MOST OF THE INPUT FOR THE LEED PROGRAMS.
C  INPUT FROM MAIN PROGRAM
C   NL=NL1*NL2= NO. OF SUBLATTICES DUE TO SUPERLATTICE
C   KNBS= NO. OF BEAM SETS TO BE READ IN
C   KNT= TOTAL NO. OF BEAMS TO BE READ IN
C   NPUN= NO. OF BEAMS FOR WHICH INTENSITIES ARE TO BE PUNCHED OUT
C   NPSI= NO. OF ENERGIES AT WHICH PHASE SHIFTS WILL BE READ IN
C   FOR FURTHER QUANTITIES SEE BELOW
C   TST IS READ IN IN F11.8 FORMAT !                                    070592
C   PHASE-SHIFTS ARE OUTPUT FOR FIRST PHASE-SHIFT ONLY                  070592
C   MNEL and MLMAX1 is required for correct dimension of PHSS           070592

      SUBROUTINE READIN(TITLE,TVA,RAR1,RAR2,TVB,
     +IDEG,NL,V,VL,JJS,KNBS,KNB,KNT,SPQF,KSYM,SPQ,TST,TSTS,NPUN,NPU,
     +THETA,FI,NPSI,ES,PHSS,L1,NEL,LMAX1,VPI,EI,EF,DE)

      INCLUDE "GLOBAL"

      Character*80 TITLE,PAR0,PAR1

      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),
     1ARB1(2),ARB2(2),RBR1(2),RBR2(2)
      DIMENSION V(NL,2),JJS(NL,IDEG),KNB(KNBS),SPQF(2,KNT),KSYM(2,KNT),
     1SPQ(2,KNT),NPU(NPUN),ES(NPSI),PHSS(NPSI,NEL,LMAX1)
      COMPLEX VL(NL,2)
C  Format parameter for phaseshift file
      INTEGER PSFMT,PSL1
C  Parameters for calculation of inner potential
	  REAL EM,C1,C2,C3,C4,C5,C6,C7,C8, WORKFN
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
      COMMON /MS/LMAX,EPS,LITER
      COMMON /ADS/ASE,VPIS,VPIO,VV
	  COMMON /INPOT/EM,C1,C2,C3,C4,C5,C6,C7,C8,WORKFN
  130 FORMAT(5F9.4)
  140 FORMAT(25H0PARAMETERS FOR INTERIOR )
  160 FORMAT(3F9.4)
  161 FORMAT(2F9.4)
  162 FORMAT(F7.4)
  170 FORMAT(10H SURF VECS,2(5X,2F9.4))
  171 FORMAT(15X,2F9.4)
  190 FORMAT(3H SO,4X,2F7.4)
C  LH Changed from 15I4 to 15I5 in order to manage more than 9999 beams
  200 FORMAT(15I5)
  201 FORMAT(26I3)
  202 FORMAT(3x,I3,2x,F6.2,1x,4F7.2,2x,F7.2,3F8.4)
  203 FORMAT(F6.2,1x,4F7.2,2x,F7.2,3F8.4)
  204 FORMAT(' PHASESHIFT FORMAT IS',I3)
  205 FORMAT(I3,1x,4F8.2)
  210 FORMAT(6H VPI = ,F7.4,10H    ASE = ,F7.4,13H    WORKFN = ,F7.4)
  220 FORMAT(3H SS,3X,2F7.4)
  230 FORMAT(24H PARAMETERS FOR SURFACE )
  250 FORMAT(1H0,9X,3HPQ1,4X,3HPQ2,8X,3HSYM)
  260 FORMAT(1H )
  270 FORMAT(2F10.5,4I3)
  280 FORMAT(1H ,1I4,F10.3,F7.3,4X,4I3)
  281 FORMAT(7H TST = ,1F11.8)
  282 FORMAT(1F11.8)
  283 FORMAT(I3,2F6.3,I3)
C  LH Changed from 15I4 to 15I5 in order to manage more than 9999 beams
  284 FORMAT(17H PUNCHED BEAMS   ,15I5)
  285 FORMAT(13H THETA  FI = ,2F7.2)
C  290 FORMAT(3H VO,4X,F7.2,2HVV,4X,F7.2)
  295 FORMAT(7H EPS = ,1F7.4,10H  LITER = ,1I3)
  305 FORMAT(7H IT1 = ,1I3,7H IT2 = ,1I3,7H IT3 = ,1I3)
  325 FORMAT(' THDB= ',1F9.4,' AM= ',1F9.4,' FPER= ',1F9.4,' FPAR= ',
     *' DR0= ',1F9.4)
  340 FORMAT(10F7.4)
  350 FORMAT(1H ,12HPHASE SHIFTS)
  351 FORMAT(1H ,30HPARAMETERS FOR INNER POTENTIAL)
  355 FORMAT(8H LMAX = ,1I3)
  359 FORMAT(29H 1st ENERGY FOR ELEMENT/SITE:,I3)
  360 FORMAT(1H ,4HE = ,1F7.4,3X,20F8.4)
  361 FORMAT(1H ,11X,13H  2ND ELEMENT,3X,10F8.4)
  362 FORMAT(1H ,11X,13H  3RD ELEMENT,3X,10F8.4)

      PI = 3.14159265359

C  consistency check

      IF (NL.ne.KNBS) then
        write(6,*)
     +  "*************************************************************"
        write(6,*)
     +  "*                                                           *"
        write(6,*)
     +  "*  Warning !! Number of beam sets, MKNBS, differs           *"
        write(6,*)
     +  "*  from number of sublattices, MNL = MNL1*MNL2 !            *"
        write(6,*)
     +  "*  Assuming incorrect input, I stop!                        *"
        write(6,*)
     +  "*                                                           *"
        write(6,*)
     +  "*************************************************************"
        STOP
      END IF

C  read title

      READ(5,'(A80)') TITLE
      WRITE(6,*) TITLE

C  read energy range here

      read(5,'(3F9.2)') EI,EF,DE

C  ARA1 AND ARA2 ARE TWO 2-D BASIS VECTORS OF THE SUBSTRATE LAYER
C  LATTICE. THEY SHOULD BE EXPRESSED IN TERMS OF THE PLANAR CARTESIAN
C  Y- AND Z-AXES (X-AXIS IS PERPENDICULAR TO SURFACE)(ANGSTROM)
      READ(5,160)(ARA1(I),I=1,2)
      READ(5,160)(ARA2(I),I=1,2)
      WRITE(6,170)(ARA1(I),I=1,2)
      WRITE(6,171)(ARA2(I),I=1,2)
      DO 460 I=1,2
      ARA1(I)=ARA1(I)/BOHR
  460 ARA2(I)=ARA2(I)/BOHR
      TVA=ABS(ARA1(1)*ARA2(2)-ARA1(2)*ARA2(1))
      ATV=2.0*PI/TVA
C  RAR1 AND RAR2 ARE THE RECIPROCAL-LATTICE BASIS VECTORS CORRESPONDING
C  TO ARA1 AND ARA2
      RAR1(1)=ARA2(2)*ATV
      RAR1(2)=-ARA2(1)*ATV
      RAR2(1)=-ARA1(2)*ATV
      RAR2(2)=ARA1(1)*ATV
      WRITE (6,230)
C  ARB1,ARB2,RBR1,RBR2 ARE EQUIVALENT TO ARA1,ARA2,RAR1,RAR2
C  BUT FOR AN OVERLAYER
      READ(5,160)(ARB1(I),I=1,2)
      READ(5,160)(ARB2(I),I=1,2)
      WRITE(6,170)(ARB1(I),I=1,2)
      WRITE(6,171)(ARB2(I),I=1,2)
      DO 467 I=1,2
      ARB1(I)=ARB1(I)/BOHR
  467 ARB2(I)=ARB2(I)/BOHR
      TVB=ABS(ARB1(1)*ARB2(2)-ARB1(2)*ARB2(1))
      ATV=2.0*PI/TVB
      RBR1(1) = ARB2(2) * ATV
      RBR1(2) =  - ARB2(1) * ATV
      RBR2(1) =  - ARB1(2) * ATV
      RBR2(2) = ARB1(1) * ATV
C  VPI IS THE IMAGINARY PART OF THE INNER POTENTIAL.
C  ASE IS THE SPACING BETWEEN THE SURFACE (WHERE MUFFIN-TIN CONSTANT
C  AND DAMPING SET IN) AND THE TOP-LAYER NUCLEI (ANGSTROM).
      READ (5,160)  VPI,ASE,WORKFN
      WRITE (6,210) VPI,ASE,WORKFN
      ASE=ASE/BOHR
      VPI=-VPI/HARTREE
      VPIS=VPI
      VPIO=VPI
C  The variables VPIS and VPIO should be removed somewhen from the code. LH 25.03.21

C  SLIND COMPARES SUBSTRATE LATTICE AND SUPERLATTICE FOR THE BENEFIT
C  OF FMAT
      CALL SLIND(V,VL,JJS,NL,IDEG,2.0E-4)                               141278
C     WRITE (6,250)
      KNTT= 0
C  READ IN KNBS SETS OF BEAMS
      DO 580 J = 1, KNBS
C     WRITE (6,260)
C  READ IN KNB(J) BEAMS IN J-TH SET
      READ (5,200)  KNB(J)
      N = KNB(J)
      DO 570 K = 1, N
      KK = K + KNTT
C  READ IN A BEAM (AS A RECIPROCAL LATTICE VECTOR IN TERMS OF RAR1,RAR2
C  I.E. IN THE FORM (0,0),(1,0),(0,1),ETC.) AND ITS SYMMETRY PROPERTY
C  (BOTH FOR THE SUPERLATTICE REGION AND THE SUBSTRATE REGION)
      READ (5,270)  SPQF(1,KK), SPQF(2,KK), (KSYM(I,KK),I=1,2)
      DO 560 I = 1, 2
  560 SPQ(I,KK)=SPQF(1,KK)*RAR1(I)+SPQF(2,KK)*RAR2(I)
C     WRITE (6,280) KK, SPQF(1,KK), SPQF(2,KK), (KSYM(I,KK),I=1,2)
  570 CONTINUE
  580 KNTT= KNTT+ KNB(J)
      IF (KNTT.ne.KNT) THEN
        write(6,*) "Warning! Number of beams in full dynamic beam",
     +             " list, ",KNTT," is inconsistent with dimension"
        write(6,*) "MKNT = ",KNT," from PARAM! Please correct!"
        STOP
      END IF
C  TST GIVES THE CRITERION FOR THE SELECTION OF BEAMS AT ANY ENERGY
C  TO BE SELECTED A BEAM MAY NOT DECAY TO LESS THAN A FRACTION TST OF
C  ITS INITIAL VALUE ON TRAVELING FROM ANY LAYER TO THE NEXT
      READ(5,282)TST
      WRITE(6,281)TST
      TSTS=TST
C  PUNCH A CARD WITH NPUN
c      WRITE(7,283)NPUN
C  NPU GIVES THOSE BEAMS FOR WHICH INTENSITIES ARE TO BE PUNCHED OUT.
C  NPU(K)=J MEANS PUNCH OUTPUT IS DESIRED FOR THE J-TH BEAM IN THE
C  INPUT LIST
      READ(5,200)(NPU(K),K=1,NPUN)
      WRITE(6,284)(NPU(K),K=1,NPUN)
C  (THETA,FI) IS THE DIRECTION OF INCIDENCE. THETA=0 AT NORMAL
C  INCIDENCE. FI=0 FOR INCIDENCE ALONG Y-AXIS IN (YZ) SURFACE PLANE
C  (X-AXIS IS PERPENDICULAR TO SURFACE)(DEGREE)
      READ (5,161)  THETA, FI
      WRITE(6,285) THETA,FI
      THETA=THETA*PI/180.0
      FI=FI*PI/180.0
C  EPS,LITER ARE CONVERGENCE CRITERIA FOR LAYER DOUBLING
      READ (5,162)  EPS
      READ (5,201)  LITER
      WRITE(6,295)EPS,LITER

C  LMAX= HIGHEST L-VALUE TO BE CONSIDERED
C  ICLEB = FLAG FUER FELDER CLM,YLM,FAC1,FAC2 - obsolete!
      READ(5,201)LMAX
      WRITE(6,355)LMAX
      IF (LMAX.ne.(LMAX1-1)) THEN
        WRITE(6,*) "Please enforce LMAX in input = MLMAX in PARAM."
        WRITE(6,*) "Inaccurate dimensions may lead to problems."
        STOP
      END IF

      WRITE (6,350)

      READ(5,201) PSFMT
	  WRITE(6,204) PSFMT
      WRITE (6,351)

      L1=LMAX+1

      IF (PSFMT.eq.1) THEN

C  NEL= NO. OF CHEMICAL ELEMENTS FOR WHICH PHASE SHIFTS ARE TO BE READ IN

      READ(5,205)NELCHECK,C4,C1,C2,C3
      IF (NELCHECK.ne.NEL) THEN
        write(6,*) "NEL from phaseshift file not consistent with MNEL!"
        write(6,*) "Please correct!"
        STOP
      ENDIF

      WRITE(6,*) "   C1      C2      C3      C4   "
	  WRITE(6,'(4F8.2)') C1,C2,C3,C4

      DO 660 IENER = 1, NPSI

C  ES= ENERGIES (HARTREES) AT WHICH PHASE SHIFTS ARE INPUT. LINEAR
C  INTERPOLATION OF THE PHASE SHIFTS WILL OCCUR FOR ACTUAL ENERGIES
C  FALLING BETWEEN THE VALUES OF ES (AND LINEAR EXTRAPOLATION ABOVE
C  THE HIGHEST ES)
C  PHSS STORES THE INPUT PHASE SHIFTS (RADIAN)

        READ (5,162)  ES(IENER)
        DO IEL=1,NEL,1
          READ(5,'(10F7.4)')(PHSS(IENER,IEL,L),L=1,L1)
          IF (LMAX .LT. 10) THEN
            READ(5,*)
          END IF
        ENDDO
  660 CONTINUE

	  ELSEIF(PSFMT.eq.2) THEN

      READ(5,201)NELCHECK
      IF (NELCHECK.ne.NEL) THEN
        write(6,*) "NEL from phaseshift file not consistent with MNEL!"
        write(6,*) "Please correct!"
        STOP
      ENDIF

      READ(5,'(1A8,1A70)') PAR0,PAR1
	  WRITE(6,'(A70)') PAR1
      READ(5,202) PSL1,EM,C1,C2,C3,C4,C5,C6,C7,C8
      WRITE(6,203) EM,C1,C2,C3,C4,C5,C6,C7,C8

      DO 661 IEL=1,NEL,1

	    READ(5,*)
        READ(5,*)

        DO IENER = 1, NPSI

          READ(5,'(F8.3,20F8.4)')ES(IENER),(PHSS(IENER,IEL,L),L=1,L1)
		  ES(IENER)=ES(IENER)/HARTREE

        ENDDO
  661 CONTINUE

	  ELSE
	    WRITE(6,*) "Phaseshift format incorrectly assigned!"
		STOP
	  ENDIF


      DO 670 IENER = 1, 1

        DO IEL = 1,NEL,1
		  WRITE(6,359)IEL
          WRITE(6,360)ES(IENER),(PHSS(IENER,IEL,L),L=1,L1)
        ENDDO

  670 CONTINUE

      WRITE(6,2000)
 2000 FORMAT(/)
      RETURN
      END
C-----------------------------------------------------------------------
C  Subroutine READTL reads input necessary for Tensor LEED calculation.
C  VB, 06/98
C  NT0  : number of Tensor LEED beams to be considered
C  PQFEX: list of output beams for which tensor should be calculated
C  IFORM: produce tensor output formatted (1) or unformatted (0)
C
C  READTL subroutine is presently only retained because superstructures
C  not known to the fd calculation can be calculated using Tensor LEED.
C  Since this feature is currently not implemented, the TLEED list PQFEX
C  is superfluous. Either include 'TLEED superstructures' later or do not
C  read in PQFEX, instead use appropriate beams NPU from fd list SPQF.

      SUBROUTINE READTL(IFORM,NT0,PQFEX)

      INTEGER IFORM,NT0
      REAL PQFEX
      DIMENSION PQFEX(2,NT0)

      READ(5,'(I4)') IFORM

      DO IBEAM = 1,NT0

        READ(5,'(2F10.5)') PQFEX(1,IBEAM),PQFEX(2,IBEAM)

      ENDDO

      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE RINT_SIMPLE COMPUTES THE REFLECTED BEAM INTENSITIES FROM THE
!  (COMPLEX) REFLECTED AMPLITUDES. INCLUDED ARE ANGULAR PREFACTORS AND A
!  SYMMETRY-RELATED PREFACTOR (SO THAT THE RESULTING INTENSITY IS
!  DIRECTLY COMPARABLE WITH EXPERIMENT). RINT PRINTS OUT NON-ZERO
!  INTENSITIES AND PUNCHES OUT (IF REQUESTED) THE INTENSITIES OF THE
!  BEAMS SPECIFIED IN NPUC.
!
!  This is a simplified version that removes support for symmetry codes
!
!   N= NO. OF BEAMS AT CURRENT ENERGY.
!   WV= INPUT REFLECTED AMPLITUDES.
!   AT= OUTPUT REFLECTED INTENSITIES.
!   PQ= LIST OF RECIPROCAL-LATTICE VECTORS G (BEAMS).
!   PQF= SAME AS PQ, BUT IN UNITS OF THE RECIPROCAL-LATTICE CONSTANTS.
!   SYM= LIST OF SYMMETRY CODES FOR ALL BEAMS.
!   VV= VACUUM LEVEL ABOVE SUBSTRATE MUFFIN-TIN CONSTANT.
!   THETA,FI= POLAR AND AZIMUTHAL ANGLES OF INCIDENT BEAM.
!   MPU= NO. OF BEAM INTENSITIES TO BE PUNCHED.
!   NPUC  INDICATES WHICH BEAM INTENSITIES ARE TO BE PUNCHED.
!   EEV= CURRENT ENERGY IN EV ABOVE VACUUM LEVEL.
!   A= STRUCTURAL PARAMETER OR OTHER IDENTIFIER TO BE PUNCHED ON CARDS.
!   NPNCH=0  NO PUNCH DESIRED.
!   NPNCH.NE.0  PUNCH DESIRED.
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY IN HARTREES ABOVE SUBSTRATE MUFFIN-TIN CONSTANT.
!   VPI= IMAGINARY PART OF CURRENT ENERGY.
!   CK2,CK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
      SUBROUTINE RINT_SIMPLE(N,WV,AT,ATP,PQ,PQF,VV,THETA,FI,
     &                       MPU,NPUC,EEV,AP,NPNCH,IFILE)
      DIMENSION  WV(N), PQ(2,N), PQF(2,N), AT(N)
      DIMENSION  NPUC(MPU),ATP(MPU)
      COMPLEX WV
      COMMON  E, CK2, CK3, VPI

   35 FORMAT(1F7.2,1F7.4,4E14.5,/,100(5E14.5,/))                         140779
      AK = SQRT(AMAX1(2.0 * E - 2.0 * VV, 0.0))
      BK2 = AK * SIN(THETA) * COS(FI)
      BK3 = AK * SIN(THETA) * SIN(FI)

      C = AK * COS(THETA)  ! == k_perp in vacuum for incident beam
      DO J = 1, N
        AK2 = PQ(1, J) + BK2
        AK3 = PQ(2, J) + BK3
        A = 2.0 * E - 2.0 * VV - AK2 * AK2 - AK3 * AK3    ! == ||k_perp||^2
        IF (A.GT.0.0) THEN  ! beam is non-evanescent --> store its intensity for printing later
          A = SQRT(A)
          WR = REAL(WV(J))
          WI = AIMAG(WV(J))
!         AT IS REFLECTED INTENSITY (FOR UNIT INCIDENT CURRENT)
          AT(J) = (WR * WR + WI * WI) * A/C
        ELSE
          AT(J) = 0.0
        ENDIF
      ENDDO

      IF (NPNCH.EQ.0) RETURN  ! nothing to print to file

      DO K = 1,MPU
        NPC = NPUC(K)
        IF (NPC.EQ.0) THEN  ! This beam was not requested
          ATP(K) = 0.0
        ELSE
          ATP(K) = AT(NPC)
        ENDIF
      ENDDO

      WRITE(IFILE, 35) EEV, AP, (ATP(K), K=1,MPU)
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE RSMF_SIMPLE GENERATES REFLECTION AND TRANSMISSION MATRICES FOR A
!  SINGLE BRAVAIS-LATTICE LAYER BY PENDRY\S FORMULA. IT USES SYMMETRY
!  AS REQUESTED AND PRODUCES MATRICES FOR UP TO FOUR ORIGIN SHIFTS
!  TOWARDS A SYMMETRY AXIS OR PLANE. THERMAL VIBRATIONS ARE TAKEN INTO
!  ACCOUNT AND MAY BE ANISOTROPIC (APPROXIMATE TREATMENT).
!
!  Created 2021-09-21, MRiva
!  As compared to RSMF support for beam symmetry and registry shifts
!  is removed.
!
!   RA,TA,...,RD,TD= REFLECTION AND TRANSMISSION MATRICES (LETTERS R AND
!    T STAND FOR REFLECTION AND TRANSMISSION, RESP.; LETTERS A,B,!,D
!    CORRESPOND TO THE FOUR POSSIBLE ORIGIN SHIFTS, S1,S2,S3,S4,RESP.).
!   N= NO. OF BEAMS IN CURRENT BEAM SET.
!   NM,NP= DIMENSIONS OF DIFFRACTION MATRICES (.GE.N).
!   AMULT,CYLM= WORKING SPACE.
!   PQ= LIST OF BEAMS G.
!   SYM= LIST OF SYMMETRY CODES FOR ALL BEAMS.
!   NT= TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
!   FLMS= LATTICE SUMS FROM SUBROUTINE FMAT.
!   FLM= WORKING SPACE.
!   V  SUBLATTICE INFORMATION FROM SUBROUTINE SLIND.
!   NL= NO. OF SUBLATTICES IN LATTICE SUM (CF. SLIND).
!   N_OFFSET= OFFSET OF CURRENT BEAM SET IN LIST PQ.
!   ID= NO. OF ORIGIN SHIFTS TO BE CONSIDERED.
!   NLL= NO. OF SUBLATTICES TO BE CONSIDERED IN CURRENT USE OF MSMF
!    (NLL=1 FOR OVERLAYER, NLL=NL FOR SUBSTRATE LAYER).
!   AF= TEMPERATURE-INDEPENDENT ATOMIC T-MATRIX ELEMENTS.
!   CAF= TEMPERATURE-DEPENDENT ATOMIC T-MATRIX ELEMENTS
!   LM= LMAX+1.
!   LX,LXI= PERMUTATIONS OF THE (L,M) SEQUENCE, FROM SUBROUTINE LXGENT.
!   LMMAX= (LMAX+1)**2.
!   KLM= (2*LMAX+1)*(2*LMAX+2)/2.
!   XEV,XOD= WORKING SPACE.
!   LEV= (LMAX+1)*(LMAX+2)/2.
!   LOD= LMAX*(LMAX+1)/2.
!   YLM,YLME,YLMO,IPLE,IPLO= WORKING SPACES.
!   CLM= CLEBSCH-GORDON COEFFICIENTS, FROM SUBROUTINE CELMG.
!   NLM= DIMENSION OF CLM (SEE MAIN PROGRAM).
!   LAY= INDEX FOR CHOICE OF SYMMETRY CODES SYM(LAY,I) AND FOR CHOICE OF
!    SET OF ORIGIN SHIFTS (NORMALLY LAY=1 FOR OVERLAYER, LAY=2 FOR
!    SUBSTRATE LAYERS).
!  IN COMMON BLOCKS
!   E,VPI= CURRENT COMPLEX ENERGY.
!   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!   TV= AREA OF UNIT CELL.
!   EMACH= MACHINE ACCURACY.
!   LMAX= LARGEST VALUE OF L.
!   EPS,LITER  NOT USED.
!   /BT/  TEMPERATURE DATA PASSED TO MFOLD.
      SUBROUTINE  RSMF_SIMPLE(RA, TA, N, AMULT, CYLM, PQ, NT, FLMS,
     &                        FLM, V, NL, N_OFFSET, NLL, AF, CAF, LM,
     &                        LX, LXI, LMMAX, KLM, XEV, XOD, LEV, LOD,
     &                        YLM, YLME, YLMO, IPLE, IPLO, CLM, NLM,
     &                        XEVST, XODST)

      INTEGER  LX, LXI, L1
      COMPLEX XODST(LOD,LOD), XEVST(LEV,LEV)
      COMPLEX  AF, XEV, XOD, RA, TA, CYLM, YLM, YLME, YLMO,
     &         AMULT, XA, YA, ST, CF, CT, CI, RU, CZ, MFOLD_RT, AM
      COMPLEX  CAF, FLMS, FLM
      COMPLEX SQRT, EXP
      DIMENSION  MFOLD_RT(2)
      DIMENSION  RA(NT, NT), TA(NT, NT), AMULT(N)
      DIMENSION  XEV(LEV, LEV), XOD(LOD, LOD), CYLM(NT, LMMAX),
     &           PQ(2, NT), YLM(LMMAX), CLM(NLM), AF(LM), YLME(LEV),
     &           YLMO(LOD)
      DIMENSION  IPLE(LEV), IPLO(LOD), V(NL, 2)
      DIMENSION  CAF(LM), FLMS(NL, KLM), FLM(KLM)
      DIMENSION  LX(LMMAX), LXI(LMMAX)
      COMMON  E, AK2, AK3, VPI, TV
      COMMON /MS/LMAX,EPS,LITER

      PI = 3.14159265
      L1 = LM
      LEV1 = LEV + 1
      CI = CMPLX(0.0,1.0)
      CZ = CMPLX(0.0,0.0)
      RU = CMPLX(1.0,0.0)

      YA = CMPLX(2.0 * E, - 2.0 * VPI + 0.000001)
      YA = SQRT(YA)
      DO IG = 1, N
        JG = IG + N_OFFSET
        BK2 = PQ(1, JG) + AK2
        BK3 = PQ(2, JG) + AK3
        C = BK2 * BK2 + BK3 * BK3
        XA = CMPLX(2.0 * E - C, - 2.0 * VPI + 0.000001)
        XA = SQRT(XA)
        B = 0.0
        CF = RU
        IF (C.GT.1.0E-7) THEN ! off-normal beam
          B = SQRT(C)
          CF = CMPLX(BK2/B, BK3/B)
        ENDIF
        CT = XA/YA
        ST = B/YA
!       PREFACTOR FOR REFLECTION AND TRANSMISSION MATRIX ELEMENTS
        AMULT(IG) =  - 8.0 * PI * PI * CI/(YA * TV * XA)
!       COMPUTE AND STORE SPHERICAL HARMONICS FOR ALL BEAM DIRECTIONS
        CALL  SPHRM_MOD(LMAX, YLM, LMMAX, CT, ST, CF)
        DO K = 1, LMMAX
          CYLM(IG, K) = YLM(K)
        ENDDO
      ENDDO

!     PERFORM SUM OVER SUBLATTICES. << Could be done by constructing XA=XA(JS) and then CGEMV with FLMS.T (or swap FLMS indices)
      DO K = 1, KLM
        FLM(K) = CZ
      ENDDO
      BK2 = PQ(1, 1 + N_OFFSET)
      BK3 = PQ(2, 1 + N_OFFSET)
      DO JS = 1, NLL
        ABR = V(JS, 1) * (BK2 * COS(V(JS, 2)) + BK3 * SIN(V(JS, 2)))
        XA = EXP(ABR * CI)
        DO I = 1, KLM
          FLM(I) = FLM(I) + FLMS(JS, I) * XA
        ENDDO
      ENDDO

!     COMPUTE INTRALAYER MULTIPLE-SCATTERING MATRIX X, SPLIT ACCORDING TO
!     L+M=EVEN (XEV) AND L+M=ODD (XOD)
      CALL  XM (FLM,XEV,XOD,LEV,LOD,AF,CAF,LM,LX,LXI,LMMAX,KLM,CLM,NLM)
      DO J = 1, LEV
        XEV(J, J) = XEV(J, J) - RU
      ENDDO
      DO J = 1, LOD
        XOD(J, J) = XOD(J, J) - RU
      ENDDO
!
!
!********************************************************************
! PERTURBATIVE LEED; STORE AWAY THE TRANSPOSE(????) OF (1-X) INTO
! XODST & XEVST.

      DO J = 1, LEV
        DO I = 1, LEV
          XEVST(I, J) = -XEV(I, J)
        ENDDO
      ENDDO

      DO J = 1, LOD
        DO I = 1, LOD
          XODST(I, J) = -XOD(I, J)
        ENDDO
      ENDDO
!*********************************************************************
!
!     PREPARE INVERSION OF 1-X (- SIGN IS PUT IN AMULT) BY GAUSSIAN
!     ELIMINATION
      CALL CGETRF(LEV, LEV, XEV, LEV, IPLE, INFO)
      CALL CGETRF(LOD, LOD, XOD, LOD, IPLO, INFO)

!     START LOOP OVER SCATTERED BEAMS
      DO JGP = 1, N
        JGPS = JGP + N_OFFSET

!       PREPARE QUANTITIES TO BE MULTIPLIED FROM LEFT BY INVERSE OF 1-X
        JLM = 1
        ST = RU
        DO L = 1, LM
          LL = L + L - 1
          DO JM = 1, LL    ! This loop is weird: JM not used in loop TODO
            YLM(JLM) = CYLM(JGP, JLM) * AF(L)/ST
            JLM = JLM + 1
          ENDDO
          ST = ST * CI
        ENDDO
!       Change ordering of (L, M) pairs
        DO JLM = 1, LEV
          JJ = LX(JLM)
          YLME(JLM) = YLM(JJ)
        ENDDO
        DO JLM = LEV1, LMMAX
          JJ = LX(JLM)
          YLMO(JLM - LEV) = YLM(JJ)
        ENDDO

!       Complete inversion and multiplication by back-substitution
        CALL CGETRS('N', LEV, 1, XEV, LEV, IPLE, YLME, LEV, INFO)
        CALL CGETRS('N', LOD, 1, XOD, LOD, IPLO, YLMO, LOD, INFO)
!       Reorder (L, M) pairs
        DO JLM = 1, LEV
          JJ = LX(JLM)
          YLM(JJ) = YLME(JLM)
        ENDDO
        DO JLM = LEV1, LMMAX
          JJ = LX(JLM)
          YLM(JJ) = YLMO(JLM - LEV)
        ENDDO

        AM = AMULT(JGP)
!       Start loop over incident beams
        DO JG = 1, N
!         Complete computation in MFOLD (include spherical harmonics)
          CALL  MFOLD_SIMPLE (JG, L1, YLM, CYLM, LMMAX, NT, MFOLD_RT)
          RA(JGPS, JG + N_OFFSET) = MFOLD_RT(1) * AM
          TA(JGPS, JG + N_OFFSET) = MFOLD_RT(2) * AM
        ENDDO
!       Add unscattered beam to transmission
        TA(JGPS, JGPS) = TA(JGPS, JGPS) + RU
      ENDDO  ! JGP, i.e., loop over scattered beams

      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE RTINV_SIMPLE computes reflection and transmission matrices for
!  an atomic layer consisting of NLAY subplanes. A Beeby-type matrix
!  inversion is used (a perturbation treatment is done in subroutine
!  MPERTI -- which does not exist any longer).
!  This is a simplified version of the original RTINV that removes
!  support for (1) Beam symmetry codes (which were never used anyway),
!  (2) registry shifts SO(I) and SS(I) (which were also never fully
!  implemented), and (3) NOPT==1.
!
!  Parameters
!  ----------
!  R_TOP, T_TOP, R_BOT, T_BOT
!      Reflection and transmission matrices for beams incident from
!      above (_TOP) or below (_BOT) the composite layer. R/T(i, j) is
!      reflection/transmission of incident beam j into beam i. The
!      complex amplitudes of both beams are taken right next to the
!      composite layer, i.e., R_TOP(i, j): beams right above;
!      T_TOP(i, j): j right above, i right below; R_BOT(i, j): both
!      beams below; T_BOT(i, j): j right below, i right above.
!
!
!  SPECIFICATIONS NOT FOUND BELOW CAN BE FOUND IN THE SPECIFICATIONS
!  OF SUBROUTINE MPERTI (van Hove/Tong 1979), WHICH HAS THE SAME STRUCTURE
!  AS MTINV.
!  NOTE  REAL AND IMAGINARY PARTS OF MUFFIN-TIN CONSTANT ARE ASSUMED
!  CONSTANT THROUGHOUT LAYER (I.E. BETWEEN OUTERMOST NUCLEAR PLANES).
!   TAUG,TAUGM HAVE DIMENSION LMT=NTAU*LMMAX HERE.
!   TH(LMNI,LMNI2) IS SQUARE HERE WITH LMNI2=LMNI=NLAY*LMMAX.
!   IPL IS ADDITIONAL WORKING SPACE.
!
!  Author of edits: Michele Riva, 2021-09-02

!  REMOVED from arguments:
!    * extra reflection/transmission for other registry shifts
!    * SYM
!    * ID
!    * LX (unused)
!    * NM, NP --> both always == NT
!  (others)
!    * /SL/ common block
      SUBROUTINE RTINV_SIMPLE(
     &                 R_TOP,T_TOP,R_BOT,T_BOT, N,
     &                 AMULT, CYLM, PQ, NT, FLMS, FLM, NL,
     &                 LXI, LT, LXM, LMMAX, KLM, XEV, LEV, LEV2,
     &                 TAU, LMT, TAUG, TAUGM, CLM, NLM, POS, POSS, MGH,
     &                 NLAY, DRL, SDRL, NUGH, NGEQ, NGOL, NLAY2, TEST,
     &                 GH, LMG, RG, RG_PERP, TS, LMN, TG, LM2N, VT,
     &                 CAA, NCAA, TH, LMNI, IPL, TSF, LPS, LPSS, NORD,
     &                 NNSUB, LMAX1, XH, HGHD, L2M, RG1, RG2, TSTORE)
      COMPLEX XEV, R_TOP, T_TOP, R_BOT, T_BOT, CYLM, AMULT, XA, YA, ST,
     +        CF, CT, CI, RU, CZ, AM, FLMS, FLM, MFOLT_RT
      COMPLEX SQRT,EXP
      COMPLEX TH(LMNI, LMNI), TSTORE(2, LMN, NT)
      COMPLEX TAU(LMT, LEV), GH(LMG, LMMAX), RG(2, NLAY, N)            ! May be worth to have 2 as last?
      COMPLEX RG_PERP(N)
      COMPLEX TS(LMN), TG(2, LM2N), VT(LM2N), TAUG(LMT), TAUGM(LMT)
      DIMENSION  FLM(KLM)  ! Used in TAUMAT as working space?            BUG? SIZE WAS 4!
      DIMENSION  R_TOP(NT, NT), T_TOP(NT, NT), R_BOT(NT, NT),
     +           T_BOT(NT, NT), AMULT(N)
      DIMENSION  XEV(LEV, LEV2), CYLM(NT, LMMAX), PQ(2, NT)            ! Swap CYLM indices?
      DIMENSION  CLM(NLM), GP(2), CAA(NCAA)
      DIMENSION  FLMS(NL, KLM)
      DIMENSION  LX(LMMAX), LXI(LMMAX), LT(LMMAX), LXM(LMMAX)
      DIMENSION  POS(NLAY, 3), POSS(NLAY, 3)
      DIMENSION  MGH(NLAY, NLAY), DRL(NLAY2, 3), SDRL(NLAY2, 3)
      DIMENSION  NUGH(NLAY2), NGEQ(NLAY2), NGOL(NLAY2)
      DIMENSION  TEST(NLAY2), IPL(LMNI), MFOLT_RT(2)
CVB
      INTEGER NNSUB, LMAX1
      COMPLEX TSF
      DIMENSION TSF(NNSUB, LMAX1)
      INTEGER LPS, LPSS, NORD
      DIMENSION LPS(NNSUB), LPSS(NNSUB), NORD(NNSUB)
      COMPLEX XH(LEV)
      INTEGER L2M
      COMPLEX HGHD(L2M)
      COMPLEX RG1(NLAY), RG2(NLAY), TEMP

      COMMON  E, AK2, AK3, VPI
      COMMON /MFB/ GP, L1
      COMMON /MS/  LMAX, EPS1, LITER
      COMMON /MPT/ NA, NS, ID, LAY, LM, NTAU, TST, TV, DCUT, NOPT, NEW ! NOPT always == 1 for Tensor LEED; ID == 1 (registry shifts), not used anymore

      PI = 3.14159265
      L1 = LM
      LMS = L2M * L2M
      LOD = LMMAX - LEV
      LEE = (LMAX/2 + 1)**2
      LOE = ((LMAX - 1)/2 + 1)**2
      CI = CMPLX(0.0, 1.0)
      CZ = CMPLX(0.0, 0.0)
      RU = CMPLX(1.0, 0.0)

      IF (LAY.EQ.1) THEN
        NLL = 1     ! Overlayer
      ELSE
        NLL = NL    ! Substrate (bulk)
      ENDIF

!     Check whether TAU needs recomputing since last call to RTINV_SIMPLE
!     NB: NEW is always == 1, so the next bit is always executed.
      IF (NEW.EQ.1) THEN
!       Compute TAU for each chemical element
        CALL TAUMAT(TAU, LMT, NTAU, XEV, LEV, LEV2, LOD, TSF, LMMAX,
     &              LMAX, FLMS, NL, KLM, LM, CLM, NLM, LXI, NT, PQ, NA,
     &              NLL, FLM, NNSUB, LMAX1, XH)
      ENDIF

!     Sort subplanes according to increasing position along +x axis,
!     i.e., from vacuum side towards deeper into the solid.
!     Figure out which old GH\S can be kept, which must be recomputed
!     and which are mutually identical
      CALL SRTLAY(POS, POSS, LPS, LPSS, MGH, NLAY, DRL, SDRL, NLAY2,
     &            NUGH, NGEQ, NGOL, NEW, NORD)

!     Compute (or copy) the GH matrices, i.e., the interlayer
!     propagators in (LM)-space.
      CALL GHMAT(GH, LMG, LMMAX, MGH, NLAY, NUGH, NGEQ, NGOL, NLAY2,
     &           TST, TEST, VT, L2M, VT, LM, TG, LMS, DRL, TV, LXM,
     &           LEV, DCUT, CAA, NCAA, LAY, HGHD, PQ(1, 1 + NA))

!     Create matrix TH to be inverted in Beeby\s scheme
      CALL THMAT(TH, LMNI, GH, LMG, LMMAX, MGH, NLAY, TAU, LMT,
     &           LEV, LPSS)

!     Prepare TH for inversion (Gaussian elimination)
      CALL CGETRF(LMNI, LMNI, TH, LMNI, IPL, INFO)

      YA = CMPLX(2.0 * E, - 2.0 * VPI + 0.000001)
      YA = SQRT(YA)
      DO IG = 1, N
!       Generate plane-wave propagators for use as prefactors for GH and TH
        JG = IG + NA
        BK2 = PQ(1, JG) + AK2
        BK3 = PQ(2, JG) + AK3
        C = BK2 * BK2 + BK3 * BK3
        XA = CMPLX(2.0 * E - C, - 2.0 * VPI + 0.000001)
        XA = SQRT(XA)
        RG_PERP(IG) = EXP(CI * XA * (POSS(NLAY, 1) - POSS(1, 1)))
        DO I = 1, NLAY
          X = BK2 * POSS(I, 2) + BK3 * POSS(I, 3)
          RG(1, I, IG) = EXP(CI * (XA * POSS(I, 1) + X))
          RG(2, I, IG) = EXP(CI * (XA * POSS(I, 1) - X))
        ENDDO
!       If NEW=1 (i.e. new energy) compute new spherical harmonics for
!       the beam directions, using SPHRM. NB: NEW is always == 1!
        IF (NEW.EQ.1) THEN
          B = 0.0
          CF = RU
          IF (C.GT.1.0E-7) THEN
            B = SQRT(C)               ! Norm of in-plane wave vector
            CF = CMPLX(BK2/B, BK3/B)  ! Complex azimuthal incidence angle
          ENDIF
          CT = XA/YA  ! Complex cosine of polar incidence angle
          ST = B/YA   ! Complex sine of polar incidence angle

!         Generate prefactor of reflection and transmission matrix elements
          AMULT(IG) =  - 16.0 * PI * PI * CI / (TV * XA)
!         Prepare and store the spherical harmonics
          CALL  SPHRM_MOD(LMAX, VT, LMMAX, CT, ST, CF)
          DO K = 1, LMMAX
            CYLM(IG, K) = VT(K)
          ENDDO
        ENDIF
      ENDDO

!     Start loop over incident beams
      DO JGP = 1, N
!       Generate quantities TAUG, TAUGM into which inverse of TH
!       will be multiplied
        CALL TAUY_SIMPLE(TAUG, TAUGM, LMT, TAU, LMT, LEV, CYLM, NT,
     &                   LMMAX, LT, NTAU, LOD, LEE, LOE, JGP)

        JGPS = JGP + NS
        GP(1) = PQ(1, JGP + NA)
        GP(2) = PQ(2, JGP + NA)
!       First consider incidence towards +x, i.e., beams propagating
!       towards the solid
        INC = 1
  285   CONTINUE  ! Will then continue from here with INC=-1

!       Include appropriate plane-wave propagating factors into TAUG/TAUGM
        DO I = 1, NLAY
          IN = (I - 1) * LMMAX
          LP = (LPSS(I) - 1) * LMMAX
          IF (INC.EQ.1) THEN
            ST = RG(1, I, JGP)
            DO K = 1, LMMAX
              TS(IN + K) = TAUG(LP + K) * ST
            ENDDO
          ELSE
            ST = RU / RG(2, I, JGP)
            DO K = 1, LMMAX
              TS(IN + K) = TAUGM(LP + K) * ST
            ENDDO
          ENDIF
        ENDDO

!       Do left multiplication with inverse of TH
        CALL CGETRS('N', LMNI, 1, TH, LMNI, IPL, TS, LMNI, INFO)

!       Include further plane-wave propagating factors, and store
!       the intra-layer scattering for computing Tensors.
        DO I = 1, NLAY
          IN = (I - 1) * LMMAX
          INSORT = (NORD(I) - 1) * LMMAX
          IF (INC.EQ.1) THEN
            ST = RU / RG(1, I, JGP)
          ELSE IF (INC.eq.-1) THEN
            ST = RG(2, I, JGP)
          END IF
!         121190 (ORIGINAL AUTHOR: P. ROUS)
!         take back sublayer reordering by SRTLAY for output via NORD here
          DO K = 1, LMMAX
            IF (INC.EQ.1) THEN
              TSTORE(1, INSORT + K, JGP) = 4.0 * PI * TS(IN + K)
            ELSE
              TSTORE(2, INSORT + K, JGP) = 4.0 * PI * TS(IN + K)
     +                                    * RG_PERP(JGP)
            ENDIF
            TS(IN + K) = TS(IN + K) * ST
          ENDDO
        ENDDO

!       Start loop over scattered beams
        DO JG = 1, N
!         Complete computation in MFOLT_SIMPLE
          CALL MFOLT_SIMPLE(JG, NA, AK2, AK3, CYLM, N, LMMAX, PQ,
     &                      NT, TS, LMN, RG, RG1, RG2, NLAY, JGP,
     &                      LXM, MFOLT_RT, INC, POSS)
          JGS = JG + NS
          AM = AMULT(JG)
!         Put final matrix elements in proper location
          IF (INC.EQ.-1) THEN
            R_BOT(JGS, JGPS) = MFOLT_RT(1) * AM
            T_BOT(JGS, JGPS) = MFOLT_RT(2) * AM
          ELSE
            R_TOP(JGS, JGPS) = MFOLT_RT(1) * AM
            T_TOP(JGS, JGPS) = MFOLT_RT(2) * AM
          ENDIF
        ENDDO

        IF (INC.EQ.1) THEN
!         Add unscattered plane wave to transmission from above
          T_TOP(JGPS, JGPS) = T_TOP(JGPS, JGPS) + RG_PERP(JGP)
!         And restart loop for beams moving from solid to vacuum
          INC = -1
          GO TO 285
        ELSE
!         Add unscattered plane wave to transmission from below
          T_BOT(JGPS, JGPS) = T_BOT(JGPS, JGPS) + RG_PERP(JGP)
        ENDIF
      ENDDO  ! END DO over incident beams (JGP)
      RETURN
      END


C-----------------------------------------------------------------------
C  SUBROUTINE SB COMPUTES SPHERICAL BESSEL FUNCTIONS.
C   X= COMPLEX ARGUMENT OF BESSEL FUNCTIONS.
C   HH= OUTPUT COMPLEX BESSEL FUNCTIONS.
C   N3= LMAX+1.
      SUBROUTINE  SB (X, HH, N3)
      COMPLEX X
      COMPLEX A, B, C, HH
      DIMENSION  HH(N3)
  330 FORMAT(1H0,11HSB INFINITE)
      A = (0.0,1.0)
      C = X
      B = A * C
      F = CABS(X)
      IF (DBLE(F)-1.0E-38)  1540, 1540, 1520
 1520 HH(1) = EXP(B) * ( - A)/C
      HH(2) = EXP(B) * (1.0/( - C) - A/C**2)
      DO 1530 J = 3, N3
      HH(J) = (2.0 * (J - 2) + 1.0)/C * HH(J - 1) - HH(J - 2)
 1530 CONTINUE
      GO TO 1550
 1540 CONTINUE
C     WRITE(6,330)
 1550 RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE SH COMPUTES SPHERICAL HARMONICS.
C   NHARM= LMAX+1.
C   Z= COMPLEX ARGUMENT COS(THETA).
C   FI= AZIMUTHAL ANGLE.
C   Y= OUTPUT COMPLEX SPHERICAL HARMONICS.
      SUBROUTINE  SH (NHARM, Z, FI, Y)
      REAL FI
      COMPLEX  A, ANOR, Q1A, YSTAR, YY, ZNW
      REAL*8     AM              ! Added 2021-10-10 MRiva. Prevent overflow to infinity when calculating factorials below. Increases NHARM range from 28 to 150.
      COMPLEX*16 BB, BB1, AA     ! Added 2021-10-10 MRiva. Prevent overflow and NaNs.
      COMPLEX  Z, Y(NHARM,NHARM) ! MRiva: still unclear whether one would also like Y to be COMPLEX*16
      AM = 1.0
      A = (0.0,1.0) * FI
      RZ = REAL(Z)
      ZNW = SQRT((1.0,0.0) - Z * Z)
      ANORA = 0.2820948
      YY = (1.0,0.0)

      DO L = 1, NHARM
        BM = FLOAT(L) - 1.0
        Q1A = EXP(BM * A)
        ANOR = ANORA * Q1A * ( - 1.0)**(L + 1)
        ANORA = ANORA * SQRT(1.0 + 0.5/(BM + 1.0))/(2.0 * BM + 1.0)
        IF((ABS(ABS(RZ) - 1.0) .GE. 1.0E-10) .OR. L .NE. 1) THEN
          YY = ZNW**INT(BM + 0.1)
        ENDIF
        AA = AM * YY
        AM = AM * (2.0 * BM + 1.0)
        BB = AM * Z * YY
        DO LL = L, NHARM
          BN = FLOAT(LL)
          Y(LL,L) = AA * ANOR
          IF (L .NE. 1) THEN
            YSTAR = Y(LL,L)/(Q1A * Q1A)
            Y(L - 1,LL) = YSTAR
          ENDIF
          BB1 = BB
          BB = BB1 * Z + (BB1 * Z - AA) * (BN + BM)/(BN - BM + 1.0)
          AA = BB1
          ANOR = ANOR*SQRT((1.0 + 1.0/(BN - 0.5))
     +                     * ((BN - BM)/(BN + BM)))
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
COMMENT SLIND SETS UP A MATRIX JJS(JS,J) CONTAINING DETAILS OF
C       HOW THE SUBLATTICES JS ARE TRANSFORMED INTO ONE
C       ANOTHER BY ROTATIONS THROUGH J*2.0*PI/IDEG.  V(JS,2)
C       CONTAINS THE ADDING VECTORS DEFINING THE SUBLATTICES
C       JS IN POLAR FORM
C  AUTHOR  PENDRY.
C   V= VECTORS POINTING TO DIFFERENT SUBLATTICES (IN POLAR COORDINATES).
C   VL= WORKING SPACE.
C   JJS= OUTPUT RELATIONSHIP OF SUBLATTICES UNDER ROTATION.
C   NL= NO. OF SUBLATTICES.
C   IDEG= DEGREE OF SYMMETRY OF LATTICE  IDEG-FOLD ROTATION AXIS.
C   EPSD= SMALL NUMBER.
C  IN COMMON BLOCKS
C   BR1,BR2= BASIS VECTORS OF SUBSTRATE LAYER LATTICE.
C   AR1,AR2,RAR1,RAR2= BASIS VECTORS OF SUPERLATTICE IN DIRECT AND
C    RECIPROCAL SPACE.
C   NL1,NL2= SUPERLATTICE CHARACTERIZATION (SEE MAIN PROGRAM).
      SUBROUTINE  SLIND (V, VL, JJS, NL, IDEG, EPSD)
      COMPLEX  VL, VLA, VLB, CI
      COMPLEX EXP
      DIMENSION  V(NL,2), JJS(NL,IDEG), VL(NL,2), AR1(2), AR2(2)
      DIMENSION  BR1(2), BR2(2), RAR1(2), RAR2(2)
      COMMON  /SL/BR1, BR2, AR1, AR2, RAR1, RAR2, NL1, NL2
C
      CI = CMPLX(0.0,1.0)
      PI = 3.14159265
C  SET UP VECTORS V DEFINING SUBLATTICES AND QUANTITIES VL FOR LATER
C  REFERENCE.
      I = 1
      S1 = 0.0
      DO 560 J = 1, NL1
      S2 = 0.0
      DO 550 K = 1, NL2
      ADR1 = S1 * BR1(1) + S2 * BR2(1)
      ADR2 = S1 * BR1(2) + S2 * BR2(2)
      V(I,1) = SQRT(ADR1 * ADR1 + ADR2 * ADR2)
      V(I,2) = 0.0
      IF (V(I,1))  540, 540, 530
  530 V(I,2) = ATAN2(ADR2,ADR1)
  540 VL(I,1) = EXP(CI * (ADR1 * RAR1(1) + ADR2 * RAR1(2)))
      VL(I,2) = EXP(CI * (ADR1 * RAR2(1) + ADR2 * RAR2(2)))
      I = I + 1
  550 S2 = S2 + 1.0
  560 S1 = S1 + 1.0
C
C  ROTATE EACH VECTOR V AND FIND TO WHICH V IT BECOMES EQUIVALENT IN
C  TERMS OF THE QUANTITIES VL. THIS EQUIVALENCE MEANS BELONGING TO THE
C  SAME SUBLATTICE.
      AINC = 2.0 * PI/FLOAT(IDEG)
      DO 590 I = 1, NL
      ADR = V(I,1)
      ANG = V(I,2)
      DO 590 K = 1, IDEG
      ANG = ANG + AINC
      CANG = COS(ANG)
      SANG = SIN(ANG)
      A = RAR1(1) * CANG + RAR1(2) * SANG
      B = RAR2(1) * CANG + RAR2(2) * SANG
      VLA = EXP(CI * ADR * A)
      VLB = EXP(CI * ADR * B)
      DO 570 J = 1, NL
      TEST = CABS(VLA - VL(J,1)) + CABS(VLB - VL(J,2))
      IF (TEST-5.0*EPSD)  580, 580, 570
  570 CONTINUE
  580 JJS(I,K) = J
  590 CONTINUE
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE SPHRM_MOD COMPUTES SPHERICAL HARMONICS.
!  AUTHOR  PENDRY.
!
!  Edited: 2021-09-22 MRiva. Get rid of GOTOs and add ENDDOs
!
!
!   LMAX= LARGEST VALUE OF L.
!   YLM= OUTPUT COMPLEX SPHERICAL HARMONICS.
!   LMMAX= (LMAX+1)**2.
!   CT= COS(THETA) (COMPLEX).
!   ST= SIN(THETA) (COMPLEX).
!   CF= EXP(I*FI).
      SUBROUTINE  SPHRM_MOD(LMAX, YLM, LMMAX, CT, ST, CF)
      COMPLEX  YLM
      COMPLEX  CT, ST, CF, SF, SA
!  (21=LMAX+1; 441=(LMAX+1)**2   FOR LMAX=20 - if you need more (what for???),
!   please introduce variable dimensions, carrying the respective quantities through
!   all calls)

      DIMENSION  FAC1(21), FAC3(21), FAC2(441), YLM(LMMAX)

      IF (LMAX.gt.20) THEN
        write(6,*) "Dimensions for auxiliary arrays FAC1,FAC2,FAC3 in"
        write(6,*) "subroutine sphrm too low for LMAX.gt.20 . Remedy:"
        write(6,*) "increase dim. to"
        write(6,*) "FAC1(LMAX+1),FAC2((LMAX+1)**2),FAC3(LMAX+1)."
        STOP
      END IF

      PI = 3.14159265
      LM = 0
      CL = 0.0
      A = 1.0
      B = 1.0
      ASG = 1.0
      LL = LMAX + 1
      DO L = 1, LL
        FAC1(L) = ASG * SQRT((2.0 * CL + 1.0) * A/(4.0 * PI * B * B))
        FAC3(L) = SQRT(2.0 * CL)
        CM =  - CL
        LN = L + L - 1
CDIR$ SHORTLOOP
        DO M = 1, LN
          LO = LM + M
          FAC2(LO) = SQRT((CL + 1.0 + CM) * (CL + 1.0 - CM)
     +                    /((2.0 * CL + 3.0) * (2.0 * CL + 1.0)))
          CM = CM + 1.0
        ENDDO
        CL = CL + 1.0
        A = A * 2.0 * CL * (2.0 * CL - 1.0)/4.0
        B = B * CL
        ASG =  - ASG
        LM = LM + LN
      ENDDO

      LM = 1
      CL = 1.0
      ASG =  - 1.0
      SF = CF
      SA = CMPLX(1.0, 0.0)
      YLM(1) = CMPLX(FAC1(1), 0.0)
CDIR$ SHORTLOOP
      DO L = 1, LMAX
        LN = LM + L + L + 1
        YLM(LN) = FAC1(L + 1) * SA * SF * ST
        YLM(LM + 1) = ASG * FAC1(L + 1) * SA * ST/SF
        YLM(LN - 1) =  - FAC3(L + 1) * FAC1(L + 1) * SA * SF * CT/CF
        YLM(LM + 2) = ASG * FAC3(L + 1) * FAC1(L + 1) * SA * CT * CF/SF
        SA = ST * SA
        SF = SF * CF
        CL = CL + 1.0
        ASG =  - ASG
        LM = LN
      ENDDO

      LM = 1
      LL = LMAX - 1
      DO L = 1, LL
        LN = L + L - 1
        LM2 = LM + LN + 4
        LM3 = LM - LN
CDIR$ SHORTLOOP
        DO M = 1, LN
          LO = LM2 + M
          LP = LM3 + M
          LQ = LM + M + 1
          YLM(LO) =  - (FAC2(LP) * YLM(LP) - CT * YLM(LQ))/FAC2(LQ)
        ENDDO
        LM = LM + L + L + 1
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE SRTLAY IS USED BY SUBROUTINES MPERTI AND MTINV TO REORDER
C  THE SUBPLANES OF A COMPOSITE LAYER ACCORDING TO INCREASING POSITION
C  ALONG THE +X AXIS (REORDERING THE CHEMICAL ELEMENT ASSIGNMENT AS
C  WELL). IF NEW=-1 IT FINDS WHICH RESULTS FOR GH WILL BE AVAILABLE BY
C  SIMPLY COPYING OLD ONES (OF A PREVIOUS CALL TO MPERTI OR MTINV),
C  WHICH ONES NEED TO BE COMPUTED AFRESH AND WHICH ONES CAN BE OBTAINED
C  BY COPYING OTHER NEW ONES.
C   POS= INPUT SUBPLANE POSITIONS (ATOMIC POSITIONS IN UNIT CELL).
C   POSS= OUTPUT REORDERED POS.
C   LPS= INPUT CHEMICAL ELEMENT ASSIGNMENT.
C   LPSS= OUTPUT REORDERED LPS.
C         dimension should be NNSUB, NLAY is sufficient.
C   MGH= OUTPUT  INDICATES ORGANIZATION OF GH IN MPERTI OR MTINV.
C    MGH(I,J)=K MEANS GH(I,J) (I,J= SUBPLANE INDICES) IS TO BE FOUND IN
C    K-TH POSITION IN COLUMNAR MATRIX GH.
C   NLAY= NO. OF SUBPLANES.
C   DRL= OUTPUT INTERPLANAR VECTORS.
C   SDRL= STORAGE FOR OLD DRL FROM PREVIOUS CALL TO MPERTI OR MTINV.
C   NLAY2= NLAY*(NLAY-1)/2.
C   NUGH= OUTPUT  NUGH(K)=1 MEANS THE GH(I,J) (I,J= SUBPLANE INDICES)
C    FOR WHICH MGH(I,J)=K MUST BE COMPUTED AFRESH. NUGH(K)=0 MEANS NO
C    COMPUTATION NECESSARY FOR THAT GH(I,J).
C   NGEQ= OUTPUT  NGEQ(K)=L MEANS THE GH(I,J) FOR WHICH MGH(I,J)=K CAN
C    BE COPIED FROM THE NEWLY CREATED GH(M,N) FOR WHICH MGH(M,N)=L.
C    NGEQ(K)=0 MEANS NO COPYING TO BE DONE FOR THAT GH(I,J).
C   NGOL= OUTPUT  SAME AS NGEQ, BUT FOR COPYING FROM OLD EXISTING VALUES
C    OF GH(M,N).
C   NORD stores sublayer order so that reordering can be traced back in TLEED)
C        (note srtlay is superfluous in itself but kept so that DRL, SDRL can
C        be handled more conveniently - should be removed in a future version)
      SUBROUTINE SRTLAY(POS,POSS,LPS,LPSS,MGH,NLAY,DRL,SDRL,NLAY2,
     1NUGH,NGEQ,NGOL,NEW,NORD)
      DIMENSION POS(NLAY,3),POSS(NLAY,3),POSA(3),DRL(NLAY2,3),
     1SDRL(NLAY2,3)
      DIMENSION LPS(NLAY),LPSS(NLAY),MGH(NLAY,NLAY),NUGH(NLAY2),
     1NGEQ(NLAY2),NGOL(NLAY2)
      INTEGER NORD,NORDA
      DIMENSION NORD(NLAY)

      DO 1 I=1,NLAY
      LPSS(I)=LPS(I)
      NORD(I)=I
      DO 1 J=1,3
1     POSS(I,J)=POS(I,J)
C  ANALYSE ORDER OF SUBPLANE POSITIONS ALONG X-AXIS AND REORDER IN
C  ASCENDING POSITION ALONG +X-AXIS (EQUALLY POSITIONED SUBPLANES ARE
C  NOT PERMUTED)
      NLAY1=NLAY-1
      DO 7 I=1,NLAY1
      KM=I
      PM=POSS(I,1)
      DO 2 K=I+1,NLAY
      IF (POSS(K,1).GE.PM) GO TO 2
      PM=POSS(K,1)
      KM=K
2     CONTINUE
      IF (KM.EQ.I) GO TO 7
      DO 3 J=1,3
3     POSA(J)=POSS(KM,J)
      LPSA=LPSS(KM)
      NORDA=NORD(KM)
      DO 5 K=KM,I+1,-1
      DO 4 J=1,3
4     POSS(K,J)=POSS(K-1,J)
      NORD(K)=NORD(K-1)
C  REORDER CHEMICAL ASSIGNMENTS CORRESPONDINGLY
5     LPSS(K)=LPSS(K-1)
      DO 6 J=1,3
6     POSS(I,J)=POSA(J)
      LPSS(I)=LPSA
      NORD(I)=NORDA
7     CONTINUE
Caws  SET ORIGIN OF COMPOSITE LAYER BRUTALLY ZERO
      DO I = 2,NLAY,1
        POSS(I,1) = POSS (I,1) - POSS (1,1)
      ENDDO
      POSS(1,1) = 0.0000
C  GENERATE INTERPLANAR VECTORS DRL
      IN=1
      DO 14 I=1,NLAY
      DO 14 J=I,NLAY
      MGH(I,J)=0
      MGH(J,I)=0
      IF (I-J) 8,14,8
8     DO 9 K=1,3
9     DRL(IN,K)=POSS(I,K)-POSS(J,K)
C  ASSIGN ORGANIZATION OF COLUMNAR MATRIX GH
      MGH(I,J)=IN
      MGH(J,I)=IN+NLAY2
      NUGH(IN)=1
      NGEQ(IN)=0
      NGOL(IN)=0
      IF (IN.EQ.1) GO TO 13
C  FIND EQUAL INTERPLANAR VECTORS, RECORD FINDINGS IN NGEQ
      IN1=IN-1
      DO 12 M=1,IN1
      DO 10 M1=1,3
10    POSA(M1)=ABS(DRL(IN,M1)-DRL(M,M1))
      IF ((POSA(1)+POSA(2)+POSA(3))-0.001) 11,11,12
11    NGEQ(IN)=M
      NUGH(IN)=0
      GO TO 13
12    CONTINUE
13    CONTINUE
      IN=IN+1
14    CONTINUE
C  IF OLD VALUES OF GH NOT RELEVANT (DIFFERENT ENERGY) SKIP NEXT
C  SECTIONS
      IF (NEW+1) 19,15,19
C  COMPARE PRESENT WITH OLD INTERPLANAR VECTORS
15    DO 18 IN=1,NLAY2
      NGOL(IN)=0
      DO 16 K=1,3
16    POSA(K)=ABS(DRL(IN,K)-SDRL(IN,K))
      IF (POSA(1)+POSA(2)+POSA(3)-0.001) 17,17,24
C  OLD GH CAN BE LEFT IN ITS PLACE
17    NUGH(IN)=0
      NGEQ(IN)=0
      NGOL(IN)=0
      GO TO 18
24    DO 23 IO=1,NLAY2
      DO 21 K=1,3
21    POSA(K)=ABS(DRL(IN,K)-SDRL(IO,K))
      IF (POSA(1)+POSA(2)+POSA(3)-0.001) 22,22,23
C  OLD GH CAN BE USED, BUT NEEDS COPYING INTO OTHER PLACE. HOWEVER, DO
C  NOT ALLOW COPYING FROM A PLACE THAT COPYING WILL ALREADY HAVE
C  OVERWRITTEN
22    IF (NGOL(IO).NE.0) GO TO 23
      NUGH(IN)=0
      NGEQ(IN)=0
      NGOL(IN)=IO
      GO TO 18
23    CONTINUE
18    CONTINUE
C  STORE PRESENT INTERPLANAR VECTORS FOR LATER COMPARISON
19    DO 20 IN=1,NLAY2
      DO 20 J=1,3
20    SDRL(IN,J)=DRL(IN,J)
      RETURN
      END

!-----------------------------------------------------------------------------
!  SUBROUTINE SUBRAS_BLAS uses layer doubling to generate reflection and
!  transmission matrices for a stack of layers of sufficient thickness
!  to achieve convergence to the semi-infinite crystal limit.
!  The crystal is build from two 'layers' with identical reflection/
!  transmission matrices. The layers may be asymmetric, i.e., the
!  reflection/transmission matrices for beams propagating toward and
!  away from the solid may be different.
!
!  Modified: MRiva 2021-09-10. Get rid of M, get rid of most GOTOs.
!            Use optimized BLAS routine CGEMM for matrix multiplication
!            in TLRTA_BLAS.
!            TODO: The only useful output of this is RA, as it is used
!            only to stack the bulk. Figure out if there is a way to
!            do this exactly with a matrix inversion!
!
!  RA,TA,RB,TB
!    Input reflection and transmission matrices for reflection and
!    transmission by 'input layers' for incidence towards the solid
!    (RA, TA) and for incidence away from it (RB, TB).
!    These matrices will be overwritten and on output represent
!    reflection and transmission by a stack of layers for incidence
!    towards the solid (RA,TA) and for incidence away from it (RB,TB).
!
!   N= total NO. OF BEAMS at current energy.
!   S1,S2,S3,S4,S5,S6,PP,XS,INV_PIVOT= WORKING SPACES.
!   PQ= LIST OF RECIPROCAL LATTICE VECTORS G (BEAMS).
!   AS= INTERLAYER VECTOR (AS(1).GT.0.0).
!  IN COMMON BLOCKS
!   EPS= CONVERGENCE CRITERION FOR LAYER DOUBLING
!   LITER= LIMIT ON NUMBER OF DOUBLINGS.
!   E,VPI  COMPLEX CURRENT ENERGY.
!   AK2,AK3  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!  NOTE  DO NOT IN GENERAL ASSIGN THE SAME STORAGE AREA FOR RA AND RB,
!  AND FOR TA AND TB, IN THE CALLING PROGRAM, EVEN THOUGH RA=RB, TA=TB
!  MAY HOLD FOR THE INDIVIDUAL LAYERS  THESE EQUALITIES WILL IN GENERAL
!  NO LONGER HOLD AFTER DOUBLING.
      SUBROUTINE SUBRAS_BLAS (RA, TA, RB, TB, N,
     &                        S1, S2, S3, S4, S5, S6,
     &                        PP, INV_PIVOT, PQ, AS)
      INTEGER  INV_PIVOT
      COMPLEX  RA, TA, RB, TB, PP, IU, XX
      COMPLEX  S1, S2, S3, S4, S5, S6
      COMPLEX  SQRT, EXP
      DIMENSION  RA(N,N), TA(N,N), RB(N,N), TB(N,N), PP(N,2)
      DIMENSION  S1(N,N), S2(N,N), S3(N,N), S4(N,N), S5(N,N), S6(N,N),
     &           PQ(2,N)
      DIMENSION  AS(3), INV_PIVOT(N)

      COMPLEX TEMP

      COMMON  /MS/  LMAX, EPS, LITER
      COMMON  E, AK2, AK3,VPI

   10 FORMAT(4H X =,E15.4)
   20 FORMAT(24H NO CONV IN SUBREF AFTER,I3,2X,9HITER, X =,E15.4)
   30 FORMAT(I4,16H  ITER IN SUBREF)

      IU = CMPLX(0.0, 1.0)
      AK = 2.0 * E
!     Generate plane-wave propagators between layers
      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)
        XX = CMPLX(AK - BK2 * BK2 - BK3 * BK3, - 2.0 * VPI + 0.000001)
        XX = SQRT(XX) * AS(1)
        X = BK2 * AS(2) + BK3 * AS(3)
        PP(I, 1) = EXP(IU * (XX + X))
        PP(I, 2) = EXP(IU * (XX - X))
      ENDDO

      L = 0
      X = 0.0                                                           280193

!     START ITERATION OVER NUMBER OF DOUBLINGS. Layers may be asymmetric,
!     but should be identical. This means that they may have differences
!     in the transmission matrices from above and from below, but both
!     have the same reflection/transmission from above, and the same
!     from below.

  800 X0 = X
      IF (N.NE.1) THEN
        CALL TLRTA_BLAS (RA, TA, RB, TB, N, S1, S2, S3, S4, S5, S6, PP,
     &                   INV_PIVOT)
      ELSE
!       If only one beam, use direct formula rather than the
!       matrix inversion of TLRTA
!       XX now is the multiple scattering effect, from 1 / (1-x) = sum x**k:
!       propagate beam down [PP(1,1)], reflect there from above [RA(1,1)],
!       propagate up again [PP(1,2)], then reflect from below [RB(1,1)]
        XX = 1.0
     +       /(CMPLX(1.0, 0.0) - PP(1,1) * RA(1,1) * PP(1,2) * RB(1,1))
        TEMP = 1.0 + TA(1,1) * PP(1,1) * TB(1,1) * PP(1,2) * XX
        RA(1,1) = RA(1,1) * TEMP
        RB(1,1) = RB(1,1) * TEMP

        TA(1,1) = TA(1,1) * PP(1,1) * TA(1,1) * XX
        TB(1,1) = TB(1,1) * PP(1,2) * TB(1,1) * XX
      ENDIF

!     PERFORM ONE DOUBLING STEP
      X = 0.0
      DO J = 1, N
        TEMP = PP(J, 1)
        DO I = 1, N
          XX = RA(I, J) * PP(I, 2) * TEMP
          X = X + ABS(REAL(XX)) + ABS(AIMAG(XX))
        ENDDO
      ENDDO
      WRITE (6, 10)  X
      L = L + 1

!     Check number of iterations
      IF (L.GE.LITER) THEN
        WRITE (6, 20)  LITER, X
        write(6, *) "Make sure that interlayer spacing is ",
     +              "large enough for plane-wave space, else ",
     +              "use composite layer instead. Stopping program."
        STOP

!     Check convergence
      ELSEIF (ABS(X-X0)/X .LT. EPS) THEN
        WRITE (6, 30)  L
      ELSE
         GO TO 800
      ENDIF

      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE TAUMAT COMPUTES THE MATRIX TAU (INTRA-SUBPLANE MULTIPLE    111181
C  SCATTERING IN (L,M)-REPRESENTATION) FOR A SINGLE BRAVAIS-LATTICE
C  LAYER FOR EACH OF SEVERAL CHEMICAL ELEMENTS. THE MATRIX ELEMENTS ARE
C  ORDERED THUS  THE (L,M) SEQUENCE IS THE \SYMMETRIZED\ ONE (CF. BOOK).
C  TAU IS SPLIT INTO BLOCK-DIAGONAL PARTS CORRESPONDING TO EVEN L+M
C  AND ODD L+M OF DIMENSIONS LEV AND LOD, RESP. IN THE MATRIX TAU THE
C  PART FOR L+M= ODD IS LEFT-JUSTIFIED AND THE MATRIX REDUCED IN WIDTH
C  FOR STORAGE EFFICIENCY. THE TAU MATRICES FOR DIFFERENT CHEMICAL
C  ELEMENTS ARE ARRANGED UNDER EACH OTHER IN COLUMNAR FASHION  AT THE
C  TOP IS THE TAU BASED ON THE ATOMIC T-MATRIX ELEMENTS TSF(1,L),
C  FOLLOWED BY THE TAU\S BASED ON TSF(2,L), TSF(3,L), ETC.
C  FOR QUANTITIES NOT EXPLAINED BELOW SEE MPERTI OR MTINV.
C   TAU= OUTPUT MATRIX (SEE ABOVE).
C   LMT= NTAU*LMMAX.
C   NTAU= NO. OF CHEMICAL ELEMENTS TO BE USED.
C   X= WORKING SPACE.
C   TSF= ATOMIC T-MATRIX ELEMENTS.
C   FLMS= LATTICE SUMS FROM SUBROUTINE FMAT.
C   NL=NO. OF SUBLATTICES CONSIDERED IN FMAT (ONLY THE FIRST SUBLATTICE
C    SUM IS USED BY TAUMAT).
C   CLM= CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CELMG.
C   NLM= DIMENSION OF CLM (SEE MAIN PROGRAM).
C   LXI= PERMUTATION OF (L,M) SEQUENCE FROM SUBROUTINE LXGENT.
C   NT,PQ,NA,NLL- SEE MPERTI,MTINV.                                     111181
C   FLM- WORKING SPACE.
C   XH - working space carried through to ensure proper dimensions                                                 111181
C  IN COMMON BLOCKS
C   E,VPI= CURRENT COMPLEX ENERGY.
      SUBROUTINE TAUMAT(TAU,LMT,NTAU,X,LEV,LEV2,LOD,TSF,
     +LMMAX,LMAX,FLMS,NL,KLM,LM,CLM,NLM,LXI,NT,PQ,NA,NLL,FLM,
     +NNSUB,LMAX1,XH)
      INTEGER NNSUB,LMAX1
      DIMENSION CLM(NLM),LXI(LMMAX),PQ(2,NT)                            111181
      DIMENSION BR1(2),BR2(2),AR1(2),AR2(2),RAR1(2),RAR2(2)             111181
      COMPLEX AK,CZ,TAU(LMT,LEV),X(LEV,LEV2),TSF(NNSUB,LMAX1),
     1FLMS(NL,KLM),DET,FLM(KLM),CI,XA                                   111181
      COMPLEX XH(LEV)
      COMMON E,AK2,AK3,VPI
      COMMON /SL/BR1,BR2,AR1,AR2,RAR1,RAR2,NL1,NL2                      111181
      CZ=(0.0,0.0)
      CI=(0.0,1.0)                                                      111181
      AK=-0.5/SQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
      DO 100 K=1,KLM                                                    111181
  100 FLM(K)=CZ                                                           .
      BK2=PQ(1,1+NA)
      BK3=PQ(2,1+NA)
      JS=1
      S1=0.
      DO 130 J=1,NL1
      S2=0.
      DO 120 K=1,NL2
      ADR1=S1*BR1(1)+S2*BR2(1)
      ADR2=S1*BR1(2)+S2*BR2(2)
      ABR=ADR1*BK2+ADR2*BK3
      XA=EXP(ABR*CI)
      DO 110 I=1,KLM
110   FLM(I)=FLM(I)+FLMS(JS,I)*XA
      IF (NLL.EQ.1) GO TO 140
      JS=JS+1
  120 S2=S2+1.
  130 S1=S1+1.                                                            .
  140 CONTINUE                                                          111181
      DO 13 IT=1,NTAU
      DO 13 IL=1,2
      DO 1 I=1,LEV
      DO 1 J=1,LEV2
1     X(I,J)=CZ
      LL=LOD
      IF (IL.EQ.2) LL=LEV
C  GENERATE MATRIX 1-X FOR L+M= ODD (IL=1), LATER FOR L+M= EVEN (IL=2)
      CALL XMT(IL,FLM,1,X,LEV,LL,TSF,NTAU,IT,LM,LXI,LMMAX,              111181
     1KLM,CLM,NLM,NST,NNSUB,LMAX1)
C
C  PREPARE QUANTITIES INTO WHICH INVERSE OF 1-X WILL BE MULTIPLIED
      IF (IL-2) 6,2,2
2     IS=0
      LD1=0
      L=0
3     LD=LD1+1
      LD1=LD+L
      DO 4 I=LD,LD1
4     X(I,I+LEV)=AK*TSF(IT,L+1)
      L=L+2
      IF (L-LMAX) 3,3,5
5     IS=IS+1
      L=1
      IF (IS-1) 3,3,10
C
6     IS=0
      LD1=0
      L=1
7     LD=LD1+1
      LD1=LD+L-1
      DO 8 I=LD,LD1
8     X(I,I+LOD)=AK*TSF(IT,L+1)
      L=L+2
      IF (L-LMAX) 7,7,9
9     IS=IS+1
      L=2
      IF (IS-1) 7,7,10
10    LL2=LL+LL
C
C  PERFORM INVERSION AND MULTIPLICATION
      CALL CXMTXT(X,LEV,LL,LL,LL2,MARK,DET,-1,XH)    ! TODO Michele: Judge if this can be replaced with CGETRF/CGETRS (likely better optimized?)
C
      LD=(IT-1)*LMMAX
      IF (IL.EQ.1) LD=LD+LEV
      DO 11 I=1,LL
      DO 11 J=1,LL
C  PUT RESULT IN TAU
11    TAU(LD+I,J)=X(I,J+LL)
13    CONTINUE
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE TAUY_SIMPLE PRODUCES THE QUANTITIES TAU*YLM(G) FOR ALL CHEMICAL
!  ELEMENTS AND, IF BEEBY-TYPE INVERSION IS CALLED FOR IN THE RSP
!  CALCULATION, TH*YLM(G) FOR ALL RELEVANT LAYERS. THESE ARE USED
!  EITHER FOR THE BEEBY-TYPE INVERSION IN MTINV OR AS INITIAL VALUES FOR
!  THE RSP ITERATION IN MPERTI.
!  FOR QUANTITIES NOT DESCRIBED BELOW SEE CALLING ROUTINES (MPERTI OR
!  MTINV).
!
!  Edits: MRiva 2021-09-06. Remove obsolete code, rewrite to remove
!         all GOTOs
!
!   TAUG,TAUGM= OUTPUT RESULTING VECTORS FOR K(G+) AND K(G-),RESP.
!   LTAUG= NTAU*LMMAX
!   LT=PERMUTATION OF (LM) ORDER, FROM LXGENT.
!   JGP= CURRENT INCIDENT BEAM.
      SUBROUTINE TAUY_SIMPLE(TAUG, TAUGM, LTAUG, TAU, LMT, LEV, CYLM,
     &                       NT, LMMAX, LT, NTAU, LOD, LEE, LOE, JGP)
      COMPLEX CZ, ST, SU, CF, CF1
      COMPLEX TAU(LMT, LEV), CYLM(NT, LMMAX)
      COMPLEX TAUG(LTAUG), TAUGM(LTAUG)
      DIMENSION LT(LMMAX)
      CZ = CMPLX(0.0, 0.0)

!
!  PERFORM MATRIX PRODUCT TAU*YLM(G+-) FOR EACH CHEMICAL ELEMENT
      DO I = 1, NTAU
        IS = (I - 1) * LMMAX
        DO JLM = 1, LEV
          ST = CZ
          JTAU = IS + JLM
          DO ILM = 1, LEV
            KLP = LT(ILM)
            IF (ILM.GT.LEE) THEN
              ST = ST - TAU(JTAU, ILM) * CYLM(JGP, KLP)
            ELSE
              ST = ST + TAU(JTAU, ILM) * CYLM(JGP, KLP)
            ENDIF
          ENDDO
          TAUG(JTAU) = ST
          TAUGM(JTAU) = ST
        ENDDO

        IS = IS + LEV
        DO JLM = 1, LOD
          ST = CZ
          SU = CZ
          JTAU = IS + JLM
          DO ILM = 1, LOD
            KLP = LT(ILM + LEV)
            CF = TAU(JTAU, ILM) * CYLM(JGP, KLP)
            IF (ILM.GT.LOE) THEN
              ST = ST - CF
              SU = SU + CF
            ELSE
              ST = ST + CF
              SU = SU - CF
            ENDIF
          ENDDO
          TAUG(JTAU) = ST
          TAUGM(JTAU) = SU
        ENDDO
      ENDDO
      END

!-----------------------------------------------------------------------
!  SUBROUTINE THMAT PRODUCES THE MATRIX THAT IS INVERTED IN THE BEEBY-
!  TYPE METHOD IN SUBROUTINE RTINV. FOR QUANTITIES NOT EXPLAINED
!  BELOW SEE RTINV.
!
!  Edited: 2021-10-07 MRiva add ENDDOs and indents. Remove LMNI2 (== LMNI)
!
!   TH= OUTPUT MATRIX, STILL TO BE INVERTED.
!   LMNI= NLAY*LMMAX.
      SUBROUTINE THMAT(TH, LMNI, GH, LMG, LMMAX, MGH, NLAY, TAU,
     1                 LMT, LEV, LPS)
      COMPLEX TH(LMNI, LMNI), GH(LMG, LMMAX), TAU(LMT, LEV)
      COMPLEX RU, CZ, ST, SU
      DIMENSION LPS(NLAY), MGH(NLAY, NLAY)

      LOD = LMMAX - LEV
      RU = (1.0, 0.0)
      CZ = (0.0, 0.0)

!     Start from identity matrix
      DO J = 1, LMNI
        DO I = 1, LMNI
          TH(I, J) = CZ
        ENDDO
        TH(J, J) = RU
      ENDDO

      IC = 0
      NLAY1 = NLAY - 1
      DO IN = 1, NLAY1
        IC = IC + 1
        ICL = (IC - 1) * LMMAX
        ID = IC
        IN1 = IN + 1
        DO INP = IN1, NLAY
          ID = ID + 1
          IDL = (ID - 1) * LMMAX
          M = (MGH(IN, INP) - 1) * LMMAX
          M1 = (MGH(INP, IN) - 1) * LMMAX
          LP = (LPS(IN) - 1) * LMMAX
          LP1 = (LPS(INP) - 1) * LMMAX

!         INSERT -TAU*GH FOR EVEN L+M (IN FIRST INDEX)
          DO I = 1, LEV
            DO J = 1, LMMAX
              ST = CZ
              SU = CZ
              DO K = 1,LEV
                ST = ST + TAU(LP + I, K) * GH(M + K, J)
                SU = SU + TAU(LP1 + I, K) * GH(M1 + K, J)
              ENDDO
              TH(ICL + I, IDL + J) = -ST
              TH(IDL + I, ICL + J) = -SU
            ENDDO
          ENDDO    ! I loop, over even L+M

          LP = LP + LEV
          LP1 = LP1 + LEV
          M = M + LEV
          M1 = M1 + LEV
          ICV = ICL + LEV
          IDV = IDL + LEV

!         INSERT -TAU*GH FOR ODD L+M (IN FIRST INDEX)
          DO I = 1, LOD
            DO J = 1, LMMAX
              ST = CZ
              SU = CZ
              DO K = 1, LOD
                ST = ST + TAU(LP + I, K) * GH(M + K, J)
                SU = SU + TAU(LP1 + I, K) * GH(M1 + K, J)
              ENDDO
              TH(ICV + I, IDL + J) = -ST
              TH(IDV + I, ICL + J) = -SU
            ENDDO
          ENDDO  ! I loop, over odd L+M
        ENDDO  ! INP loop
      ENDDO  ! IN loop
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE TLRTA_BLAS performs each layer doubling step for
!  subroutine SUBRAS_BLAS by matrix inversion. The two layers are
!  assumed to be identical, but can be asymmetric in the beam
!  propagation direction.
!
!  This modified version is optimized to use the LAPACK/BLAS routine
!  CGEMM for matrix multiplication, rather than doing explicit looping.
!
!   R1,T1,R2,T2
!    Input reflection and transmission matrices of the two layers.
!    Numbers 1 and 2 refer to incidence towards the solid and
!    away from it, respectively. They are overwritten and will
!    contain the corresponding matrices for the combination.
!
!   N= NO. OF BEAMS USED.
!   S1, S2, R1_P, R2_P, T1_P, T2_P, INV_PIVOT= WORKING SPACES.
!   PP= PLANE-WAVE PROPAGATORS, FROM SUBROUTINE SUBRAS_BLAS.

      SUBROUTINE  TLRTA_BLAS (R1, T1, R2, T2, N,
     &                        S1, S2, R1_P, R2_P, T1_P, T2_P,  ! storage only
     &                        PP, INV_PIVOT)
      COMPLEX  R1, T1, R2, T2, XX, YY, C_ZERO, PP
      COMPLEX  S1, S2, R1_P, R2_P, T1_P, T2_P
      COMPLEX  P_UP, P_DOWN, C_M_ONE, C_ONE
      DIMENSION  R1(N, N), T1(N, N), R2(N, N), T2(N, N)
      DIMENSION  R1_P(N, N), R2_P(N, N), T1_P(N, N), T2_P(N, N)
      DIMENSION  S1(N, N), S2(N, N), PP(N, 2)
      DIMENSION  INV_PIVOT(N)

      C_ZERO = CMPLX(0.0, 0.0)
      C_ONE = CMPLX(1.0, 0.0)
      C_M_ONE = CMPLX(-1.0, 0.0)

!     Store in R1_P, R2_P, T1_P, and T2_P modified versions of the
!     reflection/transmission matrices that include also propagation
!     between the layers for the incoming beam j.
!                                                        i ^
!                                          R2(i,j)         | T2(i,j)
!     --<<TOP>>--------------------------------------------------------
!                                          j ^  |        j ^
!                                            |  v i        |
!                  |  ^ i      |
!                j v  |      j v
!     --<<BOT>>--------------------------------------------------------
!                R1(i,j)       |  T1(i,j)
!                            i v
!                                                        i ^
!                                          R2_P(i,j)       |  T2_P(i,j)
!     --<<TOP>>--------------------------------------------------------
!                  |           |           j ^  |        j ^
!                  |           |             |  v i        |
!                  |  ^ i      |             |             |
!                j v  |      j v             |             |
!     --<<BOT>>--------------------------------------------------------
!                R1_P(i,j)     |  T1_P(i,j)
!                            i v
      DO J = 1, N
        P_UP = PP(J, 2)
        P_DOWN = PP(J, 1)
        DO I = 1, N
          R1_P(I, J) = R1(I, J) * P_DOWN
          R2_P(I, J) = R2(I, J) * P_UP
          T1_P(I, J) = T1(I, J) * P_DOWN
          T2_P(I, J) = T2(I, J) * P_UP
        ENDDO
      ENDDO

!     Prepare multiple-scattering matrix S1 for initial and final beams
!     propagating 'down' (i.e., toward the solid), and right below the
!     top layer. Use (Id-X)**-1 = sum_k X**k, with single-pass scattering
!     X(i,j) = sum_k R2_P(i,k) R1_P(k,j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S1 = Id - X
     &           R2_P, N, R1_P, N,
     &           C_ZERO, S1, N)
      DO I = 1, N
        S1(I, I) = S1(I, I) + C_ONE
      ENDDO

!     LU decomposition of matrix S1 (multiple scattering, beams moving down)
      CALL CGETRF(N, N, S1, N, INV_PIVOT, INFO)

!     Then multiply with transmission from above, such that T1(i, j) is:
!     beam j coming from above, transmitted through top and multiple-scattered
!     into i, which is right below the top layer and propagating down
      CALL CGETRS('N', N, N, S1, N, INV_PIVOT, T1, N, INFO)

!     Now calculate the total transmission from above into S1(i, j)
!     by taking S1(i, j) = sum_k [T1(k, j)     ! transmission through top of j into k, with multiple scattering; k right below top.
!                                 T1_P(i, k)]  ! transmission through bottom of k into i, including propagation down.
!     S1 will later be copied into T1.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T1_P, N, T1, N,
     &           C_ZERO, S1, N)

!     Go for the total reflection from above: (1) propagate-reflect
!     at the bottom layer [==multiply with R1_P(i,k)] the transmitted
!     and multiple-scattered beam k [==T1(k,j)]. Store this in S2(i,j).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, R1_P, N, T1, N,
     &           C_ZERO, S2, N)

!     Then (2) propagate-transmit at top [==multiply with T2_P(i,k)]
!     this beam k [==S2(k,j)]. Also, add to this the specular reflection
!     at the top layer R1(i,j), and store the whole thing in R1.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T2_P, N, S2, N,
     &           C_ONE, R1, N)

!     The original T1 is not needed anymore. Replace it with the total
!     transmission (temporarily in S1)
      DO J = 1, N
        DO I = 1, N
          T1(I, J) = S1(I, J)
        ENDDO
      ENDDO

!     Now repeat the whole thing for incidence from below.

!     Prepare multiple-scattering matrix S2 for initial and final beams
!     propagating 'up' (i.e., away from the solid), and right above the
!     bottom layer. Use (Id-X)**-1 = sum_k X**k, with single-pass scattering
!     X(i,j) = sum_k R1_P(i,k) R2_P(k,j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S2 = Id - X
     &           R1_P, N, R2_P, N,
     &           C_ZERO, S2, N)
      DO I = 1, N
        S2(I, I) = S2(I, I) + C_ONE
      ENDDO

!     LU decomposition of matrix S2 (multiple scattering, beams moving up)
      CALL CGETRF(N, N, S2, N, INV_PIVOT, INFO)

!     Then multiply with transmission from below, such that T2(i, j) is:
!     beam j coming from below, transmitted through bottom and multiple-
!     scattered into i, which is right above the bottom layer and
!     propagating up
      CALL CGETRS('N', N, N, S2, N, INV_PIVOT, T2, N, INFO)

!     Now calculate the total transmission from below into S2(i, j)
!     by taking S2(i, j) = sum_k [T2(k, j)     ! transmission through bottom of j into k, with multiple scattering; k right above bottom.
!                                 T2_P(i, k)]  ! transmission through top of k into i, including propagation up.
!     S2 will later be copied into T2.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T2_P, N, T2, N,
     &           C_ZERO, S2, N)

!     Go for the total reflection from below: (1) propagate-reflect
!     at the top layer [==multiply with R2_P(i,k)] the transmitted
!     and multiple-scattered beam k [==T2(k,j)]. Store this in S1(i,j).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, R2_P, N, T2, N,
     &           C_ZERO, S1, N)

!     Then (2) propagate-transmit at bottom [==multiply with T1_P(i,k)]
!     this beam k [==S1(k,j)]. Also, add to this the specular reflection
!     at the top layer R2(i,j), and store the whole thing in R2.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T1_P, N, S1, N,
     &           C_ONE, R2, N)

!     Finally, store the total transmission in T2.
      DO J = 1, N
        DO I = 1, N
          T2(I, J) = S2(I, J)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE TSCATF INTERPOLATES TABULATED PHASE SHIFTS AND PRODUCES
C  THE ATOMIC T-MATRIX ELEMENTS (OUTPUT IN AF). THESE ARE ALSO
C  CORRECTED FOR THERMAL VIBRATIONS (OUTPUT IN CAF). AF AND CAF
C  ARE MEANT TO BE STORED IN ARRAY TMAT FOR LATER USE IN RSMF, RTINV.
C
C   IEL= CHEMICAL ELEMENT TO BE TREATED NOW, IDENTIFIED BY THE INPUT
C    SEQUENCE ORDER OF THE PHASE SHIFTS (IEL=1,2 OR 3).
C   L1= LMAX+1.
C   ES= LIST OF ENERGIES AT WHICH PHASE SHIFTS ARE TABULATED.
C   PHSS= TABULATED PHASE SHIFTS.
C   NPSI= NO. OF ENERGIES AT WHICH PHASE SHIFTS ARE GIVEN.
C   EB-V= CURRENT ENERGY (V CAN BE USED TO DESCRIBE LOCAL VARIATIONS
C    OF THE MUFFIN-TIN CONSTANT).
C   PPP= CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CPPP.
C   NN1= NN2+NN3-1.
C   NN2= NO. OF OUTPUT TEMPERATURE-CORRECTED PHASE SHIFTS DESIRED.
C   NN3= NO. OF INPUT PHASE SHIFTS.
C   DR0= FOURTH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE.
C   DRPER= RMS VIBRATION AMPLITUDE PERPENDICULAR TO SURFACE.
C   DRPAR= RMS VIBRATION AMPLITUDE PARALLEL TO SURFACE.
C   T0= TEMPERATURE AT WHICH DRPER AND DRPAR HAVE BEEN COMPUTED.
C   T= CURRENT TEMPERATURE.
C   TSF0,TSF,AF,CAF  SEE ABOVE.
      SUBROUTINE TSCATF(IEL,L1,ES,PHSS,NPSI,EB,V,PPP,NN1,NN2,NN3,
     +DR0,DRPER,DRPAR,T0,T,AF,CAF,NEL,LMAX1,PHS,DEL,CTAB,SUM,BJ)

      DIMENSION PHSS(NPSI,NEL,LMAX1),PHS(LMAX1),ES(NPSI)
      DIMENSION PPP(NN1,NN2,NN3)
      COMPLEX CTAB,SUM
      DIMENSION CTAB(NN3),SUM(NN2)
      COMPLEX*16 BJ(NN1)
      COMPLEX CI,DEL(L1),CA,AF(L1),CAF(L1)

  700 FORMAT(42H TOO LOW ENERGY FOR AVAILABLE PHASE SHIFTS)
      CI=(0.0,1.0)
      E=EB-V
      IF (E.LT.ES(1)) GO TO 850
C  FIND SET OF PHASE SHIFTS APPROPRIATE TO DESIRED CHEMICAL ELEMENT
C  AND INTERPOLATE LINEARLY TO CURRENT ENERGY (OR EXTRAPOLATE TO ENERGIES
C  ABOVE THE RANGE GIVEN FOR THE PHASE SHIFTS)

  710 I = 1
  720 IF ((E-ES(I))*(E-ES(I+1)))  750, 750, 730
  730 I = I + 1
      IF (I-NPSI)  720, 740, 740
  740 I = I - 1
  750 FAC = (E - ES(I))/(ES(I + 1) - ES(I))
      DO 760 L = 1, L1
  760 PHS(L) = PHSS(I,IEL,L) + FAC * (PHSS(I + 1,IEL,L)
     +       - PHSS(I,IEL,L))
C  COMPUTE TEMPERATURE-INDEPENDENT T-MATRIX ELEMENTS
      DO 790 L = 1, L1
      A = PHS(L)
      AF(L) = A * CI
      AF(L) = EXP(AF(L))
      A = SIN(A)
      AF(L) = A * AF(L)
  790 CONTINUE
C  AVERAGE ANY ANISOTROPY OF RMS VIBRATION AMPLITUDES
      DR=SQRT((DRPER*DRPER+2.0*DRPAR*DRPAR)/3.0)
C  COMPUTE TEMPERATURE-DEPENDENT PHASE SHIFTS (DEL)
      CALL  PSTEMP (PPP,NN1,NN2,NN3,DR0,DR,T0,T,E,PHS,DEL,CTAB,SUM,BJ)
C  PRODUCE TEMPERATURE-DEPENDENT T-MATRIX ELEMENTS
      DO 840 L = 1, L1
      CA = DEL(L)
      CAF(L) = CA * CI
      CAF(L) = EXP(CAF(L))
      CA = SIN(CA)
      CAF(L) = CA * CAF(L)
  840 CONTINUE
      RETURN
  850 WRITE(6,700)
      STOP
      END

C-----------------------------------------------------------------------
C  SUBROUTINE XM PRODUCES THE INTRA-LAYER MULTIPLE-SCATTERING MATRIX X
C  FOR A SINGLE BRAVAIS-LATTICE LAYER. THE (L,M) SEQUENCE IS THE
C  'SYMMETRIZED' ONE AND X IS SPLIT INTO TWO PARTS CORRESPONDING TO
C  L+M= EVEN (XEV) AND L+M= ODD (XOD), AS A RESULT OF BLOCK-DIAGONALI-
C  ZATION.
C   FLM= LATTICE SUM FROM SUBROUTINES FMAT AND MSMF.
C   XEV,XOD= OUTPUT MATRIX X IN TWO PARTS.
C   LEV= (LMAX+1)*(LMAX+2)/2.
C   LOD= LMAX*(LMAX+1)/2.
C   AF  NOT USED.
C   CAF= TEMPERATURE-DEPENDENT ATOMIC SCATTERING T-MATRIX ELEMENTS.
C   LM= LMAX+1.
C   LX,LXI= PERMUTATIONS OF (L,M) SEQUENCE, FROM SUBROUTINE LXGENT.
C   LMMAX= (LMAX+1)**2.
C   KLM= (2*LMAX+1)*(2*LMAX+2)/2.
C   CLM= CLEBSCH-GORDON COEFFICIENTS, FROM SUBROUTINE CELMG.
C   NLM= DIMENSION OF CLM (SEE MAIN PROGRAM).
      SUBROUTINE  XM (FLM,XEV,XOD,LEV,LOD,AF,CAF,LM,LX,LXI,LMMAX,
     1 KLM,CLM,NLM)
      INTEGER  LX, LXI
      COMPLEX  XEV, XOD, FLM, AF, CZERO, ACC, CAF
      DIMENSION  CLM(NLM), XEV(LEV,LEV), XOD(LOD,LOD), FLM(KLM), AF(LM)
      DIMENSION  CAF(LM), LX(LMMAX), LXI(LMMAX)
      LMAX = LM-1
      CZERO = CMPLX(0.0,0.0)
      DO 440 I = 1, LEV
      DO 440 J = 1, LEV
  440 XEV(I,J) = CZERO
      DO 450 I = 1, LOD
      DO 450 J = 1, LOD
  450 XOD(I,J) = CZERO
      L2MAX = LMAX + LMAX
      NODD = LOD
      JSET = 1
      MM = LEV
      N = 1
C  FIRST XOD IS CREATED
  460 J = 1
      L = JSET
  470 M =  - L + JSET
      JL = L + 1
  480 K = 1
      LPP = JSET
  490 MPP =  - LPP + JSET
  500 MPA = IABS(MPP - M)
      LPA = IABS(LPP - L)
      IF (LPA-MPA)  520, 520, 510
  510 MPA = LPA
  520 MP1 = MPP - M + L2MAX + 1
      LP1 = L + LPP + 1
      ACC = CZERO
  530 JLM = (LP1 * LP1 + MP1 - L2MAX)/2
      ACC = ACC + CLM(N) * FLM(JLM)
      N = N + 1
      LP1 = LP1 - 2
      IF (LP1-1-MPA)  540, 530, 530
  540 JX = LXI(J + MM)
      KX = LXI(K + MM)
      XEV(JX,KX) = ACC * CAF(JL)
      K = K + 1
      MPP = MPP + 2
      IF (LPP-MPP)  550, 500, 500
  550 LPP = LPP + 1
      IF (LMAX-LPP)  560, 490, 490
  560 J = J + 1
      M = M + 2
      IF (L-M)  570, 480, 480
  570 L = L + 1
      IF (LMAX-L)  580, 470, 470
  580 IF (JSET)  610, 610, 590
  590 DO 600 J = 1, NODD
      DO 600 K = 1, NODD
  600 XOD(J,K) = XEV(J,K)
      JSET = 0
      MM = 0
C  NOW RETURN TO CREATE XEV
      GO TO 460
  610 RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE XMT PRODUCES THE INTRA-LAYER MULTIPLE-SCATTERING MATRIX
C  1-X FOR SUBROUTINE TAUMAT, OUTPUT IN THE TONG CONVENTION USING
C  INPUT IN THE PENDRY CONVENTION. XMT MUST BE CALLED TWICE, FIRST TO
C  PRODUCE THE VALUES FOR L+M= ODD, THEN FOR L+M= EVEN (IN THAT ORDER).
C  FOR QUANTITIES NOT EXPLAINED BELOW SEE SUBROUTINE TAUMAT.
C   IL=1 FOR L+M= ODD.
C   IL=2 FOR L+M= EVEN.
C   X= OUTPUT MATRIX 1-X, BLOCK-DIAGONALIZED, 2ND BLOCK LEFT-ADJUSTED.
C   LL= INPUT EITHER LOD OR LEV.
C   TSF= ATOMIC T-MATRIX ELEMENTS.
C   NTAU= NO. OF CHEMICAL ELEMENTS CONSIDERED.
C   IT= INDEX OF CURRENT CHEMICAL ELEMENT (.LE.NTAU).
C   N= RUNNING INDEX OF CLM  MAY NOT BE RESET BETWEEN THE TWO CALLS TO
C    XMT IN TAUMAT.
      SUBROUTINE XMT(IL,FLMS,NL,X,LEV,LL,TSF,NTAU,IT,LM,LXI,LMMAX,
     1 KLM,CLM,NLM,N,NNSUB,LMAX1)
      COMPLEX TSF(NNSUB,LMAX1)
      COMPLEX X,FLMS,CZERO,ACC,CI,ST,SU,RU
      DIMENSION CLM(NLM),X(LEV,LL),FLMS(NL,KLM),LXI(LMMAX)
      LMAX = LM-1
      RU=(1.0,0.0)
      CZERO = CMPLX(0.0,0.0)
      CI=CMPLX(0.0,1.0)
      L2MAX = LMAX + LMAX
C  IF IL=1, CONSIDER L+M= ODD ONLY
C  IF IL=2, CONSIDER L+M= EVEN ONLY
      IF (IL-2) 455,450,450
450   JSET=0
      MM=0
      GO TO 457
455   JSET = 1
      MM = LEV
      N=1
457   CONTINUE
  460 J = 1
      L = JSET
  470 M =  - L + JSET
      JL = L + 1
      ST=(-CI)**MOD(L,4)
  480 K = 1
      LPP = JSET
  490 MPP =  - LPP + JSET
      JLPP=LPP+1
      SU=((-CI)**MOD(LPP,4))/ST
  500 MPA = IABS(MPP - M)
      LPA = IABS(LPP - L)
      IF (LPA-MPA)  520, 520, 510
  510 MPA = LPA
  520 MP1 = MPP - M + L2MAX + 1
      LP1 = L + LPP + 1
      ACC = CZERO
  530 JLM = (LP1 * LP1 + MP1 - L2MAX)/2
      ACC=ACC+CLM(N)*FLMS(1,JLM)
      N = N + 1
      LP1 = LP1 - 2
      IF (LP1-1-MPA)  540, 530, 530
  540 JX = LXI(J + MM)
      KX = LXI(K + MM)
      X(KX,JX) = -ACC * TSF(IT,JLPP)*SU
      IF (J.EQ.K) X(KX,JX)=X(KX,JX)+RU
      K = K + 1
      MPP = MPP + 2
      IF (LPP-MPP)  550, 500, 500
  550 LPP = LPP + 1
      IF (LMAX-LPP)  560, 490, 490
  560 J = J + 1
      M = M + 2
      IF (L-M)  570, 480, 480
  570 L = L + 1
      IF (LMAX-L)  610, 470, 470
  610 RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE ZGE PERFORMS GAUSSIAN ELIMINATION AS THE FIRST STEP IN THE
C  SOLUTION OF A SYSTEM OF LINEAR EQUATIONS. THIS IS USED TO MULTIPLY
C  THE INVERSE OF A MATRIX INTO A VECTOR, THE MULTIPLICATION BEING DONE
C  LATER BY SUBROUTINE ZSU.
C   A= INPUT (COMPLEX) MATRIX, WILL BE OVERWRITTEN BY A TRIANGULARIZED
C    MATRIX, TO BE TRANSMITTED TO SUBROUTINE ZSU.
C   INT= STORAGE FOR PERMUTATION OF MATRIX COLUMNS (AND ROWS), TO BE
C    TRANSMITTED TO SUBROUTINE ZSU.
C   NR= FIRST DIMENSION OF A (.GE.NC).
C   NC= ORDER OF A.
C   EMACH= MACHINE ACCURACY.
      SUBROUTINE  ZGE (A, INT, NR, NC, EMACH)
      COMPLEX  A, YR, DUM
      DIMENSION  A(NR,NC), INT(NC)
      N = NC
      DO 680 II = 2, N
      I = II - 1
      YR = A(I,I)
      IN = I
      DO 600 J = II, N
      IF (CABS(YR)-CABS(A(J,I)))  590, 600, 600
  590 YR = A(J,I)
      IN = J
  600 CONTINUE
      INT(I) = IN
      IF (IN-I)  610, 630, 610
  610 DO 620 J = I, N
      DUM = A(I,J)
      A(I,J) = A(IN,J)
  620 A(IN,J) = DUM
  630 IF (CABS(YR)-EMACH)  680, 680, 640
  640 DO 670 J = II, N
      IF (CABS(A(J,I))-EMACH)  670, 670, 650
  650 A(J,I) = A(J,I)/YR
CDIR$ IVDEP
      DO 660 K = II, N
  660 A(J,K) = A(J,K) - A(I,K) * A(J,I)
  670 CONTINUE
  680 CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE ZSU TERMINATES THE SOLUTION OF A SYSTEM OF LINEAR
C  EQUATIONS, INITIATED BY SUBROUTINE ZGE, BY BACK-SUBSTITUTING THE
C  CONSTANT VECTOR.
C   A= INPUT MATRIX, PREPARED BY SUBROUTINE ZGE.
C   INT= INPUT PERMUTATION FROM SUBROUTINE ZGE.
C   X= INPUT CONSTANT VECTOR AND OUTPUT RESULTING VECTOR.
C   NR= FIRST DIMENSION OF A (.GE.NC).
C   NC= ORDER OF A.
C   EMACH= MACHINE ACCURACY.
      SUBROUTINE  ZSU (A, INT, X, NR, NC, EMACH)
      COMPLEX  A, X, DUM
      DIMENSION  A(NR,NC), X(NC), INT(NC)
      N = NC
      DO 730 II = 2, N
      I = II - 1
      IF (INT(I)-I)  690, 700, 690
  690 IN = INT(I)
      DUM = X(IN)
      X(IN) = X(I)
      X(I) = DUM
CDIR$ IVDEP
  700 DO 720 J = II, N
      IF (CABS(A(J,I))-EMACH)  720, 720, 710
  710 X(J) = X(J) - A(J,I) * X(I)
  720 CONTINUE
  730 CONTINUE
      DO 780 II = 1, N
      I = N - II + 1
      IJ = I + 1
      IF (I-N)  740, 760, 740
CDIR$ IVDEP
  740 DO 750 J = IJ, N
  750 X(I) = X(I) - A(I,J) * X(J)
  760 IF (CABS(A(I,I))-EMACH*1.0E-5)  770, 780, 780
  770 A(I,I) = EMACH * 1.0E - 5 * (1.0,1.0)
  780 X(I) = X(I)/A(I,I)
      RETURN
      END

C---------------------------------------------------------------------
      SUBROUTINE ZTU(A,INT,X,ND,N,EMACH)
C     ==============
      COMPLEX A,X,Z1,Z2
      DIMENSION A(ND,ND),INT(ND),X(ND)
C
COMMENT ZTU IS A MODIFIED BACK SUBSTITUTION SUBROUTINE
C       USING THE OUTPUT OF ZGE TO CALCULATE X TIMES
C       A INVERSE, RETURNED IN X
      X(1)=X(1)/A(1,1)
      DO 32 I=2,N
      I1=I-1
      Z1=X(I)
      DO 31 J=1,I1
   31 Z1=Z1-X(J)*A(J,I)
   32 X(I)=Z1/A(I,I)
      DO 34 NI=2,N
      I=N-NI+1
      I1=I+1
      Z1=X(I)
      DO 33 J=I1,N
   33 Z1=Z1-X(J)*A(J,I)
      IN=INT(I)
      IF(I.NE.IN) X(I)=X(IN)
   34 X(IN)=Z1
      RETURN
      END
