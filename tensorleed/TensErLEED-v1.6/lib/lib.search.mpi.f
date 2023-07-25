C  Tensor LEED subroutines for optimization algorithm 
C  v1.2 (search version v106), VB 13.04.00
C  for use with search.f v1.2 (v106)
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
C  Allocation error in Subroutine RINTAV removed via introduction 
C  of additional fields ATAV and PQ1.        Oct. 2018 L. Hammer
C
C**********************************************************************************
C
C  Please read the comment in search.f, v1.2 (v106) .
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                          C
C                              SUBROUTINES                                 C
C                                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  Subroutine r-factor computes either R2 or R-Pe and optimizes for
C  inner potential V0R using the total, integer or half-order
C  r-factor depending on flag WHICHG

      SUBROUTINE RFAKTOR(ATSMK,NSS,NBTD,NBMD,NETI,PQ,KAV,SYM,ESMK,IPR,
     +     AE,NDATT,NET,AT,ET,EINCR,XPL,YPL,NDATA,ATP,YT,VI,NBED,TSTY2,
     +     V02,V01,VINCR,ARM,ARPEM,ERANGM,RAZZM,RANNM,RAVPM,XRPEM,NST1,
     +     NST2,MITTEL,IBP,R2,RPE,EE,NEE,EET,YE,WB,BRGES,BRINS,BRHAS,
     +     BV0,NPS,V0RR,IPOP,TST,TSE,TSE2,TSEY2,BENAME,NBE,WHICHG,
     +     WHICHR,BARAV,OVLG)

CVB  Include global parameters for dimension statements etc.

      INTEGER KAV,SYM,IPR,NET,NST1,NST2,MITTEL,IBP,NEE,WHICHG,WHICHR

      DIMENSION KAV(NBTD),SYM(NBTD),NET(NBTD)
      REAL WB
      DIMENSION MITTEL(NBED),IBP(NBED),NEE(NBED),WB(NBED)
      DIMENSION NST1(NBMD),NST2(NBMD)

      REAL ATSMK,PQ,ESMK,OVLG
      DIMENSION ATSMK(NBTD,NDATT),ATAV(NBTD,NDATT)
      DIMENSION PQ(2,NBTD),PQAV(2,NBTD),ESMK(NDATT)
      REAL AT,ET,XPL,YPL,YT,ATP,TSTY2
      DIMENSION AT(NBTD,NDATA),ET(NBTD,NDATA),XPL(NDATA),YPL(NDATA)
      DIMENSION YT(NBTD,NDATA),ATP(NBTD,NDATA),TSTY2(NBTD)

CVB  for R2:

      REAL AE
      DIMENSION AE(NBED,NDATA)

      REAL TST
      DIMENSION TST(NBTD)

CVB

      REAL EINCR,VI,RPE,RAV,EE,EET,YE
      DIMENSION RPE(NBED),EE(NBED,NDATA),EET(NBED)
      DIMENSION YE(NBED,NDATA)
      REAL ARM(2),ARPEM(2),ERANGM(2),RAZZM(2),RANNM(2),RAVPM(2)
      REAL XRPEM(2)                                                    

CVB  for R2:

      REAL R2(NBED)
      REAL AR2M(2)

CVB

      REAL BV0,BRGES,BRINS,BRHAS,TSEY2
      DIMENSION BV0(NPS),BRGES(NPS),BRINS(NPS),BRHAS(NPS),TSEY2(NBED)
      DIMENSION BENAME(5,NBED)

CVB  for R2:

      REAL TSE,TSE2
      DIMENSION TSE(NBED),TSE2(NBED)

CVB

      REAL BARAV
      INTEGER IS


      OVLG=0.

C  PERFORM DOMAIN-AVERAGING
      CALL RINTAV(ATSMK,NSS,NBTD,NETI,PQ,PQAV,KAV,SYM,NBT,ESMK,IPR,ATAV)
C  CHECK FOR TOO HIGH THEOR. INTENS.             
      DO 28 IB=1,NBT                             
      CALL MAXINT(ATAV,NSS,NBTD,IB,NETI,AM,NDATT)
     
       IF (AM.LE.1.) GO TO 28                     
       WRITE(6,29)IB,AM                           
29    FORMAT(1H ,///,29H MAX. INTENS. IN THEOR. BEAM ,1I3,12H IS SUSPECT
     1(,1E13.5,29H)- ****** STOP PROGRAM ******)

      write(8,29)IB,AM

      STOP                                      
28    CONTINUE                                  

C Initialize some values

      BARAV=1.E2
      IS=1

C  INTERPOLATE THEOR. DATA ONTO WORKING GRID

      DO 40 IB=1,NBT
40    NET(IB)=NETI

      CALL STRIP2(ATAV,NSS,NBTD,NET,IS,AT,NBT,ESMK,ET,NDATT)

      CALL INTPOL(AT,1,NBTD,NET,1,NBT,ET,EINCR,IPR,XPL,YPL)

CVB   

      IF (WHICHR.eq.1) THEN

CVB

C  PRODUCE 1ST AND 2ND DERIVATIVE OF THEORY

      CALL DER(AT,NET,1,NBTD,1,NBT,ATP,EINCR)

C  PRODUCE PENDRY Y FUNCTION FOR THEORY

      CALL YPEND(AT,ATP,1,NBTD,1,NBT,NET,ET,YT,VI,IPR)

CVB

      END IF

CVB

C  PRODUCE SOME INTEGRALS OVER THEOR. DATA

      DO 50 IB=1,NBT

      IE2=NET(IB)

      IF (IE2.EQ.0) GO TO 50

CVB  for R2

      IF (WHICHR.eq.2) THEN

      CALL VARSUM(AT,AT,AT,AT,1,1,NBTD,1,1,1,IB,1,1,IE2,0,EINCR,
     10.,0.,1,TST(IB),YPL)

CVB

CVB
      ELSE IF (WHICHR.eq.1) THEN
CVB

      CALL VARSUM(YT,AT,AT,AT,1,1,NBTD,1,1,1,IB,1,1,IE2,0,EINCR,
     10.,0.,2,TSTY2(IB),YPL)

CVB
      END IF
CVB

50    CONTINUE

C  START LOOP OVER INNER POTENTIAL VALUES

55    NV0=INT((V02-V01)/VINCR+0.0001)+1
      V0=V01

      DO 160 IV=1,NV0

C      WRITE(6,60) V0
60    FORMAT(34H0THEOR. INNER POTENTIAL SHIFTED BY,1F7.2,3H EV)
                V0R=-V0RR+V0
CVB for R2:

      AR2 = 0.

CVB


      ARPE=0.
      ARAV=0.
      ERANG=0.
      OVLG=0.
      RAZZ=0.
      RANN=0.
      RAVP=0.

      DO 132 I=1,2                                                       260187

CVB  for R2:

         AR2M(I)=0.

CVB

         ERANGM(I)=0.                                                    260187
         ARPEM(I)=0.                                                     260187
         RAZZM(I)=0.                                                     260187
         RANNM(I)=0.                                                     260187
         RAVPM(I)=0.                                                     260187
132   CONTINUE                                                           260187

      DO 62 IB=1,NBTD
      NST1(IB)=0
62    NST2(IB)=0

      ICO=0

C  START LOOP OVER (EXP.) BEAMS, OR THEOR. BEAMS IF NBED=0

      IF (NBED.EQ.0) NBE=NBT
      DO 130 IBE=1,NBE
      KMIT=MITTEL(IBE)                                                  260187

C  IBT INDICATES THE THEOR. BEAM CORRESPONDING TO THE EXP. BEAM IBE

      IBT=IBP(IBE)

      IF (NBED.NE.0) GO TO 64

      NST1(IBE)=1                                                       110380
      NST2(IBE)=NET(IBT)                                                110380
      GO TO 130

64    IF (IBT.NE.0) GO TO 68

C  SKIP SOME BEAMS
      NST1(IBE)=0                                                       110380
      NST2(IBE)=0                                                       110380
      WRITE(6,66) IBE,(BENAME(I,IBE),I=1,5)
66    FORMAT(14H0EXP. BEAM NO.,1I3,10H SKIPPED (,5A4,1H))
      GO TO 130

68    CONTINUE
C69    WRITE(6,70)IBE,(BENAME(I,IBE),I=1,5),IBT,(PQAV(I,IBT),I=1,2)
C70    FORMAT(9H0BEAM NO.,1I3,10H IN EXP. (,5A4,20H), WHICH IS BEAM NO.,
C     11I3,12H IN THEORY (,2F12.4,1H))

CVB  for R2:

      R2(IBE) = 0.

CVB

      RPE(IBE)=0.

C  DETERMINE ENERGY INTERVAL COMMON TO EXP. AND THEORY. THIS INTERVAL,
C  OF LENGTH EET, IS BOUNDED BY THE GRID POINTS (NE1,NE2) AND (NT1,NT2)
C  FOR EXP. AND THEORY, RESP.

      CALL COMNEI(EE,NBED,NEE,ET,NBTD,NET,IBE,IBT,V0,EINCR,NE1,NE2,
     1NT1,NT2,EET(IBE))

      IF (NT2.GT.NT1) GO TO 71                                          040280

C     WRITE(6,705)IBE                                                     .
705   FORMAT(48H0 **** NO OVERLAP THEOR./EXP. IN (EXP.) BEAM NO.,1I3)
C                WRITE(7,4444) (BENAME(I,IBE),I=1,3),IBE,D12,DCO1,
C     *          DCO2,V0R,EMIN,EMAX,0.,0.
4444            FORMAT(3A4,I6,F8.4,4F8.2,1F8.4,8X)
      GO TO 130

71    NST1(IBE)=NT1                                                     110380
      NST2(IBE)=NT2                                                     110380

C  IF ANY IV-CURVE IS TRUNCATED, THE INTEGRALS PERFORMED BEFORE SHOULD
C  BE REDUCED ACCORDINGLY

      NE=NEE(IBE)
      NT=NET(IBT)

CVB
      IF (WHICHR.eq.2) THEN
CVB

CVB  for R2:

C    quantities for normalisation (SE) and quadratic integrals (SE2)

      SS=0.
      SU=0.

      CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,1,1,IBE,1,1,NE1,0,EINCR,
     10.,0.,1,SS,YPL)

      CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,1,1,IBE,1,NE2,NE,0,EINCR,
     10.,0.,1,SU,YPL)

      SE=TSE(IBE)-SS-SU

      SS2 = 0.
      SU2 = 0.

      CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,1,1,IBE,1,1,NE1,0,EINCR,
     10.,0.,2,SS2,YPL)

      CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,1,1,IBE,1,NE2,NE,0,EINCR,
     10.,0.,2,SU2,YPL)

      SE2=TSE2(IBE)-SS2-SU2

CVB

CVB
      ELSE IF (WHICHR.eq.1) THEN
CVB

      SSY2=0.
      SUY2=0.

      CALL VARSUM(YE,AE,AE,AE,1,1,NBED,1,1,1,IBE,1,1,NE1,0,EINCR,
     10.,0.,2,SSY2,YPL)

      CALL VARSUM(YE,AE,AE,AE,1,1,NBED,1,1,1,IBE,1,NE2,NE,0,EINCR,
     10.,0.,2,SUY2,YPL)

      SEY2=TSEY2(IBE)-SSY2-SUY2

CVB
      END IF
CVB

CVB
      IF (WHICHR.eq.2) THEN
CVB

CVB  for R2: produce quantities for normalisation

      SS=0.
      SU=0.

      CALL VARSUM(AT,AT,AT,AT,1,1,NBTD,1,1,1,IBT,1,1,NT1,0,EINCR,
     10.,0.,1,SS,YPL)

      CALL VARSUM(AT,AT,AT,AT,1,1,NBTD,1,1,1,IBT,1,NT2,NT,0,EINCR,
     10.,0.,1,SU,YPL)

      ST=TST(IBT)-SS-SU

C    normalisation factor needed in integral for R2

      C = SE/ST

CVB

CVB
      ELSE IF (WHICHR.eq.1) THEN
CVB

      SSY2=0.
      SUY2=0.

      CALL VARSUM(YT,AT,AT,AT,1,1,NBTD,1,1,1,IBT,1,1,NT1,0,EINCR,
     10.,0.,2,SSY2,YPL)

      CALL VARSUM(YT,AT,AT,AT,1,1,NBTD,1,1,1,IBT,1,NT2,NT,0,EINCR,
     10.,0.,2,SUY2,YPL)

      STY2=TSTY2(IBT)-SSY2-SUY2

CVB   
      END IF
CVB

C  PRODUCE INTEGRALS INVOLVING BOTH EXP. AND THEORY

      NV=NT1-NE1

CVB
      IF (WHICHR.eq.2) THEN
CVB

CVB  for R2:

      CALL VARSUM(AE,AT,AE,AE,1,1,NBED,NBTD,1,1,IBE,IBT,NE1,NE2,
     1NV,EINCR,0.,C,5,S2,YPL)

C  R-FACTOR BASED ON INTEGRAL OF (EXP-C*TH)**2

      R2(IBE)=S2/SE2
      AR2=AR2+WB(IBE)*EET(IBE)*R2(IBE)

CVB

CVB
      ELSE IF (WHICHR.eq.1) THEN
CVB

      CALL VARSUM(YE,YT,YE,YE,1,1,NBED,NBTD,1,1,IBE,IBT,NE1,NE2,
     1NV,EINCR,0.,1.,5,SY2,YPL)

C  R-FACTOR ACCORDING TO PENDRY (MULT. BY 0.5)

      RPE(IBE)=1.0*SY2/(SEY2+STY2)

      ARPE=ARPE+WB(IBE)*EET(IBE)*RPE(IBE)
              RAZZ=RAZZ+SY2*WB(IBE)
              RANN=RANN+SEY2*WB(IBE)+STY2*WB(IBE)

CVB
      END IF
CVB

CVB  KMIT is identification of beam group - only use these values
C  if beamgroup is valid, i.e 1 or 2

      IF ((KMIT .GT. 0).AND.(KMIT .LE. 2)) THEN                           260187

CVB  for R2:

         AR2M(KMIT)=AR2M(KMIT)+WB(IBE)*EET(IBE)*R2(IBE)         

CVB

         ARPEM(KMIT)=ARPEM(KMIT)+WB(IBE)*EET(IBE)*RPE(IBE)              260187
         RAZZM(KMIT)=RAZZM(KMIT)+SY2*WB(IBE)                            260187
         RANNM(KMIT)=RANNM(KMIT)+SEY2*WB(IBE)+STY2*WB(IBE)              260187

      ENDIF        


      ERANG=ERANG+WB(IBE)*EET(IBE)
      OVLG=OVLG+EET(IBE)

      IF ((KMIT .GT. 0).AND.(KMIT .LE. 2)) THEN                          260187
         ERANGM(KMIT)=ERANGM(KMIT)+WB(IBE)*EET(IBE)                     260187
      ENDIF                                                             260187

      OVL=EET(IBE)
C       WRITE(7,4444) (BENAME(I,IBE),I=1,3),IBE,D12,V0R,
C     * EMIN,EMAX,OVL,RPE(IBE)
      EET(IBE)=WB(IBE)*EET(IBE)

130   CONTINUE

C  END OF LOOP OVER BEAMS

CVB  for R2:

      AR2 = AR2/ERANG

CVB

      AR=ARPE/ERANG

C  now do all averaging procedures for beam groups

      DO 137 I=1,2                                                      260187

        IF (ERANGM(I) .LT. 0.0001) THEN

        ARM(I)=5.
        AR2M(I)=5.

        ELSE 

CVB  for R2:

        AR2M(I) = AR2M(I)/ERANGM(I)

CVB

        ARM(I)=ARPEM(I)/ERANGM(I)                                       260187

        END IF

137   CONTINUE                                                          260187

      OVL=ERANG
      XRPE=AR

CVB
      IF (WHICHR.eq.1) THEN
CVB

        RAVP=RAZZ/RANN

CVB   
      ELSE

        RAVP = 5.0

      END IF
CVB

      DO 142 I=1,2                                                       260187

        XRPEM(I)=ARM(I)                                                  260187

                   IF(RANNM(I) .LT. 0.0001) THEN

                   RAVPM(I)=5.0

                   ELSE

                   RAVPM(I)=RAZZM(I)/RANNM(I)                            260187

		   END IF	

142   CONTINUE                                                           260187

CVB  Now check if current R-factor is smaller than the previous ones

C  use either R2 or RPe depending on WHICHR

C  R2:
C  AR2 is total average r-factor, AR2M(1) is average r-factor
C  for integer, AR2M(2) for half-order beams

C  RPe:
C  RAVP is total average r-factor, RAVPM(1) is average r-factor
C  over integer beams, RAVPM(2) is r-factor over half-order beams

C  use variable AR as current rfactor

      IF (WHICHR.eq.1) THEN

        IF (WHICHG.eq.1) THEN
          AR=RAVPM(1)
        ELSE IF (WHICHG.eq.2) THEN
          AR=RAVPM(2)
        ELSE
          AR=RAVP
        END IF

      ELSE IF (WHICHR.eq.2) THEN

        IF (WHICHG.eq.1) THEN
          AR=AR2M(1)
        ELSE IF (WHICHG.eq.2) THEN
          AR=AR2M(2)
        ELSE
          AR=AR2
        END IF

      END IF

CVB Now use AR for search and store all other relevant information
C  accordingly (i.e. total, integer and half order rfactor and 
C  inner potential)

      IF (AR.GE.BARAV) GO TO 145

        BARAV=AR
        BV0(IPOP)=-V0RR+V0

        IF (WHICHR.eq.1) THEN  

          BRGES(IPOP)=RAVP
          BRINS(IPOP)=RAVPM(1)
          BRHAS(IPOP)=RAVPM(2)

        ELSE IF (WHICHR.eq.2) THEN

          BRGES(IPOP)=AR2
          BRINS(IPOP)=AR2M(1)
          BRHAS(IPOP)=AR2M(2)

        END IF

CVB  Continue inner potential loop

 145  V0=V0+VINCR

 160  CONTINUE
C  END OF INNER POTENTIAL LOOP

C  Beam me up to main! Now!

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SEA_RCD(NDOM,NPS,NPRMK,NSTEP,PNUM,VARST,PARIND,RPEIND,
     +                   WSK,WIDT,RMUT,NPAR,PARDEP)

C Global variables
      INTEGER NDOM,NPS,NPRMK,PNUM,NSTEP,NPAR
      INTEGER VARST,PARIND,PARDEP
      DIMENSION VARST(NPRMK),PARIND(NPRMK,NPS),PARDEP(NPRMK)
      REAL RPEIND,WSK,WIDT,RMUT
      DIMENSION WIDT(NPRMK)
      DIMENSION RPEIND(NPS),WSK(NSTEP)
      INTEGER random

C Local variables
      INTEGER MKLP1,MKLP2,MKLP5,INDEX
      REAL MKSQR,FMKSQR,HELP,MKSUM,FMKRN,MKHELP1,MKHELP2
      REAL BACKGROUND, NormG
      integer P(30), IPS, IDOM

C      write(4,*) "now in sea_rcd"
C      write(4,*) NPS,NPRMK,PNUM

      DO 1856 IPOP=1,NPS
      DO 1855 IPARAM=1,PNUM
      
C If parameter is dependent on another, we can skip everything

       IF(PARDEP(IPARAM).ne.0) THEN
         PARIND(IPARAM,IPOP) = PARIND(PARDEP(IPARAM),IPOP)
         GOTO 1855
       ENDIF

c      write(8,*) "PARAMETER",IPARAM,"IN POP",IPOP,"started"

      IF(ABS(WIDT(IPARAM)-1).LE.1E-4) THEN
      width=2.
      ELSE
      width=2.*(RMUT/(NPAR*WIDT(IPARAM)))
      ENDIF
c      WRITE(8,*) 'width=',width



      DO 1851 IPVAL=1,VARST(IPARAM)

c      write(8,*) "start probability of",IPVAL

      MKSQR=(FLOAT(IPVAL)-FLOAT(PARIND(IPARAM,IPOP)))
      FMKSQR=MKSQR*MKSQR

C  open up parameter space completely if r-factor is not good enough

      IF (RPEIND(IPOP).GT.0.8) THEN
         HELP=10000.
      ELSE
         INDEX=VARST(IPARAM)
         HELP=RPEIND(IPOP)+0.05
         HELP=HELP*HELP
      ENDIF
      WSK(IPVAL)=exp(-0.5*FMKSQR/((width*width)*HELP))

c      write(8,*)"probability of value",IPVAL," is",WSK(IPVAL)

 1851 CONTINUE

C Normalization of PB-Distribution
      MKSUM=0.
      DO 1852 IPVAL=1,VARST(IPARAM)
 1852    MKSUM=MKSUM+WSK(IPVAL)

C  MKSUM is now used to normalise distribution to an integral value of 1.
C  Instead, the gaussian distribution is now normalised to an integral value
C  of 1.-background*VARST(IPARAM), where background is a constant 
C  probability that is added to the gaussian probability of each parameter 
C  grid point.

      BACKGROUND = 0.005

      NormG = 1. - BACKGROUND * REAL(VARST(IPARAM)) 
      IF (NormG.lt.0.) THEN
        NormG = 0.
        BACKGROUND = 1/REAL(VARST(IPARAM))
      END IF

      DO 1853 IPVAL=1,VARST(IPARAM)
         WSK(IPVAL)=WSK(IPVAL)*NormG/MKSUM + BACKGROUND
c         write(8,*)"probability of value",IPVAL," is",WSK(IPVAL)
 1853 CONTINUE

c      write(8,*) "distribution normalised"

C  Determination of new random number
C  not that if name of random subroutine is changed, integer declaration of
C  random (see above) must also be changed!

      FMKRN=random()

C      write(8,*)"random",FMKRN

C Determination of new parameter

      MKHELP1=0.
      MKHELP2=0.
      DO 1854 IPVAL=1,VARST(IPARAM)
      MKHELP2=MKHELP2+WSK(IPVAL)*1000.

C      write(8,*) "random",FMKRN," between",MKHELP1," ",MKHELP2,"?"

      IF ( MKHELP1 .LE. FMKRN) THEN
      IF ( MKHELP2 .GE. FMKRN) THEN
      PARIND(IPARAM,IPOP)=IPVAL
      ENDIF
      ENDIF
      MKHELP1=MKHELP2
 1854 CONTINUE
 1855 CONTINUE
 1856 CONTINUE

      RETURN
      END

C
C ========================================================================
C

      SUBROUTINE ORDER(RPEIND,PARIND,BRGES,BRINS,BRHAS,BV0,
     +RPEHELP,BRGESHELP,BRINSHELP,BRHASHELP,BV0HELP,PARHELP,
     +PMISCH)

C Subroutine ORDER performs an ordering of the population's individua
C Excactly the same reordering has to be applied to the old, stored values

      INCLUDE "PARAM"

      INTEGER   PARIND,PARHELP
      DIMENSION PARIND(MNPRMK,MPS),PARHELP(MNPRMK,MPS)
      REAL      RPEIND,RPEHELP,BRGES,BRINS,BRHAS,BV0
      DIMENSION RPEIND(MPS),RPEHELP(MPS),BRGES(MPS),BRINS(MPS),
     +          BRHAS(MPS),BV0(MPS)
      REAL      BRGESHELP,BRINSHELP,BRHASHELP,BV0HELP
      DIMENSION BRGESHELP(MPS),BRINSHELP(MPS),BRHASHELP(MPS),
     +          BV0HELP(MPS)
      INTEGER   PMISCH
      DIMENSION PMISCH(MNDOM,MPS)
      INTEGER   PMHELP
      DIMENSION PMHELP(MNDOM,MPS)

C initialize RPEHELP

      DO 100 IMK=1,MPS
 100     RPEHELP(IMK)=5.

      DO 150 IDOM=1,MNDOM
      DO 150 IPOP=1,MPS
        PMHELP(IDOM,IPOP) = PMISCH(IDOM,IPOP)
 150  CONTINUE

C start ordering

      DO 200 IMK=1,MPS
      FLAG=0
      DO 200 IMK1=1,MPS
         IF ((RPEIND(IMK).LT.RPEHELP(IMK1)).AND.(FLAG.EQ.0)) THEN

            DO 250 IMK2=MPS,IMK1+1,-1

            RPEHELP(IMK2)=RPEHELP(IMK2-1)
            BRGESHELP(IMK2)=BRGESHELP(IMK2-1)
            BRINSHELP(IMK2)=BRINSHELP(IMK2-1)
            BRHASHELP(IMK2)=BRHASHELP(IMK2-1)
            BV0HELP(IMK2)=BV0HELP(IMK2-1)
            DO 270 JMK1=1,MNPRMK
 270        PARHELP(JMK1,IMK2)=PARHELP(JMK1,IMK2-1)
            DO 275 IDOM=1,MNDOM
 275        PMHELP(IDOM,IMK2) = PMHELP(IDOM,IMK2-1)

 250        CONTINUE

            RPEHELP(IMK1)=RPEIND(IMK)
            BRGESHELP(IMK1)=BRGES(IMK)
            BRINSHELP(IMK1)=BRINS(IMK)
            BRHASHELP(IMK1)=BRHAS(IMK)
            BV0HELP(IMK1)=BV0(IMK)
            DO 280 JMK1=1,MNPRMK
 280        PARHELP(JMK1,IMK1)=PARIND(JMK1,IMK)
            DO 285 IDOM = 1,MNDOM
 285        PMHELP(IDOM,IMK1) = PMISCH(IDOM,IMK)

            FLAG=1

         ENDIF
 200  CONTINUE

      DO 300 IDOM=1,MNDOM
      DO 300 IPOP=1,MPS

        PMISCH(IDOM,IPOP) = PMHELP(IDOM,IPOP)

 300  CONTINUE

      RETURN
      END
C
C ========================================================================
C
      SUBROUTINE OUTRF(PARIND,NDOM,NPLACES,NFILES,NPRMK,NPRAS,NPS,NFIL,
     +                 PARTYP,BRGES,BRINS,BRHAS,BV0,AVERNEW,
     +                 IGEN,WHICHG,WHICHR,PMISCH,DMISCH)

      include "PARAM"

      INTEGER NPLACES(NDOM),NFILES,NPRMK,NPS,IGEN

      INTEGER NFIL,PARTYP
      DIMENSION NFIL(NDOM,MNPLACES),PARTYP(NDOM,MNPLACES,NFILES)

      INTEGER WHICHG,WHICHR
      INTEGER PARIND
      DIMENSION PARIND(NPRMK,NPS)
      REAL AVERNEW
      REAL BRGES,BRINS,BRHAS,BV0
      DIMENSION BRGES(NPS)
      DIMENSION BRINS(NPS)
      DIMENSION BRHAS(NPS)
      DIMENSION BV0(NPS)
      character*60 text1, text2, text3, text
      integer NDOM, NPRAS(NDOM), IDOM, KDOM
      real DMISCH
      integer PMISCH(MNDOM,MPS)


      WRITE(4,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      WRITE(4,*) 'CCCCCCCCCCC    GENERATION ',IGEN,'    CCCCCCCCCCC'
      WRITE(4,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      WRITE(4,*)
      WRITE(4,*)

      text1 = 'IND ||  R-Av  |  R-In  | R-Ha   |   V0r   || Proz | '
      text2 = '                                                    '
      text3 = '                                          ||'
      text = text1

      do 100 IDOM = 1, NDOM

      WRITE(4,'(a52,500(a,2i1,a))')
     .  text,
     .  (
     .    
     .      (
     .        (
     .          'P',IPLACE,IFILE,','
     .          , IPARAM = 1, PARTYP(IDOM,IPLACE,IFILE)
     .        ), IFILE = 1, NFIL(IDOM,IPLACE)
     .      ),
     .      'C',IPLACE,IPLACE,','
     .    , IPLACE = 1,NPLACES(IDOM)
     .  )

      text = text2
  100 continue


cas      WRITE(4,5558) ' IND ||  R-Av  |  R-In  | R-Ha   |   V0r   ||',
cas     +((((('P',IPLACE,IFILE,','),IPARAM=1,PARTYP(IPLACE,IFILE)),
cas     +IFILE=1,NFIL(IPLACE)),('C',IPLACE,IPLACE,',')),IPLACE=1,MNPLACES)

 5558 FORMAT(A46,40(A,I1,I1,A1))

      WRITE(4,5557) '---------------------------------------------------
     +--------------------------------------------------'

 5557 FORMAT(A100)


      do 110 IPOP = 1, NPS

      OFFSET = 0

      DO 5569 IDOM = 1, NDOM

        if (IDOM .eq. 1) then

          WRITE(4,5559) IPOP,'|| ',BRGES(IPOP),' | ',BRINS(IPOP),
     .    ' | ',BRHAS(IPOP),' | ',BV0(IPOP),'|| ', 
     .    nint(PMISCH(IDOM,IPOP)*DMISCH*100),
     .    '% | ',
     .    (PARIND(IPARAM,IPOP),' ',IPARAM=1,NPRAS(IDOM),1)

        else

          OFFSET = OFFSET + NPRAS(IDOM - 1)

          write(4,5560) text3,
     .    nint(PMISCH(IDOM,IPOP)*DMISCH*100),
     .    '% | ',
     .    (PARIND(OFFSET+IPARAM,IPOP),' ',IPARAM=1,NPRAS(IDOM),1)

        end if

 5559 FORMAT(I3,A4,F6.4,A3,F6.4,A3,F6.4,A3,F7.3,A4,i3,a4,500(I3,A1))
 5560 format(a44,i4,a4,500(I3,A1))

 5569 CONTINUE
  110 continue

      write(4,*)

      IF (WHICHR.eq.1) THEN

        IF (WHICHG.eq.1) THEN

          write(4,'(A34,I6,A4,F10.8)') 
     +    'Average RPe-INTEGER of GENERATION ',
     +    IGEN,' :  ',AVERNEW

        ELSE IF (WHICHG.eq.2) THEN

          write(4,'(A34,I6,A4,F10.8)') 
     +    'Average RPe-HALF-O. of GENERATION ',
     +    IGEN,' :  ',AVERNEW

        ELSE

          write(4,'(A26,I6,A4,F10.8)') 
     +    'Average RPe of GENERATION ',
     +    IGEN,' :  ',AVERNEW

        ENDIF

      ELSE IF (WHICHR.eq.2) THEN

        IF (WHICHG.eq.1) THEN

          write(4,'(A33,I6,A4,F10.8)') 
     +    'Average R2-INTEGER of GENERATION ',
     +    IGEN,' :  ',AVERNEW

        ELSE IF (WHICHG.eq.2) THEN

          write(4,'(A33,I6,A4,F10.8)') 
     +    'Average R2-HALF-O. of GENERATION ',
     +    IGEN,' :  ',AVERNEW

        ELSE

          write(4,'(A25,I6,A4,F10.8)') 
     +    'Average R2 of GENERATION ',
     +    IGEN,' :  ',AVERNEW

        ENDIF

      ENDIF

      write(4,*)
      RETURN
      END
C
C ================================================================================
C 
      SUBROUTINE READSC(NDOM,NPLACES,NFILES,INFILE,NSURF,
     +                  IFORM,PNUM,VARST,NPRMK,NPRAS,PARTYP,NPS,
     +                  PARIND,STAFLA,OUTINT,FILREL,WHICHG,WHICHR,
     +                  NFIL,NCONCS,CONC,DMISCH,MAXGEN,SEANAME,
     +                  NPAR,RMUT,INIT)

      include "PARAM"

C  Dimension statements

      INTEGER NPLACES(NDOM),NFILES,NPRMK,NPS,NCONCS,NDOM
      integer NPRAS(NDOM), MAXGEN

C  WHICHG determines which beam group to use
C  WHICHR determines whether RPe or R2 are used

      INTEGER NPAR
      INTEGER PNUM,STAFLA,WHICHG,WHICHR
      INTEGER PARTYP,FILREL,OUTINT
      DIMENSION PARTYP(NDOM,MNPLACES,NFILES),FILREL(NDOM,MNPLACES)
      CHARACTER*15 INFILE(NDOM,MNPLACES,NFILES)
      INTEGER NSURF
      DIMENSION NSURF(NDOM,MNPLACES)
      INTEGER IFORM
      DIMENSION IFORM(NDOM,MNPLACES,NFILES)
      INTEGER VARST
      DIMENSION VARST(NPRMK)
      INTEGER HELP
      INTEGER PARIND
      DIMENSION PARIND(NPRMK,NPS)
      character*100 text
      integer IDOM
      character*10 SEANAME
      REAL RMUT
      INTEGER INIT

CVB  NFIL(NDOM,IPLACE) is number of files per domain NDOM and atomic site IPLACE

      INTEGER NFIL
      DIMENSION NFIL(NDOM,MNPLACES)
CVB

CVB
      REAL CONC
      DIMENSION CONC(NDOM,NCONCS,MNPLACES,NFILES)

C  (note NCONCS is MNCONCS)
CVB

C  local variable - number of conc steps for each place
C  directly saved as part of array VARST

      INTEGER NCONC

CVB

C --- DMISCH ist die Mischungsverhaeltnissschrittweite in Prozentpunkten
      real DMISCH
      integer HELPDMISCH

C  begin work

      PNUM = 0

      OPEN(21,FILE='search.steu',STATUS='UNKNOWN')

      READ(21,'(I5)') NPAR
      READ(21,'(F7.4)') RMUT
      READ(21,'(I5)') INIT
      READ(21,'(I5)') WHICHR
      READ(21,'(I5)') WHICHG
      READ(21,'(I6)') OUTINT
      READ(21,'(I6)') MAXGEN

      READ(21,'(I5)') HELPDMISCH
      DMISCH = FLOAT(HELPDMISCH) / 100.0

      READ(21,'(a10)') SEANAME

      READ(21,'(I5)') NDOM
      do 1900 IDOM = 1, NDOM
        NPRAS(IDOM) = 0
 1900 continue

      write(6,'(I5)') whichr
      write(6,'(I5)') whichg
      write(6,'(I5)') outint
      write(6,'(I8)') MAXGEN
      write(6,'(f7.4)') DMISCH
      write(6,'(a10)') SEANAME

C --- Schleife ueber alle Domaenen
      do 2000 IDOM = 1, NDOM

        read(21,'(a100)') text
        write(6,'(a100)') text

        read(21, '(i5)') NPLACES(IDOM)

        DO 1366 IPLACE = 1, NPLACES(IDOM)

          READ(21,'(a100)') text
          READ(21,'(I3)') NSURF(IDOM,IPLACE)
          READ(21,'(I3)') FILREL(IDOM,IPLACE)
          READ(21,'(I3)') NFIL(IDOM,IPLACE)

          if (NFIL(IDOM,IPLACE).gt.NFILES) then
            write(8,*) "Parameter MNFILES smaller than number of ",
     +                 "delta amplitude storage files required for"
            write(8,*) "site no. ",IPLACE," in domain no. ", IDOM,"."
            write(8,*) "Please correct."
            stop
          end if

          write(6,'(a100)') text
          write(6,'(I3)') NSURF(IDOM,IPLACE)
          write(6,'(I3)') FILREL(IDOM,IPLACE)
          write(6,'(I3)') NFIL(IDOM,IPLACE)

          DO 1368 IFILE = 1, NFIL(IDOM,IPLACE), 1

            READ(21,'(a100)') text
            READ(21,'(A15)') INFILE(IDOM,IPLACE,IFILE)
            READ(21,'(I3)')  IFORM(IDOM,IPLACE,IFILE)
            READ(21,'(I3)')  PARTYP(IDOM,IPLACE,IFILE)      

            WRITE(6,'(a100)') text
            WRITE(6,'(A15)') INFILE(IDOM,IPLACE,IFILE)
            WRITE(6,'(I3)')  IFORM(IDOM,IPLACE,IFILE)
            WRITE(6,'(I3)')  PARTYP(IDOM,IPLACE,IFILE)      

            DO 1369 IPARAM = 1, PARTYP(IDOM,IPLACE,IFILE)
              READ(21,'(2I3)') VARST(PNUM + IPARAM)
              write(6,'(2I3)') VARST(PNUM + IPARAM)
 1369       CONTINUE

            PNUM = PNUM + PARTYP(IDOM,IPLACE,IFILE)
            NPRAS(IDOM) = NPRAS(IDOM) + PARTYP(IDOM,IPLACE,IFILE)

 1368     CONTINUE

CVB  read in concentration steps now!
C    NCONC is current number of concentration steps per file
C    CONC(IDOM,ICONC,IPLACE,IFILE) is concentration of file IFILE
C    on current site IPLACE for current concentration step and domain IDOM

          READ(21,'(a100)') text

          READ(21,'(I3)') NCONC

CVB  error check for MNCONCS done here instead of CheckVal for the sake
C    of clarity

          IF (NCONC.gt.MNCONCS) THEN
            write(8,*) "Parameter MNCONCS too small!"
            STOP
          END IF

          DO 1367 ICONC = 1, NCONC, 1
            READ(21,'(5F7.4)') (CONC(IDOM,ICONC,IPLACE,IFILE),
     .                          IFILE = 1, NFIL(IDOM,IPLACE))
            write(6,'(5F7.4)') (CONC(IDOM,ICONC,IPLACE,IFILE),
     .                          IFILE = 1, NFIL(IDOM,IPLACE))
 1367     CONTINUE

C  increment number of parameters again so concentration will be
C  handled similar to all other parameters - therefore NCONC is
C  saved in array VARST directly

        PNUM = PNUM + 1
        NPRAS(IDOM) = NPRAS(IDOM) + 1
        VARST(PNUM) = NCONC

CVB

 1366   CONTINUE

 2000 continue
C --- Ende der Schleife ueber alle Domaenen

C --- Pro Domaene ein zusaetzlicher Parameter
      do 2010 IDOM = 1, NDOM
        PNUM = PNUM + 1
        VARST(PNUM) = INT(1.0/DMISCH + 1.0 + 1.0E-02)
 2010 continue

      READ(21,'(A10)') text
      READ(21,'(I3)') STAFLA

      write(6,'(A10)') text
      write(6,'(I3)') STAFLA

      IF (STAFLA.EQ.1) THEN
        DO 1370 IPOP = 1, NPS
          READ(21,'(500I3)') (PARIND(IPARAM,IPOP),IPARAM=1,NPRMK)
          write(6,'(500I3)') (PARIND(IPARAM,IPOP),IPARAM=1,NPRMK)
 1370   CONTINUE
      ENDIF

cvb      write(6,*) 'Inhalt von VARST():'
cvb      write(6,'(100i3)') (VARST(i), i = 1, PNUM)

cvb      write(6,*) 'Inhalt von NPRAS():'
cvb      write(6,'(100i3)') (NPRAS(IDOM), IDOM = 1, NDOM)

      CLOSE(21)

      RETURN
      END
C
C #######################################################################
C
      SUBROUTINE PREEXP(AE,EE,NBED,NEE,BENAME,NBEA,IPR,ISMOTH,EINCR,VI,
     +     YE,NDATA,TSE,TSE2,TSEY2,XPL,YPL,AEP,NNN,NBE)

C  Declaration of global variables

      INTEGER IPR
      INTEGER ISMOTH
      INTEGER NDATA
      INTEGER NBED
      INTEGER NBE           
      INTEGER NEE
      DIMENSION NEE(NBED)
      INTEGER NBEA 
      DIMENSION NBEA(NBED)
      DIMENSION BENAME(5,NBED)
      REAL EINCR
      REAL VI
      REAL AE
      DIMENSION AE(NBED,NDATA)
      REAL EE
      DIMENSION EE(NBED,NDATA)
      REAL YE
      DIMENSION YE(NBED,NDATA)
      REAL TSEY2
      DIMENSION TSEY2(NBED)

CVB  for R2:

      REAL TSE,TSE2
      DIMENSION TSE(NBED),TSE2(NBED)

CVB

      REAL XPL,YPL
      DIMENSION XPL(NDATA),YPL(1,NDATA)
      REAL AEP
      DIMENSION AEP(NBED,NDATA)        
      INTEGER NNN(NDATA)

C AVERAGE DATA FROM DIFFERENT EXPERIMENTS AND ORDER BY INCREASING ENERGY

      CALL EXPAV(AE,EE,NBED,NEE,BENAME,NBEA,NBE,IPR,XPL,NNN,NDATA)

C  SMOOTH EXP. DATA THE NUMBER OF TIMES DESIRED
      IF (ISMOTH.EQ.0) GO TO 15

      DO 10 I=1,ISMOTH
10    CALL SMOOTH(AE,EE,NBED,NBE,NEE,IPR)
15    CONTINUE

C  INTERPOLATE EXP. DATA TO WORKING GRID (MULTIPLES OF EINCR EV)
      CALL INTPOL(AE,1,NBED,NEE,1,NBE,EE,EINCR,IPR,XPL,YPL)

C  PRODUCE 1ST DERIVATIVE OF EXPERIMENTAL DATA
      CALL DER(AE,NEE,1,NBED,1,NBE,AEP,EINCR)

C  PRODUCE PENDRY Y FUNCTION FOR EXP. DATA
      CALL YPEND(AE,AEP,1,NBED,1,NBE,NEE,EE,YE,VI,IPR)

      DO 20 IB=1,NBE
      IE2=NEE(IB)
      IF (IE2.EQ.0) GO TO 20

CVB  for R2 :

      TSE(IB)=0.

      CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,1,1,IB,1,1,IE2,0,EINCR,
     10.,0.,1,TSE(IB),YPL)

      TSE2(IB)=0.

      CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,1,1,IB,1,1,IE2,0,EINCR,
     10.,0.,2,TSE2(IB),YPL)

CVB  for RPe:

      TSEY2(IB)=0.
      
      CALL VARSUM(YE,AE,AE,AE,1,1,NBED,1,1,1,IB,1,1,IE2,0,EINCR,
     10.,0.,2,TSEY2(IB),YPL)

20    CONTINUE

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE OUTIN(ATSMK,ESMK,NBED,NDATT)
      DIMENSION ATSMK(NBED,NDATT),ESMK(NDATT)
      DO 1919 JMK=1,NDATT
1919     WRITE(6,552) ESMK(JMK),(ATSMK(N,JMK),N=1,NBED)
  552 FORMAT (1F7.2,4E14.5,/,100(5E14.5,/))
      END
C------------------------------------------------------------------------
C  SUBROUTINE READE INPUTS EXPERIMENTAL IV-CURVES
      SUBROUTINE READE(AE,EE,NBED,NEE,NBEA,BENAME,IPR)

      INCLUDE "PARAM"

      CHARACTER*80 TEXT
      CHARACTER*80 FMT
      DIMENSION AE(NBED,MNDATA),EE(NBED,MNDATA),NEE(NBED),NBEA(NBED),
     1BENAME(5,NBED)
C  FOR CDC ONLY  ***************************
C     LEVEL 2, AE,EE
30    FORMAT(20(25I3,/))

C  READ IN EXP. BEAM AVERAGING INFORMATION. IF NBEA(I)=NBEA(J), EXP.
C  BEAMS I AND J WILL BE AVERAGED TOGETHER. THE RELATION NBEA(J).GT.
C  NBEA(I) IF J.GT.I MUST HOLD. IF NBEA(I)=0, EXP. BEAM I WILL BE SKIPPED
C  LATER ON

C  READ AND PRINT DESCRIPTION OF EXPERIMENT
      READ(12,'(A80)') TEXT
      WRITE(6,*) TEXT

      READ(12,30)(NBEA(I),I=1,NBED)
      WRITE(6,40)NBED
C     WRITE(8,40)NBED
40    FORMAT(1H0,I4,25H EXP. BEAMS TO BE READ IN)                         040481
      WRITE(6,50)(NBEA(I),I=1,NBED)
C     WRITE(8,50)(NBEA(I),I=1,NBED)
50    FORMAT(31H0AVERAGING SCHEME OF EXP. BEAMS,5(25I3,/))
C  READ INPUT FORMAT OF EXP. INTENSITIES
      READ(12,'(A80)') FMT
      DO 90 IB=1,NBED
C  READ AND PRINT NAME OF THIS BEAM
      READ(12,10)(BENAME(I,IB),I=1,5)
10    FORMAT(19A4)
      WRITE(6,70)IB,(BENAME(I,IB),I=1,5)
C     WRITE(8,70)IB,(BENAME(I,IB),I=1,5)
70    FORMAT(11H0EXP. BEAM ,1I3,2H (,5A4,1H))
C  READ IN NO. OF DATA POINTS TO BE INPUT FOR THIS BEAM AND CONSTANT
C  CORRECTION FACTOR FOR INTENSITIES. THIS FACTOR IS MEANT TO ALLOW
C  NORMALIZATION TO INTENSITIES OF THE ORDER OF 1 (NOT NECESSARY, BUT
C  SAFE), AND TO MATCH UP CURVES TO BE AVERAGED TOGETHER WHEN THEIR
C  ENERGY RANGES DIFFER (TO AVOID DISCONTINUITIES AT ENERGIES WHERE
C  THE NUMBER OF CURVES AVERAGED TOGETHER CHANGES)
      READ(12,35)NEE(IB),FAC
35    FORMAT(1I4,1E13.4)

      IF (NEE(IB).gt.MNDATA) THEN
        write(8,*) "Exp. beam no. ",IB," contains more data points",
     +             " than dimension MNDATA allows for!"
        write(8,*) "Please correct!"
        STOP
      END IF

      N=NEE(IB)
C  READ (AND MAYBE PRINT) EXP. INTENSITIES
      READ(12,FMT)(EE(IB,IE),AE(IB,IE),IE=1,N)
      IF (IPR.LT.1) GO TO 82
      WRITE(6,80)(EE(IB,IE),AE(IB,IE),IE=1,N)
80    FORMAT(31H0 EXP. ENERGIES AND INTENSITIES,
     1/,83(5(1F7.2,1E13.4,3X),/))
82    DO 86 IE=1,N
86    AE(IB,IE)=AE(IB,IE)*FAC
90    CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READRF(EMIN,EMAX,EINCR,
     1IPR,VI,V0RR,V01,V02,VINCR,ISMOTH,EOT,
     2IBP,WB,NBTD,NBED,KAV,NSS,MITTEL)

CVB  This subroutine reads information for r-factor determination from file rf.info.

      REAL EMIN,EMAX,EINCR
      REAL VI,V0RR,V01,V02,VINCR
      INTEGER IPR,ISMOTH,EOT

      INTEGER IBP(NBED),KAV(NBTD)
      INTEGER MITTEL(NBED)
      REAL WB(NBED)

C  begin readin from stdin

      READ(12,'(F7.2)') EMIN
      READ(12,'(F7.2)') EMAX
      READ(12,'(F7.2)') EINCR
      READ(12,'(I3)')   IPR
      READ(12,'(F7.4)') VI
      READ(12,'(F7.4)') V0RR
      READ(12,'(F7.4)') V01
      READ(12,'(F7.4)') V02
      READ(12,'(F7.4)') VINCR
      READ(12,'(I3)')   ISMOTH
      READ(12,'(I3)')   EOT
      READ(12,'(2I3)')  NBTD,NBED

      WRITE(6,8)NBTD,NBED                                               040481
8     FORMAT(13H0NBTD,NBED = ,2I4)                                      040481

      READ(12,'(25I3)') (KAV(J),J=1,NBTD)

      WRITE(6,5)(KAV(J),J=1,NBTD)

5     FORMAT(33H0AVERAGING SCHEME OF THEOR. BEAMS,5(25I3,/))

C Values for averaging of beams in integer and half order beams
      READ(12,*) 
      READ(12,'(25I3)') (MITTEL(I),I=1,NBED)

C beam order IBP, beam weights WR

      READ(12,*)
      READ(12,'(25I3)')   (IBP(I),I=1,NBED)
      READ(12,'(25F3.1)') (WB(I),I=1,NBED)

C  NSS is artifact present in many r-factor subroutines

      NSS = 1

      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE SMOOTH SMOOTHES A SET OF DATA (GIVEN ON A NON-UNIFORM
C  GRID, BUT CHOOSING A SIMPLIFIED FORMULA WHEN EQUAL INTERVALS ARE
C  FOUND) BY WEIGHTED THREE-POINT AVERAGING
      SUBROUTINE SMOOTH(AE,EE,NBED,NBE,NEE,IPR)

      INCLUDE "PARAM"

      DIMENSION AE(NBED,MNDATA),EE(NBED,MNDATA),NEE(NBE)
C  FOR CDC ONLY  ***************************
C     LEVEL 2, AE,EE
      DO 40 IB=1,NBE
      N=NEE(IB)-1
C  IF TOO FEW POINTS, NO SMOOTHING POSSIBLE
      IF (N.LT.3) GO TO 40
      AM=AE(IB,1)                                                       020780
      DO 30 IE=2,N
      AF=AE(IB,IE)
      E21=EE(IB,IE)-EE(IB,IE-1)
      E32=EE(IB,IE+1)-EE(IB,IE)
      IF (ABS(E21-E32).LT.0.0001) GO TO 10
      AE(IB,IE)=0.5*(AF+(E32*AM+E21*AE(IB,IE+1))/(E21+E32))             020780
      GO TO 20
10    AE(IB,IE)=0.25*(2.*AF+AM+AE(IB,IE+1))                             020780
20    CONTINUE
      AM=AF                                                             020780
30    CONTINUE
40    CONTINUE
      IF (IPR.LT.2) GO TO 70
      DO 50 IB=1,NBE
      N=NEE(IB)
      IF (N.EQ.0) GO TO 50
      WRITE(6,60)IB,(EE(IB,IE),AE(IB,IE),IE=1,N)
60    FORMAT(48H0EXP. ENERG. AND INTENS. AFTER SMOOTHING IN BEAM,1I3,
     1/,50(5(1F7.2,1E13.4,3X),/))
50    CONTINUE
70    RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE INTPOL INTERPOLATES LEED INTENSITIES ONTO A WORKING        230480
C  GRID (WITH STEPS OF EINCR EV)
      SUBROUTINE INTPOL(A,NS,NBD,NE,IS,NB,E,EINCR,IPR,X,WORYT)

      INCLUDE "PARAM"

      DIMENSION WORYT(MNDATA),X(MNDATA),A(NS,NBD,MNDATA),E(NBD,MNDATA),
     +          NE(NBD)
C  FOR CDC ONLY  ***************************
C     LEVEL 2, A,E
      ITIL=0
      ITIH=0
      DO 60 IB=1,NB
      NEM=NE(IB)
C  FIND FIRST NON-ZERO INTENSITY (FOR THEORY, WHERE NON-EMERGENCE OF
C  CURRENT BEAM CAN OCCUR)
      DO 30 IE=1,NEM
      IMIN=IE
      IF (A(IS,IB,IE).GT.1.E-6) GO TO 40
30    CONTINUE
40    CONTINUE
      IF (IMIN.EQ.NEM) GO TO 9                                           200480
      LMIN=INT((E(IB,IMIN)-0.01)/EINCR)+1
      LMIN=MAX0(LMIN,0)
      XMIN=FLOAT(LMIN)*EINCR
      NEM=NE(IB)
      LMAX=INT((E(IB,NEM)+0.01)/EINCR)
      XMAX=FLOAT(LMAX)*EINCR
C  NPTS IS NO. OF POINTS USED ON THE INTERPOLATION GRID
      NPTS=LMAX-LMIN+1
      NE(IB)=NPTS
      DO 5 I=IMIN,NEM
      X(I-IMIN+1)=E(IB,I)
5     WORYT(I-IMIN+1)=A(IS,IB,I)
      NEM=NEM-IMIN+1
      IF (NEM.GE.1) GO TO 8
C  TOO FEW ENERGY VALUES FOR INTERPOLATION. THIS BEAM WILL BE SKIPPED
C  FROM NOW ON
9     NE(IB)=0                                                           200480
C     WRITE(6,6)                                                         230480
C     WRITE(8,6)
6     FORMAT(59H0** IN PRESENT BEAM TOO FEW ENERGY VALUES FOR INTERPOLAT 230480
     1ION)
      GO TO 60
8     CONTINUE
      XVAL=XMIN-EINCR
      DO 10 I = 1, NPTS
      XVAL = XVAL+EINCR
C  INTERPOLATE (AND SET NEGATIVE INTENSITIES TO ZERO)
      A(IS,IB,I) = YVAL(XVAL,WORYT,X,NEM,ITIL,ITIH)
10    IF (A(IS,IB,I).LT.0.0) A(IS,IB,I)=0.0
      E(IB,1)=XMIN
      DO 50 IE=2,NPTS
50    E(IB,IE)=E(IB,IE-1)+EINCR
60    CONTINUE
      IF (IPR.LT.2) GO TO 90
      DO 70 IB=1,NB
      N=NE(IB)
      IF (N.EQ.0) GO TO 70
      WRITE(6,80)IB,(E(IB,IE),A(IS,IB,IE),IE=1,N)
80    FORMAT(40H0INTENSITIES AFTER INTERPOLATION IN BEAM,1I4,
     1/,50(5(1F7.2,1E13.4,3X),/))
70    CONTINUE
90    RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE DER CALCULATES 1ST DERIVATIVE (AFTER ZANAZZI-JONA)
      SUBROUTINE DER(Y,NE,NS,NBD,IS,NB,Y1,H)

      INCLUDE "PARAM"

      DIMENSION Y(NS,NBD,MNDATA),NE(NBD),Y1(NBD,MNDATA)
      DO 60 I=1,NB
      IF (NE(I).EQ.0) GO TO 60
      IF (NE(I).LT.23) GO TO 30
      JF=NE(I)-3
      DO 10 J=4,JF
10    Y1(I,J)=(Y(IS,I,J+3)-9.*Y(IS,I,J+2)+45.*Y(IS,I,J+1)
     1-45.*Y(IS,I,J-1)+9.*Y(IS,I,J-2)-Y(IS,I,J-3))/(60.*H)
      DO 20 J=1,3
      Y1(I,J)=(2.*Y(IS,I,J+3)-9.*Y(IS,I,J+2)+18.*Y(IS,I,J+1)
     1-11.*Y(IS,I,J))/(6.*H)
      M=NE(I)-J+1
20    Y1(I,M)=(11.*Y(IS,I,M)-18.*Y(IS,I,M-1)+9.*Y(IS,I,M-2)
     1-2.*Y(IS,I,M-3))/(6.*H)
      GO TO 60
30    IF (NE(I).LT.3) GO TO 50
      JF=NE(I)-1
      DO 40 J=2,JF
40    Y1(I,J)=(Y(IS,I,J+1)-Y(IS,I,J-1))/(2.*H)
      Y1(I,1)=(-3.*Y(IS,I,1)+4.*Y(IS,I,2)-Y(IS,I,3))/(2.*H)
      Y1(I,JF+1)=(-3.*Y(IS,I,JF+1)+4.*Y(IS,I,JF)-Y(IS,I,JF-1))/(2.*H)
      GO TO 60
50    Y1(I,1)=(Y(IS,I,2)-Y(IS,I,1))/H
      Y1(I,2)=Y1(I,1)
60    CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE YPEND CALCULATES THE PENDRY Y FUNCTION
C  Y = (A/AP) / ((A/AP)**2 + VI), WHERE AP/A IS THE LOGARITHMIC
C  DERIVATIVE OF THE (TABULATED) FUNCTION A
      SUBROUTINE YPEND(A,AP,NS,NBD,IS,NB,NE,E,Y,VI,IPR)

      INCLUDE "PARAM"

      DIMENSION A(NS,NBD,MNDATA),AP(NS,NBD,MNDATA),NE(NBD),E(NBD,MNDATA)
      DIMENSION Y(NBD,MNDATA)
C     LEVEL 2, A,AP,E,Y
      DO 25 IB=1,NB
      N=NE(IB)
      IF (N.EQ.0) GO TO 25
      DO 20 IE=1,N
      AF=A(IS,IB,IE)
      IF (ABS(AF).LT.1.E-7) GO TO 10
      AF=AP(IS,IB,IE)/AF
      Y(IB,IE)=AF/(1.+VI*VI*AF*AF)
      GO TO 20
10    APF=AP(IS,IB,IE)
      IF (APF.GT.1.E-7) GO TO 15
      Y(IB,IE)=0.
      GO TO 20
15    AF=AF/APF
      Y(IB,IE)=AF/(AF*AF+VI*VI)
20    CONTINUE
25    CONTINUE
      IF (IPR.LT.2) GO TO 50
      DO 30 IB=1,NB
      N=NE(IB)
      IF (N.EQ.0) GO TO 30
      WRITE(6,40)IB,(E(IB,IE),Y(IB,IE),IE=1,N)
40    FORMAT(26H PENDRY Y FUNCTION IN BEAM,1I3,/,100(5(1F7.2,1E13.4,3X),
     1/))
30    CONTINUE
50    RETURN
      END
C-----------------------------------------------------------------------
C  FUNCTION YVAL INTERPOLATES
      FUNCTION YVAL(X, WORY, WORX, LENGTH,ITIL,ITIH)
      DIMENSION WORY(LENGTH), WORX(LENGTH)
      COMMON /DATBLK/ X0,X1,X2,X3,Y0,Y1,Y2,Y3
C  IF FEWER THAN FOUR GRID POINTS AVAILABLE, USE 2ND OR 1ST ORDER
C  POLYNOMIAL INTERPOLATION
      IF (LENGTH.LT.4) GO TO 10
C  FIND REQUIRED INTERPOLATION INTERVAL
      CALL BINSRX(IL, IH, X, WORX, LENGTH)
C  SKIP NEXT STEP IF SAME INTERVAL IS FOUND AS LAST TIME
      IF(IL .EQ. ITIL .AND. IH .EQ. ITIH) GO TO 5
C  FIND FOUR NEAREST GRID POINTS AND CORRESPONDING INTENSITIES
C  FOR 3RD-ORDER POLYNOMIAL INTERPOLATION
      CALL STFPTS(IL, IH, WORX, WORY, LENGTH)
C  DO ACTUAL 3RD-ORDER POLYNOMIAL INTERPOLATION
 5    Y = XNTERP(X)
      ITIH = IH
      ITIL = IL
      YVAL = Y
      GO TO 20
10    CALL XNTRP2(X,Y,WORX,WORY,LENGTH)
      YVAL=Y
20    RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE BINSRX FINDS A REQUIRED INTERPOLATION INTERVAL
C  BY BINARY SEARCH (SUCCESSIVE HALVING OF INITIAL INTERVAL)
      SUBROUTINE BINSRX(IL, IH, X, WORX, LENGTH)
      DIMENSION WORX(LENGTH)
      I = (LENGTH + 1) / 2
      IHI =  (I + 1) / 2
      IF(X.LT. WORX(1) .OR. X .GE. WORX(LENGTH)) GO TO 100
20    CONTINUE
      IF(I .LT. 1) GO TO 40
      IF(I .GE. LENGTH) GO TO 30
      IF(X .GE.WORX(I) .AND. X .LE. WORX(I + 1)) GO TO 50
      IF(X .LT. WORX(I)) GO TO 30
40    I = IHI + I
      IHI = (IHI + 1) / 2
      GO TO 20
30    I = I - IHI
      IHI =(IHI + 1) / 2
      GO TO 20
50    IL = I
      IH = I + 1
      GO TO 110
100   IF(X .LT. WORX(1)) GO TO 105
102   IL = LENGTH - 1
      IH = LENGTH
      GO TO 110
105   IL = 1
      IH = 2
110   RETURN
      END
C------------------------------------------------------------------------
C  SUBROUTINE STFPTS FINDS, GIVEN THE INTERPOLATION INTERVAL, THE
C  FOUR NEAREST GRID POINTS AND THE CORRESPONDING ORDINATE VALUES
      SUBROUTINE STFPTS(IL, IH, WORX, WORY, LENGTH)
      DIMENSION WORX(LENGTH), WORY(LENGTH), TEMP(2, 4)
      COMMON / DATBLK / X0, X1, X2, X3, Y0, Y1, Y2, Y3
      I = IL - 1
      IF(IL .LE. 1) I = I + 1
      IF(IH .GE. LENGTH) I = I - 1
      DO 10 K = 1, 4
      N = K + I - 1
      TEMP(1,K) = WORX(N)
10    TEMP(2, K) = WORY(N)
      X0 = TEMP(1, 1)
      X1 = TEMP(1, 2)
      X2 = TEMP(1, 3)
      X3 = TEMP(1, 4)
      Y0 = TEMP(2, 1)
      Y1 = TEMP(2, 2)
      Y2 = TEMP(2, 3)
      Y3 = TEMP(2, 4)
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE XNTERP PERFORMS 3RD-ORDER POLYNOMIAL INTERPOLATION
      FUNCTION XNTERP(X)
      COMMON / DATBLK / X0, X1, X2, X3, Y0, Y1, Y2, Y3
      TERM = Y0
      FACT1 = X - X0
      FACT2 = (Y1 - Y0) / (X1 - X0)
      TERM = TERM + FACT1 * FACT2
      FACT1 = FACT1 * (X - X1)
      FACT2 = ((Y2 - Y1)/(X2 - X1) - FACT2) / (X2 - X0)
      TERM = TERM + FACT1 * FACT2
      FACT1 = FACT1 * (X - X2)
      TEMP = ((Y3 - Y2)/(X3 - X2) - (Y2 - Y1)/(X2 - X1))/(X3 - X1)
      FACT2 = (TEMP - FACT2) / (X3 - X0)
      TERM = TERM + FACT1 * FACT2
      XNTERP = TERM
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE XNTRP2 PERFORMS 2ND OR 1ST ORDER POLYNOMIAL INTERPOLATION
      SUBROUTINE XNTRP2(X,Y,XS,YS,N)
      DIMENSION XS(N),YS(N)
      IF (N.GT.2) GO TO 10
      Y=(YS(2)-YS(1))*(X-XS(1))/(XS(2)-XS(1))+YS(1)
      RETURN
10    A=(YS(1)-YS(2))/(XS(1)-XS(2))/(XS(2)-XS(3))-
     1 (YS(1)-YS(3))/(XS(1)-XS(3))/(XS(2)-XS(3))
      B=(YS(1)-YS(2))/(XS(1)-XS(2))-A*(XS(1)+XS(2))
      C=YS(1)-(A*XS(1)+B)*XS(1)
      Y=C+X*(B+A*X)
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE EXPAV AVERAGES INTENSITIES FROM DIFFERENT EXPERIMENTS.    100281
C  IT ALSO PUTS THEM IN ORDER OF INCREASING ENERGY                      100281

      SUBROUTINE EXPAV(AE,EE,NBED,NEE,BENAME,NBEA,NBE,IPR,A,NV,NDATA)

      INCLUDE "PARAM"

      DIMENSION AE(NBED,NDATA),EE(NBED,NDATA),NEE(NBED),NBEA(NBED),
     1A(NDATA), NV(NDATA),BENAME(5,NBED)

C     LEVEL 2, AE,EE
      NBE=1
      DO 90 IB=1,NBED
      IF (NBEA(IB).EQ.0) GO TO 90                                       040280
      EMIN=1.E6                                                         040280
      EMAX=0.                                                           040280

C  LOOP OVER BEAMS TO BE AVERAGED TOGETHER WITH BEAM IB TO OBTAIN       100281
C  MINIMUM AND MAXIMUM ENERGIES                                         100281

      DO 45 IB1=IB,NBED                                                 040280
      IF (NBEA(IB1).NE.NBEA(IB)) GO TO 45                               040280
      N=NEE(IB1)                                                        100281
      DO 43 IE=1,N                                                      100281
      IF (EMIN.LE.EE(IB1,IE)) GO TO 40                                  100281
      EMIN=EE(IB1,IE)                                                   100281
      IBMIN=IB1                                                         040280
40    CONTINUE                                                          100281
      IF (EMAX.LT.EE(IB1,IE)) EMAX=EE(IB1,IE)                           100281
43    CONTINUE                                                          100281
45    CONTINUE                                                          040280
C  FIND ENERGY INCREMENT                                                100281
      EMIP=1.E6                                                         100281
      N=NEE(IBMIN)                                                      100281
      DO 48 IE=1,N                                                      100281
      EP=EE(IBMIN,IE)                                                   100281
48    IF (EP.LT.EMIP.AND.EP.GT.(EMIN+.01)) EMIP=EP                      100281
      DE=EMIP-EMIN                                                      100281
      NEMAX=INT((EMAX-EMIN)/DE+0.0001)+1                                040280
      DO 50 IE=1,NEMAX                                                  040280
      NV(IE)=0
50    A(IE)=0.
      IEMAX=0
      NBEAT=NBEA(IB)                                                    040280
C  LOOP OVER SAME BEAMS AGAIN, AVERAGING OVER THESE BEAMS AND
C  REORDERING ENERGIES
      DO 70 IB1=IB,NBED                                                 040280
      IF (NBEA(IB1).NE.NBEAT) GO TO 70                                  040280
      NEMAX=NEE(IB1)
      DO 60 IE=1,NEMAX
      IEN=INT((EE(IB1,IE)-EMIN)/DE+0.0001)+1
      IF (IEN.GT.IEMAX) IEMAX=IEN
      A(IEN)=A(IEN)+AE(IB1,IE)
60    NV(IEN)=NV(IEN)+1
      NBEA(IB1)=0                                                       040280
70    CONTINUE                                                          040280
      DO 80 IE=1,IEMAX
      EE(NBE,IE)=EMIN+FLOAT(IE-1)*DE
80    AE(NBE,IE)=A(IE)/FLOAT(NV(IE))
      NEE(NBE)=IEMAX
C  KEEP NAME OF FIRST BEAM ENCOUNTERED IN SET OF BEAMS TO BE AVERAGED
      DO 85 I=1,5
85    BENAME(I,NBE)=BENAME(I,IB)
      NBE=NBE+1
90    CONTINUE
      NBE=NBE-1
      IF (IPR.LT.2) GO TO 120
      DO 100 IB=1,NBE
      N=NEE(IB)
100   WRITE(6,110)IB,(EE(IB,IE),AE(IB,IE),IE=1,N)
110   FORMAT(48H0EXP. ENERG. AND INTENS. AFTER AVERAGING IN BEAM,1I4,/,
     150(5(1F7.2,1E13.4,3X),/))
120   RETURN
      END
C------------------------------------------------------------------------
C  SUBROUTINE VARSUM INTEGRATES OVER VARIOUS COMBINATIONS OF THE INPUT
C  FUNCTIONS (TABULATED) A1,A2,B1,B2, DEPENDING ON THE VALUE OF NF.
C  NV IS A RELATIVE SHIFT OF THE X-AXIS BETWEEN FUNCTIONS. IE1,IE2 ARE
C  THE INTEGRATION LIMITS.
C  WITH THE INTEGRAND OF THE ZANAZZI-JONA R-FACTOR A 10 TIMES
C  DENSER GRID IS USED AND FIRST INTERPOLATED TO
C
C   NF     INTEGRAND
C
C    1       A1
C    2       A1**2
C    3       ABS(A1)
C    4       ABS(A1-C*A2)
C    5       (A1-C*A2)**2
C    6       ABS(B1-C*B2)*ABS(A1-C*A2)/(ABS(A1)+EPS)
C
      SUBROUTINE VARSUM(A1,A2,B1,B2,NS1,NS2,NBD1,NBD2,IS1,IS2,IB1,IB2,
     1IE1,IE2,NV,EINCR,EPS,C,NF,S,Y)

      INCLUDE "PARAM"

      DIMENSION A1(NS1,NBD1,MNDATA),A2(NS2,NBD2,MNDATA),
     1B2(NS2,NBD2,MNDATA),Y(MNDATA),Y1(MNDATA),
     2Y2(MNDATA),Y3(MNDATA),Y4(MNDATA)
      DIMENSION B1(NS1,NBD1,MNDATA), YY(4   )

C     LEVEL 2, A1,A2,B1,B2

CVB  check whether integration limits are equal

      IF (IE1.eq.IE2) THEN
    
        S=0.

        RETURN

      END IF

CVB

      N=0

C  FOR ZANAZZI-JONA R-FACTOR INTERPOLATION ONTO 10-FOLD DENSER GRID
C  IS MADE

      IF (NF.EQ.6) GO TO 100

      DO 80 IE=IE1,IE2
      N=N+1
      IES=IE+NV
      GO TO (10,20,30,40,50,60),NF
10    Y(N)=A1(IS1,IB1,IE)
      GO TO 70
20    Y(N)=A1(IS1,IB1,IE)**2
      GO TO 70
30    Y(N)=ABS(A1(IS1,IB1,IE))
      GO TO 70
40    Y(N)=ABS(A1(IS1,IB1,IE)-C*A2(IS2,IB2,IES))
      GO TO 70
50    Y(N)=(A1(IS1,IB1,IE)-C*A2(IS2,IB2,IES))**2
      GO TO 70
60    Y(N)=ABS(B1(IS1,IB1,IE)-C*B2(IS2,IB2,IES))*
     1     ABS(A1(IS1,IB1,IE)-C*A2(IS2,IB2,IES))/
     2        (ABS(A1(IS1,IB1,IE))+EPS)
70    CONTINUE
80    CONTINUE
      CALL INTSUM(Y,1,1,1,1,EINCR,1,N,S)
      RETURN

100   DO 110 IE=IE1,IE2
      N=N+1
      IES=IE+NV
      Y(N)=FLOAT(N-1)*EINCR
      Y1(N)=A1(IS1,IB1,IE)
      Y2(N)=A2(IS2,IB2,IES)
      Y3(N)=B1(IS1,IB1,IE)
110   Y4(N)=B2(IS2,IB2,IES)
      DE=0.1*EINCR
      NN=10*(N-1)+1
      DO 120 IE=1,NN
      X=FLOAT(IE-1)*DE
      ITIL=0
      ITIH=0
      AA1=YVAL(X,Y1(1),Y,N,ITIL,ITIH)
      ITIL=0
      ITIH=0
      AA2=YVAL(X,Y2(1),Y,N,ITIL,ITIH)
      ITIL=0
      ITIH=0
      AB1=YVAL(X,Y3(1),Y,N,ITIL,ITIH)
      ITIL=0
      ITIH=0
      AB2=YVAL(X,Y4(1),Y,N,ITIL,ITIH)
120   YY(IE)=ABS(AB1-C*AB2)*ABS(AA1-C*AA2)/(ABS(AA1)+EPS)
      CALL INTSUM(YY,1,1,1,1,DE,1,NN,S)
      RETURN
      END
C------------------------------------------------------------------------
C  SUBROUTINE INTSUM INTEGRATES BY THE SIMPLE TRAPEZOID RULE (AFTER
C  ZANAZZI-JONA)
      SUBROUTINE INTSUM(Y,NS,NBD,IS,IB,H,I1,I2,S)

      INCLUDE "PARAM"

      DIMENSION Y(NS,NBD,MNDATA)
      A=0.
      S=0.
      DO 10 J=I1,I2
10    A=A+Y(IS,IB,J)
      S=A-0.5*(Y(IS,IB,I1)+Y(IS,IB,I2))
      S=S*H
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE RINTAV AVERAGES LEED INTENSITIES OVER DOMAINS, USING       281279
C  WEIGHTS DETERMINED BY THE SYMMETRY PROPERTIES (SYM) USED IN THE
C  LEED PROGRAM (THIS IS NECESSARY FOR THE CASE WHEN SYMMETRY HAS
C  HAS BEEN EXPLOITED, OTHERWISE ALL WEIGHTS ARE 1.)
      SUBROUTINE RINTAV(AT,NS,NB,NE,PQ,PQ1,KAV,SYM,LAVM,ES,NPRIN,ATAV)  281279
      DIMENSION AT(NS,NB,NE),ATAV(NS,NB,NE),KAV(NB),PQ(2,NB),PQ1(2,NB)
	  DIMENSION ES(NE)
      INTEGER SYM(NB)                                                    281279
C  FOR CDC ONLY *****************************
C     LEVEL 2, AT
      LAV=1
      LAVM=0
10    CONTINUE
      DO 20 IB=1,NB
C  FIND THE FIRST BEAM OF A GROUP TO BE AVERAGED TOGETHER
      IF (LAV.NE.KAV(IB)) GO TO 20
      LAVM=LAVM+1
C  USE NAME OF FIRST BEAM ENCOUNTERED IN THIS GROUP (FOR PRINTING)
      DO 12 K=1,2
12    PQ1(K,LAV)=PQ(K,IB)
      WT=0.0
      DO 15 JB=IB,NB
C  FIND THE OTHER BEAMS OF THE GROUP TO BE AVERAGED WITH THE FIRST
C  BEAM ENCOUNTERED ABOVE
      IF (LAV.NE.KAV(JB)) GO TO 15
      LAB=SYM(JB)
C  DETERMINE WEIGHTING FACTOR
      IF (LAB-10) 45,100,60
45    IF (LAB-2) 70,50,50
50    IF (LAB-7) 80,90,90
60    IF (LAB-12) 110,65,65
65    IF (LAB-15) 120,130,130
70    WGHT=1.0
      GO TO 200
80    WGHT=2.0
      GO TO 200
90    WGHT=4.0
      GO TO 200
100   WGHT=8.0
      GO TO 200
110   WGHT=3.0
      GO TO 200
120   WGHT=6.0
      GO TO 200
130   WGHT=12.0
200   CONTINUE
      IF (WT.LT.0.001) GO TO 14
      DO 13 IS=1,NS
      DO 13 IE=1,NE
13    ATAV(IS,LAV,IE)=ATAV(IS,LAV,IE)+WGHT*AT(IS,JB,IE)
      GO TO 145
14    DO 150 IS=1,NS
      DO 150 IE=1,NE
150   ATAV(IS,LAV,IE)=WGHT*AT(IS,JB,IE)
145   CONTINUE
      WT=WT+WGHT
15    CONTINUE
      DO 17 IS=1,NS
      DO 17 IE=1,NE
17    ATAV(IS,LAV,IE)=ATAV(IS,LAV,IE)/WT
      GO TO 25
20    CONTINUE
25    LAV=LAV+1
      IF (LAV.LE.NB) GO TO 10
      IF (NPRIN.LT.2) RETURN
C  PRINT OUT AVERAGED INTENSITIES
C      DO 40 IS=1,NS
C     WRITE(6,30)IS,(PQ1(1,IB),PQ1(2,IB),IB=1,LAVM)
C  30  FORMAT(1H ,//,19H0SURFACE STRUCTURE ,I3,//,45H BEAMS AND BEAM INTE
C     1NSITIES AFTER AVERAGING  ,5(/,22X,8(1X,2A6),/,28X,8(1X,2A6)))
C      DO 32 IE=1,NE
C 32  WRITE(6,35)ES(IE),(AT(IS,IB,IE),IB=1,LAVM)
C  35  FORMAT(10H ENERGY = ,1F7.2,5H EV  ,8E13.5,5(/,28X,8E13.5,/,22X,
C     18E13.5))
C  40  CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE MAXINT FINDS THE MAXIMUM VALUE OF MATRIX ELEMENTS,         121280
C  SKIPPING UNDESIRED GEOMETRIES
      SUBROUTINE MAXINT(A,NS,NBD,IB,NE,AM,NDATT)
      DIMENSION A(NS,NBD,NDATT)
C  FOR CDC ONLY *******************
C     LEVEL 2, A
      AM=0.
      DO 30 IS=1,NS
      DO 20 IE=1,NE
20    IF (A(IS,IB,IE).GT.AM) AM=A(IS,IB,IE)
30    CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE STRIP2 COPIES ENERGIES AND INTENSITIES FOR A GIVEN
C  GEOMETRY IS
      SUBROUTINE STRIP2(YS,NS,NBD,NE,IS,Y,NB,ES,ET,NDATT)

      INCLUDE "PARAM"

      DIMENSION YS(NS,NBD,NDATT),Y(NBD,MNDATA),NE(NBD),ES(MNDATA)
      DIMENSION ET(NBD,NDATT)
C     LEVEL 2, YS,Y,ET
      DO 10 IB=1,NB
      N=NE(IB)
      DO 10 IE=1,N
      ET(IB,IE)=ES(IE)
10    Y(IB,IE)=YS(IS,IB,IE)
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE COMNEI FINDS ENERGY INTERVAL COMMON TO EXP. AND THEORY    040280
      SUBROUTINE COMNEI(EE,NBED,NEE,ET,NBTD,NET,IBE,IBT,V0,EINCR,
     1NE1,NE2,NT1,NT2,EET)
C  FOR CDC ONLY  ***************************
C     LEVEL 2, ET,EE

      INCLUDE "PARAM"

      DIMENSION EE(NBED,MNDATA),NEE(NBED),ET(NBTD,MNDATA),NET(NBTD)
      NE=NEE(IBE)
      NT=NET(IBT)
      DE1=ET(IBT,1)+V0-EE(IBE,1)
      DE2=ET(IBT,NT)+V0-EE(IBE,NE)
      IF (DE1.LT.0.) GO TO 10
      NE1=INT((DE1/EINCR)+0.0001)+1
      NT1=1
      GO TO 20
10    NE1=1
      NT1=INT((-DE1/EINCR)+0.0001)+1
20    IF (DE2.LT.0.) GO TO 30
      NE2=NE
      NT2=NT-INT((DE2/EINCR)+0.0001)
      GO TO 40
30    NE2=NE-INT((-DE2/EINCR)+0.0001)
      NT2=NT
40    CONTINUE
      IF (NE2.GT.NE1) EET=EE(IBE,NE2)-EE(IBE,NE1)                       040280
      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE READT READS IN THEORETICAL IV-CURVES TO BE USED AS 
C  PSEUDOEXPERIMENT
C  AUTHOR: R.DOELL, 06.11.92
      SUBROUTINE READT(AE,EE,NBED,NEE,NBEA,NDATA,BENAME,IPR,MAXI,IOFF)

      INCLUDE "PARAM"

      DIMENSION AE(NBED,NDATA),EE(NBED,NDATA),NEE(NBED),NBEA(NBED),
     1FMT(19),BENAME(5,NBED),EEH(MNDATA)
C
C  MAXI     MAXIMAL INTENSITY IN ONE BEAM
C  IOFF     OFFSET FOR DATA POINT NUMBER
C  CHARD    CHARACTER-DUMMY
      REAL MAXI(NBED)
      INTEGER IOFF(NBED)
      CHARACTER*8 CHARD
C
30    FORMAT(5(25I3,/))
C  READ IN BEAM AVERAGING INFORMATION. IF NBEA(I)=NBEA(J), PSEUDO-EXP.
C  BEAMS I AND J WILL BE AVERAGED TOGETHER. THE RELATION NBEA(J).GT.
C  NBEA(I) IF J.GT.I MUST HOLD. IF NBEA(I)=0, EXP. BEAM I WILL BE SKIPPED
C  LATER ON
      READ(12,30)(NBEA(I),I=1,NBED)
      write(6,30)(NBEA(I),I=1,NBED)
C
C
  101 FORMAT(I3,19A4)
  102 FORMAT(A8)
  103 FORMAT(I4,E13.4)
  104 FORMAT(A27)
  105 FORMAT(F7.2,F7.4,4E14.5,/,20(5E14.5,/))
  106 FORMAT(8F10.3)
C
C
C BEGIN READ IN PART
      READ (12,102) CHARD
      READ (12,103) NBEAMS
      write(6,102) CHARD
      write(6,103) NBEAMS
      IF (NBED .NE. NBEAMS) THEN
          WRITE(7,'(A33)') '*********************************'
          WRITE(7,'(A33)') 'SOMETHING IS WRONG WITH NUMBER OF'
          WRITE(7,'(A33)') 'PSEUDOEXPERIMENTAL BEAMS!        '
          WRITE(7,'(A33)') 'NO AVERAGING ALLOWED!            '
          WRITE(7,'(A33)') '*********************************'
          STOP 'ERROR 1 IN PSEUDOEXPERIMENT'
      ENDIF
      DO 210 IB=1,NBED   
         READ(12,101) INN,(BENAME(I,IB),I=1,5)
         write(6,101) INN,(BENAME(I,IB),I=1,5)
  210 CONTINUE
C
C  PRESET ALL ARRAYS
      DO 220 K=1,NBED   
            MAXI(K) = 1.E-30
            NEE(K) = 0
            IOFF(K) = 0
      DO 220 J=1,NDATA    
            AE(K,J) = 0.
  220 CONTINUE
C
      DO 230 J=1,NDATA
         READ (12,105,END=250) EEH(J),DUMMY,(AE(K,J),K=1,NBED)

         DO 240 K=1,NBED
            IF (MAXI(K) .LT. AE(K,J)) MAXI(K)=AE(K,J)
            IF (AE(K,J) .GT. 1.E-30) THEN
               NEE(K) = NEE(K)+1
            ELSE IF (NEE(K) .EQ. 0) THEN
               IOFF(K) = IOFF(K)+1
            ENDIF
  240    CONTINUE
  230 CONTINUE
C
C END READ IN PART
  250 CONTINUE
C
      DO 300 K=1,NBED
      DO 300 J=1,NEE(K)
         EE(K,J)=EEH(J+IOFF(K))
         AE(K,J)=AE(K,J+IOFF(K))/MAXI(K)
  300 CONTINUE
C
      RETURN
      END
*******************************************************************************
CVB
C  Subroutine HeadDoc writes header for output file

      Subroutine HeadDoc(WHICHR,WHICHG)

      INTEGER WHICHR,WHICHG

      WRITE(4,*)

      IF (WHICHR.eq.1) THEN

        IF (WHICHG.eq.1) THEN
        
          WRITE(4,*) "R-factor minimum search using",
     &               " integer beam RPe for optimization."

        ELSE IF (WHICHG.eq.2) THEN

          WRITE(4,*) "R-factor minimum search using ",
     &               "half-order beam RPe for optimization."

        ELSE

          WRITE(4,*) "R-factor minimum search using ",
     &               "RPe for all beams for optimization."

        ENDIF

      ELSE IF (WHICHR.eq.2) THEN

        IF (WHICHG.eq.1) THEN
        
          WRITE(4,*) "R-factor minimum search using ",
     &               "integer beam R2 for optimization."

        ELSE IF (WHICHG.eq.2) THEN

          WRITE(4,*) "R-factor minimum search using ",
     &               "half-order beam R2 for optimization."

        ELSE

          WRITE(4,*) "R-factor minimum search using ",
     &               "R2 for all beams for optimization."

        ENDIF

      ENDIF

      WRITE (4,*)

      RETURN
    
      END

**************************************************************************
C  Subroutine ReadFile reads in a specified delta amplitude file

      Subroutine ReadFile(CNTFIL,NDOM,IDOM,IPLACE,IFILE,NPLACES,NFILES,
     +                    NT0,NNATOMS,NSTEP,NDATT,IFORM,EMIN,EMAX,
     +                    INFILE,THETA,FI,RAR1,RAR2,PQFEX,CUNDISP,
     +                    CDISP,AID,EMK,EMK0,ESMK,VPI,VO,VV,
     +                    XISTMK,DELMK,DATTNO,XIST,DELWV,CUNDVB,CDVB)

      include "PARAM"

C  all these are dimension sizes

      INTEGER NPLACES(NDOM),NFILES,NT0,NNATOMS,NSTEP,NDATT

C  CNTFIL is file number
C  IPLACE is current atomic site
C  IFILE is file under consideration for that site

      INTEGER CNTFIL,IPLACE,IFILE,IDOM

C  IFORM is info on format of input data

      INTEGER IFORM
      DIMENSION IFORM(NDOM,MNPLACES,NFILES)

C  INFILE is array containing file names

      CHARACTER*15 INFILE(NDOM,MNPLACES,NFILES)

C  angular, geometrical and general information on current file

      REAL THETA,FI,RAR1(2),RAR2(2)

C  beamlist and displacement information
C  Note AID is not read correctly for unformatted case!

      REAL PQFEX,CUNDISP,CDISP,AID
      DIMENSION PQFEX(2,NT0)
      DIMENSION CUNDISP(NNATOMS,3,NDOM,MNPLACES,NFILES)
      DIMENSION CDISP(NSTEP,NNATOMS,3,NDOM,MNPLACES,NFILES)
      DIMENSION AID(NSTEP,NDOM,MNPLACES,NFILES)

C  energy, optical potential, inner potential
C  EMK0 to check consistency of read energy values EMK
C  ESMK to save energy values in eV

      REAL EMK,EMK0,ESMK,VPI,VO,VV
      DIMENSION EMK(NDATT),EMK0(NDATT),ESMK(NDATT),VPI(NDATT),
     +          VO(NDATT),VV(NDATT)

C  amplitudes and delta amplitudes

      COMPLEX XISTMK,DELMK
      DIMENSION XISTMK(NDOM,NDATT,NT0)
      DIMENSION DELMK(NDATT,NSTEP,NT0,NDOM,MNPLACES,NFILES)

C  limits of energy range to be considered

      REAL EMIN,EMAX

C  number of energies in delta amp files

      INTEGER DATTNO

C  XIST, DELWV, CUNDVB, CDVB are auxiliary arrays only 
C  to provide variable boundaries for unformatted readin

      COMPLEX XIST,DELWV
      DIMENSION XIST(NT0),DELWV(NSTEP,NT0)

      REAL CUNDVB,CDVB
      DIMENSION CUNDVB(NNATOMS,3),CDVB(NSTEP,NNATOMS,3)

C  local variables

      REAL TTHETA,TFI,TRAR1(2),TRAR2(2)

C  NATOMS and NCSTEP need not be known for any of the remaining program
C  - therefore local, too

      INTEGER NATOMS,NCSTEP

      INTEGER INT0

      REAL CHECK

      INTEGER IDATT

      REAL EeV

      INTEGER EndFlag

C  open current delta amplitude file

          write(6,*) 'Datei ', CNTFIL, ':',INFILE(IDOM,IPLACE,IFILE)

          IF (IFORM(IDOM,IPLACE,IFILE).eq.0) THEN

            OPEN(CNTFIL, FILE = INFILE(IDOM,IPLACE,IFILE),
     +           FORM = 'UNFORMATTED', STATUS = 'OLD')

          ELSE IF (IFORM(IDOM,IPLACE,IFILE).eq.1) THEN

            OPEN(CNTFIL, FILE = INFILE(IDOM,IPLACE,IFILE),
     +           FORM = 'FORMATTED', STATUS = 'OLD')

          ELSE

            write(8,*) 'Illegal format for file ',
     +           INFILE(IDOM,IPLACE,IFILE)
            STOP

          END IF

          write(6,*) 'Datei ', CNTFIL, 'geoeffnet.'

C  Read in headers of all files

          IF (IFORM(IDOM,IPLACE,IFILE).eq.0) THEN

            READ(CNTFIL) TTHETA,TFI,TRAR1,TRAR2,INT0,NATOMS,NCSTEP

          ELSE
            READ(CNTFIL,'(6E13.7)') TTHETA,TFI,(TRAR1(I),I=1,2),
     +                              (TRAR2(I),I=1,2)
            READ(CNTFIL,'(3I3)') INT0,NATOMS,NCSTEP

          END IF

C  Now store data or check whether files are ok          

          IF (CNTFIL.eq.11) THEN

C  first file

            THETA = TTHETA
            FI = TFI

C  NT0 is a known quantity here so it must be kept!

            IF (INT0.ne.NT0) THEN 

              write (8,*) 'Number of theoretical beams incorrect in',
     +                    INFILE(IDOM,IPLACE,IFILE),'!'
              STOP

            END IF

            DO 2224 I=1,2

              RAR1(I) = TRAR1(I) 
              RAR2(I) = TRAR2(I) 

 2224       CONTINUE

          ELSE

            CHECK = ABS(THETA - TTHETA)
            CHECK = CHECK + ABS(FI - TFI)

            DO 2225 I = 1,2

              CHECK = CHECK + ABS(RAR1(I) - TRAR1(I))
              CHECK = CHECK + ABS(RAR2(I) - TRAR2(I))

 2225       CONTINUE   

            IF ((CHECK.gt.1.0e-09).or.(INT0.ne.NT0)) THEN

              WRITE(8,*) 'Improper file combination!'
              STOP

            END IF

          END IF

C  read in remaining part of file header

          IF (IFORM(IDOM,IPLACE,IFILE).eq.0) THEN

C  need variable dimensions -> must be done in subroutine
C  note AID is not read correctly!!

            CALL ReadCD(CNTFIL,NDOM,IDOM,PQFEX,CUNDISP,CDISP,AID,NT0,
     +                  NNATOMS,NPLACES,NFILES,NSTEP,IPLACE,IFILE,
     +                  NATOMS,NCSTEP,CUNDVB,CDVB)

          ELSE

            READ(CNTFIL,'(10F10.5)',ERR=100) 
     +      ((PQFEX(I,J), I=1,2), J=1,NT0)

            GO TO 110

 100        write(8,*) "Error reading v1.2 style output beam list",
     +                 " from delta amplitude file. Trying pre-v1.2"
            write(8,*) "format instead ..."
            
            READ(CNTFIL,'(10F7.4)') 
     +      ((PQFEX(I,J), I=1,2), J=1,NT0)

            write(8,*) "Read may have been successful. Check",
     +                 " carefully, anyhow!"

 110        CONTINUE

            READ(CNTFIL,'(10F7.4)') ((CUNDISP(I,J,IDOM,IPLACE,IFILE),
     +                                 J=1,3), I=1,NATOMS)

            READ(CNTFIL,'(10F7.4)') (((CDISP(I,J,K,IDOM,IPLACE,IFILE),
     +          K=1,3), J=1,NATOMS),I=1,NCSTEP)

            READ(CNTFIL,'(10F7.4)') (AID(I,IDOM,IPLACE,IFILE),
     +          I=1,NCSTEP)

          END IF

          write(6,*) 'Header ', CNTFIL, ' komplett eingelesen.'

C  Now read in delta amplitudes in energy range EMIN, EMAX from file
C  in repeat-until-EOF loop 2226 - IDATT counts energy
C  Unfortunately, array DELWV must be kept for proper reading of 
C  unformatted files

          IDATT=0

 2226     CONTINUE

            IDATT = IDATT + 1

            IF (IFORM(IDOM,IPLACE,IFILE).eq.0) THEN

              call ReadDels(CNTFIL,NDOM,IDOM,EMK(IDATT),VPI(IDATT),
     +                VO(IDATT),VV(IDATT),
     +                XISTMK,DELMK,XIST,DELWV,IPLACE,IFILE,NDATT,NT0,
     +                NSTEP,NPLACES,NFILES,NCSTEP,IDATT,EndFlag)

              IF (EndFlag.eq.1) THEN

                GO TO 2227

              END IF

            ELSE

              READ(CNTFIL,'(6E13.7)',ERR=2227,END=2227) 
     +             EMK(IDATT),VPI(IDATT),VO(IDATT),VV(IDATT)

              READ(CNTFIL,'(6E13.7)') (XISTMK(IDOM,IDATT,I), I=1,NT0)

              READ(CNTFIL,'(6E13.7)') 
     +             ((DELMK(IDATT,I,J,IDOM,IPLACE,IFILE), 
     +             J=1,NT0), I=1,NCSTEP)

            END IF

C  check validity of current energy value
C  if energy is within range save it in array ESMK

            EeV = (EMK(IDATT)-VV(IDATT))*27.21

            IF (EeV.lt.(EMIN - 1.E-4)) THEN

              IDATT = IDATT -1

            ELSE IF (EeV.gt.(EMAX + 1.E-4)) THEN

              write(6,*) 'EeV = ', EeV, '     EMAX + 1.E-4 = ', EMAX + 1.E-4
              write(6,*) 'EeV.gt.(EMAX + 1.E-4)'
              GO TO 2227

            ELSE IF (IDATT.gt.NDATT) THEN

              write(8,*) 'Insufficient dimension MNDATT!!'

              STOP

            ELSE

              ESMK(IDATT) = EeV

C  check validity of read energy

              IF (CNTFIL.eq.11) THEN

                EMK0(IDATT) = EMK(IDATT) - VV(IDATT)

              ELSE

                CHECK = ABS(EMK(IDATT)-VV(IDATT)-EMK0(IDATT))

                IF (CHECK.gt.1.E-4) THEN

                  write(8,*) 'Illegal energy read in input file',
     +            CNTFIL

                  STOP

                END IF

              END IF

            END IF              

C  Continue energy loop until last energy has been read in

            GO TO 2226

C  jump here to leave energy loop

 2227     CONTINUE
           
C  Save total number of energies in DATTNO
C  last valid step (the one before exiting loop) must be saved

          IDATT = IDATT-1

          DATTNO = IDATT

C  close current file and proceed to the following one

          CLOSE(CNTFIL)

          write(6,*) 'Alles ', CNTFIL, ' eingelesen.'

          RETURN

          END

****************************************************************************

C  Subroutine ReadCD performs unformatted readin of CDISP etc.

      Subroutine ReadCD(CNTFIL,NDOM,IDOM,PQFEX,CUNDISP,CDISP,AID,NT0,
     +                  NNATOMS,NPLACES,NFILES,NSTEP,IPLACE,IFILE,
     +                  NATOMS,NCSTEP,CUNDVB,CDVB)

      include "PARAM"

      INTEGER CNTFIL,NT0,NNATOMS,NPLACES(NDOM),NFILES,NSTEP
      INTEGER NATOMS,NCSTEP,IDOM

      REAL PQFEX(NT0),AID(NCSTEP)

      REAL CUNDISP(NNATOMS,3,NDOM,MNPLACES,NFILES)         
      REAL CDISP(NSTEP,NNATOMS,3,NDOM,MNPLACES,NFILES)

      REAL CUNDVB(NATOMS,3),CDVB(NCSTEP,NATOMS,3)

C  now read dummy arrays

      READ(CNTFIL) PQFEX,CUNDVB,CDVB,AID

C  and save to useful arrays

      DO 2230 I=1,3,1
        DO 2230 J=1,NATOMS,1

          CUNDISP(J,I,IDOM,IPLACE,IFILE) = CUNDVB(J,I)

          DO 2230 K=1,NCSTEP,1
            CDISP(K,J,I,IDOM,IPLACE,IFILE) = CDVB(K,J,I)
 2230 CONTINUE

      RETURN
      END
*************************************************************************

C  Subroutine ReadDels reads in unformatted version of delta amplitudes
C  to allow for variable arrays cause f77 cannot handle anything intelligently.

      Subroutine ReadDels(CNTFIL,NDOM,IDOM,E,VPI,VO,VV,XISTMK,DELMK,
     +                    XIST,DELWV,IPLACE,IFILE,NDATT,NT0,
     +                    NSTEP,NPLACES,NFILES,NCSTEP,IDATT,EndFlag)

      include "PARAM"

      INTEGER CNTFIL,EndFlag

      INTEGER NDATT,NSTEP,NT0,NPLACES(NDOM),NFILES,IDOM

      REAL E,VPI,VO,VV

      COMPLEX XISTMK,DELMK
      DIMENSION XISTMK(NDOM,NDATT,NT0)
      DIMENSION DELMK(NDATT,NSTEP,NT0,NDOM,MNPLACES,NFILES)

      COMPLEX XIST(NT0),DELWV(NCSTEP,NT0)


      READ(CNTFIL,ERR=2240,END=2240) E,VPI,VO,VV,XIST,DELWV

      DO 2250 I=1,NT0
        
        XISTMK(IDOM,IDATT,I) = XIST(I)

        DO 2250 J=1,NCSTEP

          DELMK(IDATT,J,I,IDOM,IPLACE,IFILE) = DELWV(J,I)

 2250 CONTINUE

      EndFlag = 0

      RETURN

C  continue here if file ended

 2240 CONTINUE

      EndFlag = 1

      RETURN

      END

***********************************************************************************

C  Subroutine GetInt produces output intensities for given files and positions
C  IFNUM in those files. 'Concentration' of different files can be varied here.

      Subroutine GetInt(NDOM,IDOM,NPLACES,NFILES,NSTEP,NCONCS,NNATOMS,
     +                  NDATT,NT0,NSURF,IFNUM,NFIL,CONC,CDISP,NPARC,
     +                  DATTNO,EMK,VV,VO,VPI,THETA,FI,PQFEX,RAR1,RAR2,
     +                  XISTMK,DELMK,ATSAS)

      include "PARAM"

C  cast (in approximate order of appearance)

C  Quantities for dimensions

      INTEGER NPLACES(NDOM),NFILES,NSTEP,NCONCS,NNATOMS,NDATT,NT0

C  surface, current delta amp set, no of files per place

      INTEGER NSURF, IFNUM, NFIL
      DIMENSION NSURF(NDOM,MNPLACES),IFNUM(NDOM,MNPLACES,NFILES)
      DIMENSION NFIL(NDOM,MNPLACES)

C  concentration, displacement, current concentration parameter

      REAL CONC,CDISP
      DIMENSION CONC(NDOM,NCONCS,MNPLACES,NFILES)
      DIMENSION CDISP(NSTEP,NNATOMS,3,NDOM,MNPLACES,NFILES)

      INTEGER NPARC
      DIMENSION NPARC(NDOM,MNPLACES)

C  no of energies to be considered, energy values, inner potential

      INTEGER DATTNO

      REAL EMK,VV,VO,VPI
      DIMENSION EMK(NDATT),VV(NDATT),VO(NDATT),VPI(NDATT)

C  angular info

      REAL THETA,FI

C  beam list, reciprocal lattice

      REAL PQFEX(2,NT0),RAR1(2),RAR2(2)

C  undisturbed amplitudes, delta amplitudes

      COMPLEX XISTMK,DELMK
      DIMENSION XISTMK(NDOM,NDATT,NT0)
      DIMENSION DELMK(NDATT,NSTEP,NT0,NDOM,MNPLACES,NFILES)

C  resulting intensities ATSAS for each domain

      REAL ATSAS
      DIMENSION ATSAS(NDOM,NT0,NDATT)


C  local variables

C  outermost displacement 

      REAL CXDISP,XDISP

      INTEGER IAct

      REAL E

C  variables for calculation of 'propagator' between surface, vacuum,
C  and angular correction

      REAL AK,C,BK2,BK3,A,RPRE

      REAL AK2,AK3,APERP

      COMPLEX AKZ

      COMPLEX BKZ,PRE

C  delta amplitude being calculated, variables for intensity calc

      COMPLEX DelAct

      REAL AmpAbs


C  Calculate outermost displacement CXDISP of any surface atom first to later 
C  correct propagator between top layer and onset of damping above (ASE).
C  Will only work properly for flat surface in precalc since only relative
C  displacements can be considered! To avoid this, quantity CUNDISP must be used 
C  in computation of delta amplitudes.

C  NATOMS can no longer be used as it once was in superpos!

      CXDISP = 1000.

      DO 3323 IPLACE = 1, NPLACES(IDOM)

        IF (NSURF(IDOM, IPLACE).eq.1) THEN

          XDISP = 0.

          DO 3325 IFILE = 1, NFIL(IDOM, IPLACE)

            IAct = IFNUM(IDOM, IPLACE, IFILE)

            XDISP = XDISP + 
     .              CONC(IDOM, NPARC(IDOM, IPLACE), IPLACE, IFILE) *
     .              CDISP(IAct, 1, 1, IDOM, IPLACE, IFILE)

 3325     CONTINUE

          IF (XDISP .lt. CXDISP) THEN

            CXDISP = XDISP

          END IF

        END IF

 3323 CONTINUE

C  do not change onset of damping if flag surface was never set

      IF (CXDISP.eq.1000.) THEN

        CXDISP = 0.

      END IF

C  Begin loop over energies EMK

      DO 3333 IDATT = 1,DATTNO

        E = EMK(IDATT)

C  Calculate incident wave vector

        AK  = SQRT(AMAX1(2.0*E-2.0*VV(IDATT),0.))

        C   = AK*COS(THETA)
        BK2 = AK*SIN(THETA)*COS(FI)
        BK3 = AK*SIN(THETA)*SIN(FI)

C  AND PREPARE CALCULATION OF OVERLAYER PROPAGATOR
C  (phase shift and damping for beam travelling between top layer
C   and vacuum)

        BKZ=CMPLX(2.0*E-2.0*VO(IDATT)-BK2*BK2-BK3*BK3,-2.0*VPI(IDATT))
        BKZ=CSQRT(BKZ)

C  CALCULATE WAVEVECTORS OF OUTPUTBEAMS, PREPARE CALCULATION OF
C  OVERLAYER PROPAGATOR

        DO 550 IBEAM=1,NT0

          AK2 = BK2 + PQFEX(1,IBEAM)*RAR1(1) 
     +               + PQFEX(2,IBEAM)*RAR2(1)


          AK3 = BK3 + PQFEX(1,IBEAM)*RAR1(2)
     +               + PQFEX(2,IBEAM)*RAR2(2)

          AK = 2.0*E - AK2*AK2 - AK3*AK3

          AKZ = CMPLX(AK-2.0*VO(IDATT),-2.0*VPI(IDATT))
          AKZ = CSQRT(AKZ)

          APERP =  AK-2.0*VV(IDATT)

C  now calculate actual amplitude for each beam

          DelAct = XISTMK(IDOM, IDATT, IBEAM)

C  Add Delta amplitudes for current parameter set IFNUM(IPLACE,IFILE)
C  for each place and each file multiplying the atomic concentrations
C  to each amplitude.

          DO 545 IPLACE = 1, NPLACES(IDOM)

            DO 540 IFILE = 1, NFIL(IDOM, IPLACE)

              IAct = IFNUM(IDOM, IPLACE, IFILE)

              DelAct = DelAct + 
     .                 DELMK(IDATT, IAct, IBEAM, IDOM, IPLACE, IFILE)*
     .                 CONC(IDOM, NPARC(IDOM, IPLACE), IPLACE, IFILE)

 540        CONTINUE

 545      CONTINUE

C  now calculate the resulting intensity for emerging beams
C  APERP determines whether or not a beam can emerge....

          IF (APERP.ge.0.) THEN

            A = sqrt(APERP)

C  prefactor for ASE-correction of intensity

            PRE = (BKZ + AKZ) * CXDISP
            PRE = CEXP(CMPLX(0.,-1.)*PRE)

            AmpAbs = CABS(DelAct)

C  Now get intensity ATSMK(energy,beam)
C  A/C is simply angular correction for off-normal incidence.
C  cutoff for maximum intensity - caution if occurs!!

            IF (AmpAbs.gt.1.E+10) THEN

              ATSAS(IDOM, IBEAM, IDATT) = 1.E+20

            ELSE

              RPRE = CABS(PRE)
              ATSAS(IDOM, IBEAM, IDATT) = 
     .             AmpAbs*AmpAbs * RPRE*RPRE * A/C
            END IF

          ELSE

C  non-emerging beam

            ATSAS(IDOM, IBEAM, IDATT) = 0.

          END IF

c  next beam

 550    CONTINUE

c next energy

 3333 CONTINUE  

      RETURN

      END

***************************************************************************
C  Subroutine GetWid calculates the width of the gaussian distribution
C  that is used for the next generation for each parameter from results
C  of the current generation.

C  The width is determined for each 'individual' of current generation by
C  comparison with the previous results, yielding an average sensitivity
C  of each parameter wrt the r-factor which is represented in WIDT. 
C  This average is taken by determining an average change in r-factor for
C  each parameter, weighting the change in each parameter with the sum
C  of changes in all parameters, for each 'individual'. 

      Subroutine GetWid(WIDT,PARIND,PAROLD,RPEIND,RPEOLD,NormInd,
     +                  NPS,NPRMK)

C  for dimension statements etc.

      INTEGER NPS,NPRMK

C  parameter step values and r-factors for previous and current generation

      INTEGER PARIND(NPRMK,NPS),PAROLD(NPRMK,NPS)
      REAL RPEIND(NPS),RPEOLD(NPS)

C  width value to be calculated

      REAL WIDT(NPRMK)

C  local variable, but with var. dimension
C  NormInd is normalises parameter sensitivity for each individual

      INTEGER NormInd(NPS)

C  other local variables

C  DelPar is change of current parameter for current individual
C  SumSteps is sum of parameter displacements of current parameter
C  over all individuals

      INTEGER SumSteps,DelPar

C  Sens is sensitivity of current parameter in current population

      REAL Sens
      REAL SumSens


C  NormInd is sum over displacements between this and last generation for
C  each individual 

      DO 40 IPop = 1,NPS,1

        NormInd(IPOP) = 0

        DO 50 IPARAM=1,NPRMK,1

          NormInd(IPop) = NormInd(IPop)
     +                  + ABS(PARIND(IPARAM,IPop)-PAROLD(IPARAM,IPop))

 50     CONTINUE

 40   CONTINUE

C  Now calculate width for each parameter

      DO 100 IPARAM = 1,NPRMK,1

C  Norm is sum of displacements between this and last generation for each
C  parameter

C  SumSens sum over average sensitivity for current parameter in each
C  'individual'

        SumSteps = 0
        SumSens = 0.

        DO 200 IPop=1,NPS,1

          DelPar =  ABS(PARIND(IPARAM,IPop)-PAROLD(IPARAM,IPop))

          IF (DelPar.gt.0) THEN

            SumSteps = SumSteps + DelPar

            Sens = ( ABS(RPEIND(IPop) - RPEOLD(Ipop))) 
            Sens = Sens * REAL(DelPar) / REAL(NormInd(IPop)) 


            SumSens = SumSens + Sens

          END IF

 200    CONTINUE

C  Now normalise SumSens using SumSteps

        IF ((SumSens.gt.1.e-04).and.(SumSteps.gt.0)) THEN

          WIDT(IPARAM) = SumSens/REAL(SumSteps)

        ELSE

          WIDT(IPARAM) = 1.

        END IF

 100  CONTINUE

      RETURN

      END

**************************************************************************
C  Subroutine CheckVal checks consistency of data read from the various input
C  files.

      Subroutine CheckVal(NBTD,NBED,PNUM,NDOM,NPLACES,VI,VPI,V0RR,VV,
     +                    VARST,NDATT,DATTNO)

      INCLUDE "PARAM"
      INCLUDE "GLOBAL"

      INTEGER NBTD,NBED,PNUM,NPLACES(NDOM),IDOM,VARST,DATTNO,NDATT
      REAL VI,VPI,V0RR,VV
      DIMENSION VPI(NDATT),VV(NDATT)
      DIMENSION VARST(MNPRMK)

      REAL VPIAV

      IF (NBTD.ne.MNBTD) THEN

        write(8,*) 'NBTD in WEXPEL or MNBTD wrong!'
        STOP

      ELSE IF (NBED.ne.MNBED) THEN

        write(8,*) 'NBED in WEXPEL or MNBED wrong!'
        STOP

      END IF

      IF (PNUM.NE.MNPRMK) THEN

         write(8,'("PNUM:",I3,"MNPRMK:",I3)') PNUM, MNPRMK

         write(8,*) "par.no    varst"
         DO 50 IPARAM=1,PNUM,1
           write(8,*) IPARAM, VARST(IPARAM)
 50      CONTINUE

         WRITE(8,*) 'Parameter MNPRMK wrong!'
         STOP
      ENDIF

      if (NDOM .ne. MNDOM) then
        write(8,*) 'Parameter MNDOM wrong!'
        stop
      end if

      do 100 IDOM = 1, NDOM
        if (NPLACES(IDOM) .gt. MNPLACES) then
           write(8,'("NPLACES:",I3,"MNPLACES:",I3)') NPLACES,MNPLACES
           WRITE(8,*) 'Parameter MNPLACES wrong !'
           STOP
        end if
  100 continue

C  check consistency of inner potential values from READRF and ReadFile

      VPIAV = 0

      do IENER = 1,DATTNO
        VPIAV=VPIAV+VPI(IENER)
      enddo

      VPIAV = VPIAV/REAL(DATTNO)

      IF ((ABS(VI)-ABS(VPIAV)*HARTREE).gt.1.0e-04) THEN

        write(8,*)
     +  "Average optical potential value in rf.info is incorrect:"
        write(8,*)
     +  "rf.info: ",VI," eV, true average: ",ABS(VPIAV)*HARTREE," eV."
        STOP

      END IF

c      IF ((ABS(V0RR)-ABS(VV)*HARTREE).gt.1.0e-04) THEN
c
c        write(8,*) 'Inner potential value in WEXPEL file is incorrect!'
c        STOP
c
c      END IF

      RETURN

      END

********************************************************************************
C Subroutine GetDependency initializes the PARDEP array. For each parameter, if 
C it should always be equivalent to another parameter via the "Atom number"
C FILREL, the PARDEP(IPARAM) will be set to the index of that parameter. For 
C all other parameters, PARDEP will be 0, and only those parameters need to be 
C calculated by SEA_RCD.
C Added 2020-12 by Florian Kraushofer; largely a copy of GetGrid that stores 
C the information instead of having to run in every iteration of the optimization
C loop.

      Subroutine GetDependency(NDOM,NPLACES,NFILES,NPRMK,NPRAS,NPS,
     +                         NFIL,PARTYP,FILREL,PARDEP)

      include "PARAM"

C  Dimension sizes

      INTEGER NPLACES(NDOM),NFILES,NPRMK,NPS
      integer NPRAS(NDOM)

C  Global variables

C  NFIL is number of files for each place
C  PARTYP is number of parameters in current file
C  FILREL can force different atoms to be treated equally ("atom number")
C  PARDEP stores whether a parameter should be set to the same value as another one
      
      INTEGER   NFIL,FILREL
      DIMENSION NFIL(NDOM,MNPLACES),FILREL(NDOM,MNPLACES)
      INTEGER   PARTYP
      DIMENSION PARTYP(NDOM,MNPLACES,NFILES)
      INTEGER   PARDEP
      DIMENSION PARDEP(NPRMK)

C  local variables

C  CNT1 counts parameters that have already been processed in previous files
C  CNT2 has the same purpose when atom number is considered (for second atom)
C  CNTPAR is another counter for parameters within a place
C  OFFSET is used to skip the parameters PARIND of the first (LDOM - 1) domains;
C  so OFFSET is equal to 0, if LDOM = 1, and equal to NPRAS(1), if LDOM = 2,
C  and so on.
C  KDOM and LDOM are the domain number in loops over domains

      INTEGER CNT1,CNT2,CNTPAR,OFFSET
      INTEGER KDOM,LDOM

C  Set 0 for all
      PARDEP = 0

C  start loop to do this for each domain
      do 201 LDOM = 1, NDOM
      
C --- first determine the value of OFFSET corresponding to the LDOM-th domain
      OFFSET = 0
      do 200 KDOM = 1, LDOM - 1
        OFFSET = OFFSET + NPRAS(KDOM)
 200  continue

C  Check dependences

      CNT1 = OFFSET
      DO 970 IPLACE1 = 1, NPLACES(LDOM) - 1

C  Set CNT2 to CNT1, then start to compare places after IPLACE1 to IPLACE1

        CNT2 = CNT1
        DO 971 IPLACE2 = IPLACE1 + 1, NPLACES(LDOM)

C  increment CNT2 for last place (which is IPLACE1 in first run)

          DO 972 IFILE = 1, NFIL(LDOM, IPLACE2 - 1)

            CNT2 = CNT2 + PARTYP(LDOM, IPLACE2 - 1, IFILE)

 972      CONTINUE

C  never forget concentration parameter

          CNT2 = CNT2 + 1

          IF (FILREL(LDOM, IPLACE2) .eq. FILREL(LDOM, IPLACE1)) THEN

C  store dependence

            CNTPAR = 0
            DO 973 IFILE = 1, NFIL(LDOM, IPLACE2)

              DO 974 IPARAM = 1, PARTYP(LDOM, IPLACE2, IFILE)

                PARDEP(CNT2+CNTPAR+IPARAM) = CNT1+CNTPAR+IPARAM

 974          CONTINUE

              CNTPAR = CNTPAR + PARTYP(LDOM, IPLACE2, IFILE)

 973        CONTINUE

C  never forget concentration parameter

            PARDEP(CNT2+CNTPAR+1) = CNT1+CNTPAR+1

          END IF

 971    CONTINUE

C  now increment CNT1 for next place IPLACE1

        DO 975 IFILE = 1, NFIL(LDOM, IPLACE1)

          CNT1 = CNT1 + PARTYP(LDOM, IPLACE1, IFILE)

 975    CONTINUE

C  never forget concentration parameter

        CNT1 = CNT1 + 1

 970  CONTINUE

C  domain loop ends

 201  continue

      RETURN
   
      END

********************************************************************************
C  Compute current grid point in deltaamplitude-file from parameter grip points
C  for each site IFNUM and store away concentration step no. NPARC
C  grid points are computed from parameter values PARIND file by file
C  "Atom number" FILREL is taken into account
C  Handling of the PARIND is a messy affair because of it simply counting
C  parameters while IFNUM refers to a specific file. The structure of the
C  loops below is tedious and can certainly be improved using a structogram.:)

      Subroutine GetGrid(NDOM,NPLACES,NFILES,NPRMK,NPRAS,NPS,IDOM,IPOP,
     +                   NFIL,IFNUM,PARTYP,PARIND,VARST,NPARC,FILREL,
     +                   PARDEP)

      include "PARAM"

C  Dimension sizes

      INTEGER NPLACES(NDOM),NFILES,NPRMK,NPS
      integer NPRAS(NDOM)

C  Global variables

C  IPOP is current population number
C  NFIL is number of files for each place
C  IFNUM is delta amp number in delta amp file to be computed
C  PARTYP is number of parameters in current file
C  PARIND are current parameter values as determined by Sea_RCD
C  VARST is number of grid points for each parameter
C  NPARC is current concentration step number
C  FILREL can force different atoms to be treated equally ("atom number")
C  PARDEP stores whether a parameter should be set to the same value as 
C  another one (see subroutine GetDependency)

      INTEGER   IPOP
      
      INTEGER   NFIL
      DIMENSION NFIL(NDOM,MNPLACES)
      INTEGER   IFNUM,PARTYP
      DIMENSION IFNUM(NDOM,MNPLACES,NFILES),
     .          PARTYP(NDOM,MNPLACES,NFILES)
      INTEGER   PARIND,VARST,PARDEP
      DIMENSION PARIND(NPRMK,NPS),VARST(NPRMK),PARDEP(NPRMK)
      INTEGER   NPARC,FILREL
      DIMENSION NPARC(NDOM,MNPLACES),FILREL(NDOM,MNPLACES)

C  local variables

C  CNT1 counts parameters that have already been processed in previous files
C  CNT2 has the same purpose when atom number is considered (for second atom)
C  CNTPAR is another counter for parameters within a place
C  FACT is auxiliary variable
C  OFFSET is used to skip the parameters PARIND of the first (IDOM - 1) domains;
C  so OFFSET is equal to 0, if IDOM = 1, and equal to NPRAS(1), if IDOM = 2,
C  and so on.

      INTEGER CNT1,CNT2,CNTPAR,FACT,OFFSET

      INTEGER KDOM

C --- first determine the value of OFFSET corresponding to the IDOM-th domain
      OFFSET = 0
      do 200 KDOM = 1, IDOM - 1
        OFFSET = OFFSET + NPRAS(KDOM)
 200  continue

C  Take independence of atoms into account

      CNT1 = OFFSET
      DO 970 IPARAM = 1, NPRMK
        IF(PARDEP(IPARAM).ne.0) THEN
          PARIND(IPARAM,IPOP) = PARIND(PARDEP(IPARAM),IPOP)
        END IF
 970  CONTINUE

C  Begin computation of IFNUM, NPARC

      CNT1=OFFSET

      DO 981 IPLACE = 1, NPLACES(IDOM), 1

        DO 982 IFILE = 1, NFIL(IDOM, IPLACE)

          IFNUM(IDOM, IPLACE, IFILE) = 1
          FACT = 1

          DO 983 IPARAM = PARTYP(IDOM, IPLACE, IFILE), 1, -1

            IFNUM(IDOM, IPLACE, IFILE) = IFNUM(IDOM, IPLACE, IFILE)
     +           + (FACT * (PARIND(CNT1+IPARAM,IPOP) - 1))

            FACT=FACT*VARST(CNT1+IPARAM)

 983     CONTINUE

          CNT1 = CNT1 + PARTYP(IDOM, IPLACE, IFILE)

 982   CONTINUE

C  store away the concentration parameter for each place
C  the last parameter parind belonging to each place contains
C  the concentration step number

        CNT1=CNT1+1

        NPARC(IDOM, IPLACE) = PARIND(CNT1,IPOP)

 981  CONTINUE

      RETURN
   
      END

********************************************************************************
C Die Subroutine Mischung vollzieht die eigentliche Domaenenmischung.
C Die Intensitaeten der einzelnen Domaenen stehen in ATSAS, das Resultat
C wird nach ATSMK geschrieben. x ist der aus PARIND und DMISCH errechnete
C explizite Intensitaetsanteil der aktuellen (= IDOM-ten) Domaene.

CVB Mischung has been changed to include the normalisation of concentration
C steps here instead of in subroutine korrekt (now obsolete).

      subroutine Mischung(NDOM, NDATT, NT0, NPRMK, NPS, PARIND,
     +                    PMISCH,DMISCH,ATSAS, ATSMK, DATTNO, IPOP,
     +                    VARST)

      integer NDOM, NDATT, NT0,NPRMK,NPS
      integer DATTNO,IPOP,VARST

      integer PARIND,PMISCH
      DIMENSION PARIND(NPRMK, NPS),PMISCH(NDOM,NPS),VARST(NPRMK)

      real ATSAS(NDOM, NT0, NDATT), ATSMK(NT0, NDATT)
      real DMISCH

C  local variables

      integer IDOM, IDATT, IBEAM

C  PARSUM will contain the sum of PARIND values for domain fractions
C         to be used as a normalisation constant
C  NDSTEP will be the number of fraction steps for each domain,
C         counting from 0 upwards
C  NDSL1  is sum of all domain parameters PMISCH except for the last
C  Norm   will be the normalisation factor from PARIND to PMISCH
C  area   will be area fraction of current domain in percent

      integer PARSUM,NDSTEP, NDSL1
      real    NORM

C initialise step number (identical for all domains!)

      NDSTEP = VARST(NPRMK)-1

C  initialise intensity ATSMK first

      do 110 IBEAM = 1, NT0
        do 100 IDATT = 1, DATTNO
          ATSMK(IBEAM, IDATT) = 0.0
 100    continue
 110  continue

C  now compute actual area fractions (normalised to 1) 

      PARSUM = 0

      DO 200 IDOM=1,NDOM

        PMISCH(IDOM,IPOP) = PARIND(NPRMK-NDOM+IDOM,IPOP) - 1
        PARSUM = PARSUM+PMISCH(IDOM,IPOP)

 200  CONTINUE

C  wenn nicht alle Parameter 0 sind, kann man ordentlich normieren,...

      if (PARSUM.gt.0) then

        Norm = REAL(NDSTEP)/REAL(PARSUM)

      else

C  ... ansonsten wird Norm auf Dummy-Wert gesetzt; mit dieser 
C  Konvention wird im Falle aller Parameterwerte auf 0 die letzte
C  Domaene spaeter automatisch auf 100 Prozent Flaechenanteil gesetzt.

        Norm = 1.

      end if

      NDSL1 = 0

      DO 210 IDOM=1,NDOM-1

        PMISCH(IDOM,IPOP)=INT(REAL(PMISCH(IDOM,IPOP))*Norm)
        NDSL1 = NDSL1 + PMISCH(IDOM,IPOP)

 210  CONTINUE

      PMISCH(NDOM,IPOP) = NDSTEP - NDSL1

C  jetzt sind PMISCH korrekt normierte Parameterstuetzstellen,
C  zaehlend ab 0, nicht ab 1 (also anders als PARIND)

      do 140 IDOM = 1, NDOM

        area = PMISCH(IDOM,IPOP) * DMISCH

        do 130 IDATT = 1, DATTNO
          do 120 IBEAM = 1, NT0
            ATSMK(IBEAM, IDATT) = ATSMK(IBEAM, IDATT)
     .                            + area*ATSAS(IDOM, IBEAM, IDATT)
  120     continue
  130   continue
  140 continue

      return
      end

********************************************************************************
