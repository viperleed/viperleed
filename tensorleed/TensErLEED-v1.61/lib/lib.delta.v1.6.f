C  Tensor LEED subroutine library for delta amplitude calculation
C  v1.2, VB 13.04.00
C  for use with delta.f v1.2
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
C  Please read the comments in delta.f, v1.2
C
C  Subroutines are derived from original codes by P.J. Rous (Tensor LEED) and U. Loeffler
C  (inclusion of thermal vibrations), with modifications by W. Oed, R. Doell, and others.
C
C  One word of advice: If you find a bug in these subroutines, do not ignore
C  it or just patch it somehow. Look into the codes, find the bug and fix it
C  once and for all (and let others using the code know). It'll help a lot
C  of people avoid hitting the same trap over and over again.
C  Another one: Please write your comments in English!!
C
C  Subroutines are included in alphabetical order

C-----------------------------------------------------------------------------------------
C  BELMG computes Clebsh-Gordon coefficients to be used in tmatrix subroutine
       SUBROUTINE BELMG(BELM,NLMB,LMAX,NFAC,FAC)
C
       IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION BELM(NLMB)
       REAL BELM(NLMB),BLM2                                              220395
       DOUBLE PRECISION FAC(NFAC)

       PI=4.0*ATAN(1.0)

C
C      GENERATE FACTORIALS FOR ROUTINE BLM2
C
       NF=4*LMAX+2
       FAC(1)=1.0
       DO 340 IF=1,NF
  340  FAC(IF+1)=DBLE(IF)*FAC(IF)
C
C      NOW GENERATE B(LM,LPPMPP,LP-MP)*4*PI*(-1.)**M
C
       K=0
C
       DO 1 L=0,LMAX
       DO 1 LP=0,LMAX
       LL2=L+LP
       DO 1 M=-L,L
        PRE=4.0*PI*((-1.)**M)
       DO 1 MP=-LP,LP
C
       MPP=MP-M
       IMPP=IABS(MPP)
       LL1=IABS(L-LP)
C
       IF(IMPP.GT.LL1)LL1=IMPP
       IF(LL1.GT.LL2) GOTO 1
C
       DO 2 LPP=LL2,LL1,-2
       K=K+1
   2   BELM(K)=PRE*BLM2(L,M,LPP,MPP,LP,-MP,LL2,NFAC,FAC)
C
    1  CONTINUE
C
       RETURN
       END

C-----------------------------------------------------------------------------------------
COMMENT FUNCTION BLM2 PROVIDES THE INTEGRAL OF THE PRODUCT
C       OF THREE SPHERICAL HARMONICS, EACH OF WHICH CAN BE
C       EXPRESSED AS A PREFACTOR TIMES A LEGENDRE FUNCTION.
C       THE THREE PREFACTORS ARE LUMPED TOGETHER AS FACTOR
C       \C\; AND THE INTEGRAL OF THE THREE LEGENDRE FUNCTIONS
C       FOLLOWS GAUNT\S SUMMATION SCHEME SET OUT BY SLATER
C       ATOMIC STRUCTURE, VOL1, 309,310
C  AUTHOR  PENDRY
      FUNCTION  BLM2 (L1, M1, L2, M2, L3, M3, LMAX,NFAC,FAC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL BLM2
      DOUBLE PRECISION FAC
      DIMENSION FAC(NFAC)

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
  420 BLM2 = 0.0
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
  500 BLM2 = 0.0
      RETURN
C
C      CALCULATION OF FACTOR \A\
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
C      CALCULATION OF SUM \B\
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
C      CALCULATION OF FACTOR \C\
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
      BLM2 = (( - 1.0)**IC) * A * B * C
      RETURN
530   WRITE(6,40)L1,M1,L2,M2,L3,M3
      RETURN
      END

C-----------------------------------------------------------------------------------------
C      write information in long DELWV_DWG array (not very nice, but gives
C      easy compatibility to 'superpos'
      SUBROUTINE COPDELWV(DELWV_DWG,DELWV,IDEB,NCSTEP,NDEB,NT0)
C
      COMPLEX DELWV(NCSTEP,NT0), DELWV_DWG(NCSTEP*NDEB,NT0)
C
      IOFS = (IDEB-1)*NCSTEP
      DO 50 NCS=1,NCSTEP
        DO 50 INB=1,NT0
         DELWV_DWG(IOFS+NCS,INB) = DELWV(NCS,INB)
  50  CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------------------------
C     SUBROUTINE HARMONY GENERATES THE SPHERICAL HARMONICS
C     YLM FOR VECTOR C.
C
      SUBROUTINE HARMONY(YLM,C,LMAX,LMMAX,FAC1,FAC2,FAC3)
C
      INTEGER LMAX,LMMAX
      COMPLEX CT,ST,CF,YLM(LMMAX)
      DIMENSION C(3)
      REAL*8 FAC1,FAC2,FAC3
      DIMENSION FAC1(LMAX+1),FAC2(LMMAX),FAC3(LMAX+1)
C
C
      CD=C(2)**2+C(3)**2
      YA=SQRT(CD+C(1)**2)
      B=0.0
      CF=CMPLX(1.0,0.0)
C
      IF(CD-1.0E-7)80,80,70
C
  70  B=SQRT(CD)
      CF=CMPLX(C(2)/B,C(3)/B)
C
  80  CT=C(1)/YA
      ST=B/YA
C
C
      CALL SPHRM4(LMAX,YLM,LMMAX,CT,ST,CF,FAC1,FAC2,FAC3)
C
C
      RETURN
      END

C-----------------------------------------------------------------------------------------
C  SUBROUTINE INAMP READS IN ALL AMPLITUDES IN MOMENTUM SPACE
C
      SUBROUTINE INAMP(IFILE,IFORM,NT0,PSQ,AK2M,AK3M,ALM,EXLM,LPMMAX,
     1                 WORK,LMMAX)
C
      DIMENSION PSQ(2,NT0),AK2M(NT0),AK3M(NT0)
      COMPLEX ALM(LPMMAX),EXLM(LPMMAX,NT0)
      COMPLEX WORK(LMMAX)
C
C  PRESET AK2M,AK3M WITH VERY LARGE VALUE
C
      DO 11 I=1,NT0
           AK2M(I) = 1.E10
           AK3M(I) = 1.E10
   11 CONTINUE
C
   50 IF (IFORM.EQ.0) THEN
      READ (IFILE) NEXIT
      ELSE
      READ (IFILE,10) NEXIT
   10 FORMAT (I3)
      ENDIF
      IF (NEXIT.LT.0) RETURN
      IF (IFORM.EQ.0) THEN
      READ (IFILE) PQ1,PQ2,AK2,AK3,WORK
      ELSE
      READ (IFILE,20) PQ1,PQ2,AK2,AK3
      READ (IFILE,20) (WORK(I),I=1,LMMAX)
   20 FORMAT(5E12.6)
      ENDIF
      IF (NEXIT.EQ.0) THEN
      DO 100 I=1,LPMMAX
      ALM(I)=WORK(I)
  100 CONTINUE
      ELSE
      PSQ(1,NEXIT)=PQ1
      PSQ(2,NEXIT)=PQ2
      AK2M(NEXIT)=AK2
      AK3M(NEXIT)=AK3
      DO 200 I=1,LPMMAX
      EXLM(I,NEXIT)=WORK(I)
  200 CONTINUE
      ENDIF
      GOTO 50
      END

C-----------------------------------------------------------------------------------------
C  SUBROUTINE INDATA READS ALL INFORMATION BELONGING TO ONE ENERGY
C
      SUBROUTINE INDATA(IFILE,IFORM,E,L1,CAF,NT0,XIST,VPIS,VPIO,VO,VV)
C
      REAL E,VPIS,VPIO,VO,VV
      COMPLEX CAF(L1),XIST(NT0)
C
      IF (IFORM.EQ.0) THEN
      READ(IFILE,ERR=100,END=200)E,VPIS,VPIO,VO,VV,L1DAT,CAF,XIST
      ELSE
      READ(IFILE,20,ERR=100,END=200)E,VPIS,VPIO,VO,VV
      READ(IFILE,10,ERR=100,END=200)L1DAT
      READ(IFILE,20,ERR=100,END=200)(CAF(I),I=1,L1)
      READ(IFILE,20,ERR=100,END=200)(XIST(I),I=1,NT0)
   10 FORMAT (I3)
   20 FORMAT (5E16.10)
      ENDIF
      IF (L1DAT.NE.L1) THEN
      WRITE (6,*) 'L1=',L1,'L1DAT=',L1DAT,'ASSUMED WRONG INPUTFILE!'
      STOP
      ENDIF
      RETURN
C
  100 WRITE (6,*) 'ERROR OCURRED IN READING INPUTFILE'
      STOP
  200 E=-1.
      RETURN
      END

C-----------------------------------------------------------------------------------------
C      ROUTINE MATEL_DWG EVALUATES THE CHANGE IN AMPLITUDE
C      DELWV FOR EACH OF THE EXIT BEAMS FOR EACH OF THE
C      DISPLACEMENTS GIVEN THE SPH WAVE AMPLITUDES
C      CORRESPONDING TO THE INCIDENT WAVE ALM & FOR EACH
C      OF THE TIME REVERSED EXIT BEAMS EXLM.
C      It is a slightly modified version from MATEL5, to give propper 
C      variation of DebyeTemperature AND Geometry                        070493
C      A changed TMATRIX0 --> TMATRIX_DWG has been neccessary
C
C      DELWV(NCSTEP,NT0): CHANGE IN AMPLITUDE DUE TO
C                         DISPLACENT C FOR EACH DIPLACEMENT
C                         & FOR EACH EXIT BEAM
C
C      ALM(LMMAX): SPH WAVE AMPLITUDES INCIDENT AT THE
C                  ORIGIN OF THE TOP LAYER DUE TO THE
C                  INCIDENT LEED BEAM.
C
C      EXLM(NT0,LMMAX): AS ALM BUT FOR EACH TIME REVERSED
C                       EXIT BEAM.
C
C      C(3): CURRENT DISPLACEMENT, C(1)= COMPONENT ALONG
C            X INTO THE SURFACE. C(2),C(3) ALONG ARB1/ARB2
C      CSTEP(3): INCREMENT IN DISPLACEMENT.
C
C      NCSTEP: NUMBER OF DISPLACEMENTS
C *****   NOTE ALL DISPLACEMENTS IN ANGSTROMS *********
C
C      NT0: NUMBER OF EXIT BEAMS
C      NRATIO: RATIO OF AREA OF SURFACE UNIT CELL OF
C               RECONSTRUCTED SURFACE TO UNIT CELL AREA
C              OF THE UNRECONSTRUCTED SURFACE. E.G. FOR
C              P(2X2) NRATIO=4, FOR C(2X2) NRATIO=2.
C
C
      SUBROUTINE MATEL_DWG(DELWV,NCSTEP,AF,NewAF,BJ,BELM,NLMB,            070493
     1 DELTAT,GTWOC,YLM,E,VV,VPI,LMAX,LMAX1,LMMAX,NT0,EXLM,ALM,AK2M,
     1 AK3M,NRATIO,TV,LPMAX,LP1,LPMMAX,NATOMS,CDISP,CUNDISP,PSQ,
     * GTEMP,FAC1,FAC2,FAC3,LMAX21,LMMAX2)
C
      INCLUDE "GLOBAL"
C
      COMPLEX GTEMP(LPMMAX,LPMMAX)
C
      COMPLEX NewAF(LP1),AF(LP1),DELTAT(LMMAX,LMMAX)                      070493
      COMPLEX*16 BJ(2*LPMAX+1)
      COMPLEX GTWOC(LPMMAX,LPMMAX),YLM((2*LPMAX+1)**2)
      COMPLEX EXLM(LMMAX,NT0),ALM(LMMAX)
      COMPLEX DELWV(NCSTEP,NT0)
      COMPLEX AMAT,XA,CAK,CI,PK,XAC
      DIMENSION C(3),BELM(NLMB)
      DIMENSION PSQ(2,NT0), AK2M(NT0),AK3M(NT0)
      DIMENSION CDISP(NCSTEP,NATOMS,3),CUNDISP(NATOMS,3)
      INTEGER LMAX21,LMMAX2
      REAL*8 FAC1,FAC2,FAC3
      DIMENSION FAC1(LMAX21),FAC2(LMMAX2),FAC3(LMAX21)
C
C
      CI=CMPLX(0.0,1.0)
C
C     LOOP OVER MODEL STRUCTURES
C
      DO 1 NC=1,NCSTEP
C
C     SET THE CHANGE IN AMPLITUDES TO ZERO FOR EACH EXIT BEAM.
C
      DO 4 NEXIT=1,NT0
  4   DELWV(NC,NEXIT) = CMPLX(0.0,0.0)
C
C
C     LOOP OVER THE ATOMS OF THE RECONSTRUCTED UNIT CELL.
C
      DO 5 NR=1,NATOMS
C
C     SET NEW DISPLACEMENT
C     FOR THE CURRENT ATOM.
C
       CTEMP=0.0
       DO 2 J=1,3
       CTEMP=CTEMP+ABS(CDISP(NC,NR,J))
   2   C(J)=CDISP(NC,NR,J)/BOHR
C
C
C      THE VECTOR C MUST BE EXPRESSED W.R.T A RIGHT HANDED SET
C      OF AXES. CDISP() & CUNDISP() ARE INPUT W.R.T. A LEFT HANDED
C      SET OF AXES. TO CONVERT C() TO A R.H. SET OF AXES THE
C      TRANSFORMATION  Y > -Y = C(2) > -C(2) IS USED.
C
       C(3)=-C(3)
C
C      EVALUATE DELTAT MATRIX FOR CURRENT DISPLACEMENT.
C
       CALL TMATRIX_DWG(AF,NewAF,BJ,C,DELTAT,GTWOC,YLM,BELM,NLMB,E,
     1 VPI,LPMAX,LP1,LPMMAX,LMAX,LMAX1,LMMAX,GTEMP,.TRUE.,FAC1,FAC2,
     2 FAC3,LMAX21,LMMAX2)
C
C      LOOP OVER EXIT BEAMS
C
         DO 11 NEXIT=1,NT0
C
C      EVALUATE MATRIX ELEMENT.
C
         EMERGE=2.0*(E-VV)-AK2M(NEXIT)**2-AK3M(NEXIT)**2
         IF(EMERGE.LT.0.0) GOTO 11
       AMAT=CMPLX(0.0,0.0)
C
          DO 3 L=0,LMAX
          DO 3 M=-L,L
C
          AM= (-1.0)**M
          I=L+1
          I=I*I-L+M
          IM=I-2*M
C
           DO 3 LP=0,LMAX
           DO 3 MP=-LP,LP
           IP= LP+1
           IP= IP*IP-LP+MP
C
   3   AMAT=AMAT+AM*EXLM(IM,NEXIT)*DELTAT(I,IP)*ALM(IP)
C
C
C     EVALUATE PREFACTOR
C
       D2 = AK2M(NEXIT)
       D3 = AK3M(NEXIT)
       D= D2*D2 + D3*D3
       CAK=CMPLX(2.0*E,-2.0*VPI+0.0000001)
       CAK=CSQRT(CAK)
C
C
       IF(D.GE.2.0*(E)) GOTO 11
C
C
C      XA IS EVALUATED RELATIVE TO THE MUFFIN TIN ZERO I.E. IT
C      USES ENERGY= INCIDENT ELECTRON ENERGY + INNER POTENTIAL.
C
       XA= CMPLX( 2.0*E-D,-2.0*VPI+0.0000001)
C
       XA=CSQRT(XA)
C
C
C     PK=CEXP(CI*DELTAK*UNDISLACED POSN IN UNIT CELL)
C
       DELTK=PSQ(1,NEXIT)*CUNDISP(NR,2)+PSQ(2,NEXIT)*CUNDISP(NR,3)
       DELTK=DELTK/BOHR
       PK=CMPLX(0.0,DELTK)
       PK=CEXP(PK)
C
C
       AMAT=AMAT*PK/(2.0*CAK*TV*XA*FLOAT(NRATIO))
       DELWV(NC,NEXIT)=DELWV(NC,NEXIT)+AMAT
C
  11   CONTINUE
C
  5    CONTINUE
C
  1    CONTINUE
C
       RETURN
       END
C-----------------------------------------------------------------------------------------
C  SUBROUTINE OPENIN OPENS THE INPUTFILE FOR THE INFORMATION PRODUCED IN
C  PROGRAM 'PRETUNG' (historic first reference calculation by Rous)
C
      SUBROUTINE OPENIN(IFILE,FILENAM,IFORM)
C
      CHARACTER*(*) FILENAM
C
      IF (IFORM.EQ.0) THEN
C  UNFORMATTED INPUT
      OPEN (IFILE,FILE=FILENAM,FORM='UNFORMATTED',STATUS='OLD')
      ELSE
C  FORMATTED INPUT
      OPEN (IFILE,FILE=FILENAM,FORM='FORMATTED',STATUS='OLD')
      ENDIF
      REWIND(IFILE)
      RETURN
      END

C-----------------------------------------------------------------------------------------
C  SUBROUTINE OUTDELT_DWG OPENS THE OUTPUTFILE FOR ALL DELTA AMP'S PRODUCED BY
C  PROGRAM ..NACH AND WRITES ALL ADDITIONAL REQUIRED INFORMATION.
C  Parameter for Debye and geoms are written, AID is used for debye information;
C  original AID gets lost !!!!
C  for compatibility to SUPERPOS, formats not changed.
C
      SUBROUTINE OUTDELT_DWG(IFILE,FILENAM,IFORM,NT0,PQFEX,THETA,FI,
     1               RAR1,RAR2,NATOMS,CUNDISP,NCSTEP,CDISP,CDISP_DWG,     281093
     2               AID,NDEB,DR0_A,DRPER_A,DRPAR_A)
C
      CHARACTER*(*) FILENAM
      DIMENSION PQFEX(2,NT0),RAR1(2),RAR2(2)
      DIMENSION CDISP(NCSTEP,NATOMS,3)
      DIMENSION CDISP_DWG(NCSTEP*NDEB,NATOMS,3)
      DIMENSION CUNDISP(NATOMS,3),AID(NCSTEP*NDEB)
      DIMENSION DR0_A(NDEB),DRPER_A(NDEB),DRPAR_A(NDEB)
C
      NTOT = NCSTEP*NDEB
      DO 30 ID=1,NDEB
       DO 30 IA=1,NATOMS
         DO 30 IC=1,NCSTEP
          DO 25 I=1,3
            CDISP_DWG((ID-1)*NCSTEP+IC,IA,I) = CDISP(IC,IA,I)
   25     CONTINUE
          IAID = (ID-1)*NCSTEP+IC                                        230793
          AID(IAID) = 0.                                                 230793
          IF (IC .EQ. 1) THEN
            AID(IAID) = DRPER_A(ID)                                      230793
          ENDIF
   30 CONTINUE
C
      IF (IFORM.EQ.0) THEN
C  UNFORMATTED OUTPUT
      OPEN (IFILE,FILE=FILENAM,FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE (IFILE) THETA,FI,RAR1,RAR2,NT0,NATOMS,NCSTEP*NDEB
      WRITE (IFILE) PQFEX,CUNDISP,CDISP_DWG,AID
      ELSE
C  FORMATTED OUTPUT
      OPEN (IFILE,FILE=FILENAM,FORM='FORMATTED',STATUS='UNKNOWN')
      WRITE (IFILE,100) THETA,FI,(RAR1(I),I=1,2),(RAR2(I),I=1,2)
      WRITE (IFILE,110) NT0,NATOMS,NCSTEP*NDEB
      WRITE(6,*) ' OUTDELT_DWG: NT0, NATOMS,NCSTEP*NDEB '
      WRITE (6,110) NT0,NATOMS,NCSTEP*NDEB
      WRITE (IFILE,'(10F10.5)')((PQFEX(I,J),I=1,2),J=1,NT0)
      WRITE (IFILE,120)((CUNDISP(I,J),J=1,3),I=1,NATOMS)
      WRITE (IFILE,120)(((CDISP_DWG(I,J,K),K=1,3),J=1,NATOMS),I=1,
     1 NCSTEP*NDEB)
      WRITE (IFILE,120)(AID(I),I=1,NCSTEP*NDEB)
C
  100 FORMAT (6E13.7)
  110 FORMAT (3I3)
  120 FORMAT (10F7.4)
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------------------------
C  SUBROUTINE OUTRINT WRITES ALL DELTA AMP'S TO FILE# IFILE
C
      SUBROUTINE OUTRINT(IFILE,IFORM,E,VV,VO,VPI,NT0,NCSTEP,
     1                   XIST,DELWV)
C
      COMPLEX XIST(NT0),DELWV(NCSTEP,NT0)
C
      IF (IFORM.EQ.0) THEN
      WRITE (IFILE) E,VPI,VO,VV,XIST,DELWV
      ELSE
      WRITE (IFILE,100) E,VPI,VO,VV
      WRITE (IFILE,100) (XIST(I),I=1,NT0)
      WRITE (IFILE,100) ((DELWV(I,J),J=1,NT0),I=1,NCSTEP)
  100 FORMAT (6E13.7)
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------------------------
C  Subroutine ReadBas reads basic information on preceding reference calculation as well
C  as information regarding desired perturbations of present atomic site from stdin, unit 5.

      subroutine ReadBas(EI,EF,RAR1,RAR2,TV,THETA,FI,FORMIN,PQFEX,NT0,
     +                   ES,PHSS,NPSI,NEL,L1,LMAX,FORMOUT,IEL,CUNDISP,
     +                   CDISP,NATOMS,NCSTEP,DR0_A,DRPER_A,DRPAR_A,
     +                   NDEB)

      INCLUDE "GLOBAL"

C  variable dimensions

      integer NT0,NPSI,NEL,NATOMS,NCSTEP,NDEB

C  global variables

      real EI,EF
      real RAR1,RAR2
      dimension RAR1(2),RAR2(2)
      real TV,THETA,FI
      integer FORMIN
      real PQFEX
      dimension PQFEX(2,NT0)
      real ES,PHSS
      dimension ES(NPSI),PHSS(NPSI,NEL,L1)
      integer FORMOUT
      real CUNDISP,CDISP
      dimension CUNDISP(NATOMS,3),CDISP(NCSTEP,NATOMS,3)
      integer IEL
      real DRPER_A,DRPAR_A,DR0_A
      dimension DRPER_A(NDEB),DRPAR_A(NDEB),DR0_A(NDEB)

      integer LMAX

C  local variables

      real PI
      Character*80 TITLE
      real ARA1,ARA2,TVA,ATV
      dimension ARA1(2),ARA2(2)
      real ARB1,ARB2,TVB
      dimension ARB1(2),ARB2(2)
      integer NELCHECK,NCCHECK,NDCHECK

C  set PI

      PI = 3.14159265

C  begin readin

C  first, read information on the reference calculation that was already given there!

C  title

      read(5,'(A80)') TITLE
      write(6,*) TITLE

C  energy range - do not convert to Hartrees

      read(5,'(3F7.2)') EI,EF

C  info on lateral periodicity: substrate lattice
C  ARA1 AND ARA2 ARE TWO 2-D BASIS VECTORS OF THE SUBSTRATE LAYER
C  LATTICE. THEY SHOULD BE EXPRESSED IN TERMS OF THE PLANAR CARTESIAN
C  Y- AND Z-AXES (X-AXIS IS PERPENDICULAR TO SURFACE)(ANGSTROM)

      READ(5,'(2F7.4)')(ARA1(I),I=1,2)
      READ(5,'(2F7.4)')(ARA2(I),I=1,2)

c      WRITE(6,'(2F7.4)')(ARA1(I),I=1,2)
c      WRITE(6,'2F7.4)')(ARA2(I),I=1,2)

      DO I=1,2
        ARA1(I)=ARA1(I)/BOHR
        ARA2(I)=ARA2(I)/BOHR
      ENDDO

      TVA=ABS(ARA1(1)*ARA2(2)-ARA1(2)*ARA2(1))
      ATV=2.0*PI/TVA

C  RAR1 AND RAR2 ARE THE RECIPROCAL-LATTICE BASIS VECTORS CORRESPONDING
C  TO ARA1 AND ARA2

      RAR1(1)=ARA2(2)*ATV
      RAR1(2)=-ARA2(1)*ATV
      RAR2(1)=-ARA1(2)*ATV
      RAR2(2)=ARA1(1)*ATV

C  info on lateral periodicity: superlattice 

      READ(5,'(2F7.4)')(ARB1(I),I=1,2)
      READ(5,'(2F7.4)')(ARB2(I),I=1,2)
c      WRITE(6,'(2F7.4)')(ARB1(I),I=1,2)
c      WRITE(6,'(2F7.4)')(ARB2(I),I=1,2)
      DO 467 I=1,2
      ARB1(I)=ARB1(I)/BOHR
  467 ARB2(I)=ARB2(I)/BOHR
      TVB=ABS(ARB1(1)*ARB2(2)-ARB1(2)*ARB2(1))

C  pass on TV = TVB to main (note that this holds only for current handling of TLEED)

      TV = TVB

C  incidence angles theta, fi

      READ (5,'(2F7.2)')  THETA, FI

C s.w. conversion in radians, 02.02.00
      THETA=THETA*PI/180.0
      FI=FI*PI/180.0

C  FORMIN determines whether tensor input file AMP was written formatted (1) or
C         unformatted (0) 

      read(5,'(I4)') FORMIN

C  read beam list for which TLEED is to be performed, PQFEX, in units of RAR1,RAR2

      DO IBEAM = 1,NT0

        READ(5,'(2F10.5)') PQFEX(1,IBEAM),PQFEX(2,IBEAM)

      ENDDO

C  read phase shifts

C  NEL= NO. OF CHEMICAL ELEMENTS FOR WHICH PHASE SHIFTS ARE TO BE READ
C  IN

      READ(5,'(I3)') NELCHECK

      IF (NELCHECK.ne.NEL) THEN
        write(6,*) "NEL from phaseshift file not consistent with MNEL!"
        write(6,*) "Please correct!"
        STOP
      ENDIF
 
      DO 660 IENER = 1, NPSI

C  ES= ENERGIES (HARTREES) AT WHICH PHASE SHIFTS ARE INPUT. LINEAR
C  INTERPOLATION OF THE PHASE SHIFTS WILL OCCUR FOR ACTUAL ENERGIES
C  FALLING BETWEEN THE VALUES OF ES (AND LINEAR EXTRAPOLATION ABOVE
C  THE HIGHEST ES)
C  PHSS STORES THE INPUT PHASE SHIFTS (RADIAN)

        READ (5,'(F7.4)')  ES(IENER)

        DO IEL=1,NEL,1
          READ(5,'(10F7.4)')(PHSS(IENER,IEL,L),L=1,L1)
          IF (LMAX .LT. 10) THEN
            READ(5,*)
          END IF
        ENDDO

  660 CONTINUE

C  next, read original input for delta amplitude calculation

C  OUTPUT: formatted or unformatted

      read (5,'(I4)') FORMOUT

      read(5,*)
      read(5,*)
      read(5,*)

C  IEL is number of the element that is assumed in this calculation - 
C      preparation of chemical TLEED

      read (5,'(I4)') IEL

C  info on atomic positions - note NATOMS doesn't work, at least not with thermal TLEED

      do IATOM = 1,NATOMS

        read (5,*)
        read (5,*)
        read (5,*)
        read (5,'(3F7.4)') (CUNDISP(IATOM,IPOS), IPOS=1,3)

        do IPOS = 1,3
          CUNDISP(IATOM,IPOS) = CUNDISP(IATOM,IPOS)
        enddo

        read (5,*)
        read (5,*)
        read (5,*)
        read (5,'(I4)') NCCHECK

        IF (NCCHECK.ne.NCSTEP) THEN
          write(6,*) "NCSTEP from input not consistent with MNCSTEP!"
          write(6,*) "Please correct!"
          STOP
        ENDIF
 
        do ICSTEP = 1,NCSTEP

          read (5,'(3F7.4)') (CDISP(ICSTEP,IATOM,IPOS), IPOS=1,3)

          do IPOS = 1,3

            CDISP(ICSTEP,IATOM,IPOS) = CDISP(ICSTEP,IATOM,IPOS)
 
            CDISP(ICSTEP,IATOM,IPOS) = CDISP(ICSTEP,IATOM,IPOS) 
     +                               - CUNDISP(IATOM,IPOS)

          enddo

        enddo

      enddo

      read(5,*)
      read(5,*)
      read(5,*)

C  number of vibrational amplitude steps to be performed

      read(5,'(I4)') NDCHECK

      IF (NDCHECK.ne.NDEB) THEN
        write(6,*) "NDEB from input not consistent with MNDEB!"
        write(6,*) "Please correct!"
        STOP
      ENDIF
 
C  read vibrational amplitudes DRPER_A, set DRPAR_A, DR0_A accordingly

      DO IDEB = 1,NDEB

        read(5,'(F7.4)') DRPER_A(IDEB)

        DRPER_A(IDEB) = DRPER_A(IDEB)/BOHR
        DRPAR_A(IDEB) = DRPER_A(IDEB)
        DR0_A(IDEB)   = 0.

      ENDDO

      end

C-----------------------------------------------------------------------------------------
C  SUBROUTINE SPHRM COMPUTES SPHERICAL HARMONICS.
C  AUTHOR  PENDRY.
C   LMAX= LARGEST VALUE OF L.
C   YLM= OUTPUT COMPLEX SPHERICAL HARMONICS.
C   LMMAX= (LMAX+1)**2.
C   CT= COS(THETA) (COMPLEX).
C   ST= SIN(THETA) (COMPLEX).
C   CF= CEXP(I*FI).
      SUBROUTINE  SPHRM4(LMAX,YLM,LMMAX,CT, ST, CF, FAC1,FAC2,FAC3)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX  YLM
      COMPLEX  CT, ST, CF, SF, SA

C   Note that the lmax used here is 2*lmax+1 as defined in the input. 

      REAL*8 FAC1,FAC2,FAC3
      DIMENSION  FAC1(LMAX+1),FAC3(LMAX+1),FAC2(LMMAX),YLM(LMMAX)                 200792

      PI = 3.14159265
      LM = 0
      CL = 0.0
      A = 1.0
      B = 1.0
      ASG = 1.0
      LL = LMAX + 1
      DO 550 L = 1, LL
      FAC1(L) = ASG * SQRT((2.0 * CL + 1.0) * A/(4.0 * PI * B * B))
      FAC3(L) = SQRT(2.0 * CL)
      CM =  - CL
      LN = L + L - 1
      DO 540 M = 1, LN
      LO = LM + M
      FAC2(LO) = SQRT((CL + 1.0 + CM) * (CL + 1.0 - CM)/((2.0 * CL + 3.
     10) * (2.0 * CL + 1.0)))
  540 CM = CM + 1.0
      CL = CL + 1.0
      A = A * 2.0 * CL * (2.0 * CL - 1.0)/4.0
      B = B * CL
      ASG =  - ASG
  550 LM = LM + LN
      LM = 1
      CL = 1.0
      ASG =  - 1.0
      SF = CF
      SA = CMPLX(1.0,0.0)
      YLM(1) = CMPLX(FAC1(1),DBLE(0.0))
      DO 560 L = 1, LMAX
      LN = LM + L + L + 1
      YLM(LN) = FAC1(L + 1) * SA * SF * ST
      YLM(LM + 1) = ASG * FAC1(L + 1) * SA * ST/SF
      YLM(LN - 1) =  - FAC3(L + 1) * FAC1(L + 1) * SA * SF * CT/CF
      YLM(LM + 2) = ASG * FAC3(L + 1) * FAC1(L + 1) * SA * CT * CF/SF
      SA = ST * SA
      SF = SF * CF
      CL = CL + 1.0
      ASG =  - ASG
  560 LM = LN
      LM = 1
      LL = LMAX - 1
      DO 580 L = 1, LL
      LN = L + L - 1
      LM2 = LM + LN + 4
      LM3 = LM - LN
      DO 570 M = 1, LN
      LO = LM2 + M
      LP = LM3 + M
      LQ = LM + M + 1
      YLM(LO) =  - (FAC2(LP) * YLM(LP) - CT * YLM(LQ))/FAC2(LQ)
  570 CONTINUE
  580 LM = LM + L + L + 1
      RETURN
      END

C-----------------------------------------------------------------------------------------
C     SUBROUTINE TMATRIX_DWG GENERATES THE TMATRIX(L,L') MATRIX
C     FOR A GIVEN ENERGY & DISPLACEMENT VECTOR.
C     This modified Version of TMATRIX0 enables simultaneous displacement
C     and a Change of the Debye Temperature.                              070493
C
C     E,VPI: CURRENT ENERGY (REAL,IMAGINARY)
C     C(3) : DISPLACEMENT VECTOR;
C            C(1)= COMPONENT ALONG X AXIS INTO THE SURFACE
C            C(2)= COMPONENT ALONG Y AXIS (ARB1)
C            C(3)= COMPONENT ALONG Z AXIS (ARB2)
C     DELTAT(LMMAX,LMMAX): CHANGE IN T MATRIX CAUSED BY
C                          THE DISPLACEMENT.
C     AF(LMAX1): EXP(CI*PHS(L))*SIN(PHS(L)). NOTE THAT
C                ATOMIC T MATRIX IS CI*AF
C     BJ(LMAX1): BESSEL FUNCTIONS JN(CAPPA*C) FOR EACH L
C     YLM(LMMAX): SPHERICAL HARMONICS OF VECTOR C
C     GTWOC(LMMAX,LMMAX): PROPAGATOR FROM ORIGIN TO C
C     LMAX1=LMAX+1
C
C
      SUBROUTINE TMATRIX_DWG(AF,NewAF,BJ,C,DELTAT,GTWOC,YLM,BELM,NLMB,    070493
     1E,VPI,LMAX,LMAX1,LMMAX,LSMAX,LS1,LSMMAX,GTEMP,AFLAG,FAC1,FAC2,
     2FAC3,LMAX21,LMMAX2)
C
      COMPLEX AF(LMAX1),NewAF(LMAX1),
     1DELTAT(LSMMAX,LSMMAX)                                               070493
      COMPLEX*16 BJ(2*LMAX+1)                                             281093
      COMPLEX GTEMP(LMMAX,LMMAX)
      COMPLEX GTWOC(LMMAX,LMMAX),YLM((2*LMAX+1)**2)
      DIMENSION C(3),BELM(NLMB)
      LOGICAL AFLAG
      COMPLEX CAPPA,Z,PRE,CSUM,CI,CZ,TL1
      INTEGER LMAX21,LMMAX2
      REAL*8 FAC1,FAC2,FAC3
      DIMENSION FAC1(LMAX21),FAC2(LMMAX2),FAC3(LMAX21)
C
C
C     SET INITIAL VARIABLES

      PI=4.0*ATAN(1.0)
      CI=CMPLX(0.0,1.0)
      CZ=CMPLX(0.0,0.0)
C
C
      LMAX2=2*LMAX
ctest
      if (LMAX21.ne.LMAX2+1) then
        write(6,*) "Dimension error in LMAX21:"
        write(6,*) "LMAX21 = MN : ",LMAX21
        write(6,*) "LMAX2 + 1   : ",LMAX2+1
        stop
      end if
ctest - end

ctest
      if (LMMAX2.ne.LMAX21*LMAX21) then
        write(6,*) "Dimension error in LMMAX2: "
        write(6,*) "LMMAX2 = MNN : ",LMMAX2
        write(6,*) "LMAX21*LMAX21: ",LMAX21*LMAX21
        stop
      end if
ctest - end
C
C     if displacement = 0, calculate deltaT and jump to END
      CL=SQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))
      IF (CL .le. 1.0e-7) THEN
C        set deltat to zero
         DO 26 I=1,LSMMAX
         DO 26  J=1,LSMMAX
  26     DELTAT(I,J)=CZ
C        calculate new deltat
         DO 2 L=0,LSMAX
         DO 2 M=-L,L
            I = L+1
            I= I*I-L+M
   2     DELTAT(I,I)=CI*(NewAF(L+1)-AF(L+1))
         GOTO 33
      ENDIF
C
      IF (AFLAG) THEN
C
C
C     GENERATE BESSEL FUNCTIONS FOR Z=CAPPA*MOD(C)
C
      CL=SQRT(C(1)*C(1)+C(2)*C(2)+C(3)*C(3))
      CAPPA=CMPLX(2.0*E,-2.0*VPI)
      Z=CSQRT(CAPPA)*CL
C
      CALL BESSEL(BJ,Z,LMAX2,LMAX21)
C
C     GENERATE SPHERICAL HARMONICS YLM FOR VECTOR C
C
      CALL HARMONY(YLM,C,LMAX2,LMMAX2,FAC1,FAC2,FAC3)
C
C     GENERATE THE PROPOGATOR GTWOC
C
C     DO 21 I=1,LMMAX
C     DO 21 J=1,LMMAX
C 21  GTWOC(I,J)=CZ
C
C
      K=0
      DO 1 L=0,LMAX
      IS=(L+1)*(L+1)-L
      DO 1 LP=0,LMAX
      ISP=(LP+1)*(LP+1)-LP
      LL2=LP+L
      LL1S=IABS(L-LP)
      PRE=CI**(L+LP)
      DO 1 M=-L,L
      DO 1 MP=-LP,LP
C
C     I=(L+1)**2-L+M
C     IP=(LP+1)**2-LP+MP
      I=IS+M
      IP=ISP+MP
      GTWOC(I,IP)=CZ
C     GTWOC(I,IP)=CZ
C
      MPP=MP-M
      IMPOS=IABS(MPP)
C     LL1=IABS(L-LP)
C     IF(IMPOS.GT.LL1) LL1=IMPOS
      LL1=MAX(LL1S,IMPOS)
C     IF(LL1.GT.LL2) GOTO 1
      DO 4 LPP=LL2,LL1,-2
      K=K+1
C
C     IPPM=(LPP+1)**2-LPP-MPP
      IPPM=LPP*LPP+LPP+1-MPP
      CSUM=BJ(LPP+1)*YLM(IPPM)*BELM(K)*CI**(-LPP)
C
   4  GTWOC(I,IP)=GTWOC(I,IP)+PRE*CSUM
C
   1  CONTINUE
C

C     GENERATE T MATRIX USING
C     GTWCO=(-1.0)**(L+LP)
C
C
      DO 23 I=1,LMMAX
      I1=0
      DO 23 L=0,LMAX
      DO 23 M=-L,L
      I1=I1+1
 23   GTEMP(I,I1)=GTWOC(I,I1)*CI*NewAF(L+1)                              070493
C
C
      DO 210 I=1,LSMMAX
      DO 25  J=1,LSMMAX
  25  DELTAT(I,J)=CZ
      DO 210 L=1,LSMMAX
      DO 210 J=1,LSMMAX
      DELTAT(I,J)=DELTAT(I,J)+GTEMP(I,L)*GTWOC(L,J)
 210  CONTINUE                                                           260689
C
C
C
C     TEST CONSERVATION OF CROSSECTION
C

C     CALL FLUX(DELTAT,NewAF,LSMAX,LSMMAX,LS1,E)                         070493
C
C     EVALUATE DELTAT BY SUBTRACTING ATOMIC T MATRIX
C
      ELSE
C
      DO 725 I=1,LSMMAX
      DO 725  J=1,LSMMAX
 725  DELTAT(I,J)=CZ
C
      ENDIF
C
      DO 3 L=0,LSMAX
      DO 3 M=-L,L
C
      I = L+1
      I= I*I-L+M
C
   3  DELTAT(I,I)=DELTAT(I,I)-CI*AF(L+1)
  33  CONTINUE
C
C
      RETURN
      END

C-----------------------------------------------------------------------------------------
