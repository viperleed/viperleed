C  Utility for the generation of VHT style beam lists
C  for use with TensErLEED ref-calc.f v1.7 .
C  VB, 14.04.00, modified wrt field dimensions and input formats by LH, 11.03.2020
C  removed input of symmetry-code (set fix to SSYM=1),
C  because not usable in TensorLEED routines, LH 11.03.2020
C
C  Program includes possible
C  symmetry usage. Symmetry, however, is currently only applicable for
C  mere full dynamic calculation, not Tensor LEED run. Possible future
C  improvement are inclusion of beamgen into ref-calc.f and the use of
C  symmetry in its full-dynamic part, expanding the symmetrised beam list
C  into the full one for all calculations that involve off-normal incidence
C  later.
C

C  PROGRAM BEAMGEN_U
C
      PROGRAM BEAMGEN_U
C  THIS PROGRAM GENERATES A LIST OF BEAMS G FOR INPUT INTO A LEED       220179
C  PROGRAM. THE BEAMS ARE ORDERED BY INCREASING ABS(G) AND GROUPED
C  BY SUBSETS INDUCED BY A SUPERLATTICE. THE LIST IS SHORTENED IN
C  THE PRESENCE OF SYMMETRY, THE TYPE OF SYMMETRY OF THE SURFACE
C  BEING INPUT THROUGH THE CODE NUMBER SSYM.
C
C  Updated Version 14.02.91
C  (c) 1991 U.Loeffler
C  Since the old Version produces not ALL necessary beams, when there
C  is not a rectangular Coordinatesystem, this version is designed to
C  deal with this deficiency.
C  Any angle between the lattice vectors is allowed and this program
C  should work properly now.
C  But it is not tested under ANY circumstances, so perfect performance
C  is NOT granted.
C  The following changes are made:
C  - SALPH=sin(alpha) is calculated, where alpha is the angle between
C    the reciprocal lattice vectors RBR1 and RBR2.
C  - NIT, the cutoff criterion in RBR2-direction was calculated as
C    Max(NI2) at NI1=0. This is changed to NIT=SQRT(GMAX2)/sin(alpha)*
C    |RBR2|. Since NIT must be an integer and the program uses beams
C    including NIT+1, NIT will be int(NIT+1) - 1.
C  - NIT2 is introduced to match the old NIT in order to compare NIT
C    with NIT2
C  - WRITE(13, ) statements are introduced to have some additional debug
C    aid, but they are as a comment now.
C  - All changes are documented as '140291'
C  - The file 'DATA_U' contains test data-files and a short documentaion
C    about the performed tests.
C
C  Modified version 21.04.92
C  R.Doell
C  - Version that can be used also for larger superlattices (DIM KNBS(999))
C  - Error occured with SSYM=10 in DELBM because of reordering of beams
C    in UNFOLD
C  - WRITE (9, ) statements are introduced to have some additional debug
C    aid, but they are as a comment now.
C  - All changes are documented as '210492'
C
C  SSYM IS DEFINED AS FOLLOWS (Y AND Z ARE THE COORDINATES IN THE
C  SURFACE PLANE, ROTATION AXES ARE PERPENDICULAR TO THE SURFACE) -
C   SSYM       SYMMETRY
C     1      NONE (1-FOLD ROT. AXIS)
C     2      2-FOLD ROT. AXIS
C     3      MIRROR PLANE Z=-Y
C     4      MIRROR PLANE Z=Y
C     5      MIRROR PLANE Z=0
C     6      MIRROR PLANE Y=0
C     7      4-FOLD ROT. AXIS
C     8      2 MIRROR PLANES Z=-Y AND Z=Y
C     9      2 MIRROR PLANES Z=0 AND Y=0
C    10      4 MIRROR PLANES Z=-Y, Z=Y, Z=0 AND Y=0
C    11      3-FOLD ROT. AXIS
C    12      3-FOLD ROT. AXIS AND MIRROR PLANE Z=0
C    13      3-FOLD ROT. AXIS AND MIRROR PLANE Y=0
C    14      6-FOLD ROT. AXIS
C    15      6-FOLD ROT. AXIS AND 2 MIRROR PLANES Z=0 AND Y=0
C
      INTEGER SSYM,LATMAT(2,2),KSYM(2,99999),NBSET(99999),KNB(999)          210492 LH 110320
      DIMENSION ARA1(2),ARA2(2),SPQF(2,99999),G(2,99999),GS1(2),GS2(2)
      DIMENSION RAR1(2),RAR2(2)
C  ARA1,ARA2 - BASIC SUBSTRATE LATTICE VECTORS (ANGSTROM)
      OPEN(5,FILE='DATA',STATUS='UNKNOWN')
      OPEN(6,FILE='PROT',STATUS='UNKNOWN')
      OPEN(7,FILE='BELIST',STATUS='UNKNOWN')
      READ(5,160)ARA1
      READ(5,160)ARA2
160   FORMAT(2F9.4)
      CALL RMPRT(ARA1,1,2,4HARA1)
      CALL RMPRT(ARA2,1,2,4HARA2)
C  LATMAT - MATRIX TRANSFORMING BASIC SUBSTRATE LATTICE INTO
C  SUPERLATTICE
      DO 10 I=1,2
10    READ(5,200)(LATMAT(I,J),J=1,2)
      CALL IMPRT(LATMAT,2,2,6HLATMAT)
200   FORMAT(10I4)  ! Also change error condition and message at L136 and following!
C  SSYM - CODE NUMBER DESCRIBING THE OVERALL SURFACE SYMMETRY
      SSYM=1
      WRITE(6,201)SSYM
201   FORMAT(5H SSYM,I5)
C  EMAX - MAXIMUM ENERGY OF LEED CALCULATION (EV).
C  DMIN - MINIMUM INTERLAYER SPACING OF LEED CALCULATION (ANGSTROM).
C  TST - CRITERION FOR INCLUSION OF EVANESCENT BEAMS (SEE LEED
C         PROGRAM).
      READ(5,169)EMAX,DMIN
169   FORMAT(F6.1,F7.4)
      READ(5,170)TST
170   FORMAT(1F8.6)
      WRITE(6,171)EMAX,DMIN,TST
171   FORMAT(5H EMAX,F9.4,5H DMIN, F9.4,4H TST,F9.4)
C  KNBMAX - NO. OF BEAMS NOT TO BE EXCEEDED BY THIS PROGRAM.
C    THIS NO. APPLIES BEFORE SYMMETRIZATION
      READ(5,206)KNBMAX
206   FORMAT(I5)
      WRITE(6,172)KNBMAX
172   FORMAT(7H KNBMAX,I5)
      CALL BEAMG(ARA1,ARA2,LATMAT,SSYM,EMAX,DMIN,TST,KNBMAX,
     1SPQF,G,KNT,KSYM,NBSET,KNBS,KNB,RAR1,RAR2)
      WRITE(6,202)SSYM
202   FORMAT(5H SSYM,I5)
      WRITE(6,203)KNT
203   FORMAT(4H KNT,I5)
      WRITE(6,204)KNBS
204   FORMAT(5H KNBS,I5)
      CALL IMPRT(KNB,1,KNBS,3HKNB)
205   FORMAT(1H ,10I3)
      WRITE(6,164)
164   FORMAT(6X,4HSPQF,14X,1HG,9X,6HABS(G),2X,4HKSYM,3X,5HNBSET,
     19H BEAM NO.)
      DO 20 II=1,KNT
      AG=SQRT(G(1,II)*G(1,II)+G(2,II)*G(2,II))
20    WRITE(6,165)(SPQF(I,II),I=1,2),(G(I,II),I=1,2),AG,
     1(KSYM(I,II),I=1,2),NBSET(II),II
165   FORMAT(1H ,2F7.4,2X,2F7.4,2X,F7.4,2I3,2I6)
      IOFF=0
      DO 50 IBS=1,KNBS
      IM=KNB(IBS)
	  if (IM.gt.9999) then  ! Write to stderr
          write(0, *) "### ERROR ### Fortran beamgen.v.1.7 "
          write(0, *) "# Too many beams for format I4: KNB(IBS)= ", IM
          write(0, *) "# If you want to have this many, change format "
          write(0, *) "# specifier at label 200"
          write(0, *) "#############"
          STOP
	  endif
      WRITE(7,200)IM
      DO 40 II=1,IM
      IJ=II+IOFF
      GS1(1)=SPQF(1,IJ)*RAR1(1)
      GS1(2)=SPQF(1,IJ)*RAR1(2)
      GS2(1)=SPQF(2,IJ)*RAR2(1)
      GS2(2)=SPQF(2,IJ)*RAR2(2)
      GS=(GS1(1)+GS2(1))**2 + (GS1(2)+GS2(2))**2
      ES=0.5*GS*27.21
      WRITE(7,180)(SPQF(I,IJ),I=1,2),(KSYM(I,IJ),I=1,2),ES
180   FORMAT(2F10.5,2I3,10X,'E = ',F10.4,' EV')
40    CONTINUE
      IOFF=IOFF+IM
50    CONTINUE
CHL SUBROUTINE NUMBER ERZEUGT EINE AUFSTEIGEND NUMMERIERTE BEAM-LISTE AUF
CHL DEM FILE 'NBLIST'.
      CALL NUMBER
      STOP
      END
C---------------------------------------------------------------------------
      SUBROUTINE BEAMG(ARA1,ARA2,LATMAT,SSYM,EMAX,DMIN,TST,KNBMAX,
     1SPQF,G,KNT,KSYM,NBSET,KNBS,KNB,RAR1,RAR2)
      INTEGER SSYM,LATMAT(2,2),KSYM(2,99999),NBSET(99999),KNB(400)         210492
      INTEGER SSYMA(2,99999),LSYMSA(2,99999)                               210492
      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2),RBR1(2),RBR2(2)
      DIMENSION GDUMMY(2,99999)                                            210492
      DIMENSION ALMR(2,2),SPQF(2,99999),G(2,99999),ST(4),GU(2,12)
      DIMENSION DG(2),GS(2,12)
      INTEGER LSYMN(15),LSYMT(10)
      DATA LSYMN/1,1,1,1,1,1,6,4,4,9,1,3,3,5,10/
      PI=3.1415926535
C  SET UP POSSIBLE SYMMETRY CODE NUMBERS FOR BEAM SUBSETS
      LSYMT(1)=1
      GO TO (150,150,150,150,150,150,20,40,60,70,150,90,100,
     1110,120),SSYM
20    DO 30 I=2,6
30    LSYMT(I)=I
      GO TO 150
40    DO 50 I=2,4
50    LSYMT(I)=I
      GO TO 150
60    LSYMT(2)=2
      LSYMT(3)=5
      LSYMT(4)=6
      GO TO 150
70    DO 80 I=2,9
80    LSYMT(I)=I
      GO TO 150
90    LSYMT(2)=5
      LSYMT(3)=11
      GO TO 150
100   LSYMT(2)=6
      LSYMT(3)=11
      GO TO 150
110   LSYMT(2)=2
      LSYMT(3)=5
      LSYMT(4)=6
      LSYMT(5)=11
      GO TO 150
120   DO 130 I=2,6
130   LSYMT(I)=I
      DO 140 I=7,10
140   LSYMT(I)=I+4
150   CONTINUE
C  ATOMIC UNITS USED
      DO 160 I=1,2
      ARA1(I)=ARA1(I)/0.529
160   ARA2(I)=ARA2(I)/0.529
      TVA=ABS(ARA1(1)*ARA2(2)-ARA1(2)*ARA2(1))
      ATV=2.0*PI/TVA
C  RAR1 AND RAR2 ARE THE RECIPROCAL-LATTICE BASIS VECTORS CORRESPONDING
C  TO ARA1 AND ARA2
      RAR1(1)=ARA2(2)*ATV
      RAR1(2)=-ARA2(1)*ATV
      RAR2(1)=-ARA1(2)*ATV
      RAR2(2)=ARA1(1)*ATV
C  ALMR IS MATRIX RELATING SUBSTRATE TO OVERLAYER LATTICES
      DET=LATMAT(1,1)*LATMAT(2,2)-LATMAT(1,2)*LATMAT(2,1)
      IF (ABS(DET).GE.1.E-5)GO TO 166
      WRITE(6,165)
165   FORMAT(28H LATMAT HAS ZERO DETERMINANT)
      STOP
166   ALMR(1,1)=FLOAT(LATMAT(2,2))/DET
      ALMR(2,2)=FLOAT(LATMAT(1,1))/DET
      ALMR(1,2)=-FLOAT(LATMAT(2,1))/DET
      ALMR(2,1)=-FLOAT(LATMAT(1,2))/DET
C  RBR1,RBR2 IS BASIS OF RECIPROCAL LATTICE OF SUPERLATTICE
      RBR1(1)=ALMR(1,1)*RAR1(1)+ALMR(1,2)*RAR2(1)
      RBR1(2)=ALMR(1,1)*RAR1(2)+ALMR(1,2)*RAR2(2)
      RBR2(1)=ALMR(2,1)*RAR1(1)+ALMR(2,2)*RAR2(1)
      RBR2(2)=ALMR(2,1)*RAR1(2)+ALMR(2,2)*RAR2(2)
C  COS(ALPHA) = RBR1*RBR2 / |RBR1| * |RBR2|                             140291
      CALPH=RBR1(1)*RBR2(1) + RBR1(2)*RBR2(2)                           140291
      CALPH=CALPH/(SQRT(RBR1(1)*RBR1(1) + RBR1(2)*RBR1(2))*             140291
     1SQRT(RBR2(1)*RBR2(1) + RBR2(2)*RBR2(2)))                          140291
C  SIN(ALPHA) = SQRT(1-COS^2(ALPHA)                                     140291
      SALPH=SQRT(1-CALPH*CALPH)                                         140291
      CALL RMPRT(ALMR,2,2,4HALMR)
      CALL RMPRT(RAR1,1,2,4HRAR1)
      CALL RMPRT(RAR2,1,2,4HRAR2)
      CALL RMPRT(RBR1,1,2,4HRBR1)
      CALL RMPRT(RBR2,1,2,4HRBR2)
      GMAX2=2.*EMAX/27.21+(ALOG(TST)/DMIN*0.529)**2
C  GENERATE ALL BEAMS WITHIN A BEAM CIRCLE OF RADIUS SQRT(GMAX2),
C  LIMITED TO KNBMAX BEAMS
      NIT=INT(SQRT(GMAX2)/(ABS(SALPH)*                                  140291
     1SQRT(RBR2(1)*RBR2(1) + RBR2(2)*RBR2(2)))+1)                       140291
C  SINCE NIT+1 IS EFFECTLY USED, IT IS SUFFICIENT TO USE NIT=NIT-1      140291
      NIT=NIT-1                                                         140291
      WRITE(6,3) NIT                                                    140291
    3 FORMAT( 'NIT = ' , I5)                                            140291
C     WRITE(13,5) GMAX2                                                 140291
C   5 FORMAT( 'GMAX2= ', F5.2)                                          140291
      KNT=0
170   NI1=-1
180   NI1=NI1+1
      NOP=0                                                             280679
      NI2=-1
190   NI2=NI2+1                                                         280679
      DO 280 K=1,4
      GO TO (200,210,220,230),K
200   II1=NI1
      II2=NI2
      GO TO 240
210   IF (NI1.EQ.0.AND.NI2.EQ.0) GO TO 290
      IF (NI1.EQ.0) GO TO 280
      II1=-NI1
      II2=NI2
      GO TO 240
220   IF (NI1.EQ.0.OR.NI2.EQ.0) GO TO 280
      II1=-NI1
      II2=-NI2
      GO TO 240
230   IF (NI2.EQ.0) GO TO 280
      II1=NI1
      II2=-NI2
240   CONTINUE
      GT1=FLOAT(II1)*RBR1(1)+FLOAT(II2)*RBR2(1)
      GT2=FLOAT(II1)*RBR1(2)+FLOAT(II2)*RBR2(2)
C     WRITE(13,245) II1,II2,(GT1*GT1+GT2*GT2)                           140291
C 245 FORMAT(' BEAM ',I3,I3,' LEN = ',F5.2)                             140291
      IF ((GT1*GT1+GT2*GT2).GT.GMAX2) GO TO 280
      KNT=KNT+1
      IF (KNT.GE.KNBMAX) GO TO 260
      NOP=1
      IF (NI1.EQ.0) NIT2=NI2                                            140291
250   G(1,KNT)=GT1
      G(2,KNT)=GT2
      SPQF(1,KNT)=ALMR(1,1)*II1+ALMR(2,1)*II2
      SPQF(2,KNT)=ALMR(1,2)*II1+ALMR(2,2)*II2
280   CONTINUE
290   IF (NI2.LE.NIT) GO TO 190                                         280679
      IF (NOP.EQ.1) GO TO 180                                           280679
      GO TO 295
260   WRITE(6,270)
270   FORMAT(23H MORE BEAMS THAN KNBMAX)
      STOP                                                              280679
295   CONTINUE
C     WRITE(13,7) NIT                                                   140291
C     WRITE(13,8) NIT2                                                  140291
C   7 FORMAT( 'NIT = ', I3)                                             140291
C   8 FORMAT( 'NIT2= ', I3)                                             140291
C  ORDER BEAMS BY INCREASING ABS(G)
      KNT1=KNT-1
      DO 320 I=1,KNT1
      AM=G(1,I)*G(1,I)+G(2,I)*G(2,I)
      I1=I+1
      DO 310 J=I1,KNT
      AM1=G(1,J)*G(1,J)+G(2,J)*G(2,J)
      IF (AM1.GE.AM) GO TO 310
      ST(1)=G(1,J)
      ST(2)=G(2,J)
      ST(3)=SPQF(1,J)
      ST(4)=SPQF(2,J)
      DO 300 KK=I1,J
      K=J+I1-KK
      SPQF(1,K)=SPQF(1,K-1)
      SPQF(2,K)=SPQF(2,K-1)
      G(1,K)=G(1,K-1)
300   G(2,K)=G(2,K-1)
      G(1,I)=ST(1)
      G(2,I)=ST(2)
      SPQF(1,I)=ST(3)
      SPQF(2,I)=ST(4)
      AM=AM1
310   CONTINUE
320   CONTINUE
C
C  ORDER BEAMS BY BEAM SET
      TWPI=2.*3.1415926535
      I=1
      KNBS=1
330   NBSET(I)=KNBS
      KNB(KNBS)=1
      JL=I+1
      IF (JL.GT.KNT) GO TO 370
      DO 360 J=JL,KNT
      DG(1)=G(1,J)-G(1,I)
      DG(2)=G(2,J)-G(2,I)
C  TEST WHETHER BEAMS I AND J BELONG TO SAME SUBSET. IF THEY DO,
C  BRING THEM TOGETHER IN THE LIST
      B=ABS(DG(1)*ARA1(1)+DG(2)*ARA1(2))+0.001
      B=AMOD(B,TWPI)/TWPI
      IF (ABS(B).GE.0.002) GO TO 360
      B=ABS(DG(1)*ARA2(1)+DG(2)*ARA2(2))+0.001
      B=AMOD(B,TWPI)/TWPI
      IF (ABS(B).GE.0.002) GO TO 360
C  IF J=I+1 NO REORDERING NEEDED
      IF (J.EQ.I+1) GO TO 350
      ST(1)=G(1,J)
      ST(2)=G(2,J)
      ST(3)=SPQF(1,J)
      ST(4)=SPQF(2,J)
      I2=I+2
      DO 340 KK=I2,J
      K=J+I2-KK
      SPQF(1,K)=SPQF(1,K-1)
      SPQF(2,K)=SPQF(2,K-1)
      G(1,K)=G(1,K-1)
340   G(2,K)=G(2,K-1)
      G(1,I+1)=ST(1)
      G(2,I+1)=ST(2)
      SPQF(1,I+1)=ST(3)
      SPQF(2,I+1)=ST(4)
350   I=I+1
      KNB(KNBS)=KNB(KNBS)+1
      NBSET(I)=KNBS
      IF (I.EQ.KNT) GO TO 370
360   CONTINUE
      I=I+1
      IF (I.EQ.KNT) GO TO 370
      KNBS=KNBS+1
      GO TO 330
370   CONTINUE
      CALL RMPRT(G,2,KNT,4HGMAT)
      CALL RMPRT(SPQF,2,KNT,4HSPQF)
      CALL IMPRT(KNB,1,KNBS,3HKNB)
      CALL IMPRT(NBSET,1,KNT,5HNBSET)
C  SYMMETRIZE, I.E. DETERMINE SYMMETRY CODE NUMBERS KSYM
C  FOR EACH BEAM, REMOVING THOSE BEAMS SYMMETRICAL TO OTHERS.
C  KSYM(1,I) APPLIES TO SUPERLATTICE REGION OF SURFACE, KSYM(2,I)
C  TO SUBSTRATE LAYERS (SEE LEED PROGRAM). (0,0) BEAM IS ASSUMED
C  IN THE FIRST PLACE IN THE LIST
      KSYM(1,1)=1
      KSYM(2,1)=1
      DO 380 II=2,KNT
      KSYM(1,II)=SSYM
380   KSYM(2,II)=SSYM
      IF (SSYM.EQ.1) GO TO 445
      DO 440 II=2,KNT
      IF (KSYM(2,II).EQ.0) GO TO 440
      IF (KSYM(1,II).EQ.1) GO TO 440
      DO 382 JD=II,KNT                                                  210492
         JDI=JD-II+1                                                    210492
         GDUMMY(1,JDI)=G(1,JD)                                          210492
         GDUMMY(2,JDI)=G(2,JD)                                          210492
382   CONTINUE                                                          210492
C     CALL UNFOLD(1,GU,NG,1,G(1,II),SSYM,KNT)                           210492
      SSYMA(1,1)=SSYM                                                   210492
      CALL UNFOLD(1,GU,NG,1,GDUMMY,SSYMA,KNT)                           210492
C  MARK SYMMETRICAL BEAMS FOR DELETION
      CALL DELBM(G,KNT,II,GU,NG,KSYM,NBSET,NSSET,SSYM,NBR)
      IF (NSSET.EQ.0) GO TO 432
C  DETERMINE KSYM FOR REDUCED SETS (I.E. WHEN SYMMETRICAL
C  BEAMS ARE FOUND IN DIFFERENT SUBSETS)
      ILST=LSYMN(SSYM)
      DO 420 ILS=1,ILST
      LSYMS=LSYMT(ILS)
      DO 384 JD=II,KNT                                                  210492
         JDI=JD-II+1                                                    210492
         GDUMMY(1,JDI)=G(1,JD)                                          210492
         GDUMMY(2,JDI)=G(2,JD)                                          210492
384   CONTINUE                                                          210492
C     CALL UNFOLD(1,GS,NGS,1,G(1,II),LSYMS,KNT)                         210492
      LSYMSA(1,1)=LSYMS                                                 210492
      CALL UNFOLD(1,GS,NGS,1,GDUMMY,LSYMSA,KNT)                         210492
      IF (NGS.NE.NBR) GO TO 420
      ILSN=0
      I100=II
      DO 410 IISS=1,NGS
      DO 390 IISS1=1,NGS
      IF ((ABS(G(1,I100)-GS(1,IISS1))+ABS(G(2,I100)-GS(2,IISS1)))
     1.GE.0.002) GO TO 390                                              280679
      ILSN=ILSN+1
      GO TO 400
390   CONTINUE
400   I100=I100+1
      IF (I100.GT.KNT) GO TO 410
      IF (KSYM(2,I100).NE.100) GO TO 400
410   CONTINUE
      IF (ILSN.EQ.NGS) GO TO 430
420   CONTINUE
      WRITE(6,425)
425   FORMAT(31H NO MATCH IN SYMMETRY OF SUBSET)
      STOP
430   KSYM(2,II)=LSYMS
432   DO 435 J=II,KNT
435   IF (KSYM(2,J).EQ.100) KSYM(2,J)=0
440   CONTINUE
445   CONTINUE
      CALL IMPRT(KSYM,2,KNT,4HKSYM)
C  COMPRESS THE LISTS G, SPQF, KSYM AND COUNT BEAMS PER BEAM SET
      DO 447 I=1,KNBS
447   KNB(I)=0
      KNTS=1
      DO 460 II=1,KNT
      IF (KSYM(2,II).EQ.0) GO TO 460
      NBSET(KNTS)=NBSET(II)
      N=NBSET(KNTS)
      KNB(N)=KNB(N)+1
      DO 450 I=1,2
      G(I,KNTS)=G(I,II)
      SPQF(I,KNTS)=SPQF(I,II)
450   KSYM(I,KNTS)=KSYM(I,II)
      KNTS=KNTS+1
460   CONTINUE
C  COMPRESS KNB
      KNBT=KNBS
      KNBS=1
      DO 470 I=1,KNBT
      IF (KNB(I).EQ.0) GO TO 470
      KNB(KNBS)=KNB(I)
      KNBS=KNBS+1
470   CONTINUE
      KNBS=KNBS-1
      KNT=KNTS-1
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE DELBM(G,KNT,II,GU,NG,KSYM,NBSET,NSSET,SSYM,NBR)
      INTEGER SSYM
      DIMENSION G(2,KNT),GU(2,NG),KSYM(2,KNT),NBSET(KNT)
C  TEST FOR DEGENERACIES. IF A BEAM LIES ON A SYMMETRY PLANE
C  ITS SYMMETRY CODE NUMBER SHOULD OFTEN BE DIFFERENT FROM SSYM,
C  BECAUSE THE SET OF BEAMS SYMMETRICAL TO THAT BEAM MAY CONTAIN
C  DEGENERACIES
      KSYMT=SSYM
      GO TO (160,10,20,30,40,50,10,60,70,80,10,90,120,10,140),SSYM
10    IF ((ABS(G(1,II))+ABS(G(2,II))).GT.0.001) GO TO 160
      KSYMT=1
      GO TO 160
20    IF (ABS(G(1,II)+G(2,II)).GE.0.001) GO TO 160
      KSYMT=1
      GO TO 160
30    IF (ABS(G(1,II)-G(2,II)).GE.0.001) GO TO 160
      KSYMT=1
      GO TO 160
40    IF (ABS(G(2,II)).GE.0.001) GO TO 160
      KSYMT=1
      GO TO 160
50    IF (ABS(G(1,II)).GE.0.001) GO TO 160
      KSYMT=1
      GO TO 160
60    IF (ABS(ABS(G(2,II))-ABS(G(1,II))).GE.0.001) GO TO 160
      KSYMT=2
      GO TO 160
70    IF (ABS(G(1,II)).GE.0.001.AND.ABS(G(2,II)).GE.0.001) GO TO 160
      KSYMT=2
      GO TO 160
80    IF (ABS(ABS(G(2,II))-ABS(G(1,II))).GE.0.001.AND.
     1ABS(G(1,II)).GE.0.001.AND.ABS(G(2,II)).GE.0.001) GO TO 160
      KSYMT=7
      GO TO 160
90    DO 110 I=1,NG
      IF (ABS(GU(2,I)).GE.0.001) GO TO 110
      KSYMT=11
110   CONTINUE
      GO TO 160
120   DO 130 I=1,NG
      IF (ABS(GU(1,I)).GE.0.001) GO TO 130
      KSYMT=11
130   CONTINUE
      GO TO 160
140   DO 150 I=1,NG
      IF (ABS(GU(1,I)).GE.0.001.AND.ABS(GU(2,I)).GE.0.001) GO TO 150
      KSYMT=14
150   CONTINUE
160   CONTINUE
      KSYM(1,II)=KSYMT
      KSYM(2,II)=KSYMT
      NSSET=0
      NBR=1
  225 FORMAT(A9,/,A5,I3,/,A9)                                            210492
C     WRITE (9,225) '*********',' II= ',II,'*********'                   210492
      II1=II+1
      IF (NG.EQ.1.OR.II1.GT.KNT) RETURN
      DO 180 I=2,NG
      DO 170 J=II1,KNT
      IF (KSYM(2,J).EQ.0.OR.KSYM(2,J).EQ.100) GO TO 170
      RDHELP=ABS(G(1,J)-GU(1,I))+ABS(G(2,J)-GU(2,I))                     210492
  222 FORMAT(A3,I3,A5,I3)                                                210492
  223 FORMAT(A11,I3)                                                     210492
  224 FORMAT(A11,F8.5,A11,F8.5)                                          210492
C     IF (II.EQ.15) THEN                                                 210492
C     WRITE (9,222)                                                      210492
C     WRITE (9,222) 'I= ',I,'  J= ',J                                    210492
C     WRITE (9,224) '   G(1,J)= ',G(1,J),'   G(2,J)= ',G(2,J)            210492
C     WRITE (9,224) '  GU(1,I)= ',GU(1,I),'  GU(2,I)= ',GU(2,I)          210492
C     WRITE (9,223) 'KSYM(2,J)= ',KSYM(2,J)                              210492
C     WRITE (9,224) '   RDHELP= ',RDHELP                                 210492
C     ENDIF                                                              210492
C     IF ((ABS(G(1,J)-GU(1,I))+ABS(G(2,J)-GU(2,I))).GE.0.001)            210492
C    1GO TO 170                                                          210492
      IF (RDHELP.GE.0.001) GO TO 170                                     210492
      NBR=NBR+1
C  KSYM(2,J)=100 IDENTIFIES THOSE SYMMETRICAL BEAMS THAT ARE IN THE
C  SAME SUBSET AS G(K,II)
      KSYM(2,J)=100
      IF (NBSET(J).EQ.NBSET(II)) GO TO 180
C  SYMMETRICAL BEAMS IN DIFFERENT SUBSETS ARE ENCOUNTERED
      NSSET=1
      NBR=NBR-1
      KSYM(2,J)=0
      GO TO 180
170   CONTINUE
180   CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE UNFOLD DETERMINES WHICH BEAMS ARE RELATED TO A GIVEN BEAM
C  BY A TYPE OF SYMMETRY SPECIFIED IN SYM.
C   IG= SERIAL NO. OF THE BEAM WHOSE SYMMETRICAL COUNTERPARTS ARE SOUGHT
C   G= OUTPUT LIST OF SYMMETRICAL BEAMS SOUGHT.
C   NG= OUTPUT NUMBER OF SYMMETRICAL BEAMS IN G.
C   LAY= INDEX FOR CHOICE OF SYMMETRY CODES SYM(LAY,I) (NORMALLY LAY=1
C    FOR OVERLAYER, LAY=2 FOR SUBSTRATE LAYERS).
C   PQ= LIST OF INPUT BEAMS G.
C   SYM= LIST OF SYMMETRY CODES FOR ALL BEAMS.
C   NT= TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
      SUBROUTINE  UNFOLD (IG, G, NG, LAY, PQ, SYM, NT)
      INTEGER  SYM
      DIMENSION  G(2,12), PQ(2,NT), SYM(2,NT)
      ACU=1.7320508
      DO 340 I = 1, 2
  340 G(I,1) = PQ(I,IG)
      LAB = SYM(LAY,IG)
C     DELETION OF ERROR                                                  210492
C     IF ((LAB .EQ. 10) .AND. (G(2,1) .LT. 0.0))  G(2,1) =  - G(2,1)     210492
      GO TO (350,360,380,400,420,430,440,460,480,500,
     1511,512,513,514,515), LAB
  350 NG = 1
      GO TO 520
  360 NG = 2
      DO 370 I= 1, 2
  370 G(I,2) =  - G(I,1)
      GO TO 520
  380 NG = 2
      DO 390 I = 1, 2
  390 G(I,2) =  - G(3 - I,1)
      GO TO 520
  400 NG = 2
      DO 410 I = 1, 2
  410 G(I,2) = G(3 - I,1)
      GO TO 520
  420 NG = 2
      G(1,2) = G(1,1)
      G(2,2) =  - G(2,1)
      GO TO 520
  430 NG = 2
      G(1,2) =  - G(1,1)
      G(2,2) = G(2,1)
      GO TO 520
  440 NG = 4
      G(1,2) =  - G(2,1)
      G(2,2) = G(1,1)
      DO 450 I = 1, 2
      G(I,3) =  - G(I,1)
  450 G(I,4) =  - G(I,2)
      GO TO 520
  460 NG = 4
      DO 470 I = 1, 2
      G(I,2) = G(3 - I,1)
      G(I,3) =  - G(I,1)
  470 G(I,4) =  - G(I,2)
      GO TO 520
  480 NG = 4
      G(1,4) = G(1,1)
      G(2,4) =  - G(2,1)
      DO 490 I = 1, 2
      G(I,2) =  - G(I,4)
  490 G(I,3) =  - G(I,1)
      GO TO 520
  500 NG = 8
      G(1,8) = G(1,1)
      G(2,8) =  - G(2,1)
      DO 510 I = 1, 2
      G(I,2) = G(3 - I,1)
      G(I,3) = G(3 - I,8)
      G(I,4) =  - G(I,8)
      G(I,5) =  - G(I,1)
      G(I,6) =  - G(I,2)
  510 G(I,7) =  - G(I,3)
      GO TO 520
  511 NG=3
      G(1,12)=-ACU*G(2,1)
      G(2,12)=ACU*G(1,1)
      DO 521 I=1,2
      G(I,2)=0.5*(-G(I,1)+G(I,12))
  521 G(I,3)=0.5*(-G(I,1)-G(I,12))
      GO TO 520
  512 NG=6
      G(1,12)=-ACU*G(2,1)
      G(2,12)=ACU*G(1,1)
      DO 522 I=1,2
      G(I,3)=0.5*(-G(I,1)+G(I,12))
  522 G(I,5)=0.5*(-G(I,1)-G(I,12))
      DO 532 I=2,6,2
      G(1,I)=G(1,7-I)
  532 G(2,I)=-G(2,7-I)
      GO TO 520
  513 NG=6
      G(1,12)=-ACU*G(2,1)
      G(2,12)=ACU*G(1,1)
      DO 523 I=1,2
      G(I,3)=0.5*(-G(I,1)+G(I,12))
  523 G(I,5)=0.5*(-G(I,1)-G(I,12))
      G(1,2)=-G(1,1)
      G(2,2)=G(2,1)
      DO 533 I=4,6,2
      G(1,I)=-G(1,9-I)
  533 G(2,I)=G(2,9-I)
      GO TO 520
  514 NG=6
      G(1,12)=-ACU*G(2,1)
      G(2,12)=ACU*G(1,1)
      DO 524 I=1,2
      G(I,2)=0.5*(G(I,1)+G(I,12))
      G(I,3)=0.5*(-G(I,1)+G(I,12))
      G(I,4)=-G(I,1)
      G(I,5)=-G(I,2)
  524 G(I,6)=-G(I,3)
      GO TO 520
  515 NG=12
      G(1,12)=-ACU*G(2,1)
      G(2,12)=ACU*G(1,1)
      DO 525 I=1,2
      G(I,3)=0.5*(G(I,1)+G(I,12))
      G(I,5)=0.5*(-G(I,1)+G(I,12))
      G(I,7)=-G(I,1)
      G(I,9)=-G(I,3)
  525 G(I,11)=-G(I,5)
      DO 535 I=2,12,2
      G(1,I)=G(1,13-I)
  535 G(2,I)=-G(2,13-I)
  520 RETURN
      END
C-------------------------------------------------------------------
C  SUBROUTINE RMPRT PRINTS A REAL MATRIX R(N1,N2), WITH TITLE (GIVEN    121278
C  AS A HOLLERITH STRING) AND ROW AND COLUMN LABELING.
C  FOR NON-REAL MATRICES SEE SUBROUTINES CMPRT AND IMPRT.
      SUBROUTINE RMPRT(R,N1,N2,TITLE)
      DIMENSION R(N1,N2)
      DIMENSION IND(99999)
      N=MAX0(N1,N2)
      DO 5 J=1,N
    5 IND(J)=J
      WRITE(6,10)TITLE
   10 FORMAT(1H0,1A8)
      J2=0
   15 J1=J2+1
      IF (J1.GT.N2) GO TO 30
      J2=J1+9
      J2=MIN0(J2,N2)
      WRITE(6,17)(IND(J),J=J1,J2)
   17 FORMAT(8X,I3,9(9X,I3))
      DO 20 K=1,N1
   20 WRITE(6,25)IND(K),(R(K,J),J=J1,J2)
   25 FORMAT(1H ,I3,10E12.5)
      GO TO 15
   30 RETURN
      END
C-----------------------------------------------------------------------------
C  SUBROUTINE IMPRT PRINTS A INTEGER MATRIX I(N1,N2), WITH TITLE (GIVEN 121278
C  AS A HOLLERITH STRING) AND ROW AND COLUMN LABELING.
C  FOR NON-INTEGER MATRICES SEE SUBROUTINES CMPRT AND RMPRT.
      SUBROUTINE IMPRT(I,N1,N2,TITLE)
      DIMENSION I(N1,N2)
      DIMENSION IND(99999)
      N=MAX0(N1,N2)
      DO 5 J=1,N
    5 IND(J)=J
      WRITE(6,10)TITLE
   10 FORMAT(1H0,1A8)
      J2=0
   15 J1=J2+1
      IF (J1.GT.N2) GO TO 30
      J2=J1+24
      J2=MIN0(J2,N2)
      WRITE(6,17)(IND(J),J=J1,J2)
   17 FORMAT(4X,I3,24(2X,I3))
      DO 20 K=1,N1
   20 WRITE(6,25)IND(K),(I(K,J),J=J1,J2)
   25 FORMAT(1H ,I3,25I5)
      GO TO 15
   30 RETURN
      END
CHL
C*****************************************************************************
C  SUBROUTINE    N U M B E R
C  SUBROUTINE NUMBER NUMERIERT DIE BEAMS IN AUFSTEIGENDER REIHENFOLGE UND
C  ERZEUGT EINE NUMERIERTE BEAMLISTE AUF DEM FILE NBLIST
C
C
      SUBROUTINE NUMBER
C  DEKLARATIONEN
      CHARACTER TEXT*50
C  DEFINITION DER EIN-/AUSGABE-FILES
      CLOSE(7)
      OPEN(8,FILE='BELIST',STATUS='UNKNOWN')
      OPEN(9,FILE='NBLIST',STATUS='UNKNOWN')
      REWIND(8)
C  ERZEUGUNG DER NUMERIERTEN BEAM-LISTE
      NR=1
1     CONTINUE
         READ(8,10,END=3) NBEAM
         WRITE(9,10) NBEAM
         DO 2 IBEAM=1,NBEAM
            READ(8,20,END=3) TEXT
            WRITE(9,30) TEXT,NR
            NR=NR+1
2        CONTINUE
      GOTO1
10    FORMAT(I5)
20    FORMAT(A50)
30    FORMAT(A50,2X,'NR.',I5)
3     CONTINUE
      RETURN
      END
