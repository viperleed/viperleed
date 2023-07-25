C  Tensor LEED subroutines for intensity calculation from given delta amplitudes
C  v1.2, VB 13.04.00
C  for use with superpos.f v1.2
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
C  Please read the comment in superpos.f, v1.2 .
C
C-----------------------------------------------------------------------
C  SUBROUTINE INDELT READS THE SECOND PART OF FILEHEADER
C
       SUBROUTINE INDELT(IFILE,IFORM,NT0,NATOMS,NCSTEP,PQFEX,
     1                   CUNDISP,CDISP,AID)
C
      DIMENSION PQFEX(2,NT0)
      DIMENSION CDISP(NCSTEP,NATOMS,3),CUNDISP(NATOMS,3),AID(NCSTEP)
C
      IF (IFORM.EQ.0) THEN
         READ (IFILE) PQFEX,CUNDISP,CDISP,AID
      ELSE
            READ(IFILE,'(10F10.5)',ERR=100) 
     +      ((PQFEX(I,J), I=1,2), J=1,NT0)

            GO TO 110

 100        write(3,*) "Error reading v1.2 style output beam list",
     +                 " from delta amplitude file. Trying pre-v1.2"
            write(3,*) "format instead ..."
            
            READ(IFILE,'(10F7.4)') 
     +      ((PQFEX(I,J), I=1,2), J=1,NT0)

            write(3,*) "Read may have been successful. Check",
     +                 " carefully, anyhow!"

 110        CONTINUE

         READ (IFILE,120) ((CUNDISP(I,J),J=1,3),I=1,NATOMS)
         READ (IFILE,120) (((CDISP(I,J,K),K=1,3),J=1,NATOMS),I=1,NCSTEP)
         READ (IFILE,120) (AID(I),I=1,NCSTEP)
  120    FORMAT (10F7.4)
      ENDIF
      RETURN
      END
C------------------------------------------------------------------------
C  SUBROUTINE INRINT READS IN ALL DELTA AMP'S OF ONE FILE
C
      SUBROUTINE INRINT(IFILE,IFORM,E,VV,VO,VPI,NT0,NCSTEP,XIST,DELWV)
C
      COMPLEX XIST(NT0), DELWV(NCSTEP,NT0)                              11/11/90
C
      IF (IFORM.EQ.0) THEN
         READ (IFILE,ERR=200,END=200) E,VPI,VO,VV,XIST,DELWV
      ELSE
         READ (IFILE,100,ERR=200,END=200) E,VPI,VO,VV
         READ (IFILE,100) (XIST(I),I=1,NT0)
         READ (IFILE,100) ((DELWV(I,J),J=1,NT0),I=1,NCSTEP)
  100    FORMAT (6E13.7)
      ENDIF
      RETURN
C  SET CONDITION FOR TERMINATION OF CURRENT RUN:
  200 E = -1.
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE ADDWV ADDS THE ACTUAL DELTA AMP' TO PREVIOUS AMPLITUDE
C
      SUBROUTINE ADDWV(DELACT,DELWV,NT0,NCSTEP,N,IACT,CONC)
C
      COMPLEX DELACT,DELWV(NCSTEP,NT0)
C
      DELACT = DELACT + DELWV(IACT,N)*CONC
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE XMIN SEARCHES FOR THE MINIMAL X-DISPLACEMENT OF ALL ATOMS
C
      SUBROUTINE XMIN(CXDISP,CDISP,IFVAR,NATOMS,IACT)
C
      DIMENSION CDISP(IFVAR,NATOMS)
C
      DO 5 IATOM=1,NATOMS
           IF (CDISP(IACT,IATOM).LT.CXDISP) THEN
              CXDISP = CDISP(IACT,IATOM)
           ENDIF
    5 CONTINUE
C
      RETURN
      END
