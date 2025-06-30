module tenserleed_intpol
      use, intrinsic :: iso_fortran_env
      integer, parameter :: tdp = REAL64
    contains
!C-------------------------------------------------------------------------------
!C  Subroutine INTPOL interpolates LEED intensities onto a working grid (with
!C  steps of EINCR eV).

      SUBROUTINE INTPOL(A,NE,E,NGP,EINCR,IPR,X,WORYT,IB)

        REAL(tdp) :: A(NGP), E(NGP), EINCR, X(NGP), WORYT(NGP), XMIN, XVAL
        INTEGER NE, IPR, IB, ITIL, ITIH, IMIN, LMIN, LMAX
  
   10   FORMAT(54H0 IN THIS BEAM TOO FEW ENERGY VALUES FOR INTERPOLATION)
   20   FORMAT(40H*INTENSITIES AFTER INTERPOLATION IN BEAM,1I4,&
        
   100(5(1F7.2,1E13.4,3X),/))
  
        ITIL = 0
        ITIH = 0
  
  !C  find first non-zero intensity and skip beam, if only one or less are found
  
        DO IE = 1,NE
          IMIN = IE
          IF (A(IE).GT.1.E-8) GO TO 100
        ENDDO
   100  CONTINUE
  
        IF (IMIN.EQ.NE) THEN                                               !200480
          NE = 0                                                           !200480
          WRITE(6,10)                                                     !230480
          RETURN
        END IF
  
  !C  calculate energy range in units of EINCR
  
        LMIN = INT((E(IMIN)-0.0001) / EINCR) + 1
        LMIN = MAX0(LMIN,0)
        XMIN = FLOAT(LMIN) * EINCR
        LMAX = INT((E(NE)+0.0001) / EINCR)
  
        DO I = IMIN,NE
          X(I-IMIN+1) = E(I)
          WORYT(I-IMIN+1) = A(I)
        ENDDO
        NEM = NE - IMIN + 1
  
  !C  NPTS is no. of working grid points
  
        NPTS = LMAX - LMIN + 1
        NE = NPTS
  
        XVAL = XMIN - EINCR

  
  !C  interpolate intensities (negative intensities are set to zero)
  
        DO I = 1,NPTS
          XVAL = XVAL + EINCR
          A(I) = YVAL(XVAL,WORYT,X,NEM,ITIL,ITIH)
          IF (A(I).LT.0.0) THEN
            A(I)=0.0
          ENDIF
        ENDDO
  
  !C  calculate energy grid
  
        E(1) = XMIN
        DO IE = 2,NPTS
          E(IE) = E(IE-1) + EINCR
        ENDDO
  
  !C  print working grid and interpolated intensities, if desired
  
        IF (IPR.EQ.2.AND.NE.GT.0) THEN
          WRITE(6,20) IB, (E(IE),A(IE),IE=1,NE)
        ENDIF
  
        RETURN
        END
  !C-----------------------------------------------------------------------
  
  
!  C-----------------------------------------------------------------------
!  C  FUNCTION YVAL INTERPOLATES
      real(tdp) FUNCTION YVAL(X, WORY, WORX, LENGTH,ITIL,ITIH)
      real(tdp) :: WORY(LENGTH), WORX(LENGTH)
      real(tdp) :: X, Y
      COMMON /DATBLK/ X0,X1,X2,X3,Y0,Y1,Y2,Y3
!  C  IF FEWER THAN FOUR GRID POINTS AVAILABLE, USE 2ND OR 1ST ORDER
!  C  POLYNOMIAL INTERPOLATION
      IF (LENGTH.LT.4) GO TO 10
!  C  FIND REQUIRED INTERPOLATION INTERVAL
      CALL BINSRX(IL, IH, X, WORX, LENGTH)
!  C  SKIP NEXT STEP IF SAME INTERVAL IS FOUND AS LAST TIME
      IF(IL .EQ. ITIL .AND. IH .EQ. ITIH) GO TO 5
!  C  FIND FOUR NEAREST GRID POINTS AND CORRESPONDING INTENSITIES
!  C  FOR 3RD-ORDER POLYNOMIAL INTERPOLATION
      CALL STFPTS(IL, IH, WORX, WORY, LENGTH)
!  C  DO ACTUAL 3RD-ORDER POLYNOMIAL INTERPOLATION
  5    Y = XNTERP(X)
      ITIH = IH
      ITIL = IL
      YVAL = Y
      GO TO 20
  10    CALL XNTRP2(X,Y,WORX,WORY,LENGTH)
      YVAL=Y
  20    RETURN
      END
  
!  C-----------------------------------------------------------------------
!  C  SUBROUTINE BINSRX FINDS A REQUIRED INTERPOLATION INTERVAL
!  C  BY BINARY SEARCH (SUCCESSIVE HALVING OF INITIAL INTERVAL)
      SUBROUTINE BINSRX(IL, IH, X, WORX, LENGTH)
      real(tdp) :: WORX(LENGTH), X
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
!  C-----------------------------------------------------------------------
  
!  C------------------------------------------------------------------------
!  C  SUBROUTINE STFPTS FINDS, GIVEN THE INTERPOLATION INTERVAL, THE
!  C  FOUR NEAREST GRID POINTS AND THE CORRESPONDING ORDINATE VALUES
        SUBROUTINE STFPTS(IL, IH, WORX, WORY, LENGTH)
        real(tdp) :: WORX(LENGTH), WORY(LENGTH), TEMP(2, 4)
        COMMON / DATBLK/ X0, X1, X2, X3, Y0, Y1, Y2, Y3
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

  !C-----------------------------------------------------------------------
!C  SUBROUTINE XNTERP PERFORMS 3RD-ORDER POLYNOMIAL INTERPOLATION
      FUNCTION XNTERP(X)
            real(tdp) :: X
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
!C-----------------------------------------------------------------------
!C  SUBROUTINE XNTRP2 PERFORMS 2ND OR 1ST ORDER POLYNOMIAL INTERPOLATION
      SUBROUTINE XNTRP2(X,Y,XS,YS,N)
      real(tdp) :: XS(N),YS(N), X, Y
      IF (N.GT.2) GO TO 15
      Y=(YS(2)-YS(1))*(X-XS(1))/(XS(2)-XS(1))+YS(1)
      RETURN
    15    A=(YS(1)-YS(2))/(XS(1)-XS(2))/(XS(2)-XS(3))- &
     (YS(1)-YS(3))/(XS(1)-XS(3))/(XS(2)-XS(3))
      B=(YS(1)-YS(2))/(XS(1)-XS(2))-A*(XS(1)+XS(2))
      C=YS(1)-(A*XS(1)+B)*XS(1)
      Y=C+X*(B+A*X)
      RETURN
      END SUBROUTINE XNTRP2

end module tenserleed_intpol
