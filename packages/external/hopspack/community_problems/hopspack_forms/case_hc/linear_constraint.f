
      SUBROUTINE linear_constraint(DIM,X,Y,Z,WELRATE,NLAY,
     & NROW, NCOL,FLAG)
      IMPLICIT NONE
C****************************************************************
! Checks the bound on the net pumping rate
! Sets FLAG = 1 if violated
! description of variables:
! MINRATE is the minimum allowable pumping rate before we
!         consider a well shut off
! QTMIN is the constraint on the net  rate 
! QT is the net rate
! TEMP(2-5) are used in computing the penalty terms
!**************************************************************
      INTEGER DIM,NLAY,NROW,NCOL,I,FLAG
      INTEGER X(DIM), Y(DIM), Z(DIM)
      DOUBLE PRECISION WELRATE(DIM),TEMP,MINRATE
      REAL D, QT, QTMAX


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Define constant coefficients, etc...
      QT = 0.0
      QTMAX = -0.031998
      MINRATE = 0.000001


C Check the constraint on the net pumping rate
      DO 20 , I=1,DIM
C       WRITE(*,*) X(I),Y(I),WELRATE(I)
       IF(ABS(WELRATE(I)) .GT. MINRATE)THEN
          QT = QT + WELRATE(I)
       ENDIF
 20   CONTINUE

        IF((QTMAX-QT).GT. 0.0)THEN
           WRITE(*,*) 'Net rate violated',QTMAX, QT 
           FLAG = 1
        ENDIF
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


     
      END




