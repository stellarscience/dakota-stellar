
      SUBROUTINE costao(DIM,X,Y,Z,WELRATE,HEAD,NLAY,NROW,NCOL,COSTO,
     & FLAG)
      IMPLICIT NONE
C****************************************************************
! description of variables:
! ZGS is the ground surface elevation
! COSTO is the operational cost (output)
! MINRATE is the minimum allowable pumping rate before we
!         consider a well shut off
! QTMIN is the constraint on the net  rate 
! QT is the net rate
! TEMP(2-5) are used in computing the penalty terms
! QM is the design pumping rate (this will change later)
! HMIN/HMAX are the minimum and maximum allowable heads
!**************************************************************
      INTEGER DIM,NLAY,NROW,NCOL,I,FLAG
      INTEGER X(DIM), Y(DIM), Z(DIM)
      DOUBLE PRECISION WELRATE(DIM),ZGS,TEMP,COSTO,MINRATE
      REAL QTMIN,QT
      REAL TEMP2,TEMP3,TEMP4,TEMP5,TF,QM(DIM)
      REAL*8 HEAD(NROW,NCOL,NLAY), HMIN, HMAX,DRWDOWN(DIM)

C Coefficients 
      REAL*8 C2,C3

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Define constant coefficients, etc...
      C2 = 0.00029
      C3 = 0.000145
      ZGS = 30.0
      TEMP = 0.0
      TEMP2 = 0.0
      TEMP3 = 0.0
      TEMP4 = 0.0
      QTMIN = -0.031998
      QT = 0.0
      HMIN= 10.0
      HMAX = ZGS
      MINRATE = 0.0001

C IMPORTANT: Set the final time
      TF = 1.56e8


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Compute the cost  

      DO 10 , I=1,DIM
      IF(ABS(WELRATE(I)).GT. MINRATE) THEN
C If the well is extracting....
        IF(WELRATE(I) .LT. 0.0)THEN
        TEMP = TEMP +  C2*WELRATE(I)*(HEAD(X(I),Y(I),Z(I))-ZGS)
C or if the well is injecting...
          ELSEIF(WELRATE(I) .GT. 0.0) THEN 
          TEMP = TEMP + C3*WELRATE(I)
      END IF
      ENDIF
 10   CONTINUE

C Multiply by the final time 
 
      TEMP = TF*TEMP

      COSTO = TEMP
     
      END




