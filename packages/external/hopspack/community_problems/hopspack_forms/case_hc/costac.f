

      SUBROUTINE costac(DIM,X,Y,Z,WELRATE,NLAY,NROW,NCOL,COSTC,FLAG)
      IMPLICIT NONE

      INTEGER DIM,NLAY,NROW,NCOL,I,FLAG
      INTEGER X(DIM), Y(DIM), Z(DIM)
      DOUBLE PRECISION WELRATE(DIM),ZGS,HMIN, TEMP,COSTC
      DOUBLE PRECISION TEMP2 
      DOUBLE PRECISION QM(DIM), DEPTH(DIM)
      REAL*8 DRWDOWN(DIM), MINRATE

C Exponents
      REAL*8 B0,B1,B2
C Coefficients 
      REAL*8 C0,C1
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C define constant coefficients, exponents, etc...
     
      B0 = 0.3
      B1 = 0.45
      B2 = 0.64
      C0 = 5500.0
      C1 = 5750.0
      ZGS = 30.0 
      hmin = 10.0 
      TEMP = 0.0
      COSTC = 0.0
      TEMP2 = (ZGS-HMIN)**B2
      MINRATE = 0.000001
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C define the depth of each well
      DO 10 , I=1,DIM
      DEPTH(I) = ZGS
      QM(I) = 1.5*WELRATE(I)
 10   CONTINUE
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
C Compute cost to drill the wells
      DO 20 , I=1,DIM
         IF(ABS(WELRATE(I)) .GT. MINRATE)THEN
C         WRITE(*,*) 'Well', I, 'is ON'
         TEMP = TEMP + C0*DEPTH(I)**B0
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C If the well is extracting, add on the cost of a pump
         IF(WELRATE(I) .LT. 0.0)THEN
         TEMP = TEMP + C1*(ABS(QM(I))**B1)*TEMP2
C         WRITE(*,*) '-',I
         ENDIF
         ENDIF
 20      CONTINUE
      COSTC = TEMP 
      END

