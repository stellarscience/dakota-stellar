
      SUBROUTINE head_constraint(DIM,X,Y,Z,WELRATE,HEAD,NLAY,NROW,
     & NCOL,COSTO,FLAG)
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
      REAL D
      REAL TEMP2,TEMP3,TEMP4,TEMP5,TF,QM(DIM)
      REAL*8 HEAD(NROW,NCOL,NLAY), HMIN, HMAX,DRWDOWN(DIM)


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Define constant coefficients, etc...
      D = 0.0001
      ZGS = 30.0
      TEMP = 0.0
      TEMP2 = 0.0
      TEMP3 = 0.0
      TEMP4 = 0.0

      HMIN= 10.0
      HMAX = ZGS
      MINRATE = 0.000001


C Check the constraints
      DO 20 , I=1,DIM
      IF(ABS(WELRATE(I)) .GT. MINRATE)THEN
C
C Constraints on head values in wells*********************************
         TEMP5 = HEAD(X(I),Y(I),Z(I))-HMIN
         IF(TEMP5.LT.0.0)THEN
         WRITE(*,*) 'Min head violated',HEAD(X(I),Y(I),Z(I)),X(I),Y(I),Z(I),I
         FLAG = 2
         ENDIF
C      

         TEMP5 = HMAX-HEAD(X(I),Y(I),Z(I))
         IF(TEMP5.LT.0.0)THEN
         FLAG = 2
         WRITE(*,*) 'Max head violated', HEAD(X(I),Y(I),Z(I)), X(I),Y(I), I
         ENDIF

        ENDIF
 20     CONTINUE

*********************************************************************
C Constraints on head difference
        TEMP5 = HEAD(27,18,4)-HEAD(28,18,4)
C        WRITE(*,*) HEAD(27,18,4),HEAD(28,18,4)
C        WRITE(*,*) '1',TEMP5
        IF(TEMP5 .LT. D)THEN
        FLAG = 2
        WRITE(*,*) 'gradient constraint 1 violated', TEMP5
        ENDIF

        TEMP5 = HEAD(23,24,4)-HEAD(24,24,4)
C        WRITE(*,*) '2',TEMP5
        IF(TEMP5 .LT. D)THEN
        FLAG = 2
        WRITE(*,*) 'gradient constraint 2 violated', TEMP5
        ENDIF

        TEMP5 = HEAD(26,33,4)-HEAD(27,33,4)
C        WRITE(*,*) '3',TEMP5
        IF(TEMP5 .LT. D)THEN
        FLAG = 2
        WRITE(*,*) 'gradient constraint 3 violated', TEMP5
        ENDIF

        TEMP5 = HEAD(35,39,4)-HEAD(35,38,4)
C        WRITE(*,*) '4',TEMP5
        IF(TEMP5 .LT. D)THEN
        FLAG = 2
        WRITE(*,*) 'gradient constraint 4 violated', TEMP5
        ENDIF

        TEMP5 = HEAD(46,38,4)-HEAD(46,37,4)
C        WRITE(*,*) '5',TEMP5
        IF(TEMP5 .LT. D)THEN
        FLAG = 2
        WRITE(*,*) 'gradient constraint 5 violated', TEMP5
        ENDIF
**********************************************************

     
      END




