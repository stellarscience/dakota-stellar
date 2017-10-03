
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
      ZGS = 60.0
      TEMP = 0.0
      TEMP2 = 0.0
      TEMP3 = 0.0
      TEMP4 = 0.0

      HMIN= 40.0
      HMAX = ZGS
      MINRATE = 0.0001


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

     
      END




