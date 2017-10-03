
      SUBROUTINE location_constraint(DIM,X,Y,FLAG)
      IMPLICIT NONE
C****************************************************************
! Checks that the wells don't occupy the same location
! Sets FLAG = 1 if violated
! description of variables:
!**************************************************************
      INTEGER DIM,NLAY,FLAG,I,J
      INTEGER X(DIM), Y(DIM)


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Check to make sure 2 wells do not occupy the same location
      DO 669 , I=1,DIM
        DO 670 , J=I+1,DIM
           IF((X(I) .EQ. X(J)) .AND. (Y(I) .EQ. Y(J)))THEN 
           WRITE(*,*) 'Location violated', X(I), Y(I), X(J), Y(J)
           FLAG = 1
           ENDIF
 670       CONTINUE
 669       CONTINUE


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


     
      END




