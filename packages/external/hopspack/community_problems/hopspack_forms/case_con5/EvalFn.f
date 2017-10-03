      PROGRAM EvalFn

C************************************************************************* 
!     SUBROUTINE TO EVALUATE THE OBJECTIVE FUNCTION FOR APPS
!       
!     THIS PROGRAM WILL:
!       1. READ IN THE NUMBER OF PARAMETERS AND CURRENT ITERATE
!       2. CALL THE MODIFIED feval.f ROUTINE  
!       3. RETURN A FILE CONTAINING THE FUNCTION VALUE
C************************************************************************
!     Program created by Genetha Gray 08/04/03
C************************************************************************
C Description of variables:
!
C**********************************************************************
      CHARACTER*80 inputfile, outputfile
      CHARACTER*80 mfnfile,welfile,hedfile    
      CHARACTER*2  reqtype

      INTEGER DIM,FLAG
      INTEGER I

      DOUBLE PRECISION W(100),FV

C***** Get the command line arguements
      CALL GETARG(1,inputfile) 
      CALL GETARG(2,mfnfile) 
      CALL GETARG(3,welfile)
      CALL GETARG(4,hedfile) 
      CALL GETARG(5,outputfile) 

C***** Read in the dimension number and current iterate
C***** For HOPSPACK, also read the request type (which should be 'F')
      OPEN(UNIT=4, FILE=inputfile, STATUS='OLD')
      READ(4,*) reqtype
      READ(4,*) DIM
      DO I = 1,DIM
        READ(4,*) W(I)
      ENDDO
      CLOSE(4)

C***** Call the function evaluation routine
      call feval(DIM,W,welfile,mfnfile,hedfile,FV,FLAG)

C***** Check the return flag and write function value to a file (if needed)
C***** For HOPSPACK, also write the length of the vector
      IF(FLAG .EQ. 0)THEN      
        OPEN(UNIT=5, FILE=outputfile, STATUS='UNKNOWN')
        WRITE(5,*) 1
        WRITE(5,*) FV
        CLOSE(5)
      ENDIF

      IF(FLAG .EQ. 1)THEN      
        OPEN(UNIT=5, FILE=outputfile, STATUS='UNKNOWN')
        WRITE(5,*) 'C1/C2 Viol'
        CLOSE(5)
      ENDIF

      IF(FLAG .EQ. 2)THEN      
        OPEN(UNIT=5, FILE=outputfile, STATUS='UNKNOWN')
        WRITE(5,*) 'C3 Viol'
        CLOSE(5)
      ENDIF

      END
