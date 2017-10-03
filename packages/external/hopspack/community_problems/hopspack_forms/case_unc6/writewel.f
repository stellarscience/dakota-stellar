      SUBROUTINE writewel(DIM,X,Y,Z,WELRATE,welfile)
C*****************************************************************
! WARNING: You must hardwire in the appropriate well file name
! FILE = 'foo.wel'
! and know the format of the input file
C***************************************************************
C***********************************************************
C  Modification by GG 08/04/03 to get stand alone executable
C    name of the .wel file will be predetermined and passed into the
C    the write wel subroutine instead of passing in mype
C***************************************************************
      IMPLICIT NONE 
      
      INTEGER DIM,I
      DOUBLE PRECISION WELRATE(DIM)
      INTEGER X(DIM),Y(DIM),Z(DIM)
      CHARACTER*80 welfile


C      if(mype<10)then
C      write(welfile,'(a2,i1,a4)') 'A1',mype,'.wel' 
C      endif
C      if(mype>9)then
C      write(welfile,'(a2,i2,a4)') 'A1',mype,'.wel' 
C      endif

      OPEN(UNIT = 13, FILE = welfile, STATUS = 'UNKNOWN')
      WRITE(13,*) DIM,0
      WRITE(13,*) DIM
      DO 10 , I=1,DIM 
      WRITE(13,100)  Z(I),X(I),Y(I),WELRATE(I) 
 10   CONTINUE
C FORMATS
 100  FORMAT(I3,I3,I3,F12.7)
 101  FORMAT(I3)
      
      END
