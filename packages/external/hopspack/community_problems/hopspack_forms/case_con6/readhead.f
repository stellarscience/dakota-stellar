

      SUBROUTINE readhead(NROW,NCOL,NLAY,HEAD,HEDFILE)
C***************************************************************
! WARNING: You must know how the head file is formatted to be 
!          read in properly
!*************************************************************
      IMPLICIT NONE
  
      INTEGER I,J,K,NROW,NCOL,NLAY
      REAL*8 HEAD(NROW,NCOL,NLAY)
      CHARACTER*80 HEDFILE

C***********************************************************
C  Modification by GG 08/04/03 to get stand alone executable
C    name of .hed file passed in instead of using mype
C     
!      if(mype<10)then
!      write(HEADFILE,'(a2,i1,a4)') 'A1',mype,'.hed'
!      endif
!      if(mype>9)then
!      write(HEADFILE,'(a2,i2,a4)') 'A1',mype,'.hed'
!      endif
C***********************************************************

      OPEN(UNIT = 30, FILE = HEDFILE, STATUS = 'OLD')
      DO 10 , K=1,NLAY
      DO 20 , I=1,NROW
      READ(30,203) (HEAD(I,J,K), J=1,10)
      READ(30,203) (HEAD(I,J,K), J=11,20)
      READ(30,203) (HEAD(I,J,K), J=21,30)
      READ(30,203) (HEAD(I,J,K), J=31,40)
      READ(30,203) (HEAD(I,J,K), J=41,50)

 20   CONTINUE
 10   CONTINUE
 203  FORMAT(10(F10.3))
      CLOSE(30)
      END

