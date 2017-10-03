

      SUBROUTINE readhead(NROW,NCOL,NLAY,HEAD,HEADFILE)
C***************************************************************
! WARNING: You must hardwire in the name of the head file
! 'blah.hed'
! and know how it is formatted to be read in properly
!*************************************************************
      IMPLICIT NONE
  
      INTEGER I,J,K,NROW,NCOL,NLAY,mype
      REAL*8 HEAD(NROW,NCOL,NLAY)
      CHARACTER*80 HEADFILE
      
C      if(mype<10)then
C      write(HEADFILE,'(a2,i1,a4)') 'A3',mype,'.hed'
C      endif
C      if(mype>9)then
C      write(HEADFILE,'(a2,i2,a4)') 'A3',mype,'.hed'
C      endif

      OPEN(UNIT = 30, FILE = HEADFILE, STATUS = 'OLD')
      DO 10 , K=1,NLAY
      DO 20 , I=1,NROW
      READ(30,203) (HEAD(I,J,K), J=1,10)
C      WRITE(*,203) (HEAD(I,J,K), J=1,10)
      READ(30,203) (HEAD(I,J,K), J=11,20)
C      WRITE(*,203) (HEAD(I,J,K), J=11,20)
      READ(30,203) (HEAD(I,J,K), J=21,30)
C      WRITE(*,203) (HEAD(I,J,K), J=21,30)
      READ(30,203) (HEAD(I,J,K), J=31,40)
C      WRITE(*,203) (HEAD(I,J,K), J=31,40)
      READ(30,203) (HEAD(I,J,K), J=41,50)
C      WRITE(*,203) (HEAD(I,J,K), J=41,50)
      READ(30,203) (HEAD(I,J,K), J=51,60)
C      WRITE(*,203) (HEAD(I,J,K), J=51,60)
      READ(30,203) (HEAD(I,J,K), J=61,70)
C      WRITE(*,203) (HEAD(I,J,K), J=61,70)
      READ(30,203) (HEAD(I,J,K), J=71,80)
C      WRITE(*,203) (HEAD(I,J,K), J=71,80)
      READ(30,203) (HEAD(I,J,K), J=81,90)
C      WRITE(*,203) (HEAD(I,J,K), J=81,90)
      READ(30,203) (HEAD(I,J,K), J=91,100)
C      WRITE(*,203) (HEAD(I,J,K), J=91,100)





 20   CONTINUE
 10   CONTINUE
 203  FORMAT(10(F10.3))
      CLOSE(30)
      END

