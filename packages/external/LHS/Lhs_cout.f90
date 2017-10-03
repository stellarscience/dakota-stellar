C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Jun 101   10:58 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_COROUT
      SUBROUTINE LHS_COROUT(MXNUMVAR,IERROR,C1MTRX,C2MTRX,
     xNumCor,NumVar,IPosDef)
c     LHS_COROUT returns four correlation matrices to the user:
c     two triangular matrices each in C1MTRX and C2MTRX
c
c     Input Data:
c     MXNUMCR = the maximum number of correlations that can be requested
c     MXNUMVAR = the maximum number of variables allowed in a LHS run
c     Returned Data:
c     C1MTRX - original correlation matrix in lower triangle;
c              adjusted (to positive definite) matrix in upper triangle
c     C2MTRX - raw data correlation matrix created by sample data in lower
c              triangle; rank data correlation matrix in upper triangle
c     NumCor - the actual number of correlations requested in LHS run
c     NUMVAR- the number of variables for which observations
c             are generated (with no duplication for same as)
c     IPosDef - IPosDef indicator set to one if adjusted were made to matrix
c     IError - integer error flag
c           If an error is found, IError = 1 is returned
c
      USE KILLFILE
      USE PARMS
c     PARAMS provides:   NCVAR,NVAR
      USE CPARAM
c     CPARAMS provides:  ICM,NV
      USE InByCall
cc    InByCall provides: LINT,LPREP,LFILES,LDIST,LRUN,NNames,LPOSDEF
cc                       and vectors LCMSav, VCTR1 and VCTR2:
cc                       VCTR1(NVAR*(NVAR+1)),VCTR2(NVAR*(NVAR+1))
      USE CCMATR
cc    CCMATR provides:  NCM, and arrays LCM and CORR
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C1MTRX(MXNUMVAR,MXNUMVAR),C2MTRX(MXNUMVAR,MXNUMVAR)
c
c
c check that LRUN has been executed
      IF (LRUN /= 1) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9041)
         WRITE (99,9041)
         WRITE(4,9041)
         Return
      END IF
      IF (MXNUMVAR > NVAR) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9043) NVAR
         WRITE (99,9043) NVAR
         WRITE(4,9043) NVAR
         Return
      END IF
c check that NV /> MXNUMVAR
      IF (NV > MXNUMVAR) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9044) NV
         WRITE (99,9044) NV
         WRITE(4,9044) NV
         Return
      END IF

c check that ICM=1 (from CPARAM module)
      IF (ICM /= 1) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9046)
         WRITE (99,9046)
         WRITE(4,9046)
         Return
      END IF
c
c     store identity matrix in C1MTRX
      DO I = 1,NV
         do j = 1,NV
            IF (I .eq. j) THEN
               C1MTRX(I,j) = 1.0
            Else
               C1MTRX(I,j) = 0.0
            END IF
         end do
      END DO
c     VCTR1 stored in LHS_PREP.for
      NumCor = NCM
c     store original correlation matrix in lower triangle of C1MTRX
      K = 0
      DO I = 1,NCM
         DO J = 1,I
            K = K + 1
            if (I .ne. j) then
              C1MTRX(LCMSav(I),LCMSav(j)) = VCTR1(K)
            end if            
         END DO
      END DO
c     store adjusted (to positive definite) matrix in C1MTRX
c               in upper triangle (i.e. reverse order)
      DO I = 1,NCM
         do j = 1,I
            K = K + 1
cccc            C1MTRX(J,I) = VCTR1(K)
            if (I .ne. j) then
              C1MTRX(LCMSav(j),LCMSav(i)) = VCTR1(K)
            end if
         end do
      END DO
c
c     VCTR2 stored in CorOut.for
      NumVar = NV
c     store raw data correlation matrix created by sample data in lower
c              triangle of
      K = 0
      DO I = 1,NV
         DO J = 1,I
            K = K + 1
            C2MTRX(I,J) = VCTR2(K)
         END DO
      END DO
c     store rank data correlation matrix in upper triangle of C2MTRX
      DO I = 1,NV
         do j = 1,I
            K = K + 1
            C2MTRX(J,I) = VCTR2(K)
         end do
      END DO
c
c     store positive definite indicator for return to user
      IPosDef = LPOSDEF
c
      Return
c
 9041 FORMAT(//,5X,'LHS_RUN must be called prior to calling LHS_COROUT'
     x,//)
 9043 FORMAT(//,5X,'Maximum number of variables in call list',
     x1X,'parameter exceeds LHS dimension NVAR =',I6,/,5X,
     x'Change LHS dimension using SIPRA.INI')
 9044 FORMAT(//,5X,'Number of LHS variables exceeds the',1X,
     x'maximum number of variables provided in call list',
     x1X,'number of LHS variables =',I6)
 9046 FORMAT(//,5X,'LHS_CORR and LHS_RUN must be called prior to',1X
     x,'calling LHS_COROUT',//)
c
      END SUBROUTINE
