C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  20 Jun 101   10:02 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_RtvSEED
      SUBROUTINE LHS_RtvSEED(IError,LastSeed)
c     This routine returns the last random number seed used by LHS
c     LHS_RSEED returns:
c        LastSeed - last random number used by LHS
c        IError - integer error flag
c           If an error is found, IError = 1 is returned
c
      USE KILLFILE
c
      USE CPARAM
c     CPARAMS provides:  ISeed,N,NV,NREP,IVarNm,PValue,IDIST,List
      USE InByCall
cc    InByCall provides subroutine flags: NNames,LINIT,LPREP,LFILES,
cc                                        IScrh1,IScrh6,LRUN,LPOSDEF
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c check that LRUN has been executed
      IF (LRUN /= 1) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9041)
         WRITE (99,9041)
         WRITE(4,9041)
         Return
      END IF
c
      LastSeed = ISeed
      Return
c
 9041 FORMAT(//,5X,'LHS_RUN must be called prior to calling LHS_RtvSEED'
     x,//)
c
      END SUBROUTINE
