C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Jun 101    9:05 am
      MODULE InByCall
c     This module was added to provide information needed by the various
c     LHS Input-By-Call Subroutines
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::VCTR1
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::VCTR2
      DOUBLE PRECISION, ALLOCATABLE :: VCTR1(:),VCTR2(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LCMSav
      INTEGER, ALLOCATABLE :: LCMSav(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LINIT
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LPREP
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LRUN
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LFILES
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LDIST
      INTEGER :: LINIT,LPREP,LRUN,LFILES,LDIST
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NNames
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IScrh1
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IScrh6
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LPOSDEF
      INTEGER :: NNames,IScrh1,IScrh6,LPOSDEF
c     NNames is number of variable names (was local to RDPAR2)
c     IScrh1=1 for scratch file,IScrh1=2 for user file (sample output file)
c     IScrh6=1 for scratch file,IScrh6=2 for user file (message output file)
c     Vectors VCTR1 and VCTR2 are used to store correlation matrix information
C
C      Here is the initialization subroutine for this module
C
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::InByCall_INIT
      SUBROUTINE InByCall_INIT()
C
      USE PARMS
c     initialize subroutine called tracking switches:
        LINIT = 0
        LPREP = 0
        LRUN = 0
        LFILES = 0
        LDIST = 0
c     other flag indicators
        IScrh1 = 0
        IScrh6 = 0
        LPOSDEF= 0
C
        ALLOCATE( VCTR1(NVAR*(NVAR+1)), VCTR2(NVAR*(NVAR+1)) )
        VCTR1 = 0.0
        VCTR2 = 0.0

        ALLOCATE( LCMSav(NVAR) )
        LCMSAV(1) = 0
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::InByCall_CLOSE
      SUBROUTINE InByCall_CLOSE()
C
        DEALLOCATE( VCTR1, VCTR2 )
C
        DEALLOCATE( LCMSav )
        RETURN
C
      END SUBROUTINE
C
      END MODULE

