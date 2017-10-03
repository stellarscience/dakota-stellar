C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  21 Mar 101   10:19 am
C     -- INCLUDE FILE FOR THE COMMON BLOCK CMATR
C     COMMON /CMATR/ CORR((NVAR*(NVAR+1))/2), LCM(NVAR), NCM
      MODULE CCMATR
cc    only 2001 sld changes were comments                               sld01
C
C
C       These are the elements from the old common block
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CORR
        DOUBLE PRECISION, ALLOCATABLE :: CORR(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LCM
        INTEGER, ALLOCATABLE :: LCM(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NCM
        INTEGER NCM
C
C       Now here is the initialization subroutine for this module
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CCMATR_INIT
      SUBROUTINE CCMATR_INIT()
C
cc      PARMS provides NVAR                                             sld01
        USE PARMS
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        ALLOCATE( CORR( (NVAR*(NVAR+1))/2 ) )
        CORR = 0.0
C
        ALLOCATE ( LCM( NVAR ) )
        LCM = 0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CCMATR_CLOSE
      SUBROUTINE CCMATR_CLOSE()
C
        DEALLOCATE( CORR )
C
        DEALLOCATE ( LCM )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
