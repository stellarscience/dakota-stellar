C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  21 Mar 101   10:16 am
C     -- INCLUDE BLOCK FOR THE COMMON BLOCK WORKC
C      COMMON/WORKC/Q((NVAR*(NVAR+1))/2),S((NVAR*(NVAR+1))/2)
C
C==============================================================
C
      MODULE CWORKC
cc    only 2001 sld changes were comments                               sld01
C
C
C       Here are the elements from the old common block
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::Q
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::S
      DOUBLE PRECISION, ALLOCATABLE :: Q(:), S(:)
C
C       Now here is the initialization routine for this module
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CWORKC_INIT
      SUBROUTINE CWORKC_INIT()
C
cc      PARMS provides NVAR                                             sld01
        USE PARMS
C
        ALLOCATE( Q((NVAR*(NVAR+1))/2) )
        Q = 0.0
C
        ALLOCATE( S((NVAR*(NVAR+1))/2) )
        S = 0.0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CWORKC_CLOSE
      SUBROUTINE CWORKC_CLOSE()
C
        DEALLOCATE( Q )
C
        DEALLOCATE( S )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
