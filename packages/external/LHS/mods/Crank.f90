C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  21 Mar 101   10:10 am
C     -- INCLUDE FILE FOR THE COMMON BLOCK RANK
C      COMMON /RANK/ XV(NMAX), RXV(NMAX), IWK(NMAX)
C
c===============================================================
C
      MODULE CRANK
cc    only 2001 sld changes were comments                               sld01
C
C
C       Here are the elements from the old common block
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::XV
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::RXV
        DOUBLE PRECISION, ALLOCATABLE :: XV(:), RXV(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IWK
        INTEGER, ALLOCATABLE :: IWK(:)
C
C       Now here is the initialization routine for this module
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CRANK_INIT
      SUBROUTINE CRANK_INIT()
C
cc       PARMS provides NMAX
        USE PARMS
C
        ALLOCATE( XV(NMAX) )
        XV = 0.0
C
        ALLOCATE( RXV(NMAX) ) 
        RXV = 0.0
C
        ALLOCATE( IWK(NMAX) )
        IWK = 0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CRANK_CLOSE
      SUBROUTINE CRANK_CLOSE()
C
        DEALLOCATE( XV )
C
        DEALLOCATE( RXV )
C
        DEALLOCATE( IWK )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
