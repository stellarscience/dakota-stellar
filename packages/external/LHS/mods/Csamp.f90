C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  21 Mar 101   10:12 am
C     -- INCLUDE FILE FOR THE COMMON BLOCK SAMP
C      COMMON /SAMP/ X(MAXNNV), XSAVE(MAXNNV)
C
c===============================================================
C
      MODULE CSAMP
cc    only 2001 sld changes were comments                               sld01
C
C
C       Here are the elements from the old common block
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::X
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::XSAVE
      DOUBLE PRECISION, ALLOCATABLE :: X(:), XSAVE(:)
C
C       Now here is the initialization routine for this module
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CSAMP_INIT
      SUBROUTINE CSAMP_INIT()
C
cc      PARMS provides MAXNNV                                           sld01
        USE PARMS
C
        ALLOCATE( X(MAXNNV) )
        X = 0.0
C
        ALLOCATE( XSAVE(MAXNNV) )
        XSAVE = 0.0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CSAMP_CLOSE
      SUBROUTINE CSAMP_CLOSE()
C
        DEALLOCATE( X )
C
        DEALLOCATE( XSAVE )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
