C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  21 Mar 101   10:17 am
C     -- INCLUDE BLOCK FOR THE COMMON BLOCK WORKX
C      COMMON/WORKX/XX(2*MAXTB)
C
C===============================================================
C
      MODULE CWORKX
cc    only 2001 sld changes were comments                               sld01
C
C
C       Here are the elements from the old common block
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::XX
        DOUBLE PRECISION, ALLOCATABLE :: XX(:)
C
C       Since there are equivalence statements in the code, and 
C       equivalence can't be used with allocatable arrays, we will
C       define additional arrays here to eliminate the equivalence
C       statements.  The XVLZ array will be handled as an alias of XX.
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::XTABLE
        DOUBLE PRECISION, ALLOCATABLE :: XTABLE(:,:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PRBZ
        DOUBLE PRECISION, ALLOCATABLE :: PRBZ(:)
C
C       Now here is the initialization routine for this module
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CWORKX_INIT
      SUBROUTINE CWORKX_INIT()
C
cc      PARMS provides MAXTB
        USE PARMS
C
        ALLOCATE( XX(2*MAXTB) )
        XX = 0.0
C
        ALLOCATE( XTABLE(MAXTB,2) )
        XTABLE = 0.0
C
        ALLOCATE( PRBZ(MAXTB) )
        PRBZ = 0.0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CWORKX_CLOSE
      SUBROUTINE CWORKX_CLOSE()
C
        DEALLOCATE( XX )
C
        DEALLOCATE( XTABLE )
C
        DEALLOCATE( PRBZ )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
