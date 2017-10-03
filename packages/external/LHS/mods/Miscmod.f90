C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  24 May 101   10:35 am
C
C=======================================================================
C
C  This file contains several miscellaneous modules for the LHS code.
C  They are derived from common blocks that were not separated out
C  as their own include files.  The module names are the same as the
C  names of the common blocks that they replace, and are seen in the
C  following USE statements, which should be substituted for the COMMON
C  statements within the program units.
C
C      USE STAR
C      USE UICORR
C      USE CHRCRD
C      USE OBSTR
C      USE PDMAT
C      USE FIRSTS
C
C  In addition, there is a module to declare any local arrays that
C  were previously dimensioned static.
C
C      USE LOCALVARS
C
C  Each module has an initialization subroutine to allocate and initialize
C  arrays.  The name of the initialization subroutine is always the name
C  of the module plus the characters "_INIT", so the initialization subroutine
C  for STAR is STAR_INIT.  There are never any arguments to these routines.
C
C=======================================================================
C
      MODULE STAR
C
C      This module replaces the common block:
C      COMMON/STAR/NSUBOB(NINTMX),SUBINT(NINTMX+1),NINT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NSUBOB
      INTEGER, ALLOCATABLE :: NSUBOB(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::SUBINT
      DOUBLE PRECISION, ALLOCATABLE :: SUBINT(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NINT
      INTEGER :: NINT = 0
C
C      Here is the initialization subroutine for this module
C
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::STAR_INIT
      SUBROUTINE STAR_INIT()
C
        USE PARMS
C
        ALLOCATE( NSUBOB(NINTMX), SUBINT(NINTMX+1) )
        NSUBOB = 0
        SUBINT = 0.0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::STAR_CLOSE
      SUBROUTINE STAR_CLOSE()
C
        DEALLOCATE( NSUBOB, SUBINT )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
C=======================================================================
C
      MODULE UICORR
C
C      This module replaces the common block:
C      COMMON/UICORR/ICVAR(NCVAR),JCVAR(NCVAR),CVAR(NCVAR),NCV
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ICVAR
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::JCVAR
      INTEGER, ALLOCATABLE :: ICVAR(:), JCVAR(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CVAR
      DOUBLE PRECISION, ALLOCATABLE :: CVAR(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NCV
      INTEGER :: NCV = 0
C
C      Here is the initialization subroutine for this module
C
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::UICORR_INIT
      SUBROUTINE UICORR_INIT()
C
        USE PARMS
C
        ALLOCATE( ICVAR(NCVAR), JCVAR(NCVAR), CVAR(NCVAR) )
        ICVAR = 0
        JCVAR = 0
        CVAR = 0.0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::UICORR_CLOSE
      SUBROUTINE UICORR_CLOSE()
C
        DEALLOCATE( ICVAR, JCVAR, CVAR )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
C=======================================================================
C
      MODULE CHRCRD
C
C      This module replaces the common block:
C      CHARACTER*(LENC) CRDSTR(NVAR)
C      COMMON/CHRCRD/CRDSTR
C
C      Must use parameters defined in the PARMS module
        USE PARMS
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CRDSTR
        CHARACTER (LEN=LENC), ALLOCATABLE :: CRDSTR(:)
C
C      Here is the initialization subroutine for this module
C
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CHRCRD_INIT
      SUBROUTINE CHRCRD_INIT()
C
        USE PARMS
C
        ALLOCATE( CRDSTR(NVAR) )
        CRDSTR = ' '
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CHRCRD_CLOSE
      SUBROUTINE CHRCRD_CLOSE()
C
        DEALLOCATE( CRDSTR )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
C=======================================================================
C
      MODULE OBSTR
C
C      This module replaces the common block:
C      COMMON/OBSTR/NSTR,NOBSTR(NVAR)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NSTR
        INTEGER :: NSTR = 0
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NOBSTR
        INTEGER, ALLOCATABLE :: NOBSTR(:)
C
C      Here is the initialization subroutine for this module
C
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::OBSTR_INIT
      SUBROUTINE OBSTR_INIT()
C
        USE PARMS
C
        ALLOCATE( NOBSTR(NVAR) )
        NOBSTR = 0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::OBSTR_CLOSE
      SUBROUTINE OBSTR_CLOSE()
C
        DEALLOCATE( NOBSTR )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
C=======================================================================
C
      MODULE PDMAT
C
C      This module replaces the common block:
C      COMMON/PDMAT/Z(NVAR,NVAR),D(NVAR)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::Z
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::D
        DOUBLE PRECISION, ALLOCATABLE :: Z(:,:),D(:)
C
C      Here is the initialization subroutine for this module
C
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PDMAT_INIT
      SUBROUTINE PDMAT_INIT()
C
        USE PARMS
C
        ALLOCATE( Z(NVAR,NVAR), D(NVAR) )
        Z = 0.0
        D = 0.0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PDMAT_CLOSE
      SUBROUTINE PDMAT_CLOSE()
C
        DEALLOCATE( Z, D )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
C
C=======================================================================
C
      MODULE FIRSTS
C
C     These are all of the variables that have to be reinitialized each
C     time a new execution of the code is started.
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C       For routine RMCNPI
C         DATA ISARG / 0 /
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ISARG
        INTEGER :: ISARG
ccc     Added variable for second random number routine RMCNPI2         SLD
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::JSARG
        INTEGER :: JSARG                                                SLD
ccc
C
C       For Routine ERSTGT
C         DATA LNF/-1/,LNT/0/
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LNF
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LNT
        INTEGER :: LNF, LNT
C
      CONTAINS 
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::FIRSTS_INIT
      SUBROUTINE FIRSTS_INIT()
C
cc    Also initalized in GAMMA and IGAUS                                SLD
        JSARG = 0                                                       SLD
cc                                                                      SLD
        ISARG = 0
        LNF = -1
        LNT = 0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::FIRSTS_CLOSE
      SUBROUTINE FIRSTS_CLOSE()
C
C       Nothing to deallocate
C
        Return
C
      END SUBROUTINE
C
      END MODULE
C
C=======================================================================
C
      MODULE LOCALVARS
c
c       This module contains arrays that were previously dimensioned
c       locally in various subroutines.  They were moved here to make sure
c       that the allocation and deallocation of these arrays during
c       the program does not wipe out data tat the original programmer
c       assumed would be saved from one execution to another.
c
        USE PARMS
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c       From CORCAL:
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::XM
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::SSQ
        DOUBLE PRECISION, ALLOCATABLE :: XM(:), SSQ(:)
c
c       From CMCRD:
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::RIJ
        DOUBLE PRECISION, ALLOCATABLE :: RIJ(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IJCVAR
        INTEGER, ALLOCATABLE :: IJCVAR(:)
c
c       From POSDEF:
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::WK
        DOUBLE PRECISION, ALLOCATABLE :: WK(:)
c
c
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LOCALVARS_INIT
      SUBROUTINE LOCALVARS_INIT()
C
        USE PARMS
c
c       From CORCAL:
        ALLOCATE( XM(NVAR), SSQ(NVAR) )
c
c       From CMCRD:
        ALLOCATE( RIJ(NCVAR*2),IJCVAR(2*NCVAR) )
        RIJ = 0.0
        IJCVAR = 0
c
c       From POSDEF:
        ALLOCATE( WK((NVAR*(NVAR+1))/2) )
        WK = 0.0
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LOCALVARS_CLOSE
      SUBROUTINE LOCALVARS_CLOSE()
c
c       From CORCAL:
        DEALLOCATE( XM, SSQ )
c
c       From CMCRD:
        DEALLOCATE( RIJ, IJCVAR )
c
c       From POSDEF:
        DEALLOCATE( WK )
C
        Return
C
      END SUBROUTINE
C
c
c
      END MODULE
