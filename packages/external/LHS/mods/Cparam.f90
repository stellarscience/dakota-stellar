C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 May 101   10:25 am
C
C      COMMON /CHAR1/ TITLE
Cc      CHARACTER TITLE*LENT
C      CHARACTER (len = LENT) TITLE
C      Character*(NamLen) List
C      Character*60 TreeFl, SFile, MFile, CmdLin
C      COMMON /PARAM/ ISEED, N, NV, IRS, ICM, NREP, IDATA, IHIST,
C     1               ICORR, IDIST(NVAR), IRP, IV1, IRSET,
C     2               List(NVar), IVarNm(NVar), PValue(NVar),
C     3               NamOut, IPtVal, TreeFl, SFile, MFile, CmdLin
C
c===============================================================
C
      MODULE CPARAM
C
C       Here are the elements from the old common block
C
cc      PARMS provides NVAR,LENT,NamLen                                 sld01
        USE PARMS
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::TITLE
        CHARACTER (LEN = LENT) TITLE
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::TreeFl
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::SFile
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::MFile
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CmdLin
        CHARACTER (LEN = 256) :: TreeFl, SFile, MFile, CmdLin
cc        INTEGER :: ISEED, N, NV, IRS, ICM, NREP, IDATA, IHIST         SLD
cc        INTEGER :: ICORR                                              SLD
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ISEED
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::N
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NV
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IRS
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ICM
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NREP
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IDATA
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IHIST
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ICORR
        INTEGER :: ISEED, N, NV, IRS, ICM, NREP, IDATA, IHIST, ICORR    SLD
ccc     Added variables for second random number routine                SLD
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ISEEDSV
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::JSEED
        INTEGER :: ISEEDSV, JSEED                                       SLD
ccc
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IDIST
        INTEGER, ALLOCATABLE :: IDIST(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IRP
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IV1
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IRSET
        INTEGER :: IRP, IV1, IRSET
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::List
        Character*(NamLen), ALLOCATABLE :: List (:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IVarNm
        INTEGER, ALLOCATABLE :: IVarNm(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PValue
        DOUBLE PRECISION, ALLOCATABLE :: PValue(:)
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NamOut
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::I1Col
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IPtVal
        INTEGER NamOut, I1Col, IPtVal
C
C       Now here is the initialization routine for this module
      CONTAINS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CPARAM_INIT
      SUBROUTINE CPARAM_INIT()
C
        USE PARMS
C
        ALLOCATE( IDIST(NVAR) )
        IDIST = 0 
C
        ALLOCATE( IVarNm(NVAR) )
        IVarNm = 0
C
        ALLOCATE( PValue(NVar) )
        PValue = 0.0
C
        ALLOCATE( List(NVar) ) 
        List = " "
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CPARAM_CLOSE
      SUBROUTINE CPARAM_CLOSE()
C
        DEALLOCATE( IDIST )
C
        DEALLOCATE( IVarNm )
C
        DEALLOCATE( PValue )
C
        DEALLOCATE( List )
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
