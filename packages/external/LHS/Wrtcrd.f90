C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    8:36 am
C****************************************************************
C SUBROUTINE WRTCRD PROCESSES THE DISTRIBUTION PARAMETER
C STATEMENTS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::WRTCRD
      SUBROUTINE WRTCRD(ITYPE,LABEL)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE   not needed as WRTCRD has no error conditions	sld01
cc    WRTCRD is called from routines:  RDPAR,RDPAR2                     sld01
cc    WRTCRD does not call any other routines                           sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS
cc    PARAMS provides: LENC						sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides: NV, IDIST
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER LABEL*(*),LABVAR*(LENC),BLANK
      PARAMETER (BLANK=' ')
C
      NV=NV+1
      IDIST(NV)=ITYPE
      NC=LEN(LABEL)
      IC=0
      LABVAR=BLANK
  100 CONTINUE
      IC=IC+1
      IF(IC.GT.NC)GO TO 200
      IF(LABEL(IC:IC).EQ.BLANK)GO TO 100
      LABVAR=LABEL(IC:NC)
  200 CONTINUE
C Comment out because conflicts with NPSOL
C      WRITE(9)LABVAR
      RETURN
      END
