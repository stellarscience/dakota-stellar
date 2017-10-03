C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:32 am
C****************************************************************
C SUBROUTINE VIF COMPUTES THE VARIANCE INFLATION FACTOR OF A
C CORRELATION MATRIX
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::VIF
      SUBROUTINE VIF
cc    VIF is called from routine:  COROUT                               sld01
cc    VIF calls routine:  DSINV                                         sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  NV                                              sld01
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                        
cc    CCMATR provides:  CORR array                                      sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    Statement Function:                                               sld01
      LOC1(I,J)=J+(I*I-I)/2
C
      CALL DSINV(NV)
      If(KLLERR) Return
      CRKMX=0.0
      DO 630 I=1,NV
        CRK=CORR(LOC1(I,I))
        IF(CRK.GT.CRKMX)CRKMX=CRK
  630 CONTINUE
      WRITE(4,9001)CRKMX
      RETURN
 9001 FORMAT('0','THE VARIANCE INFLATION FACTOR FOR THIS MATRIX IS',
     1       F6.2)
      END
