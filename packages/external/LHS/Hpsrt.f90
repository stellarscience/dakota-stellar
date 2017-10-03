C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:30 am
C****************************************************************
C SUBROUTINE HPSRT IS USED IN THE RANKING OF THE SAMPLED DATA
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::HPSRT
      SUBROUTINE HPSRT
cc    HPSRT is called from routine:  RANKER                             sld01
cc    HPSRT does not call any other external routines                   sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE        -- not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N                                               sld01
C     INCLUDE 'CRANK.INC'                                               GDW-96  
      USE CRANK                         
cc    CRANK provides:  XV and RXV arrays                                sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER R
C
      L=N/2+1
      R=N
   10 IF (L.LE.1) GO TO 70
      L=L-1
      XHOLD=XV(L)
      YHOLD=RXV(L)
   20 J=L
   30 I=J
      J=2*J
      IF (J-R) 40,50,60
   40 IF (XV(J).LT.XV(J+1)) J=J+1
   50 IF (XHOLD.GE.XV(J)) GO TO 60
      XV(I)=XV(J)
      RXV(I)=RXV(J)
      GO TO 30
   60 XV(I)=XHOLD
      RXV(I)=YHOLD
      GO TO 10
   70 XHOLD=XV(R)
      YHOLD=RXV(R)
      XV(R)=XV(1)
      RXV(R)=RXV(1)
      R=R-1
      IF (R.GT.1) GO TO 20
      XV(1)=XHOLD
      RXV(1)=YHOLD
      RETURN
      END
