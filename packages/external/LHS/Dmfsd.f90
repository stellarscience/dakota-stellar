C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:29 am
C****************************************************************
C SUBROUTINE DMFSD IS USED IN INVERTING A CORRELATION MATRIX
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::DMFSD
      SUBROUTINE DMFSD (N,IPARM)
cc    DMFSD is called from routine:  DSINV                              sld01
cc    DMFSD does not call any other external routines                   sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                        
cc    CCMATR provides:  CORR array                                      sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      KPIV=0
      DO 60 K=1,N
      KPIV=KPIV+K
      IND=KPIV
      LEND=K-1
      TOL=ABS(.01*(CORR(KPIV)))
      DO 60 I=K,N
      DSUM=0.0
      IF (LEND.EQ.0) GO TO 20
      DO 10 L=1,LEND
      LANF=KPIV-L
      LIND=IND-L
   10 DSUM=DSUM+CORR(LANF)*CORR(LIND)
   20 DSUM=CORR(IND)-DSUM
      IF (I.NE.K) GO TO 50
      IF (DSUM-TOL) 30,30,40
   30 IF (DSUM.LE.0.0) GO TO 70
      KT=K-1
      WRITE(4,80)KT
   40 DPIV=SQRT(DSUM)
      CORR(KPIV)=DPIV
      DPIV=1.0/DPIV
      GO TO 60
   50 CORR(IND)=DSUM*DPIV
   60 IND=IND+I
      RETURN
   70 WRITE(4,90)K
      IPARM=-K
      RETURN
C
C
   80 FORMAT(20X,'ROUNDING ERROR IN ROW ',I2)
   90 FORMAT(20X,'MATRIX IS SINGULAR AT ROW ',I2)
      END
