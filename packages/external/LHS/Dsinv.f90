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
C SUBROUTINE DSINV IS USED IN INVERTING A CORRELATION MATRIX
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::DSINV
      SUBROUTINE DSINV (N)
cc    DSINV is called from routine:  VIF                                sld01
cc    DSINV calls routine:  DMFSD                                       sld01
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
      IPARM=0
      CALL DMFSD (N,IPARM)
      If(KLLERR) Return
cc      IF (IPARM.LT.0) RETURN                                          sld01
      IF (IPARM.LT.0) then                                              sld01
       KLLERR = .True.                                                  sld01
      END If                                                            sld01
      IPIV=N*(N+1)/2
      IND=IPIV
      DO 40 I=1,N
      DIN=1.0/CORR(IPIV)
      CORR(IPIV)=DIN
      MIN=N
      KEND=I-1
      LANF=N-KEND
      IF (KEND.LE.0) GO TO 30
      J=IND
      DO 20 K=1,KEND
      WORK=0.0
      MIN=MIN-1
      LHOR=IPIV
      LVER=J
      DO 10 L=LANF,MIN
      LVER=LVER+1
      LHOR=LHOR+L
   10 WORK=WORK+CORR(LVER)*CORR(LHOR)
      CORR(J)=-WORK*DIN
   20 J=J-MIN
   30 IPIV=IPIV-MIN
   40 IND=IND-1
      DO 60 I=1,N
      IPIV=IPIV+I
      J=IPIV
      DO 60 K=I,N
      WORK=0.0
      LHOR=J
      DO 50 L=K,N
      LVER=LHOR+K-I
      WORK=WORK+CORR(LHOR)*CORR(LVER)
   50 LHOR=LHOR+L
      CORR(J)=WORK
   60 J=J+K
      RETURN
      END
