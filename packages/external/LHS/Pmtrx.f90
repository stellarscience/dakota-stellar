C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:21 am
C****************************************************************
C SUBROUTINE PMTRX PRINTS OUT A CORRELATION MATRIX
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PMTRX
      SUBROUTINE PMTRX(NC,IFLAG)
cc    PMTRX is called from:  LHS,COROUT                                 sld01
cc    PMTRX does not call any other external routines                   sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE      -- not needed					sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  IRS,TITLE
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                        
cc    CCMATR provides:  CORR and LCM arrays                             sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NROW=28)
      PARAMETER (NCOL=15)
C
      LOC1(I,J)=J+(I*I-I)/2
C
      I=1
      I1=1
      I2=NROW
      NPROW=NC/NROW
      IF (NPROW*NROW.NE.NC) NPROW=NPROW+1
      NPCOL=NC/NCOL
      IF (NPCOL*NCOL.NE.NC) NPCOL=NPCOL+1
      DO 100 NR=1,NPROW
      J1=1
      J2=NCOL
      DO 90 IC=1,NPCOL
      IF (J1.GT.I2) GO TO 90
      IF (IFLAG.EQ.2) GO TO 10
      IF (IFLAG.EQ.3) GO TO 20
      IF(IFLAG.EQ.4)GO TO 20
      WRITE(4,180)TITLE
      IF(IRS.EQ.0)WRITE(4,110)I
      IF(IRS.NE.0)WRITE(4,115)I
      GO TO 30
   10 WRITE(4,180)TITLE
      IF(IRS.EQ.0)WRITE(4,120)I
      IF(IRS.NE.0)WRITE(4,125)I
      GO TO 30
   20 WRITE(4,180)TITLE
      IF(IFLAG.EQ.3)WRITE(4,130)I
      IF(IFLAG.EQ.4)WRITE(4,190)I
   30 IF (NC.LT.I2) I2=NC
      IF (NC.LT.J2) J2=NC
      DO 60 K=I1,I2
      K1=K
      IF (K.GE.J1) GO TO 40
      WRITE(4,140)
      GO TO 60
   40 IF (K.GT.J2) GO TO 50
      WRITE(4,150)LCM(K),(CORR(LOC1(K,J)),J=J1,K)
      GO TO 60
   50 WRITE(4,150)LCM(K),(CORR(LOC1(K,J)),J=J1,J2)
   60 CONTINUE
      IF (K1.LT.J2) GO TO 70
      WRITE(4,160)(LCM(II),II=J1,J2)
      WRITE(4,170)
      GO TO 80
   70 WRITE(4,160)(LCM(II),II=J1,K1)
      WRITE(4,170)
   80 CONTINUE
      J1=J1+NCOL
      J2=J2+NCOL
      I=I+1
   90 CONTINUE
      I1=I1+NROW
      I2=I2+NROW
  100 CONTINUE
      RETURN
  110 FORMAT('0','CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE ',
     1       'LATIN HYPERCUBE SAMPLE FOR RAW DATA',30X,'PAGE ',I3)
  115 FORMAT('0','CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE ',
     1       'RANDOM SAMPLE FOR RAW DATA',40X,'PAGE ',I3)
  120 FORMAT('0','CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE ',
     1       'LATIN HYPERCUBE SAMPLE FOR RANK DATA',30X,'PAGE ',I3)
  125 FORMAT('0','CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE ',
     1       'RANDOM SAMPLE FOR RANK DATA',40X,'PAGE ',I3)
  130 FORMAT('0','INPUT RANK CORRELATION MATRIX',95X,'PAGE ',I3)
  140 FORMAT('0')
  150 FORMAT('0',I5,15(F8.4))
  160 FORMAT('0',5X,15I8)
  170 FORMAT('0','VARIABLES')
  180 FORMAT('1',A)
  190 FORMAT('0','ADJUSTED RANK CORRELATION MATRIX',92X,'PAGE ',I3)
      END
