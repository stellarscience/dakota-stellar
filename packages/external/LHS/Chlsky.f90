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
C SUBROUTINE CHLSKY COMPUTES THE CHOLESKY FACTORIZATION OF A
C CORRELATION MATRIX STORED IN LOWER TRIANGULAR SYMMETRIC FORM
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CHLSKY
      SUBROUTINE CHLSKY
cc    CHLSKY is called from routines:  LHS,MIX                          sld01
cc    CHLSKY does not call any other external routines                  sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE             not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides: NV                                               sld01
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                        
cc    CCMATR provides:  CORR array                                      sld01
C     INCLUDE 'CWORKC.INC'                                              GDW-96  
      USE CWORKC                        
cc    CWORKC provides:  Q and S arrays                                  sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    Statement Function:                                               sld01
      LOC1(I,J)=J+(I*I-I)/2
C
      NVX=(NV*(NV+1))/2
      DO 70 I1=1,NVX
         Q(I1)=0.0
   70 CONTINUE
      I1=1
      DO 80 I2=1,NV
         Q(LOC1(I2,I1))=CORR(LOC1(I2,I1))
   80 CONTINUE
C
      DO 200 I1=2, NV
         I1MIN=I1-1
         DO 100 I2=1,I1MIN
            Q(LOC1(I1,I1))=Q(LOC1(I1,I1))+Q(LOC1(I1,I2))**2
  100    CONTINUE
         Q(LOC1(I1,I1))=SQRT(1.0-Q(LOC1(I1,I1)))
         IF (I1.LT.NV) THEN
            IPLUS=I1+1
            DO 120 I2=IPLUS,NV
               DO 110 K=1,I1MIN
                  Q(LOC1(I2,I1))=Q(LOC1(I2,I1))+
     1                    Q(LOC1(I2,K))*Q(LOC1(I1,K))
  110          CONTINUE
               Q(LOC1(I2,I1))=(CORR(LOC1(I2,I1))-Q(LOC1(I2,I1)))/
     1                    Q(LOC1(I1,I1))
  120       CONTINUE
         END IF
  200 CONTINUE
C
      RETURN
      END
