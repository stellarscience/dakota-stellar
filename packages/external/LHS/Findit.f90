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
C SUBROUTINE FINDIT IS USED IN THE POSITIVE DEFINITE CHECK
C OF THE CORRELATION MATRIX
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::FINDIT
      SUBROUTINE FINDIT(NP,M,EIG,ICONV)
cc    FINDIT is called from:  POSDEF                                    sld01
cc    FINDIT does not call any other external routines                  sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE      --  not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N                                               sld01
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                        
cc    CCMATR provides:  CORR                                            sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X array                                          sld01
C
C     These statements removed to make modules work - GDW-96
C     COMMON/PDMAT/Z(NVAR,NVAR),D(NVAR)
      USE PDMAT
cc    PDMAT provides: Z and D arrays                                    sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOC(I,J)=(J-1)*N+I
      NEV = 0
      DO 9 I = 1,NP
         IF (D(I) .LT. 0.0) NEV=NEV+1
    9 CONTINUE
C
      IF (NEV .EQ. 0) THEN
C
         ICONV=1
C
      ELSE
C
         DO 11 I = 1,NEV
            D(I)=EIG
   11    CONTINUE
C
         L1=NEV+1
         L2=NEV+NEV
         DO 12 I = L1,L2
            IF(D(I).LT.EIG)D(I)=EIG
   12    CONTINUE
C
         DO 4 I = 1,M
         DO 4 J = 1,M
            X(LOC(I,J)) = 0.0
    4    CONTINUE
C
         DO 5 I = 1,NP
         DO 5 J = 1,NP
         DO 5 K = 1,NP
            X(LOC(I,J)) = X(LOC(I,J)) + Z(I,K)*D(K)*Z(J,K)
    5    CONTINUE
C
         DO 30 I = 1,NP
            X(LOC(I,I)) = 1.0
   30    CONTINUE
C
         KI = 0
         DO 10 I = 1,NP
         DO 10 J = 1,I
            KI = KI + 1
            CORR(KI) = X(LOC(I,J))
   10    CONTINUE
C
      ENDIF
C
      RETURN
      END
