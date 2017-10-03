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
C SUBROUTINE EXPON GENERATES EXPONENTIAL DISTRIBUTIONS
C WITH PARAMETER LAMBDA
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::EXPON
      SUBROUTINE EXPON(J,IDT)
cc    EXPON is called from routine:  LHS                                sld01
cc    EXPON calls routine:  RNUMLHS1                                       sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE       -- not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS                                           sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X array                                          sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    RNUMLHS1 is an external function                                     sld01
cc    statement function:                                               sld01
      LOC(I,J)=(J-1)*N+I
C
      IF (IDT .EQ. 25) THEN
C        -- TRUNCATED EXPONENTIAL
         READ (8) RLAMBD, A, B
      ELSE IF (IDT .EQ. 26) THEN
C        -- BOUNDED EXPONENTIAL - COMPUTE QUANTILES FOR A AND B
         READ (8) RLAMBD, A, B
         A = 1.0 - EXP(-RLAMBD*A)
         B = 1.0 - EXP(-RLAMBD*B)
      ELSE
C        -- UNBOUNDED EXPONENTIAL - QUANTILES ARE 0.0 AND 1.0
         READ (8) RLAMBD
         A = 0.0
         B = 1.0
      END IF
C
      PROBINC = (B - A) / FLOAT(N)
      IF (IRS .EQ. 1) PROBINC = B - A
      STRTPT = A
      DO I = 1, N
         R=PROBINC*RNUMLHS1()+STRTPT
         X(LOC(I,J))=-LOG(1.0-R)/RLAMBD
         IF (IRS == 0) STRTPT = DBLE(I) / DBLE(N)
      END DO
C
      RETURN
      END
