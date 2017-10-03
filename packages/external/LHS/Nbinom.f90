C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101   10:44 am
C**********************************************************
C SUBROUTINE NBINOM GENERATES A NEGATIVE BINOMIAL DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NBINOM
      SUBROUTINE NBINOM(J)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
cc    NBINOM is called from routine:  LHS                               sld01
cc    NBINOM calls routines:  RNUMLHS1,INTRPD,FACTOR,FACTR2                sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARMS provides:  MAXTB                                            sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  IRS,N                                           sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X                                                sld01
C     INCLUDE 'CWORKX.INC'                                              GDW-96  
      USE CWORKX                        
cc    CWORKX provides:  XTABLE                                          sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      This statement removed to make modules work - gdw-96
C      DIMENSION XTABLE(MAXTB,2)
C      EQUIVALENCE (XTABLE(1,1),XX(1))
      DOUBLE PRECISION TOTAL						SLD01
      DOUBLE PRECISION PDT,PINPUT,PINPL,PCOMPL,STOTAL
      DOUBLE PRECISION FACTOR,FACTR2
cc    RNUMLHS1, FACTOR and FACTR2 are External Functions
C
      LOC(I,J) = (J-1)*N+I
C
C     --- READ THE INPUT PARAMETER AND GENERATE A DISCRETE CUMULATIVE
C     --- DISTRIBUTION BASED ON THAT PARAMETER
C
      READ (8) AAA, NR
C
      II=0
      TOTAL=0.0
      STOTAL = 0.0D0
      DELTA=1.0/(MAXTB-1)
      TLIMIT=1.0-DELTA-DELTA
      NQ=NR-1
      PINPUT=AAA
      PINPL=LOG(PINPUT)
      PCOMPL=LOG(1.0D0-PINPUT)
C
      DO 100 I=0, 999999
         ITMP = I + NQ
         PDT= FACTOR(NQ,ITMP) - FACTR2(1,I) + NR*PINPL + I*PCOMPL
         PDT= EXP(PDT)
         STOTAL = STOTAL + PDT
         IF (STOTAL .GE. DELTA) THEN
            TOTAL=TOTAL+STOTAL
            II=II+1
            XTABLE(II,1)=I
            XTABLE(II,2)=TOTAL
            IF (TOTAL .GT. TLIMIT) GO TO 101
            STOTAL=0.0D0
         END IF
  100 CONTINUE
C
C     --- IF PROGRAM GETS HERE, THE CREATION OF THE CUMULATIVE FUNCTION
C     --- WAS UNSUCCESSFUL, SO WRITE A MESSAGE AND STOP.
      WRITE (4,999) 'CREATION OF A NEGATIVE BINOMIAL DISTRIBUTION ',
     1 'WAS NOT SUCCESSFUL.  TRY MORE REASONABLE INPUT PARAMETERS.'
      WRITE (99,999) 'CREATION OF A NEGATIVE BINOMIAL DISTRIBUTION ',
     1 'WAS NOT SUCCESSFUL.  TRY MORE REASONABLE INPUT PARAMETERS.'
  999 FORMAT ('1',5X,A,A)
      KLLERR = .TRUE.
      RETURN
C
  101 XTABLE(II,2)=1.0
C
      PROBINC=1.0/N
      IF (IRS .NE. 0) PROBINC=1.0
      STRTPT = 0.0
      IMIN = 1
C
      DO I=1, N
         PROB = PROBINC * RNUMLHS1() + STRTPT
         CALL INTRPD(PROB, BX, XTABLE, MAXTB, IMIN, II)
cc         If(KLLERR) Return -- INTRPD has no error conditions
         X(LOC(I,J)) = BX
         IF (IRS .EQ. 0) THEN
            STRTPT = DBLE(I) / DBLE(N)
         ELSE
            IMIN = 1
         END IF
      END DO
C
      RETURN
C
      END
