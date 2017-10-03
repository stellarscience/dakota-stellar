C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   2 Apr 101   10:51 am
C**********************************************************
C SUBROUTINE BINOM GENERATES A BINOMIAL DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::BINOM
      SUBROUTINE BINOM(J)
cc    BINOM is called from routine:  LHS                                sld01
cc    BINOM calls routines:  FACTOR,FACTR2,RNUMLHS1,INTRPD                 sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE            not needed				sld01
cc      FACTOR,FACTR2,RNUMLHS1,INTRPD have no error conditions             sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARAMS provides:  MAXTB                                           sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS                                           sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X array                                          sld01
C     INCLUDE 'CWORKX.INC'                                              GDW-96  
      USE CWORKX                        
cc    CWORKX provides:  XTABLE array                                    sld01
C
C      These statements removed to make modules work - GDW-96
C      DIMENSION XTABLE(MAXTB,2)
C      EQUIVALENCE (XTABLE(1,1),XX(1))
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PDT,PINPUT,PINPL,PCOMPL,STOTAL
      DOUBLE PRECISION FACTOR,FACTR2
      DOUBLE PRECISION TOTAL						SLD01
C
      LOC(I,J) = (J-1)*N+I
C
C     --- READ THE INPUT PARAMETER AND GENERATE A DISCRETE CUMULATIVE
C     --- DISTRIBUTION BASED ON THAT PARAMETER
C
      READ (8) AAA, NTOT
C
      II=0
      TOTAL=0.0
      STOTAL = 0.0D0
      DELTA=1.0/(MAXTB-1)
      TLIMIT=1.0-DELTA-DELTA
      PINPUT=AAA
      PINPL=LOG(PINPUT)
      PCOMPL=LOG(1.0D0-PINPUT)
C
      DO 100 I=0, NTOT
         ITMP=NTOT-I
         PDT= FACTOR(I,NTOT) - FACTR2(1,ITMP) + I*PINPL + ITMP*PCOMPL
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
  101 XTABLE(II,2)=1.0
C
      PROBINC=1.0/N
      IF (IRS .NE. 0) PROBINC=1.0
      STRTPT = 0.0
      IMIN = 1
C
      DO I = 1, N
         PROB = PROBINC * RNUMLHS1() + STRTPT
         CALL INTRPD(PROB, BX, XTABLE, MAXTB, IMIN, II)
cc         If(KLLERR) Return  INTRPD has no error conditions            sld01
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
