C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 Mar 101   11:01 am
C**********************************************************
C SUBROUTINE GEOM GENERATES A GEOMETRIC DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::GEOM
      SUBROUTINE GEOM(J)
cc    GEOM is called from:  LHS                                         sld01
cc    GEOM calls routines:  INTRPD,RNUMLHS1                                sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARMS provides:  MAXTB                                            sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N                                               sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X array                                          sld01
C     INCLUDE 'CWORKX.INC'                                              GDW-96  
      USE CWORKX                        
cc    CWORKX provides:  XTABLE array                                    sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      This statement removed to make modules work - gdw-96
C      DIMENSION XTABLE(MAXTB,2)
C      EQUIVALENCE (XTABLE(1,1),XX(1))
      DOUBLE PRECISION PDT,PINPUT,PINPL,PCOMPL,STOTAL
      DOUBLE PRECISION TOTAL						SLD01
C
cc    Statement Function:                                               sld01
      LOC(I,J) = (J-1)*N+I
C
C     --- READ THE INPUT PARAMETER AND GENERATE A DISCRETE CUMULATIVE
C     --- DISTRIBUTION BASED ON THAT PARAMETER
C
      READ (8) AAA
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
      DO 100 I=0, 9999999
         PDT= I*PCOMPL + PINPL
         PDT= EXP(PDT)
         STOTAL = STOTAL + PDT
         IF (STOTAL .GE. DELTA) THEN
            TOTAL=TOTAL+STOTAL
            II=II+1
cc          test for dimension limit                                    sld01
            IF (II > MAXTB) THEN                                        sld01
               WRITE(4,998) MAXTB                                       sld01
               WRITE(99,998) MAXTB                                      sld01
               KLLERR = .True.                                          sld01
               return                                                   sld01
            END IF                                                      sld01
cc                                                                      sld01
            XTABLE(II,1)=I
            XTABLE(II,2)=TOTAL
            IF (TOTAL .GT. TLIMIT) GO TO 101
            STOTAL=0.0D0
         END IF
  100 CONTINUE
C
C     --- IF PROGRAM GETS HERE, THE CREATION OF THE CUMULATIVE FUNCTION
C     --- WAS UNSUCCESSFUL, SO WRITE A MESSAGE AND STOP.
      WRITE (4,999) 'CREATION OF A GEOMETRIC DISTRIBUTION WAS ',
     1 'NOT SUCCESSFUL.  THE INPUT PARAMETER WAS TOO SMALL.'
      WRITE (99,999) 'CREATION OF A GEOMETRIC DISTRIBUTION WAS ',
     1 'NOT SUCCESSFUL.  THE INPUT PARAMETER WAS TOO SMALL.'
  999 FORMAT ('1',5X,A,A)
  998 FORMAT (' MAXTB limit exceed in routine GEOM, MAXTB = ',I5)       sld01
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
      DO I = 1, N
         PROB = PROBINC * RNUMLHS1() + STRTPT
         CALL INTRPD(PROB, BX, XTABLE, MAXTB, IMIN, II)
cc         If(KLLERR) Return     INTRPD has no error conditions         sld01
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
