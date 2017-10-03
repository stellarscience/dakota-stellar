C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  21 Jun 101    8:08 am
C****************************************************************
C SUBROUTINE GAMMA GENERATES GAMMA DISTRIBUTIONS WITH PARAMETERS
C ALPHA AND BETA.  THE DISTRIBUTION IS EXACT IF ALPHA = 1.
C OTHERWISE, AN ACCEPTANCE-REJECTION SCHEME IS USED TO GENERATE
C 10,000 RANDOM OBSERVATIONS FROM A GAMMA DISTRIBUTION WITH
C PARAMETERS ALPHA AND BETA AND THE RESULTING EDF IS SAMPLED TO
C GENERATE THE DESIRED DISTRIBUTION. ALPHA < 1 THE ALGORITHM OF
C BEST (1983) IS USED.  FOR ALPHA > 1 THE ALGORITHM OF MINH (1988)
C IS USED.
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::GAMMA
      SUBROUTINE GAMMA(J)
cc    GAMMA is called from routine:  LHS                                sld01
cc    GAMMA calls routines:  GAMMAM,GAMMAB,SIFT,RNUMLHS1                   sld01
cc    Note:  RNUMLHS2 is called by GAMMAM and GAMMAB                       SLD
cc    RNUMLHS2 was added for use in Acceptance-Rejection scheme            SLD
cc          found in GAMMA & IGAUS                                      SLD
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE      -- no longer needed                           sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARMS provides:  MAXTB                                            sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS,ISEEDSV,JSEED                             sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides: X array                                           sld01
C     INCLUDE 'CWORKX.INC'                                              GDW-96  
      USE CWORKX                        
cc    CWORKX provides:  XX array                                        sld01
cc
      USE FIRSTS, ONLY: JSARG
cc
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION VALUE,Z,B
C
cc    Statement Function:                                               sld01
      LOC(I,J)=(J-1)*N+I
C
      PROBINC = 1./FLOAT(N)
      IF(IRS .EQ. 1) PROBINC = 1.0
      STRTPT = 0.
      READ(8) ALPHA, BETA
C
C     CHECK FOR ALPHA = 1, WHICH CAN BE DONE BY DIRECT INVERSION
C
      IF (ALPHA .EQ. 1.) THEN
         DO 1 I = 1,N
            R = PROBINC*RNUMLHS1() + STRTPT
            RES = -LOG(1. - R)/BETA
            X(LOC(I,J)) = DMAX(RES,1.D-10)
            IF(IRS .EQ. 0) STRTPT = STRTPT + PROBINC
    1    CONTINUE
         RETURN
      END IF
C
C     GENERATE "MAXTB*2" RANDOM DEVIATES FROM THE GAMMA DISTRIBTION WITH
C     PARAMETERS ALPHA AND BETA
C
cc    ReSeed RNUMLHS2                                                      SLD
      JSARG = 0                                                         SLD
      JSEED = ISEEDSV                                                   SLD
cc                                                                      SLD
      IF (ALPHA .LT. 1.) THEN
         Z = .07 + .75*(1. - ALPHA)**.5
         B = 1. + EXP(-Z)*ALPHA/Z
      ENDIF
      NOBS=2*MAXTB
      DO 3 I = 1,NOBS
         IF (ALPHA .GT. 1.)THEN
C           -- CALL ROUTINE BY MINH
            CALL GAMMAM(ALPHA,VALUE)
cc            If(KLLERR) Return    GAMMAM has no error conditions       sld01
         ELSE
C           -- CALL ROUTINE BY BEST
            CALL GAMMAB(ALPHA,VALUE,Z,B)
cc            If(KLLERR) Return   GAMMAB has no error conditions        sld01
         ENDIF
         XX(I) = VALUE/BETA
    3 CONTINUE
C
C     SORT THE RANDOM DEVIATES TO FORM THE EDF
C
      CALL SIFT(XX,NOBS)
cc      If(KLLERR) Return      SIFT has no error conditions             sld01
C
C     GENERATE THE DESIRED SAMPLE BY SAMPLING FROM THE EDF
C
      DO I = 1, N
         R = PROBINC*RNUMLHS1() + STRTPT
cc         NR = MAX(1 , NINT(R*NOBS))                                   SLD
cc         RES = XX(NR)                                                 SLD
         FNR = R*(NOBS-1) + 1                                           SLD
         NR = INT(FNR)                                                  SLD
         NR1 = NR + 1                                                   SLD
cccc                                                                    SLD
         FRAC = FNR-FLOAT(NR)                                           SLD
         RES = FRAC*(XX(NR1)-XX(NR)) + XX(NR)                           SLD
cccc     RES = FRACTION(FNR)*(XX(NR1)-XX(NR)) + XX(NR)
cccc                                                                    SLD
         X(LOC(I,J)) = DMAX(RES,1.D-10)
         IF (IRS == 0) STRTPT = DBLE(I) / DBLE(N)
      END DO
C
      RETURN
      END
