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
C SUBROUTINE IGAUS GENERATES INVERSE GAUSSIAN DISTRIBUTIONS WITH
C PARAMETERS MU AND LAMBDA.  AN ACCEPTANCE-REJECTION SCHEME DUE TO
C ATKINSON IS USED TO GENERATE 10,000 RANDOM OBSERVATIONS FROM
C AN INVERSE GAUSSIAN DISTRIBUTION WITH PARAMETERS MU AND LAMBDA
C AND THE RESULTING EDF IS SAMPLED TO GENERATE THE DESIRED DISTRIBUTION.
C THE MEAN OF THIS DISTRIBUTION IS MU AND THE VARIANCE IS MU**3/LAMBDA.
C
C REFERENCE: A. C. ATKINSON (1982) "THE SIMULATION OF GENERALIZED INVERSE
C            GAUSSIAN AND HYPERBOLIC RANDOM VARIABLES,"
C            SIAM J. SCI. STAT. COMPUT.
C            VOL. 3, NO. 4, DEC. 1982, 502-515.
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IGAUS
      SUBROUTINE IGAUS(J)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE -- not needed, as IGAUS and the routines it calls  sld01
cc                      have no error conditions                        sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96
      USE PARMS
cc    PARMS provides: MAXTB                                             sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides: N,IRS,ISEEDSV,JSEED
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides: X array                                           SLD01
C     INCLUDE 'CWORKX.INC'                                              GDW-96  
      USE CWORKX                        
cc    CWORKX provides:  XX array					sld01
cc
      USE FIRSTS, ONLY: JSARG
cc
cc    IGAUS is called from routine:  LHS                                sld01
cc    IGAUS calls routines:  IGAUS1, IGAUSF, SIFT, RNUMLHS1                sld01
cc    Note:  RNUMLHS2 is called by IGAUSF                                  SLD
cc    RNUMLHS2 was added for use in Acceptance-Rejection scheme            SLD
cc          found in GAMMA & IGAUS                                      SLD
C
cc    IGAUSF is an External Function                                    sld01
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MU,LAMBDA,IGAUSF
C
      LOC(I,J)=(J-1)*N+I
C
cc    ReSeed RNUMLHS2                                                      SLD
      JSARG = 0                                                         SLD
      JSEED = ISEEDSV                                                   SLD
cc                                                                      SLD
      PROBINC=1./FLOAT(N)
      IF (IRS .EQ .1) PROBINC=1.0
      STRTPT=0.
      READ (8) MU,LAMBDA
C
C     CONVERT MU AND LAMBDA TO THE CHI AND PSI NOTATION USED BY
C     ATKINSON AS FOLLOWS:
C
      CHI = LAMBDA
      PSI = LAMBDA/MU**2
C
C     TO GENERATE THE INVERSE GAUSSIAN WITH ATKINSON'S ALGORITHM IT IS
C     NECESSARY TO USE FL = -.5.  HOWEVER, FOR VALUES OF FL < 0 THE
C     FOLLOWING RELATIONSHIP SHOULD BE USED:
C
C          IF     X ~ INVERSE GAUSS(FL,PSI,CHI)
C          THEN 1/X ~ INVERSE GAUSS(-FL,CHI,PSI)
C
C     THAT IS, INTERCHANGE CHI AND PSI IN THE CALLING STATEMENT AND
C     TAKE 1/X AS THE GENERATED VALUE RATHER THAN X.
C     SET THE CONSTANTS REQUIRED IN FUNCTION GIGLT1
C
      FL = .5
      CALL IGAUS1(FL,PSI,CHI)
cc      If(KLLERR) Return   IGAUS1 has no error conditons               sld01
C
C     GENERATE "NSIM*10" RANDOM DEVIATES FROM THE INVERSE GAUSSIAN
C     DISTRIBUTION WITH PARAMETERS MU AND LAMBDA
C
      NOBS = 2*MAXTB
      DO 3 I = 1,NOBS
         XX(I) = 1./IGAUSF()
    3 CONTINUE
C
C     SORT THE RANDOM DEVIATES TO FORM THE EDF
C
      CALL SIFT(XX,NOBS)
cc      If(KLLERR) Return  SIFT has no error conditons                  sld01
C
C     GENERATE THE DESIRED SAMPLE BY SAMPLING FROM THE EDF
C
      DO I = 1, N
         R = PROBINC * RNUMLHS1() + STRTPT
cc         NR = MAX( 1, NINT(R*NOBS))                                   SLD
cc         RES = XX(NR)                                                 SLD
         FNR = R*(NOBS-1) + 1                                           SLD
         NR = INT(FNR)                                                  SLD
         NR1 = NR + 1                                                   SLD
cccc                                                                    SLD                                                         SLD
         FRAC = FNR - FLOAT(NR)                                         SLD
         RES = FRAC*(XX(NR1)-XX(NR)) + XX(NR)                           SLD
cccc         RES = FRACTION(FNR)*(XX(NR1)-XX(NR)) + XX(NR)
cccc                                                                    SLD
         X(LOC(I,J)) = DMAX(RES,1.0D-10)
         IF (IRS == 0) STRTPT = DBLE(I) / DBLE(N)
      END DO
C
      RETURN
      END
