C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 May 101    7:40 am
C****************************************************************
C FUNCTION TO GENERATE RANDOM INVERSE GAUSSIAN DEVIATES
C
C   WRITTEN BY:
C           PROF. A. C. ATKINSON
C           THE LONDON SCHOOOL OF ECONOMICS AND POLITICAL SCIENCE
C           DEPARTMENT OF STATISTICAL AND MATHEMATICAL SCIENCES
C           HOUGHTON STREET,
C           LONDON WC2A 2AE
C           PH: (071) 405-7686
C
C REFERENCE: A. C. ATKINSON (1982) "THE SIMULATION OF GENERALIZED INVERSE
C            GAUSSIAN AND HYPERBOLIC RANDOM VARIABLES,"
C            SIAM J. SCI. STAT. COMPUT.
C            VOL. 3, NO. 4, DEC. 1982, 502-515.
C
C     THIS FUNCTION NAMED GIGLT1 IN RON'S PROGRAM
C
      DOUBLE PRECISION FUNCTION IGAUSF()
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION FLM1,PV,SV,QV,R,V,CON1,CON2,CON3,HC1,HC2,HC3,
     1       FLGS1,FLGS2,FLGS3,SIMULT1
      COMMON/IGAUSC/FLM1,PV,SV,QV,R,V,CON1,CON2,CON3,HC1,HC2,HC3,
     1              FLGS1,FLGS2,FLGS3,SIMULT1
C
C     ALGORITHM GIGLT1: SIMULATES THE GENERALIZED INVERSE GAUSSIAN
C     DISTRIBUTION WITH 0 <= LAMBDA < 1 BY REJECTION WITH A THREE
C     PART EXPONENTIAL ENVELOPE.  CONSTANTS ARE SET IN SUBROUTINE
C     IGAUS1.  NOTE IN THE CALLING PROGRAM FOR IGAUS1, LAMBDA IS
C     REPRESENTED AS FL.
C
cc    COMMON IGAUSC is shared with routine:  IGAUS1                     sld01
cc    IGAUSF is called from routine:  IGAUS                             sld01
cc    IGAUSF calls routine:  RNUMLHS2 (changed from RNUMLHS1)                 SLD
ccc                                                                     SLD
ccc    1 UVAR = RNUMLHS1()                                                 SLD
ccc      USTAR = RNUMLHS1()                                                SLD
    1 UVAR = RNUMLHS2()                                                    SLD
      USTAR = RNUMLHS2()                                                   SLD
ccc                                                                     SLD
      IF (UVAR .LE. R) THEN
         GIGLT1 = SV*LOG(1. + CON1*UVAR)
         IF (LOG(USTAR) .GT. FLM1*LOG(GIGLT1) - HC1/GIGLT1 - HC2*GIGLT1
     1     - FLGS1) GO TO 1
      ELSE IF (UVAR .LE. V) THEN
         GIGLT1 = -PV*LOG(CON2 - SIMULT1*UVAR)
         IF (LOG(USTAR) .GT. FLM1*LOG(GIGLT1) - HC1/GIGLT1 - HC3*GIGLT1
     1     - FLGS2) GO TO 1
      ELSE
         GIGLT1 = CON3 - QV*LOG(1. - UVAR)
         IF(LOG(USTAR) .GT. FLM1*LOG(GIGLT1) - HC1/GIGLT1 - FLGS3)
     1      GO TO 1
      END IF
C
      IGAUSF = GIGLT1
C
      RETURN
      END
