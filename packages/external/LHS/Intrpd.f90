C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101   10:47 am
C****************************************************************
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::INTRPD
      SUBROUTINE INTRPD (Y, X, XTABLE, ISMAX, IMIN, ITOP)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE         not needed					sld01
cc    INTRPD is called from:  BINOM,CUMULD,GEOM,HGEOM,NBINOM,POISON     sld01
cc    INTRPD does not call any other external routines                  sld01
cc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTABLE(ISMAX,2)
C
C     -- 'INTERPOLATE' ON A TABLE OF CUMULATIVE PROBABILITIES
C     -- FOR A DISCRETE PROBABILITY DISTRIBUTION (REAL)
C
C     XTABLE(I,1) IS X
C     XTABLE(I,2) IS Y
C
C     -- SET UP BOUNDS ON TABLE FOR INTERPOLATION
C
      IMAX = ITOP
C     -- This line changed to correct out of bounds condition 3/5/96
C     IF ( IMIN+1 .EQ. IMAX ) IMIN=IMAX-2
      IF ( IMIN+1 .GE. IMAX ) IMIN = MAX( IMAX-2, 1 )
C
C     -- SEARCH TABLE USING BISECTION
C
      IF (Y .LE. XTABLE(IMIN,2)) THEN
         IVAL=IMIN
      ELSE
  100    I = IMIN + (IMAX-IMIN)/2
         IF ( Y .LE. XTABLE(I,2) ) THEN
            IMAX=I
         ELSE
            IMIN=I
         END IF
         IF (IMIN+1 .LT. IMAX) GO TO 100
         IVAL=IMAX
      END IF
C
C     -- SELECT THE PROPER VALUE AND RETURN
C
      X=XTABLE(IVAL,1)
      IMIN=IVAL
C
      RETURN
      END
