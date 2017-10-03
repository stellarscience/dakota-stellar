C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  28 Mar 101    7:51 am
C****************************************************************
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::INTERP
      SUBROUTINE INTERP (Y, X, XTABLE, ISMAX, IMIN, IMX1, INTFL)
cc    INTERP is called from routines:  BETA,CUMULC                      sld01
cc    INTERP does not call any other external routines                  sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE     -- not needed					sld01
cc                                                                      sld01
cc    X is returned to calling program                                  sld01
cc                                                                      sld01
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTABLE(ISMAX,2)
C
C     INTERPOLATE ON Y TO RETURN AN X VALUE
C
C     IF INTFL = 0  USE LINEAR INTERPOLATION.  OTHERWISE USE LOGARITHMIC.
C
C     XTABLE(I,1) IS X
C     XTABLE(I,2) IS Y
C
      IMAX=IMX1
C     -- This line changed to correct out of bounds condition 3/5/96
C     IF ( IMIN+1 .EQ. IMAX ) IMIN=IMAX-2
      IF ( IMIN+1 .GE. IMAX ) IMIN = MAX( IMAX-2, 1 )
C
  100 I = IMIN + (IMAX-IMIN)/2
      IF ( Y .LT. XTABLE(I,2) ) THEN
         IMAX=I
      ELSE
         IMIN=I
      END IF
      IF (IMIN+1 .LT. IMAX) GO TO 100
C
C     --- NOW INTERPOLATE
C
      IF (INTFL .EQ. 0) THEN
C
C        --- LINEAR INTERPOLATION:
         X = XTABLE(IMIN,1) + (XTABLE(IMAX,1)-XTABLE(IMIN,1))*
     1      (Y-XTABLE(IMIN,2))/(XTABLE(IMAX,2)-XTABLE(IMIN,2))
C
      ELSE
C
C        --- FOR LOGARITHMOC INTERPOLATION:
         XMX=LOG(XTABLE(IMAX,1))
         XMN=LOG(XTABLE(IMIN,1))
         X = XMN + (XMX-XMN) * (Y-XTABLE(IMIN,2)) /
     1               (XTABLE(IMAX,2)-XTABLE(IMIN,2))
         X=EXP(X)
C
      END IF
C
      RETURN
      END
