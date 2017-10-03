C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 May 101    8:53 am
C****************************************************************
C SUBROUTINE TO GENERATE GAMMA DEVIATES WHEN THE SHAPE PARAMETER
C ALPHA < 1
C REFERENCE: D. J. BEST (1983). "A NOTE OF GAMMA VARIATE GENERATORS
C            WITH SHAPE PARAMETER LESS THAN UNITY," COMPUTING, 30, 185-188.
C
C FORMERLY NAMED BEST IN RON'S PROGRAM
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::GAMMAB
      SUBROUTINE GAMMAB(ALPHA,X,Z,B)
cc    GAMMAB is called from routine:  GAMMA                             sld01
cc    GAMMAB calls routine:  RNUMLHS2 (instead of RNUMLHS1)                   SLD
cc      because of addition of second random number generator           SLD
cc      for use in Acceptance-Rejection scheme found in GAMMA & IGAUS   SLD
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE        -- not needed				sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X,Z,B,P,Y
C
cc    1 P = B*RNUMLHS1()                                                   SLD
cc      US = RNUMLHS1()                                                    SLD
    1 P = B*RNUMLHS2()                                                     SLD
      US = RNUMLHS2()                                                      SLD
C
      IF (P.LE.1.) THEN
cc   Note:  ALPHA was tested in calling program, ALPHA > 1              sld01
         X = Z*P**(1./ALPHA)
         IF (US .LE. (2.-X)/(2.+X)) RETURN
         IF (US .GT. EXP(-X)) GO TO 1
      ELSE
         X = -LOG(Z*(B-P)/ALPHA)
         Y = X/Z
         IF (US*(ALPHA + Y - ALPHA*Y) .LT. 1.) RETURN
         IF (US .GT. Y**(ALPHA-1.)) GO TO 1
      END IF
C
      RETURN
      END
