C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C*******************************************************************
C FUNCTIONS DMIN() AND DMAX() ARE DOUBLE PRECISION VERSIONS OF MIN()
C AND MAX() THAT AVOID TYPE COERSION AND LOSS OF PRECISION WHEN USED
C WITH DOUBLE PRECISION ARGUMENTS
C*******************************************************************
      DOUBLE PRECISION FUNCTION DMIN(D1,D2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION D1, D2
C
      IF (D1.LE.D2) THEN
        DMIN = D1
      ELSE
        DMIN = D2
      ENDIF
C
      RETURN
      END
C
C*******************************************************************
C
      DOUBLE PRECISION FUNCTION DMAX(D1,D2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION D1, D2
C
      IF (D1.GE.D2) THEN
        DMAX = D1
      ELSE
        DMAX = D2
      ENDIF
C
      RETURN
      END
