C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  28 Mar 101   10:53 am
C****************************************************************
C SUBROUTINE FINVNOR IS USED IN GENERATING THE NORMAL AND
C LOGNORMAL DISTRIBUTIONS
C
      DOUBLE PRECISION FUNCTION FINVNOR (X)
      USE KILLFILE                                                      sld01
cc    FINVNOR is called from:  MIX,NORMAL                               sld01
cc    FINVNOR calls routine:  RIERFC1                                   sld01
cc    RIERFC1 is an External Function                                   sld01
cc                                                                      sld01
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (X-0.5) 20,10,30
   10 FINVNOR=0.0
      RETURN
   20 F=-1.0
      Y=X
      GO TO 40
   30 F=1.0
      Y=1.0-X
   40 FINVNOR=SQRT(2.0)*F*RIERFC1(2.*Y)
cc    Functon RIERFC1 has an error condition                            sld01
      IF (KLLERR) Return                                                sld01
      RETURN
      END
