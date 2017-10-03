C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 May 101    8:27 am
C FUNCTION RMCNP2 RETURNS THE NEXT PSEUDO-RANDOM NUMBER.
C IT IS PART OF THE MCNP RANDOM NUMBER GENERATOR USED IN THE
C ACCEPTANCE-REJECTION SCHEME IN ROUTINES GAMMA AND IGAUS
C
      DOUBLE PRECISION FUNCTION RMCNP2(JSEED)
cc    RMCNP2 is called from:  RNUMLHS2
cc    RMCNP2 does not call any other external routines
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MCNPRN2/RANS2,RANB2,RANI2,RANJ2,NRN2
cc    Common MCNPRN2 is shared with routine RMCNPI2
C
C ----------------------------------------------------------------------
C
      PARAMETER (P=16777216.,Q=1./P,R=Q*Q,GB=1136868.,GS=6328637.)
C
C     -- THIS ROUTINE IS PORTABLE ONLY IF GB+GS<2**24
C
C
      A=GS*RANS2
      B=GB*RANS2+GS*RANB2+AINT(A*Q)
      RANS2=A-INT(A*Q)*P
      RANB2=B-INT(B*Q)*P
      RMCNP2=(RANB2*P+RANS2)*R
      NRN2=NRN2+1
C     -- PUT NRN2 IN JSEED
      JSEED = NRN2
C
      RETURN
C
      END
