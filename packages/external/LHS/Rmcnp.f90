C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 Mar 101    7:48 am
C****************************************************************
C FUNCTION RMCNP RETURNS THE NEXT PSEUDO-RANDOM NUMBER.
C IT IS PART OF THE MCNP RANDOM NUMBER GENERATOR.
C
      DOUBLE PRECISION FUNCTION RMCNP(ISEED)
cc    only 2001 sld changes were comments                               sld01
cc    RMCNP is called from:  RNUMLHS1                                      sld01
cc    RMCNP does not call any other external routines                   sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MCNPRN/RANS,RANB,RANI,RANJ,NRN
cc    Common MCNPRN is shared with routine RMCNPI                       sld01
C
C ----------------------------------------------------------------------
C
      PARAMETER (P=16777216.,Q=1./P,R=Q*Q,GB=1136868.,GS=6328637.)
C
C     -- THIS ROUTINE IS PORTABLE ONLY IF GB+GS<2**24
C
      A=GS*RANS
      B=GB*RANS+GS*RANB+AINT(A*Q)
      RANS=A-INT(A*Q)*P
      RANB=B-INT(B*Q)*P
      RMCNP=(RANB*P+RANS)*R
      NRN=NRN+1
C     -- PUT NRN IN ISEED
      ISEED = NRN
C
      RETURN
C
      END
