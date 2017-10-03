C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 May 101    9:51 am
C SUBROUTINE RMCNPI2 PROVIDES INITIAL CONTROL FOR THE PSEUDO-RANDOM
C NUMBER SEQUENCE FOR THE MCNP RANDOM NUMBER GENERATOR USED IN THE
C ACCEPTANCE-REJECTION SCHEME USED IN ROUTINES GAMMA AND IGAUS
C
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::RMCNPI2
      SUBROUTINE RMCNPI2(JSEED)
cc    RMCNPI2 is called from:  RNUMLHS2
cc    RMCNPI2 does not call any other external routines
C
C     TAKEN FROM MCNP, AUGUST 1990
C
C     Changes for DLL compatibility
      USE FIRSTS, ONLY: JSARG
cc    Module FIRSTS initializes JSARG = 0
cc    IGAUS and GAMMA also initialize JSARG = 0
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MCNPRN2/RANS2,RANB2,RANI2,RANJ2,NRN2
cc    Common MCNPRN2 is shared with routine RMCNP2
C
      PARAMETER (P=16777216.,Q=1./P)
C
C
      IF(JSARG.GT.0) RETURN
      JSARG = 1
      NRN2=0
C     -- SET COUNTER TO JSEED
      NRN2 = JSEED
c
      R = DBLE(JSEED)
      RANI2=AINT(R*Q)
      RANJ2=R-RANI2*P
      RANB2=RANI2
      RANS2=RANJ2
      RETURN
C
      END



