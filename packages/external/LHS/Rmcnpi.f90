C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  24 May 101    7:26 am
C****************************************************************
C SUBROUTINE RMCNPI PROVIDED INITIAL CONTROL FOR THE PSEUDO-RANDOM
C NUMBER SEQUENCE FOR THE MCNP RANDOM NUMBER GENERATOR.
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::RMCNPI
      SUBROUTINE RMCNPI(ISEED)
cc    only 2001 sld changes were comments                               sld01
cc    RMCNPI is called from:  RNUMLHS1                                     sld01
cc    RMCNPI does not call any other external routines                  sld01
C
C     TAKEN FROM MCNP, AUGUST 1990
C
C     Changes for DLL compatibility
      USE FIRSTS, ONLY: ISARG
cc    FIRSTS initializes ISARG = 0                                      sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MCNPRN/RANS,RANB,RANI,RANJ,NRN
cc    Common MCNPRN is shared with routine RMCNP                       sld01
C
      PARAMETER (P=16777216.,Q=1./P,FB=4867484.,FS=10256733.)
C
C     Changes for DLL compatibility - now uses module here.
C     DATA ISARG / 0 /
C
      IF(ISARG.GT.0) RETURN
      ISARG = 1
      NRN=0
C     -- SET COUNTER TO ISEED
      NRN = ISEED
c      R = DFLOAT(ISEED)
      R = DBLE(ISEED)
      RANI=AINT(R*Q)
      RANJ=R-RANI*P
      RANB=RANI
      RANS=RANJ
      RETURN
C
      END
