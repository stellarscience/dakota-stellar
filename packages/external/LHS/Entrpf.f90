C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101   11:00 am
C****************************************************************
C FUNCTION ENTRPF IS USED IN THE GENERATION OF THE MAXIMUM ENTROPY
C DISTRIBUTION.  IT CODIFIES THE NONLINEAR FUNCTION THAT IS SOLVED
C USING BISECTION IN THE SUBROUTINE ENTRPY.  THIS FUNCTION IS:
C
C                B*EXP(B*BETA) - A*EXP(A*BETA)              1
C        AMU = ---------------------------------    -    ------
C                  EXP(B*BETA) - EXP(A*BETA)              BETA
C
C AND IS REWRITTEN FOR THIS ROUTINE AS:
C
C                B - A*EXP((A-B)*BETA)             1
C        AMU = -------------------------    -    ------
C                  1 - EXP((A-B)*BETA)            BETA
C
C
      DOUBLE PRECISION FUNCTION ENTRPF(BETA,A,AMU,B)
cc    only 2001 sld changes were comments                               sld01
cc    ENTRPF is called from:  ENTRPY                                    sld01
cc    ENTRPF does not call any other external routines                  sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION EAB,TERM,BETA1
C
      BETA1=BETA
      EAB=EXP(BETA1*(A-B))
      TERM=(B-A*EAB)/(1.0-EAB)
      ENTRPF=TERM-AMU-1.0/BETA1
C
      RETURN
      END
