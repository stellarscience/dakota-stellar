C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:21 am
C****************************************************************
C SUBROUTINE PARETO GENERATES PARETO DISTRIBUTIONS
C WITH PARAMETERS ALPHA > 2 AND BETA > 0
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PARETO
      SUBROUTINE PARETO(J)
cc    PARETO is called from routine:  LHS                               sld01
cc    PARETO calls routine:  RNUMLHS1                                      sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE    -- not needed					sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS                                           sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X array                                          sld01
cc    RNUMLHS1 is an external function                                     sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    statement function:                                               sld01
      LOC(I,J)=(J-1)*N+I
C
      PROBINC=1./FLOAT(N)
      IF (IRS .EQ. 1) PROBINC=1.0
      READ (8) ALPHA,BETA
      STRTPT=0.
      DO I = 1, N
         R = PROBINC*RNUMLHS1()+STRTPT
         RES = BETA/(1.-R)**(1./ALPHA)
         X(LOC(I,J)) = DMAX(RES, 1.0D-10)
         IF (IRS == 0) STRTPT = DBLE(I) / DBLE(N)
      END DO
C
      RETURN
      END
