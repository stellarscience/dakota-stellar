C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    9:29 am
C****************************************************************
C SUBROUTINE CUMULD IS USED TO GENERATE USER-SPECIFIED CUMULATIVE
C DISCRETE DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CUMULD
      SUBROUTINE CUMULD(J)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE    not needed					sld01
cc    CUMULD is called from routine:  LHS                               sld01
cc    CUMULD calls routine:  INTRPD                                     sld01

C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
C     INCLUDE 'CWORKX.INC'                                              GDW-96  
      USE CWORKX                        
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      These statements removed to make modules work - GDW-96
C      DIMENSION XTABLE(MAXTB,2)
C      EQUIVALENCE (XTABLE(1,1), XX(1))
C
      LOC(I,J) = (J-1)*N+I
C
      PROBINC = 1./FLOAT(N)
      IF (IRS .NE. 0) PROBINC=1.0
      READ (8) NP
      READ (8) (XTABLE(III,1),XTABLE(III,2),III=1,NP)
      STRTPT = 0.
      IMIN=1
      DO I = 1, N
         PROB = PROBINC*RNUMLHS1()+ STRTPT
         CALL INTRPD(PROB,BX,XTABLE,MAXTB,IMIN,NP)
cc         If(KLLERR) Return -- INTRPD has no error conditons           sld01
         X(LOC(I,J)) = BX
         IF(IRS.EQ.0) THEN
            STRTPT = DBLE(I) / DBLE(N)
         ELSE
            IMIN=1
         END IF
      END DO
C
      RETURN
      END
