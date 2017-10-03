C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    9:26 am
C****************************************************************
C SUBROUTINE CUMULC IS USED TO GENERATE USER-SPECIFIED CUMULATIVE
C CONTINUOUS DISTRIBUTION USING LINEAR OR LOGARITHMIC INTERPOLATION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CUMULC
      SUBROUTINE CUMULC(J,IDT)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE   not neeed					sld01
cc    CUMULC is called from routine:  LHS                               sld01
cc    CUMULC calls routine:  INTERP                                     sld01
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
C      This statement removed to make modules work - gdw-96
C      DIMENSION XTABLE(MAXTB,2)
C      EQUIVALENCE (XTABLE(1,1), XX(1))
C
      LOC(I,J) = (J-1)*N+I
C
C     -- INTFL IS SET TO 0 FOR LINEAR INTERPOLATION
      INTFL=0
      IF (IDT .EQ. 10) INTFL=1
      PROBINC = 1./FLOAT(N)
      IF (IRS .NE. 0) PROBINC=1.0
      READ (8) NP
      READ (8) (XTABLE(III,1),XTABLE(III,2),III=1,NP)
      STRTPT = 0.
      IMIN=1
      DO I = 1, N
         PROB = PROBINC*RNUMLHS1()+ STRTPT
         CALL INTERP(PROB,BX,XTABLE,MAXTB,IMIN,NP,INTFL)
cc         If(KLLERR) Return -- INTERP has no error condition		sld01
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
