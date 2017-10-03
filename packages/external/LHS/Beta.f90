C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  GDW  25 Apr 101   12:36 pm
C****************************************************************
C SUBROUTINE BETA IS USED TO GENERATE A BETA DISTRIBUTION ON THE
C INTERVAL (A,B) AND WITH PARAMETERS P AND Q
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::BETA
      SUBROUTINE BETA(J)
cc    BETA is called from routine:  LHS                                 sld01
cc    BETA calls routines:  ERXSET,INTERP,TABLE with function BETAFN    sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARMS provides:  MAXTB                                            sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS                                           sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X and XSAVE arrays                               sld01
C     INCLUDE 'CWORKX.INC'                                              GDW-96  
      USE CWORKX                        
cc    CWORKX provides:  XTABLE array                                    sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    Labeled Common PQ shared with routine BETAFN                      sld01
      COMMON /PQ/P,Q,NZ
C      This statement removed to make modules work - gdw-96
C      DIMENSION XTABLE(MAXTB,2)
C      EQUIVALENCE (XTABLE(1,1), XX(1))
      EXTERNAL BETAFN
C
cc    Statement Function:                                               sld01
      LOC(I,J) = (J-1)*N+I
C
      CALL ERXSET(10,0)
      If(KLLERR) Return
      PROBINC = 1./FLOAT(N)
      IF(IRS.NE.0)PROBINC=1.0
      READ (8)A,B,P,Q
      STRTPT = 0.
      ISIZE=250
      CALL TABLE(BETAFN,XTABLE,MAXTB,ISIZE)
      If(KLLERR) Return
      IMIN=1
      DO I = 1, N
         PROB = PROBINC*RNUMLHS1()+ STRTPT
         CALL INTERP(PROB,BX,XTABLE,MAXTB,IMIN,ISIZE,0)
cc         If(KLLERR) Return  -- INTERP has no error conditions         sld01
         X(LOC(I,J)) = A + (B-A)*BX
         IF (IRS .EQ. 0) THEN
            STRTPT = DBLE(I) / DBLE(N)
         ELSE
            IMIN=1
         END IF
      END DO
C
      RETURN
      END
