C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:32 am
C****************************************************************
C SUBROUTINE DATOUT WILL OUTPUT THE SAMPLE AND ITS RANKS IF
C REQUESTED
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::DATOUT
      SUBROUTINE DATOUT
cc    DATOUT is called from routine:  LHS                               sld01
cc    DATOUT calls routines:  OUTDAT,RANKER                             sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE           not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,NV                                            sld01
C     INCLUDE 'CRANK.INC'                                               GDW-96  
      USE CRANK                         
cc    CRANK provides:  XV and RXV arrays                                sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X and XSAVE arrays                               sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    Statement Function:                                               sld01
      LOC(I,J)=(J-1)*N+I
C
      NNV=N*NV
      DO 10 I=1, NNV
         X(I)=XSAVE(I)
   10 CONTINUE
      CALL OUTDAT(0)
cc      If(KLLERR) Return   OUTDAT has no error conditions              sld01
      DO 620 J=1,NV
        DO 610 I=1,N
  610   XV(I)=X(LOC(I,J))
        CALL RANKER
cc        If(KLLERR) Return  RANKER has no error conditions             sld01
        DO 620 I=1,N
        X(LOC(I,J))=RXV(I)
  620 CONTINUE
      CALL OUTDAT(1)
cc      If(KLLERR) Return    OUTDAT has no error conditions             sld01
      RETURN
      END
