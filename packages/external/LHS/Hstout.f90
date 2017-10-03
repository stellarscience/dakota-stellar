C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:30 am
C****************************************************************
C SUBROUTINE HSTOUT IS USED TO GENERATE HISTOGRAMS OF THE
C VARIABLES
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::HSTOUT
      SUBROUTINE HSTOUT
cc    HSTOUT is called from routine:  LHS                               sld01
cc    HSTOUT calls routine:  HISTO                                      sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE          not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,NV,IDIST,TITLE                                sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X and XSAVE arrays                               sld01
C     INCLUDE 'CRANK.INC'                                               GDW-96  
      USE CRANK                         
cc    CRANK provides:  XV array                                         sld01
C     INCLUDE 'DISTNM.INC'                                              GDW-96  
      USE DISTNM                        
cc    DISTNM provides:  DIST,IDSEND,IDSST arrays                   	sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    Statement Function:                                               sld01
      LOC(I,J)=(J-1)*N+I
C
      NNV=N*NV
      DO 10 I=1, NNV
         X(I)=XSAVE(I)
   10 CONTINUE
      DO 590 I=1,NV
         IDT=IDIST(I)
         WRITE(4,9001)TITLE
         WRITE(4,9002)I,DIST(IDSST(IDT):IDSEND(IDT))
         DO 530 J=1,N
            XV(J)=X(LOC(J,I))
  530    CONTINUE
         CALL HISTO
cc         If(KLLERR) Return HISTO has no error conditions              sld01
  590 CONTINUE
C
      RETURN
C
 9001 FORMAT('1',3X,A)
 9002 FORMAT('0','  HISTOGRAM FOR VARIABLE NO.',I3,5X,
     1        A,'DISTRIBUTION')
C
      END
