C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  14 May 101    9:32 am
C****************************************************************
C SUBROUTINE RANKER IS USED TO FIND THE RANKS OF A VECTOR OF DATA
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::RANKER
      SUBROUTINE RANKER
cc    RANKER is called from routines:  COROUT,DATOUT,MIX                sld01
cc    RANKER calls routine:  HPSRT                                      sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE          -- not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                        sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N                                               sld01
C     INCLUDE 'CRANK.INC'                                               GDW-96  
      USE CRANK                         
cc    CRANK provides:  XV,RXV,IWK arrays                                sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DO 10 I=1,N
         RXV(I)=FLOAT(I)
   10 CONTINUE
C
      CALL HPSRT
cc      If(KLLERR) Return    HPSRT has no error conditions              sld01
      DO 20 I=1,N
cc         IWK(I)=IFIX(RXV(I))                                          sld01
         IWK(I)=INT(RXV(I))
cc                                                          14May01     sld01
         RXV(I)=FLOAT(I)
   20 CONTINUE
C
C     FIND TIES
C
      I=0
   30 I=I+1
      IF (I.GE.N) GO TO 80
      IF (XV(I).NE.XV(I+1)) GO TO 30
C
C     COUNT TIES
C
      NTIES=2
   40 II=I+NTIES
      IF (II .LE. N) THEN
         IF (XV(I) .EQ. XV(II)) THEN
            NTIES=NTIES+1
            GO TO 40
         END IF
      END IF
C
C     AVERAGE TIED RANKS
C
      AVG=0.0
      DO 60 J=1,NTIES
         AVG=AVG+RXV(I+J-1)
   60 CONTINUE
      AVG=AVG/FLOAT(NTIES)
      I=I-1
      DO 70 J=1,NTIES
         I=I+1
         RXV(I)=AVG
   70 CONTINUE
      GO TO 30
C
   80 CONTINUE
C
C     REORDER
C
      DO 200 I=1, N-1
   90    K=IWK(I)
         IF (K.NE.I) THEN
            XHOLD=XV(I)
            RHOLD=RXV(I)
            XV(I)=XV(K)
            RXV(I)=RXV(K)
            XV(K)=XHOLD
            RXV(K)=RHOLD
            IWK(I)=IWK(K)
            IWK(K)=K
            GO TO 90
         END IF
  200 CONTINUE
C
      RETURN
      END
