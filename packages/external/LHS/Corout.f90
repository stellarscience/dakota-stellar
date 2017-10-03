C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Jun 101   10:27 am
C****************************************************************
C SUBROUTINE COROUT WILL OUTPUT THE RAW AND RANK CORRELATION
C MATRICES OF THE SAMPLE IF REQUESTED
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::COROUT
      SUBROUTINE COROUT
cc    COROUT is called from routine:  LHS                               sld01
cc    COROUT calls routines:  CORCAL,PMTRX,VIF,RANKER                   sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE   	      not needed				sld01
cc      Routines: CORCAL,PMTRX,VIF,RANKER have no error conditons       sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,ICM,NV                                        sld01
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                        
cc    CCMATR provides:  LCM array                                       sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X and XSAVE arrays                               sld01
C     INCLUDE 'CRANK.INC'                                               GDW-96  
      USE CRANK                         
cc    CRANK provides:  XV and RXV arrays                                sld01
      USE InByCall                                                      SLD
cc    InByCall provides: LINT,LPREP,LFILES,LDIST,LRUN,NNames            SLD
cc                       and vectors VCTR1 and VCTR2:                   SLD
cc                       VCTR1(NVAR*(NVAR+1)),VCTR2(NVAR*(NVAR+1))            SLD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    Statement Function:                                               sld01
      LOC(I,J)=(J-1)*N+I
C
      NNV=N*NV
      DO 10 I=1, NNV
         X(I)=XSAVE(I)
   10 CONTINUE
      DO 600 I=1,NV
  600 LCM(I)=I
      CALL CORCAL
cc      If(KLLERR) Return   CORCAL has no error condtions               sld01

cc    store raw correlation data for retrieval by LHS_COROUT            SLD
      kkk = (NV*(NV+1))/2                                               SLD
      do jjj = 1,kkk                                                    SLD
         VCTR2(jjj) = CORR(jjj)                                         SLD
      end do                                                            SLD
cc
      CALL PMTRX(NV,1)
cc      If(KLLERR) Return  PMTRX has no error conditions                sld01
c      IF(N.GT.NV.AND.ICM.EQ.0)CALL VIF
      IF(N.GT.NV.AND.ICM.EQ.0) THEN
         CALL VIF
cc         If(KLLERR) Return   VIF has no error conditions              sld01
      ENDIF
      DO 620 J=1,NV
        DO 610 I=1,N
  610   XV(I)=X(LOC(I,J))
        CALL RANKER
cc        If(KLLERR) Return   RANKER has no error conditions            sld01
        DO 620 I=1,N
        X(LOC(I,J))=RXV(I)
  620 CONTINUE
      CALL CORCAL
cc      If(KLLERR) Return   CORCAL has no error condtions               sld01

cc    store rank correlation data for retrieval by LHS_COROUT           SLD
      do jjj = 1,kkk                                                    SLD
         VCTR2(kkk+jjj) = CORR(jjj)                                     SLD
      end do                                                            SLD
cc                                                                      SLD
      CALL PMTRX(NV,2)
cc      If(KLLERR) Return    PMTRX has no error conditions              sld01
      IF(N.GT.NV.AND.ICM.EQ.0)CALL VIF
      RETURN
      END
