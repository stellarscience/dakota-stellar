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
C SUBROUTINE CORCAL COMPUTES A CORRELATION MATRIX FROM THE SAMPLE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CORCAL
      SUBROUTINE CORCAL
cc    CORCAL is called from:  COROUT,MIX                                sld01
cc    CORCAL does not call any other external routines                  sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE         -- not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  NV,N                                            sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X array                                          sld01
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                   
cc   CCMATR provides:  CORR array                                       sld01
c
      USE LOCALVARS, ONLY: XM, SSQ
cc    LOCALVARS provides:  XM and SSQ arrays                            sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     To make compatible with modules - gdw-96
c     Moved to the LOCALVARS module
C      DIMENSION XM(NVAR),SSQ(NVAR)
C      DOUBLE PRECISION, ALLOCATABLE :: XM(:), SSQ(:)
C
cc    Statement Functions:                                              sld01
      LOC(I,J)=(J-1)*N+I
      LOC1(I,J)=J+(I*I-I)/2
C
c     Added for modules
c     Moved to the LOCALVARS module
C      ALLOCATE( XM(NVAR), SSQ(NVAR) )
C
      IF(NV.EQ.1)THEN
        CORR(1)=1.0
c       Added for modules
c     Moved to the LOCALVARS module
C        DEALLOCATE( XM, SSQ )
        RETURN
      ENDIF
      DO 10 I=1,NV
        XM(I)=0.0
        SSQ(I)=0.0
   10 CONTINUE
      K=(NV*(NV+1))/2
      DO 20 I=1,K
         CORR(I)=0.0
   20 CONTINUE
      FN=N
C
C     COMPUTE THE MEAN
C
      DO 30 J=1,NV
         DO 30 K=1,N
            XM(J)=XM(J)+X(LOC(K,J))
   30 CONTINUE
      DO 40 J=1,NV
         XM(J)=XM(J)/FN
   40 CONTINUE
C
C     SUBTRACT THE MEAN FROM ALL OF THE ELEMENTS OF THE MATRIX
C     AND COMPUTE THE SUM OF SQUARES
C
      DO 50 J=1,NV
         DO 50 K=1,N
           X(LOC(K,J))=X(LOC(K,J))-XM(J)
           SSQ(J)=SSQ(J)+X(LOC(K,J))*X(LOC(K,J))
   50 CONTINUE
C
C     COMPUTE THE CORRELATION
C
      DO 70 I=2,NV
         IL1=I-1
         DO 70 J=1,IL1
            DO 60 K=1,N
               CORR(LOC1(I,J))=CORR(LOC1(I,J))+X(LOC(K,I))*X(LOC(K,J))
   60       CONTINUE
   70 CONTINUE
C
      DO 80 I=2,NV
         IL1=I-1
         DO 80 J=1,IL1
            IF (CORR(LOC1(I,J)) .NE. 0.0)
     1         CORR(LOC1(I,J))=CORR(LOC1(I,J))/SQRT(SSQ(I)*SSQ(J))
   80 CONTINUE
C
      DO 90 I=1,NV
   90    CORR(LOC1(I,I))=1.0
C
c     Added for modules
c     Moved to the LOCALVARS module
c      DEALLOCATE( XM, SSQ )
C      
      RETURN
      END
