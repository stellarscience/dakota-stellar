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
C SUBROUTINE UNIFRM GENERATES THE UNIFORM, LOGUNIFORM, UNIFORM*,
C OR THE LOGUNIFORM* DISTRIBUTIONS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::UNIFRM
      SUBROUTINE UNIFRM(J,IDT)
cc    UNIFRM is called from routine:  LHS                               sld01
cc    UNIFRM calls routine:  RNUMLHS1                                      sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE       -- not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS                                           sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X array                                          sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    Statement Function:                                               sld01
      LOC(I,J)=(J-1)*N+I
C
      I=0
      IF(IDT.EQ.6.OR.IDT.EQ.7)THEN
        READ(8)NINT
      ELSE
        NINT=1
        NPT=N
      ENDIF
      DO 40 K=1,NINT
        IF(IDT.EQ.4.OR.IDT.EQ.5)THEN
           READ(8)A,B
        ELSE
           READ(8)NPT,A,B
           IF(NPT.EQ.0)GO TO 40
        ENDIF
        IF(IDT.EQ.5.OR.IDT.EQ.7)THEN
          A=LOG10(A)
          B=LOG10(B)
        ENDIF
        PROBINC=(B-A)/FLOAT(NPT)
        IF(IRS.NE.0)PROBINC=B-A
        DO 30 ID=1,NPT
          I=I+1
          IF(IRS.EQ.0)THEN
            STRTPT=(ID-1)*PROBINC
          ELSE
            STRTPT=0.0
          ENDIF
          X(LOC(I,J))=A+STRTPT+PROBINC*RNUMLHS1()
          IF(IDT.EQ.5.OR.IDT.EQ.7)X(LOC(I,J))=10.**(X(LOC(I,J)))
   30   CONTINUE
   40 CONTINUE
      RETURN
      END
