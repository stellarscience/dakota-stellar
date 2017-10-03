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
C SUBROUTINE TRIANG GENERATES THE TRIANGULAR DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::TRIANG
      SUBROUTINE TRIANG(J)
cc    TRIANG is called from routine:  LHS                               sld01
cc    TRIANG calls routine:  RNUMLHS1                                      sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE             not needed				sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                        sld01
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
      PROBINC=1./FLOAT(N)
      IF(IRS.EQ.1)PROBINC=1.0
      READ(8)A,B,C
      C1=C-A
      C2=(B-A)/C1
      STRTPT=0.
      DO I = 1, N
        R=PROBINC*RNUMLHS1()+STRTPT
        IF(R.LE.C2)THEN
          X(LOC(I,J))=A+SQRT(R*C1*(B-A))
        ELSE
          X(LOC(I,J))=C-SQRT((1.-R)*C1*(C-B))
        ENDIF
        IF (IRS == 0) STRTPT = DBLE(I) / DBLE(N)
      END DO
C
      RETURN
      END
