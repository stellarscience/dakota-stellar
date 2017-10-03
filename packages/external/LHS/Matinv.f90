C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:20 am
C****************************************************************
C SUBROUTINE MATINV IS USED TO INVERT A LOWER TRIANGULAR MATRIX
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::MATINV
      SUBROUTINE MATINV
cc    MATINV is called from routine:  MIX                               sld01
cc    MATINV does not call any other external routines                  sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE     -- not needed					sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  NV                                              sld01
C     INCLUDE 'CWORKC.INC'                                              GDW-96  
      USE CWORKC                        
cc    CWORKC provides:  Q array                                         sld01
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOC1(I,J)=J+(I*I-I)/2
C
      DO 10 I=1,NV
         Q(LOC1(I,I))=1.0/Q(LOC1(I,I))
   10 CONTINUE
C
      DO 100 K=NV, 2, -1
         DO 90 J=K-1, 1, -1
            TEMP=0.0
            JPLUS=J+1
            DO 40 I=JPLUS,K
               TEMP=TEMP+Q(LOC1(K,I))*Q(LOC1(I,J))
   40       CONTINUE
            Q(LOC1(K,J))=-TEMP*Q(LOC1(J,J))
   90    CONTINUE
  100 CONTINUE
C
      RETURN
      END
