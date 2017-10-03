C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:25 am
C****************************************************************
C SUBROUTINE GUMBEL GENERATES GUMBEL DISTRIBUTIONS
C WITH PARAMETERS ALPHA AND BETA
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::GUMBEL
      SUBROUTINE GUMBEL(J)
cc    GUMBEL is called from:  LHS                                       sld01
cc    GUMBEL calls routine:  RNUMLHS1                                   sld01
cc    RNUMLHS1 is an External Function                                  sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE           -- not needed				sld01
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
      PROBINC=1.0D0/FLOAT(N)
      IF (IRS .EQ. 1) PROBINC=1.0D0
      READ (8) ALPHA,BETA
cc    CHKDAT routine determines that ALPHA,BETA are >=.001 & <=1.E+7    sld01
      A = 1.0D0/ALPHA
      STRTPT=0.0D0
      DO I = 1, N
         R=PROBINC*RNUMLHS1()+STRTPT
         RES =  BETA-A*(LOG(-1.0D0*LOG(R)))
         X(LOC(I,J)) = RES
C         PRINT *, 'R= ', R, 'RES = ', RES, 'X = ', X(LOC(I,J))
         IF (IRS == 0) STRTPT = DBLE(I) / DBLE(N)
      END DO
C
      RETURN
      END
