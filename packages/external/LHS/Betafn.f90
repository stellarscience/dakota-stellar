C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   3 Apr 101    9:04 am
C****************************************************************
C SUBROUTINE BETAFN IS USED IN GENERATING A BETA DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::BETAFN
      SUBROUTINE BETAFN(X,FOFX)
cc    BETAFN is called from:  TABLE but passed in call from BETA        sld01
cc    BETAFN calls routine:  BETAIC                                     sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE
cc    KILLFILE provides:  KLLERR                                        sld01
cc                                                                      sld01
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(1)
      COMMON /PQ/P,Q,NZ
      IF(X.LT.0.)X=0.
      IF(X.GT.1.)X=1.
      OMX = 1.-X
      N=1
      CALL BETAIC(X,OMX,P,Q,N,Y,NZ)
cc    Test NZ flag for error condition:                                 sld01
      IF (NZ .ne. 0) THEN                                               sld01
         WRITE(4,9001)                                                  sld01
         WRITE(99,9001)                                                 sld01
         KLLERR = .True.                                                sld01
         Return                                                         sld01
      END IF                                                            sld01
cc      If(KLLERR) Return                                               sld01

      FOFX = Y(1)
      RETURN
 9001 FORMAT(' Error condition returned from BETAIC, flag NZ /= 0')     sld01
      END
