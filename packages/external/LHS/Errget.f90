C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   2 Apr 101    9:38 am
C****************************************************************
C SUBROUTINE ERRGET IS AN ERROR CHECKING ROUTINE USED IN
C GENERATING A BETA DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ERRGET
      SUBROUTINE ERRGET(NFATAL,NTRACE)
cc    ERRGET is called from:  BETAIC,BETALN,HYPGEO                      sld01
cc    ERRGET calls routine:  ERSTGT                                     sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C
C     ABSTRACT
C         ERRGET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK.
C         ERRGET ASSIGNS TO NFATAL AND NTRACE RESPECTIVELY THE VALUES
C         OF NF AND NT IN COMMON BLOCK MLBLK0 THEREBY ASCERTAINING THE
C         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK.
C
C     DESCRIPTION OF ARGUMENTS
C         BOTH ARGUMENTS ARE OUTPUT ARGUMENTS OF DATA TYPE INTEGER.
C         NFATAL - CURRENT VALUE OF NF (SEE DESCRIPTION OF ERXSET.)
C         NTRACE - CURRENT VALUE OF NT (SEE DESCRIPTION OF ERXSET.)
C
cc    call to retrieve error flags:
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CALL ERSTGT(1,NFATAL,NTRACE)
cc      If(KLLERR) Return    ERSTGT has no error conditions             sld01
      IF (NFATAL < 0) THEN                                              sld01
         KLLERR = .True.                                                sld01
      END IF                                                            sld01
      RETURN                                                            sld01
      END
