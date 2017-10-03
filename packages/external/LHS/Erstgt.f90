C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   2 Apr 101    9:32 am
C****************************************************************
C SUBROUTINE ERSTGT IS AN ERROR CHECKING ROUTINE USED IN
C GENERATING A BETA DISTRIBUTION
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ERSTGT
      SUBROUTINE ERSTGT(K,NFATAL,NTRACE)
cc    ERSTGT is called from from:  ERRGET,ERRSET                        sld01
cc    ERSTGT does not call any other external routines                  sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE                 not needed                         sld01
      USE FIRSTS, ONLY: LNF, LNT
cc    FIRSTS initialized with LNF=-1, LNT=0                             sld01
C
C     THIS ROUTINE IS A SLAVE TO ERRGET AND ERRSET WHICH KEEPS
C     THE FLAGS AS LOCAL VARIABLES.
C
C     *** IF LOCAL VARIABLES ARE NOT NORMALLY RETAINED BETWEEN
C     CALLS ON THIS SYSTEM, THE VARIABLES LNF AND LNT CAN BE
C     PLACED IN A COMMON BLOCK AND PRESET TO THE FOLLOWING
C     VALUES IN THE MAIN PROGRAM.
C
C       Changes for DLL compatibility - use FIRSTS module now.
C     DATA LNF/-1/,LNT/0/
C
cc    Note:  Values are "returned" if called with K>0,                  sld01
cc           Values are "set" if called with K<=0                       sld01
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (K.LE.0) THEN
         LNF = NFATAL
         LNT = NTRACE
      ELSE
         NFATAL = LNF
         NTRACE = LNT
      END IF
      RETURN
      END
