C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  S     4 Dec 96   12:19 pm
C
C      COMMON /FTNFILE/ CmdLin
C      Character*128 CmdLin
C
      MODULE FLNAME
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER (LEN=128) :: CmdLin = ' '
C
      END MODULE
