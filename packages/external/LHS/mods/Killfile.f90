C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  S     4 Dec 96   12:18 pm
C
C        COMMON /KERRCHK/ KLLERR, IKLERR
C        LOGICAL KLLERR
C
      MODULE KILLFILE
C
        INTEGER :: IKLERR = 0
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::KLLERR
        LOGICAL :: KLLERR = .FALSE.
C
      END MODULE
