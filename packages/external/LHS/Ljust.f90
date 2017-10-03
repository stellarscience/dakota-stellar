C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    9:42 am
c
c  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LJust
      Subroutine LJust(Card)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE   -- LJust has no error conditions			sld01
cc    LJust is called from:  READ,RDPAR2                                sld01
cc    LJust does not call any other external routines                   sld01
c
c     Left justify the character variable
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Character*(*) Card
c
c     -- Left justify the card
c
      Maxi = Len(Card)
      Do i=1, Maxi
         If ( Card(i:i) /= ' ' ) Exit
      End Do
c
      If (i > Maxi) i = maxi
      Card = Card(i:Maxi)
c
      Return
      End
