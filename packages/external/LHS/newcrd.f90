C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   9 May 101   12:14 pm
c
c  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NewCrd
      Subroutine NewCrd(Card,IUnit,IEnd)
cc    IEnd is an "end" flag                                             sld01
cc       IEnd = 0 for normal return                                     sld01
cc       IEnd = 1 indicates end of file reached                         sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE        not needed					sld01
cc    NewCrd is called from:  LREAD,READ,RDPAR2                         sld01
cc    NewCrd does not call any other external routines                  sld01
c
c     Read a new card, left justify it, strip off trailing comments,
c     and convert it to upper case and convert tabs and commas to blanks
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Character*(*) Card
c
c     -- Read a card and make sure it's not blank
c
      IEnd = 0
 100  Read (IUnit,101,End=999) Card
 101  Format (A256)
      If (Card == ' ') Go To 100
c
c     -- Now strip off trailing comments
c
      IComt = Index(Card,'$')
      If (IComt == 1) Then
c        -- record contains only a comment - reject it and read another
         Go To 100
      Else If (IComt > 1) Then
c        -- trailing comment found - strip it off
         Card = Card(1:IComt-1)
      End If
c
c     -- Now convert the card to upper case and convert tabs (ASCII 9)
c     -- and commas to blank spaces
c
      Maxi = Len(Card)                                                  sld01
      Do i=1, Maxi
         If ( IChar(Card(i:i)) > 96  .AND.  IChar(Card(i:i)) < 123 )
     1         Card(i:i) = Char( IChar(Card(i:i)) - 32 )
         If ( IChar(Card(i:i)) == 9 )  Card(i:i) = ' '
         If ( Card(i:i) == ',' )       Card(i:i) = ' '
      End Do
c
c     -- Left justify the card
c
cc      Maxi = Len(Card)                                                sld01
      Do i=1, Maxi
         If ( Card(i:i) /= ' ' ) Exit
      End Do
c
      Card = Card(i:Maxi)
c
      Return
c
c     -- Transfer to here and set a flag if the end of the file was reached
c
 999  IEnd = 1
      Return
c
      End
