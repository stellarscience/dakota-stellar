C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 Mar 101    8:19 am
C****************************************************************
C SUBROUTINE SIFT WILL SORT A VECTOR OF DATA IN INCREASING ORDER
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::SIFT
      SUBROUTINE SIFT (XV,N)
cc    only 2001 sld changes were comments                               sld01
cc    SIFT is called from routines:  IGAUS,CMCRD,GAMMA,HISTO,LHS	sld01
cc    SIFT calls no other external routines                             sld01
c      INCLUDE 'KILLFILE.INC'  --  include not needed 10-96
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XV(N)
      M=N
   10 M=M/2
      IF (M) 30,20,30
   20 RETURN
   30 K=N-M
      J=1
   40 I=J
   50 L=I+M
      IF (XV(I)-XV(L)) 70,70,60
   60 A=XV(I)
      XV(I)=XV(L)
      XV(L)=A
      I=I-M
      IF (I) 70,70,50
   70 J=J+1
      IF (J-K) 40,40,10
      END
