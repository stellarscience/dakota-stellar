C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 May 101    9:01 am
C****************************************************************
C SUBROUTINE TO GENERATE GAMMA DEVIATES WHEN THE SHAPE PARAMETER
C ALPHA > 1
C   WRITTEN BY DO LE MINH
C           DEPARTMENT OF MANAGEMENT SCIENCE
C           CALIFORNIA STATE UNIVERSITY, FULLERTON CA 92634
C           PH: (714) 773-2221
C REFERENCE: DO LE MINH (1988). "GENERATING GAMMA VAIRATES,"
C            ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C            VOL. 14, NO. 3, SEPT. 1988, 261-266.
C
C FORMERLY CALLED SUBROUTINE MINH IN RON'S PROGRAM
c
c changes made Dec 1995 by S. L. Daniel
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::GAMMAM
      SUBROUTINE GAMMAM(ALPHA,X)
cc    GAMMAM is called from routine:  GAMMA                             sld01
cc    GAMMAM is calls routine:  RNUMLHS2 (instead of RNUMLHS1)                SLD
cc      because of addition of second random number generator           SLD
cc      for use in Acceptance-Rejection scheme found in GAMMA & IGAUS   SLD
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE           -- not needed                            sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A,D,DL,X1,X2,X4,X5,XLL,XLR,F1,F2,F4,F5,
     1                 P1,P2,P3,P4,U,W,XP
      DOUBLE PRECISION X
      DATA ASAVE/-1./
      SAVE
C
cc    changed back to original code, calling routine determines that    sld01
cc    alpha is greater than 1                                           sld01
c      IF (ALPHA .NE. ASAVE) THEN                                       sld01
cc       if (ALPHA .eq. ASAVE .and. ASAVE .eq. -1.) then                sld01
cc         write (99,*) '***** fatal error in GAMMAM *****',            sld01
cc     x      ' ALPHA and ASAVE  = -1 causes division by 0'             sld01
cc         KLLERR = .TRUE.                                              sld01
cc         RETURN                                                       sld01
cc       else                                                           sld01
      IF (ALPHA .NE. ASAVE) THEN                                        sld01
         ASAVE = ALPHA
         A = ALPHA -1.
         D = SQRT(A)
         IF (ALPHA.LE.2) THEN
            DL = A/2
            X1 = 0
            X2 = DL
            XLL = -1
            F1 = 0
         ELSE
            DL = D-0.5
            X2 = A - DL
            X1 = X2 - DL
            XLL = 1. - A/X1
c            F1 = EXP(A*ALOG(X1/A) + 2*DL)
            F1 = EXP(A*LOG(X1/A) + 2*DL)
         END IF
c         F2 = EXP(A*ALOG(X2/A) + DL)
         F2 = EXP(A*LOG(X2/A) + DL)
         X4 = A + D
         X5 = X4 + D
         XLR = 1.0 - A/X5
c         F4 = EXP(A*ALOG(X4/A) - D)
c         F5 = EXP(A*ALOG(X5/A) - 2*D)
         F4 = EXP(A*LOG(X4/A) - D)
         F5 = EXP(A*LOG(X5/A) - 2*D)

         P1 = 2*F4*D
         P2 = 2*F2*DL + P1
         P3 = F5/XLR + P2
         P4 = -F1/XLL + P3
      END IF
C
cc  100 U = RNUMLHS1() * P4                                                SLD
  100 U = RNUMLHS2() * P4                                                  SLD
C
      IF (U.LE.P1) THEN
         W = U/D - F4
         IF (W.LE.0) THEN
            X = A + U/F4
            RETURN
         END IF
         IF (W.LE.F5) THEN
            X = X4 + (W*D)/F5
            RETURN
         END IF
cc         X = X4 + RNUMLHS1()*D                                           SLD
         X = X4 + RNUMLHS2()*D                                             SLD
         XP = 2.*X4 - X
         IF (W.GE.(F4+ ((F4-1.)/(X4-A))*(X-X4))) THEN
            X = XP
            RETURN
         END IF
         IF (W.LE.(A/X4-1)*F4*(X-X4)+F4) RETURN
         IF (W.LT.2*F4-1) GO TO 500
c         IF (W.LT.(2*F4-EXP(A*ALOG(XP/A)+A-XP)))  GO TO 500
         IF (W.LT.(2*F4-EXP(A*LOG(XP/A)+A-XP)))  GO TO 500
         X = XP
         RETURN
      END IF
C
      IF (U.LE.P2) THEN
         W  = (U - P1)/DL - F2
         IF (W.LE.0) THEN
            X = A - (U-P1)/F2
            RETURN
         END IF
         IF (W.LE.F1) THEN
            X = X1 + (W*DL)/F1
            RETURN
         END IF                                                         SLD
cc         X = X1 + RNUMLHS1()*DL                                          SLD
         X = X1 + RNUMLHS2()*DL
         XP = 2.*X2 - X
         IF (W.GE.(F2+ ((F2-1.)/(X2-A))*(X-X2))) THEN
            X = XP
            RETURN
         END IF
         IF (W.LE.F2*(X-X1)/DL) RETURN
         IF (W.LT.2*F2-1) GO TO 500
c         IF (W.LT.(2*F2-EXP(A*ALOG(XP/A)+A-XP)))  GO TO 500
         IF (W.LT.(2*F2-EXP(A*LOG(XP/A)+A-XP)))  GO TO 500
         X = XP
         RETURN
      END IF
C
C     THE TWO EXPONENTIAL REGIONS.
C
cc      W = RNUMLHS1()                                                     SLD
      W = RNUMLHS2()                                                       SLD
      IF (U.LE.P3) THEN
         U = (P3-U) / (P3-P2)
c         X = X5 - ALOG(U) / XLR
         X = X5 - LOG(U) / XLR
         IF (W.LT. (XLR*(X5-X)+1.)/U) RETURN
         W = W*F5*U
      ELSE
         U = (P4-U)/(P4-P3)
c         X = X1 - ALOG(U)/XLL
         X = X1 - LOG(U)/XLL
         IF (X.LE.0) GO TO 100
         IF (W.LT. (XLL*(X1-X)+1.)/U) RETURN
         W = W*F1*U
c         X = X1 - ALOG(U)/XLL
         X = X1 - LOG(U)/XLL
         IF (X.LE.0) GO TO 100
         IF (W.LT. (XLL*(X1-X)+1.)/U) RETURN
         W = W*F1*U
      END IF
C
c  500 IF (ALOG(W).GT.A*ALOG(X/A)+A-X) GO TO 100
  500 IF (LOG(W).GT.A*LOG(X/A)+A-X) GO TO 100
C
      RETURN
      END
