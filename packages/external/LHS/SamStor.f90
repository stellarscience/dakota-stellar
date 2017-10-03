C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   11:02 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::SamStor
      SUBROUTINE SamStor(MAXVAR,MAXOBS,MAXNAM,
     x LSTDNAM,INDXNAM,PTVALST,SMATX)
c     This routine stores information to returned in LHS_RUN
c     call list
c      SamStor calls routines:
c
      USE CPARAM
c     CPARAMS provides:  ISeed,N,NV,IREP
      USE InByCall
cc    InByCall provides: LINT,LPREP,LFILES,LDIST,NNames
      USE CSAMP
c     CSAMP provides: X and XSAVE arrays
      USE CRANK
c     CRANK provides: XV array
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION :: SMATX(MAXVAR,MAXOBS)
      DOUBLE PRECISION :: PTVALST(MAXNAM)
      INTEGER :: INDXNAM(MAXNAM)
      CHARACTER(LEN=16) :: LSTDNAM(MAXNAM)
c
C     -- STATEMENT FUNCTION
      LOC(I,J)=(J-1)*N+I
c
c
c     store observations in matrix form
      DO I = 1,N
         do J = 1,NV
            SMATX(J,I) = X(LOC(I,J))
         end do
      END DO
c
cccc  Next section modified from SAMOUT routine
c     -- Now loop over variable names and
c     -- determine their point values
c

c     i = 0
c     Do While ( IVarNm(i+1) /= 0 )
c        i = i + 1
cccc MSE/MPA: modify to traverse arrays to their known ends (MAXNAM).
cccc The original logic requires hard-coded array lengths larger than
cccc the actual data with trailing 0's from unused memory.
      Do 200 i = 1,MAXNAM
         LSTDNAM(i) = List(i)

         If ( IVarNm(i) == -9999999 )  Then
            INDXNAM(i) = 0
         ELSE IF (IVarNm(i) > 0 ) THEN
c           -- do nothing -- it already has the right value
            INDXNAM(i) = IVarNm(i)
         Else
            IVabNm = ABS(IVarNm(i))
            INDXNAM(i) = IVarNm(IVabNm)
         END If
ccc
c
c        -- Store them in the order encountered in the input file
c
c        -- If it is a "Same As" distribution, then do nothing and cycle
         If ( IVarNm(i) < 0  .AND.  IVarNm(i) /= -9999999 ) Cycle
c
c        -- Otherwise, calculate the point value
c
         If ( IVarNm(i) == -9999999 ) Then
c           -- If it is a constant, just store it
            Value = PValue(i)
c
         Else If (IPtVal == 0) Then
c           -- store the optional point values
            Value = PValue(i)
c
         Else If (IPtVal == 1 .AND. IVarNm(i) /= 0 ) Then
c
c           -- calculate the mean
            iFirst = Loc(1,IVarNm(i))
            iLast = Loc(N,IVarNm(i))
C$$$        Value = 0.0 ! Modified by slbrow, 2-19-2004
            Value = 0.0D0
            Do j=iFirst, iLast
               Value = Value + X(j)
            End Do
            Value = Value / N
c
         Else If ( iVarNm(i) /= 0 ) Then
c
c           -- calculate the median - Note that XSave is already sorted
            iFirst = Loc(1,IVarNm(i))
            iLast = Loc(N,IVarNm(i))
            i2 = iFirst + (N/2)
            If ( 2 * ((iLast-iFirst) / 2) == (iLast-iFirst) ) Then
c              -- Even number of samples
               i1 = i2 - 1
               Value = 0.5 * ( XSave(i1) + XSave(i2) )
            Else
c              -- Odd number of samples
               Value = XSave(i2)
            End If
c
         End If
ccc
         PTVALST(i) = Value
ccc
c        -- Check for other distributions that are the same as this one
cccc         Line = List(i)
cccc         NextP = 21

c        j = 0
c        Do While ( IVarNm(j+1) /= 0 )
c           j = j + 1
cccc MSE/MPA: modify to traverse arrays to their known ends (MAXNAM).
cccc The original logic requires hard-coded array lengths larger than
cccc the actual data with trailing 0's from unused memory.
         Do 100 j = 1,MAXNAM
            If ( IVarNm(j) == -i ) Then
cccc               If ( NextP > 50 ) Then
cccc                  Write (1,807) Line(1:NextP)
cccc                  Line = ' '
cccc                  NextP = 4
cccc               End If
cccc               Line(NextP:) = List(j)
               PTVALST(j) = Value
cccc               NextP = NextP + 17
            End If
c        End Do
 100     Continue
c
c
c        -- Write out the variable name and its point value
cc         Write (1,804) Line(1:NextP), Value

c        -- process the next variable
c     End Do
 200  Continue

c
      Return

      END SUBROUTINE
