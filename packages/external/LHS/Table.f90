C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  GDW  25 Apr 101    1:23 pm
C****************************************************************
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::TABLE
      SUBROUTINE TABLE(FUNCT, XTABLE, ISMAX, ISIZE)
cc    only 2001 sld changes were comments                               sld01
cc    TABLE is called from routine:  BETA with BETAFN as the external   sld01
cc                                   function passed in call list       sld01
cc    TABLE calls routine:  FUNCT which is BETAFN for LHS code          sld01
c
c     Note: this routine ASSUMES that the function being evaluated
c     accepts values 0 to 1 inclusive, and returns vlaues 0 to 1 inclusive.
c
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTABLE(ISMAX,2)
      EXTERNAL FUNCT
      SAVE
C
      DeltaMx = 1.0 / ISIZE
      DeltaMn = DeltaMx / 10000.
      PosDir = 1.0
c
c     -- If the table gets too large, return to this point and
c     -- start the table generation process over.
c
 10   Continue
      I = 1
      XTABLE(I,1) = 0.0
      XTABLE(I,2) = 0.0
      X = 0.0
      Q = 0.0
      FOFX=0.0
      DELTAX = DeltaMn
      DELTAQ = DeltaMn
c
      Do While (Q < 1.0  .AND.  X < 1.0)
c
c        -- Select a step size DeltaQ so that the maximum step size is
c        -- taken near the center of the distribution (DeltaQ = 1/ISIZE)
c        -- and progressively smaller steps are taken in the tails of the
c        -- distribution.  In no case, let the step size become less than
c        -- the value implied by ISIZE, called DELTAMN.
c
         DeltaQ = MIN( DeltaMx, Q/2.0, (1.0-Q)/4.0 )
c        -- make sure DeltaQ is bigger than minimum step size
         DeltaQ = MAX( DeltaQ, DeltaMn )
c
c        -- Q is the maximum value of FOFX that will be accepted
c        -- make sure DeltaQ doesn't push Q past 1.0
         IF ( Q + DeltaQ > 1.0 ) THEN
            DeltaQ = 1.0 - Q
            Q = 1.0
         ELSE
            Q = Q + DELTAQ
         END IF
c
c        -- XNEW is the estimate for a value of X that will produce
c        -- an acceptable Q
         XNEW = X + DELTAX
         IF ( XNew > 1.0 ) THEN
            XNEW = 1.0
            DELTAX = XNEW - X
         END IF
C
         CALL FUNCT(XNEW,FOFX)
         If(KLLERR) Return
c
         IF ( FofX > Q ) THEN
c
c           -- Cut step size so that FofX does not exceed Q
c
c           -- May have to repeat to allow multiple step size reductions
            DO WHILE ( FofX > Q )
               DELTAX = 0.6666667 * DeltaX
               XNEW = X + DeltaX
               IF ( XNew > NEAREST(X, PosDir) ) THEN
c                 -- DeltaX not losing significance
                  CALL FUNCT(XNEW,FOFX)
                  If(KLLERR) Return
               Else
c                 -- DeltaX is losing significance.  Ensure that XNew
c                 -- is greater than X by at least 1 bit, then evaluate
c                 -- the function and stop cutting the time step.  We
c                 -- have to accept this value for DeltaX even if it
c                 -- produces too big a change in Q because we are only
c                 -- changing X by one bit.
                  XNew = NEAREST(X, PosDir)
                  DeltaX = XNew - X
                  CALL FUNCT(XNEW,FOFX)
                  If(KLLERR) Return
                  Exit
               END IF
            END DO
c
         ELSE IF ( (FofX + 0.5 * DeltaQ) < Q ) THEN
c
c           -- We are taking a much smaller step than necessary.
c           -- We will accept this answer, but specify a larger step
c           -- for the next entry in the table.
            DeltaX = 1.5 * DeltaX
c
         Else
c
c           -- Step size is appropriate -- don't change anything
c
         END IF
c
c
c        -- Check that adding this entry to the table won't overflow it
         I=I+1
         IF ( I >= ISMAX-2 ) THEN
c           -- Table is too big, so adjust the minimum and maximum step
c           -- size to be larger (to make a smaller table) and go back
c           -- to the top and start over.
            DeltaMx = DeltaMx * 2.0
            DeltaMn = DeltaMn * 2.0
            GOTO 10
         END IF
c
c        -- Now add the value to the table and prepare for another step
         XTABLE(I,1)=XNEW
         XTABLE(I,2)=FOFX
         X = XNEW
         Q = FOFX
c
      End Do
c
c     -- Ensure that the last value in the table is identically (1.0, 1.0)
c     -- to make sure that no numerical roundoff errors occur in the last step
c     -- (note that the Do loop cannot exit unless one of these values is
c     -- already 1.0, so this makes sure that the table is complete).
c
      IF ( XTABLE(I,1) /= 1.0  .OR.  XTABLE(I,2) /= 1.0 ) THEN
         I = I + 1
         XTABLE(I,1)=1.0
         XTABLE(I,2)=1.0
      END IF

      ISIZE=I
C
      RETURN
      END
