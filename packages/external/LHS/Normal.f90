C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:18 am
C****************************************************************
C SUBROUTINE NORMAL GENERATES NORMAL AND LOGNORMAL DISTRIBUTIONS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NORMAL
      SUBROUTINE NORMAL(J,IDT)
cc    NORMAL is called from routine:  LHS                               sld01
cc    NORMAL calls routines:  FINVNOR,RNUMLHS1                             sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS                                           sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides: X,XSAVE                                           sld01

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOC(I,J)=(J-1)*N+I
C
      ILOG=0
      ATOP=1.0
      ABOT=0.0
      ITopFlg = 0
      IBotFlg = 0
C
      IF (IV1 .EQ. 1  .OR.  IDT .EQ. 27  .OR.  IDT .EQ. 28) THEN
C        --- VERSION 1 INPUT COMPATIBILITY AND (LOG)NORMAL-B DISTRIBUTIONS
C        --- A IS THE LOWER .001 QUANTILE OF THE NORMAL DISTRIBUTION.
C        --- B IS THE UPPER .001 QUANTILE OF THE NORMAL DISTRIBUTION.
         READ(8)A,B
         IF(IDT .EQ. 3  .OR.  IDT .EQ. 28)THEN
            A=LOG(A)
            B=LOG(B)
            ILOG=1
         ENDIF
         PMU=(A+B)/2.
         SIG=(B-PMU)/3.0902323062
         ATOP=0.999
         ABOT=0.001
C
      ELSE IF (IDT .EQ. 2  .OR.  IDT .EQ. 29) THEN
C        -- READ PARAMETERS OF THE UNBOUNDED NORMAL DISTRIBUTION DIRECTLY
         READ (8) PMU,SIG
         IF (IDT .EQ. 29) ILOG=1
C
      ELSE IF (IDT .EQ. 3) THEN
C        -- READ MEAN AND ERROR FACTOR FOR UNBOUNDED LOGNORMAL DISTRIBUTION
         READ (8) PMEAN,ERF
         ILOG=1
         SIG=LOG(ERF)/1.645
         PMU=LOG(PMEAN)-0.5*SIG*SIG
C
      ELSE IF (IDT .EQ. 31  .OR.  IDT .EQ. 35) THEN
C        -- READ PARAMETERS FOR TRUNCATED (LOG)NORMAL(-N) DISTRIBUTION
         READ (8) PMU,SIG,ABOT,ATOP
         IF (IDT .EQ. 35) ILOG=1
C
      ELSE IF (IDT .EQ. 33) THEN
C        -- READ PARAMETERS FOR TRUNCATED LOGNORMAL DISTRIBUTION
         READ (8) PMEAN,ERF,ABOT,ATOP
         ILOG=1
         SIG=LOG(ERF)/1.645
         PMU=LOG(PMEAN)-0.5*SIG*SIG
C
      ELSE
C
C        --- ALL BOUNDED DISTRIBUTIONS ARE HERE
C
         IF (IDT .EQ. 30  .OR.  IDT .EQ. 34) THEN
C           -- READ PARAMETERS FOR BOUNDED (LOG)NORMAL(-N) DISTRIBUTION
            READ (8) PMU,SIG,ALOBND,AHIBND
            IF (IDT .EQ. 34) THEN
               ILOG=1
               ALOBND=LOG(ALOBND)
               AHIBND=LOG(AHIBND)
            END IF
C
         ELSE
C           -- READ PARAMETERS FOR BOUNDED LOGNORMAL DISTRIBUTION (IDT=32)
            READ (8) PMEAN,ERF,ALOBND,AHIBND
            ILOG=1
            SIG=LOG(ERF)/1.645
            PMU=LOG(PMEAN)-0.5*SIG*SIG
            ALOBND=LOG(ALOBND)
            AHIBND=LOG(AHIBND)
         END IF
C
C        --- CONVERT BOUNDS TO QUANTILES BY BISECTION
C
C        -- FIRST, UPPER BOUND
C
         NITER=0
         BHI=0.999999
         BLO=0.000001
         BL1 = 0.001
         RESHI=FINVNOR(BHI)*SIG+PMU-AHIBND
         IF (KLLERR) Return                                             sld01
cc       FINVNOR has no error conditions but Function RIERFC1 called	sld01
cc               by it does                                             sld01
         RESLO=FINVNOR(BL1)*SIG+PMU-AHIBND
         IF (KLLERR) Return                                             sld01
CC       -- err gdw 7/31/00   ITopFlg = 0
         If ( ResHi < 0 ) Then
            ATop = 1.0
            ITopFlg = 1
         Else If ( ResLo > 0 ) Then
C           --ERROR IN BISECTION PARAMETERS
            WRITE (4,9998) J
            WRITE (99,9998) J
            KLLERR = .TRUE.
            RETURN
         Else
C
C           -- Perform the bisection to locate ATop
C
 100        NITER=NITER+1
            IF (NITER .GT. 1000) THEN
C              - PROBABLY AN INFINITE LOOP
               WRITE (4,9997)
               WRITE (99,9997)
               KLLERR = .TRUE.
               RETURN
            END IF
            BMD=0.50*(BHI+BLO)
            RESMD=FINVNOR(BMD)*SIG+PMU-AHIBND
         IF (KLLERR) Return                                             sld01
            IF (RESLO*RESMD .GT. 0.0) THEN
               BLO=BMD
            ELSE
               BHI=BMD
            END IF
            IF (BHI/BLO .GT. 1.00001) GO TO 100
C
            ATOP=0.50*(BHI+BLO)
C
         End If
C
C        -- NOW, LOWER BOUND
C
         NITER=0
         BHI=0.999999
         BH1 = 0.999
         BLO=0.000001
         RESHI=FINVNOR(BH1)*SIG+PMU-ALOBND
         IF (KLLERR) Return                                             sld01
         RESLO=FINVNOR(BLO)*SIG+PMU-ALOBND
         IF (KLLERR) Return                                             sld01
CC       -- err gdw 7/31/00     IBotFlg = 0
         If ( ResHi < 0 ) Then
C           --ERROR IN BISECTION PARAMETERS
            WRITE (4,9999) J
            WRITE (99,9999) J
            KLLERR = .TRUE.
            RETURN
         Else If ( ResLo > 0 ) Then
            ABot = 0.0
            IBotFlg = 1
         Else
C
C           -- Perform the bisection to locate ABot
C
 200        NITER=NITER+1
            IF (NITER .GT. 1000) THEN
C              - PROBABLY AN INFINITE LOOP
               WRITE (4,9997)
               WRITE (99,9997)
               KLLERR = .TRUE.
               RETURN
            END IF
            BMD=0.50*(BHI+BLO)
            RESMD=FINVNOR(BMD)*SIG+PMU-ALOBND
         IF (KLLERR) Return                                             sld01
            IF (RESLO*RESMD .GT. 0.0) THEN
               BLO=BMD
            ELSE
               BHI=BMD
            END IF
            IF (BHI/BLO .GT. 1.00001) GO TO 200
C
            ABOT=0.50*(BHI+BLO)
C
         End If
C
C        --- END OF BISECTION SECTION
C
      END IF
C
      PROBINC = (ATOP - ABOT) / DBLE(N)
      IF(IRS.EQ.1)PROBINC=ATOP-ABOT
      STRTPT=ABOT
      DO I = 1, N
        R=PROBINC*RNUMLHS1()+STRTPT
        R=DMIN(DMAX(R,1.0D-10),9.9999997D-01)
        X(LOC(I,J))=FINVNOR(R)*SIG+PMU
         IF (KLLERR) Return                                             sld01
C
C       --- Handle cases where bounded distribution bounds are outside
C       --- of the range 0.000001 to 0.999999
C
        If ( ITopFlg == 1 ) Then
           X(Loc(I,J)) = Min( AHiBnd, X(Loc(I,J)) )
        End If
        If ( IBotFlg == 1 ) Then
           X(Loc(I,J)) = Max( ALoBnd, X(Loc(I,J)) )
        End If
C
        IF(ILOG.EQ.1)X(LOC(I,J))=EXP(X(LOC(I,J)))
        IF (IRS==0) STRTPT = DBLE(I) * (ATOP - ABOT) / DBLE(N) + ABOT
      END DO
C
      RETURN
C
 9997 FORMAT('1',10x,'Program error: Normal Routine Bisection ',
     1        'did not converge.')
 9998 FORMAT('1',10X,'Upper bound of a bounded normal or lognormal',
     1  /,10X,'distribution must be greater than the 0.001 quantile.',
     2  /,10X,'Found in Distribution #', I5)
 9999 FORMAT('1',10X,'Lower bound of a bounded normal or lognormal',
     1  /,10X,'distribution must be less than the 0.999 quantile.'
     2  /,10X,'Found in Distribution #', I5)
C
      END
