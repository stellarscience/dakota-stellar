C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Jun 101   10:25 am
C****************************************************************
C LATIN HYPERCUBE OR RANDOM SAMPLING PROGRAM
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C              ISSUED BY SANDIA NATIONAL LABORATORIES
C                   A PRIME CONTRACTOR TO THE
C                UNITED STATES DEPARTMENT OF ENERGY
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C UNITED STATES GOVERNMENT. NEITHER THE UNITED STATES NOR THE
C UNITED STATES DEPARTMENT OF ENERGY
C NOR ANY OF THEIR EMPLOYEES, NOR ANY OF THEIR CONTRACTORS,
C SUBCONTRACTORS, OR THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS
C OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY
C FOR THE ACCURACY, COMPLETENESS OR USEFULNESS OF ANY INFORMATION,
C APPARATUS, PRODUCT OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS
C USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C AUTHORS   SANDIA NATIONAL LABORATORIES, ALBUQUERQUE, NM 87185
C R. L. IMAN                PHONE (505)844-8834   ORG. 6415
C M. J. SHORTENCARIER       PHONE (505)846-1662   ORG. 7223
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C      PROGRAM LHS
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS
      SUBROUTINE LHS
cc    LHS is called from:  LHSDLL or LHSDRV (for stand alone)           sld01
cc    LHS calls routines:  FILEOC,READ,BANNER,SAMOUT,WRTPAR,PMTRX,	sld01
cc                         POSDEF,CHLSKY,BETA,NORMAL,UNIFRM,TRIANG,     sld01
cc                         CUMULC,CUMULD,POISON,GEOM,BINOM,NBINOM,      sld01
cc                         HGEOM,EXPON,WEIBUL,PARETO,GAMMA,IGAUS,	sld01
cc                         ENTRPY,SIFT,DATOUT,HSTOUT,COROUT             sld01
cc                         GUMBEL,FRECHET                               lps08
 
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,NV,NREP,IRP,ICM,ICORR,IDATA,IHIST,IDIST       sld01
C     INCLUDE 'CSAMP.INC'                                               GDW-96  
      USE CSAMP                         
cc    CSAMP provides:  X and XSAVE arrays                               sld01
C     INCLUDE 'CRANK.INC'                                               GDW-96  
      USE CRANK                         
cc    CRANK provides:  XV array                                         sld01
C     INCLUDE 'CCMATR.INC'                                              GDW-96  
      USE CCMATR                        
cc    CCMATR provides:  CORR array                                      sld01
C     INCLUDE 'CWORKC.INC'                                              GDW-96  
      USE CWORKC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc    CWORKC provides:  Q and S arrays                                  sld01
C====================
      CHARACTER*11 TIMRES
C====================
c     next lines added 1-17-96 to convert date/time routine to LF90
      integer :: dt(8)
      character (len=5) :: zone
      character * 8 DATRES
      character * 2 Chr2
C
C     -- STATEMENT FUNCTIONS
C
      LOC(I,J)=(J-1)*N+I
      LOC1(I,J)=J+(I*I-I)/2
C
C     -- HERE ROUTINE FILEOC OPENS THE INPUT, OUTPUT AND WORKING FILES
C
C     -- ROUTINE READ READS IN THE PARAMETER STATEMENTS AND DEFINES
C     -- THE VARIABLES IN COMMON /PARAM/
C
C     -- ROUTINES BANNER AND WRTPAR ECHO A DESCRIPTION OF THE SAMPLE
C     -- AS INPUT BY THE USER AND A DESCRIPTION OF THE DISTRIBUTIONS
C     -- USED WITH THE INPUT VARIABLES
C
ccc***The following information is duplicated in LHS_INIT
      CALL FILEOC(1)
cc      If(KLLERR) Return   FILEOC has no error conditions              sld01
ccc***End of information duplicated in LHS_INIT
      CALL READ
C====================
      If(KLLERR) Return
ccc***The following information is duplicated in LHS_PREP
      Call date_and_time (DATRES,TIMRES,zone,dt)
      Chr2 = TIMRES(8:9)
C
      If (IPrint > 0 ) Write (*,9010) 'START TIME:            ',
     1                                dt(5), dt(6), dt(7), Chr2
C====================
C
      CALL BANNER(1)
cc      If(KLLERR) Return     BANNER has no error conditions            sld01
ccc***End of information duplicated in LHS_PREP
C
ccc***The following information is duplicated in LHS_RUN
C     -- If only constant and same as distribution types were found,
C     -- write the output file here and terminate gracefully without
C     -- sampling.
C
      If (NV == 0) Then
C        -- Write the sample to Unit 1
         Call SamOut(1)
cc         If(KLLERR) Return SamOut does not have any error conditions  sld01
         Write (6,*) ' '
         Write (6,*) '   Only Constant and Same As distribution ',
     1      'types were found in the LHS input.'
         Write (6,*) '   Thus, no sampling is being performed ',
     1      'by this LHS program run.'
         Write (6,*) '   The requested constants were written into ',
     1      'the LHS output file'
         Write (6,*) '   for processing by other codes.'
         Write (6,*)
ccc***The following information is duplicated in LHS_CLOSE
         Call FileOC(0)
ccc***End of information duplicated in LHS_CLOSE
c
cc         If(KLLERR) Return     FILEOC has no error conditions         sld01
cccc  Note: This is no longer treated as an error   May 22,2001         SLD
cccc         KLLERR = .TRUE.                                            SLD
cccc                                                                    SLD
         RETURN
      End If
ccc***End of information duplicated in LHS_RUN
c
c     -- Start to set up the sampling process.
c
ccc***The following information is duplicated in LHS_PREP
      CALL WRTPAR
cc      If(KLLERR) Return   WRTPAR has no error conditions              sld01
C
C     -- IF THE USER HAS SPECIFIED A CORRELATION STRUCTURE (ICM=1) THEN
C     -- THE CORRELATION MATRIX IS ECHOED AND CHECKED TO MAKE SURE THAT
C     -- IT IS POSITIVE DEFINITE. THE CHOLESKY FACTORIZATION IS COMPUTED
C     -- TO BE USED LATER AS PART OF THE PROCESS FOR INDUCING THE DESIRED
C     -- CORRELTION STRUCTURE
C
      IF(ICM.NE.0)THEN
         CALL PMTRX(NCM,3)
cc         If(KLLERR) Return  PMTRX has no error conditions             sld01
         IF(IRP.EQ.1)THEN
            IRP=0
            WRITE(6,9003)
         ENDIF
         IF(N.LE.NV)THEN
            ICM=0
            WRITE(6,9001)
            GO TO 10
         ENDIF
         CALL POSDEF(ITER)
         If(KLLERR) Return
         IF(ITER.GT.1)THEN
            WRITE(6,9002)
            CALL PMTRX(NCM,4)
cc            If(KLLERR) Return   PMTRX has no error conditions         sld01
         ENDIF
         NVX=(NV*(NV+1))/2
         DO 130 I=1,NVX
            CORR(I)=0.0
  130    CONTINUE
         DO 140 I=1,NV
            CORR(LOC1(I,I))=1.0
  140    CONTINUE
         LB=1
         REWIND 3
         READ(3)S
         DO 170 I=1,NCM
            DO 160 K=LB,NCM
               CORR(LOC1(LCM(K),LCM(I)))=S(LOC1(K,I))
  160       CONTINUE
            LB=LB+1
  170    CONTINUE
         CALL CHLSKY
cc         If(KLLERR) Return    CHLSKY has no error conditions          sld01
         REWIND 3
         WRITE(3)Q
C
C     -- THE FOLLOWING LINES ARE INCLUDED TO HANDLE THE
C     -- PATHOLOGICAL CASE N=3, NV=2 AND ICM=0
C
      ELSE IF(N.EQ.3.AND.NV.EQ.2)THEN
         CORR(1)=1.0
         CORR(2)=0.0001
         CORR(3)=1.0
         CALL CHLSKY
cc         If(KLLERR) Return     CHLSKY has no error conditions         sld01
         REWIND 3
         WRITE(3)Q
      ENDIF
   10 CONTINUE
ccc***End of information duplicated in LHS_PREP
C
ccc***The following information is duplicated in LHS_RUN
C     -- THIS LOOP IS EXECUTED ONCE FOR EACH REPETITION REQUESTED.
C     -- UNITS 7, 8, AND 9 HAVE BEEN DEFINED IN SUBROUTINE RDPAR
C
      DO 1000 IREP=1,NREP
         REWIND 7
         REWIND 8
         REWIND 9
c         IF (IREP .GT. 1) CALL BANNER(IREP)
         IF (IREP .GT. 1) THEN
            CALL BANNER(IREP)
cc            If(KLLERR) Return    BANNER has no error conditions       sld01
         End If
C
C        -- GENERATE THE DISTRIBUTION REQUESTED FOR VARIABLE J
         DO 100 J=1,NV
            IDT=IDIST(J)
            IF(IDT.EQ.1)THEN
               CALL BETA(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 2  .OR.  IDT .EQ. 3  .OR.
     1               (IDT .GE. 27  .AND.  IDT .LE. 35) ) THEN
               CALL NORMAL(J,IDT)
               If(KLLERR) Return
            ELSE IF(IDT.GE.4.AND.IDT.LE.7)THEN
               CALL UNIFRM(J,IDT)
               If(KLLERR) Return
            ELSE IF(IDT.EQ.8)THEN
               CALL TRIANG(J)
               If(KLLERR) Return
            ELSE IF(IDT .GE. 9  .AND.  IDT .LE. 11) THEN
               CALL CUMULC(J,IDT)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 12  .OR.  IDT .EQ. 13) THEN
               CALL CUMULD(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 14) THEN
               CALL POISON(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 15) THEN
               CALL GEOM(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 16) THEN
               CALL BINOM(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 17) THEN
               CALL NBINOM(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 18) THEN
               CALL HGEOM(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 19  .OR.  IDT .EQ. 25  .OR.
     1               IDT .EQ. 26) THEN
               CALL EXPON(J,IDT)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 20) THEN
               CALL WEIBUL(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 21) THEN
               CALL PARETO(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 22) THEN
               CALL GAMMA(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 23) THEN
               CALL IGAUS(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 24) THEN
               CALL ENTRPY(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 36) THEN
               CALL GUMBEL(J)
               If(KLLERR) Return
            ELSE IF (IDT .EQ. 37) THEN
               CALL FRECHET(J)
               If(KLLERR) Return
            ENDIF
  100    CONTINUE
C
C        -- IF A RANDOM SAMPLE HAS BEEN GENERATED IT MUST BE SORTED FROM
C        -- FROM SMALLEST TO LARGEST ON EACH VARIABLE AS PART OF THE
C        -- STRUCTURING TO REDUCE THE POSSIBILITY OF SPURIOUS CORRELATIONS
C====================
      Call date_and_time (DATRES,TIMRES,zone,dt)
      Chr2 = TIMRES(8:9)
C
      If (IPrint > 0 ) Write (*,9010) 'END OF SAMPLING TIME:  ',
     1                                dt(5), dt(6), dt(7), Chr2
C====================
         IF (N .NE. 1) THEN
            IF (IRS .EQ. 1) THEN
               DO 250 I=1,NV
                  DO 240 J=1,N
                     XV(J)=X(LOC(J,I))
  240             CONTINUE
                  CALL SIFT (XV,N)
cc                 If(KLLERR) Return      SIFT has no error conditions  sld01
                  DO 245 J=1,N
                     X(LOC(J,I))=XV(J)
  245             CONTINUE
  250          CONTINUE
            END IF
            NNV=N*NV
            DO 260 III=1, NNV
               XSAVE(III)=X(III)
  260       CONTINUE
C
C           -- SUBROUTINE MIX PAIRS THE SAMPLE OBSERVATIONS TO MATCH THE
C           -- DESIRED CORRELATION STRUCTURE
C====================
            Call date_and_time (DATRES,TIMRES,zone,dt)
            Chr2 = TIMRES(8:9)
C
            If (IPrint > 0 ) Write (*,9010) 'END OF SIFT TIME:      ',
     1                                dt(5), dt(6), dt(7), Chr2
            CALL MIX
            If(KLLERR) Return
C
            Call date_and_time (DATRES,TIMRES,zone,dt)
            Chr2 = TIMRES(8:9)
C
            If (IPrint > 0 ) Write (*,9010) 'END OF MIX TIME:       ',
     1                                dt(5), dt(6), dt(7), Chr2
C====================
         END IF
C
C        -- THE SAMPLE IS WRITTEN OUT TO UNIT 1
         Call SamOut(IRep)
         If(KLLERR) Return
c
c        -- Prepare the XSave matrix to reflect the post-MIX values
         NNV=N*NV
         DO 300 III=1, NNV
            XSAVE(III)=X(III)
  300    CONTINUE
C
C        -- SUBROUTINE DATOUT OUTPUTS THE LATEST SAMPLE AND CORRESPONDING
C        -- RANKS.  SUBROUTINE HSTOUT OUTPUTS HISTOGRAMS FOR THE CURRENT
C        -- SAMPLE.  SUBROUTINE COROUT OUTPUTS RAW AND RANK CORRELATIONS
C        -- FOR THE CURRENT SAMPLE.
c         IF (IDATA .NE. 0) CALL DATOUT
         IF (IDATA .NE. 0) THEN
            CALL DATOUT
            If(KLLERR) Return
         End If
c         IF (IHIST .NE. 0) CALL HSTOUT
         IF (IHIST .NE. 0) THEN
            CALL HSTOUT
            If(KLLERR) Return
         End If
c         IF (ICORR .NE. 0) CALL COROUT
         IF (ICORR .NE. 0) THEN
            CALL COROUT
            If(KLLERR) Return
         End IF
C====================
      Call date_and_time (DATRES,TIMRES,zone,dt)
      Chr2 = TIMRES(8:9)
c
      If (IPrint > 0 ) Write (*,9010) 'END TIME:              ',
     1                                dt(5), dt(6), dt(7), Chr2
C====================
C
 1000 CONTINUE
ccc***End of information duplicated in LHS_RUN
C
C     -- CLOSE ALL FILES AND STOP
C
ccc***The following information is duplicated in LHS_CLOSE
      CALL FILEOC(0)
ccc***End of information duplicated in LHS_CLOSE
c     Stop
C     normal termination from LHS
      RETURN
C
 9001 FORMAT('0',3(9('*'),' CAUTION USER PLEASE NOTE '),10('*'),//,1X,
     1       'SINCE THE SAMPLE SIZE IS LESS THAN OR EQUAL TO THE ',
     2       'NUMBER OF VARIABLES',/,' THIS IS NOT A FULL RANK CASE ',
     3       'SO THE REQUESTED',/,' CORRELATION STRUCTURE, SHOWN ',
     4       'ABOVE, CANNOT BE GENERATED',/,' THEREFORE THE INPUT ',
     5       'MATRIX WILL BE GENERATED AS IF THE',/,' INPUT ',
     6       'VARIABLES WERE INDEPENDENT',//,1X,115('*'))
 9002 FORMAT('0',3(9('*'),' CAUTION USER PLEASE NOTE '),10('*'),//,1X,
     1       'THE INPUT RANK CORRELATION MATRIX IS NOT POSITIVE ',
     2       'DEFINITE',/,' AN ITERATIVE PROCEDURE HAS BEEN USED TO ',
     3       'PRODUCE A SUBSTITUTE RANK CORRELATION MATRIX',/,' THIS ',
     3       'ADJUSTED RANK CORRELATION MATRIX APPEARS ON THE NEXT ',
     4       'PAGE',/,' THE USER SHOULD EXAMINE THIS MATRIX TO MAKE ',
     5       'SURE THAT THE CORRELATION REQUIREMENTS ARE STILL ',
     6        'SATISFIED',//,1X,115('*'))
 9003 FORMAT('0',3(9('*'),' CAUTION USER PLEASE NOTE '),10('*'),//,1X,
     1       'RANDOM PAIRING CANNOT BE USED WHEN THE USER SPECIFIES ',
     2       'A CORRELATION STRUCTURE',/,' THUS, THE RANDOM PAIRING ',
     3       'PARAMETER HAS BEEN IGNORED',//,1X,115('*'))
C9004 FORMAT('0',9('*'),' LAST VALUE OF ISEED  = ',I15)
 9010 FORMAT (5x,A,I2,":",I2,":",I2,".",A)
      END
