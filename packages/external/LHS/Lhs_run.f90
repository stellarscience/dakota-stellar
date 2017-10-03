C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   11:15 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_RUN
      SUBROUTINE LHS_RUN(MAXVAR,MAXOBS,MAXNAM,IError,
     x LSTNAME,INDXNAM,PTVALST,NUMNAM,SMATX,NUMVAR,RMATX,RFLAG)
c     LHS_RUN runs the LHS sample engine
c     LHS_RUN inputs:
c        MAXVAR - Maximum number of variables, integer
c        MAXOBS - Maximum number of observation, integer
c        MAXNAM - Maximum number of variable names including "same as"
c        RFLAG - flag to indicate whether the rank matrix is being passed in 
c     LHS_RUN returns:
c        LSTDNAM - a list of distribution names (type character)
c        INDXNAM - integer array containing index number(position)
c                  of names in sample data
c        PTVALST - an array of the associated point values (type real)
c        SMATX - the sample data matrix dimension (MAXVAR,MAXOBS)
c        RFLAG - flag to indicate if no rank data needs to be passed (RFLAG=0),
c                if rank data will be set (RFLAG = 1),
c                if rank data will be obtained (RFLAG = 2), 
c                or if rank data will both be "set" and "get" (RFLAG=3) 
c        RMATX - the rank data matrix dimension (MAXVAR,MAXOBS)
c        NUMNAM - the integer number of variable names including "same as"
c        NUMVAR- the number of variables for which observations
c                are generated (with no duplication for same as)
c        IError - integer error flag
c           If an error is found, IError = 1 is returned
c
      USE KILLFILE
      USE PARMS
c     PARAMS provides:   NMAX
      USE CPARAM
c     CPARAMS provides:  ISeed,N,NV,NREP,IVarNm,PValue,IDIST,List
      USE InByCall
cc    InByCall provides: LINT,LPREP,LFILES,LDIST,NNames
      USE CSAMP
c     CSAMP provides: X and XSAVE arrays
      USE CRANK
c     CRANK provides: XV array

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION :: SMATX(MAXVAR,MAXOBS)
      DOUBLE PRECISION :: RMATX(MAXVAR,MAXOBS)
      DOUBLE PRECISION :: PTVALST(MAXNAM)
      INTEGER :: INDXNAM(MAXNAM)
      CHARACTER(LEN=16) :: LSTNAME(MAXNAM) 
      INTEGER :: RFLAG

C====================
      CHARACTER*11 TIMRES
C====================
c     next lines added 1-17-96 to convert date/time routine to LF90
      integer :: dt(8)
      character (len=5) :: zone
      character * 8 DATRES
      character * 2 Chr2
C
C     -- STATEMENT FUNCTION
      LOC(I,J)=(J-1)*N+I
      NUMNAM = NNames
      NUMVAR = NV
c
c     Check that LHS_INIT or LHS_INIT_MEM, at least one LHS_xDIST,
c        and LHS_PREP have been called prior calling LHS_RUN
      I3 = LINIT + LDIST + LPREP
      IF (I3 /= 3) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9022)
         WRITE(4,9022)
         WRITE(99,9022)
         Return
      END IF

c
c     Following section from Standard LHS.for routine
C
C     -- If only constant and same as distribution types were found,
C     -- write the output file here and terminate gracefully without
C     -- sampling.
C
      If (NV == 0) Then
C        -- Write the sample to Unit 1
         Call SamOut(1)
cc         If(KLLERR) Return SamOut does not have any error conditions  sld01
         Write (4,*) ' '
         Write (4,*) '   Only Constant and Same As distribution ',
     1      'types were found in the LHS input.'
         Write (4,*) '   Thus, no sampling is being performed ',
     1      'by this LHS program run.'
         Write (4,*) '   The requested constants were written into ',
     1      'the LHS output file'
         Write (4,*) '   for processing by other codes.'
         Write (4,*)
CCCC         Call FileOC(0)       sld: this is done in LHS_CLOSE
cc         If(KLLERR) Return     FILEOC has no error conditions         sld01
ccc         KLLERR = .TRUE.
ccc
ccc      IError = 1
csld asked Greg if this really should be treated as an error
ccc---sld, Answer was no.
         RETURN
      End If
c
c     End of section from standard LHS routine
c
c
c
c     Following section from standard LHS routine
C
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

cc          This is used in most cases when the ranks are not pre-specified    
            IF ((RFLAG .EQ. 0).OR.(RFLAG .EQ. 2)) THEN
            CALL MIX
            END IF
cc          This is used in the case where one is sampling incrementally 
cc          and has predefined ranks            
            IF ((RFLAG .EQ. 1).OR.(RFLAG .EQ. 3)) THEN
cc            WRITE (*,*) 'We are in the loop'
 
              DO 280 J=1,NV
                DO 270 I=1,N
                  RXV(I)=RMATX(J,I)
                  X(LOC(I,J))=RXV(I)
cc                WRITE (*,*) 'RXV(I) ', RXV(I), 'I ', I, 'J ', J
  270           CONTINUE
  280         CONTINUE
              DO 290 J=1,NV
                DO 285 I=1,N
                  K=INT(X(LOC(I,J))+0.01)
                  X(LOC(I,J))=XSAVE(LOC(K,J))
  285           CONTINUE
  290          CONTINUE 
            END IF

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
ccc
ccc    Following section modified from LHS.for
ccc    (i.e. SAMSTOR call added and SamOut only called
ccc          for user defined output file)
ccc    Call SAMSTOR to store matrix, point values, variable names and
ccc    positions in matrix for return to user
       CALL SamStor(MAXVAR,MAXOBS,MAXNAM,LSTNAME,INDXNAM,
     x              PTVALST,SMATX)
ccc
ccc only call SamOUt if unit 1 is user defined output file
      IF (IScrh1 == 2) THEN
C        -- THE SAMPLE IS WRITTEN OUT TO UNIT 1
         Call SamOut(IRep)
         If(KLLERR) Return
      END IF
ccc   End of modifications
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
      
         IF ((RFLAG .EQ. 2).OR.(RFLAG.EQ.3)) THEN
            DO 620 J=1,NV
               DO 610 I=1,N
 610              XV(I)=X(LOC(I,J))
                  CALL RANKER
               DO 620 I=1,N
                  RMATX(J,I)=RXV(I)
 620        CONTINUE
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
cccc  End of standard LHS section
c
c     re-initialize subroutine called tracking switches
c     except for LRUN in case LHS_COROUT is called
      LINIT = 0
      LPREP = 0
      LRUN = 1
      LFILES = 0
      LDIST = 0
      IScrh1 = 0
      IScrh6 = 0
c
      Return
c
 9010 FORMAT (5x,A,I2,":",I2,":",I2,".",A)

 9022 FORMAT('1',5X,' LHS_INIT or LHS_INIT_MEM,'/,5X,
     x' at least one distribution call ',/,8X,
     x'(i.e. LHS_DIST, LHS_UDIST or LHS_SDIST) and',/,5X,
     x'LHS_PREP must be called prior to calling LHS_RUN')
      END SUBROUTINE
