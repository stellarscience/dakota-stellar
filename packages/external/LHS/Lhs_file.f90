C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   11:24 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_FILES
      SUBROUTINE LHS_FILES(LHSOUT,LHSMSG,LHSTITL,FLOPTS,IError)
c     LHS_FILES defines output and message file; title for LHS
c     run; and file options.
c
c     LHS_FILES inputs:
c        LHSOUT - character string containing the LHS output file name
c                 (this file contains the sample information)
c        LHSMSG - character string containing the LHS message file name
c        LHSTITL - character string containing the Tile of the LHS run
c        FLOPTS - character string containing any requested options
c        FLOPTS available are:  CORR HIST DATA LHNONAM LHSSCOL LHSWCOL
c        as described in manual for keywords:  LHSRPTS, LHNONAM,
c                                              LHSWCOL and LHSSCOL
c     LHS_FILES returns:
c        IError - integer error flag
c           If an error is found, IError = 1 is returned
c
c
      USE KILLFILE
      USE PARMS
c     PARAMS provides:  ISamW,LENC
      USE CPARAM
c     CPARAMS provides:  TITLE,NamOut,I1Col,IDATA,IHIST,ICORR,
c                        Mfile,Sfile,TreeFl,Cmdlin
      USE InByCall
cc    InByCall provides subroutine flags: LINIT,LPREP,LFILES,
cc                                        IScrh1,IScrh6
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(LEN=*) :: LHSOUT
      CHARACTER(LEN=*) :: LHSMSG
      CHARACTER(LEN=*) :: LHSTITL
      CHARACTER(LEN=*) :: FLOPTS
c
      CHARACTER(LEN=LENC) :: LCard
      CHARACTER(LEN=LENC) :: KCard
      CHARACTER PDATA*4,PHIST*4,PCORR*4,PSCOL*7,PWCOL*7,PNONAM*8
      PARAMETER (PCORR='CORR',PHIST='HIST',PDATA='DATA',
     x           PNONAM='LHNONAM',PSCOL='LHSSCOL',PWCOL='LHSWCOL')


c     Open the LHS message file (Unit=6) (LHS keyword LHSMSG)
      IL = LEN_TRIM(LHSMSG)
      IF (IL == 0) THEN
c        Open scratch file, S6
            OPEN(4, FILE='S4', Form='FORMATTED')
c Modified by SFW 3/18/2002
c    1            Carriage Control='FORTRAN')
         IScrh6 = 1
      ELSE
c        Open User Message file
         LCard = TRIM(LHSMSG)
         CALL LJust(LCard)
         Open(4, File=LCard, Form='FORMATTED')
C    MODIFIED BY SFW 3/18/2002
c    1            Carriage Control='FORTRAN')
         IScrh6 = 2
         Mfile=LCard
      END IF
c
c     Now that message file is open check other flags
c     check for second call to LHS_Files without re-initializing
      IF (LFILES /= 0) THEN
         WRITE(4,9002)
         WRITE(99,9002)
         WRITE(*,9002)
         KLLERR = .TRUE.
         IError = 1
         Return
      END IF
cc      LFILES = 1   set before return to calling program
c     Test that LHS_INIT has been called
      IF (LINIT /= 1) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9006)
         WRITE (99,9006)
         WRITE(4,9006)
         Return
      END IF
c     Test that LHS_PREP has not been called prematurely
      IF (LPREP /= 0) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9001)
         WRITE (99,9001)
         WRITE(4,9001)
         Return
      END IF
c
c     Make sure all FLOPTS are in upper case for testing
      KCard = FLOPTS
c     -- Now convert Kcard to upper case
      Maxi = Len(KCard)
      Do i=1, Maxi
         If ( IChar(KCard(i:i)) > 96  .AND.  IChar(KCard(i:i)) < 123 )
     1         KCard(i:i) = Char( IChar(KCard(i:i)) - 32 )
      End Do

c     open unit 1 as scratch or user indicated with designated width
c     determine if wide column width indicator LHSWCOL was among keywords
c     LHS keyword LHSWCOL
      IF ( INDEX(KCard,PWCOL) /= 0) THEN
         ISamW = 2
      END IF
c     Open the LHS output file (Unit=1) (LHS keyword LHSOUT)
c     If string is blanks open default
      IL = LEN_TRIM(LHSOUT)
      IF (IL /= 0) THEN
         LCard = TRIM(LHSOUT)
         CALL LJust(LCard)
         IF (ISamW == 2) THEN
            Open(1, File=LCard, RECL=32000)
         Else
            Open(1, File=LCard)
         END IF
         IScrh1 = 2
         Sfile = LCard
      Else
c        open scratch file with default format, ISamW = 1
         OPEN(1, FILE='S1')
         IScrh1 = 1
      END IF

c  Following information found in READ.for (but not identical)
c     Store the Title for this LHS run (LHS keyword LHSTITL)
      Title(1:97) = LHSTITL

c     Determine options:
c     Use Index function to determine which options are requested
c     FLOPTS available are:  CORR HIST DATA LHSNONAM LHSSCOL
c     as described in manual for keywords:  LHSRPTS LHSNONAM and LHSSCOL
c     LHS keyword LHNONAM
      IF ( INDEX(KCard,PNONAM) /= 0) THEN
         NamOut = 0
      END IF
c     LHS keyword LHSSCOL
      IF ( INDEX(KCard,PSCOL) /= 0) THEN
         I1Col = 1
      END IF
c     Regular LHS locates 'HIST', 'CORR' and 'DATA' in call to routine OUTCRD
c     which is called from routine READ on keyword 'LHSRPTS'
c     LHS keyword LHSRPTS
      IF ( INDEX(KCard,PDATA) /= 0) THEN
         IDATA = 1
      END IF
      IF ( INDEX(KCard,PHIST) /= 0) THEN
         IHIST = 1
      END IF
      IF ( INDEX(KCard,PCORR) /= 0) THEN
         ICORR = 1
      END IF
c  End of information found in READ.for and OUTCRD.for
c     Set flag to indicate that LHS_FILES has been called
      LFILES = 1
c
      Return
c
 9001 FORMAT('1',5X,'LHS_PREP has been called prematurely ',/,5X,
     x'Call LHS_PREP just before call to LHS_RUN')
 9002 FORMAT('1',5X,'LHS_FILES has been called twice without ',
     x're-initializing'/,5X,'LHS_INIT must be called to re-start')
 9006 FORMAT(//,5x,'LHS_INIT or LHS_INIT_MEM must be called before ',
     x'any other LHS Input-By-Call Subroutines')
c
      END SUBROUTINE
