C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   10:36 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_INIT
      SUBROUTINE LHS_INIT(LHSOBS,LHSSEED,IError)
c     This routine provides the initial setup for using the LHS
c     Input-by-Call Version of LHS
c     LHS_INIT calls routines:  SetDef,FILEOC,ChkDIM
c     LHS_INIT initializes all module calls
c
c     LHS_INIT inputs:
c        LHSOBS - integer number of observations requested
c        LHSSEED - integer random number seed
c     LHS_INIT returns:
c        IError - integer error flag
c           If an error is found, IError = 1 is returned
c
      USE KILLFILE
      USE PARMS                         
cc    PARMS provides:  ISamW,IPtVal
      USE CPARAM                        
cc    CPARAM provides: N,ISeed,I1Col,NamOut
      USE InByCall
cc    InByCall provides subroutine flags: LINT,LPREP,LRUN,LFILES,LDIST
C
      USE DISTNM
      USE CSAMP
      USE CWORKC
      USE CWORKX
      USE CRANK
      USE CCMATR
      USE STAR
      USE UICORR
      USE CHRCRD
      USE OBSTR
      USE PDMAT
      USE FIRSTS
      USE LOCALVARS
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: LHSOBS,LHSSEED
c
      CHARACTER PNOBS*5,LOBS*12
      PARAMETER (PNOBS='NOBS ',LOBS='OBSERVATIONS')
      LOGICAL :: MsgOpn
c
c
C       Here are the calls that need to be made to initialize the modules.
C       Note that PARMS must be initialized first since it contains the
C       constants that are to be used to initialize the other modules.
C       Note that it must also be compiled first!
C
      CALL PRAMS_INIT()
C
      CALL DISTNM_INIT()
      CALL CPARAM_INIT()
      CALL InByCall_INIT()
      CALL CSAMP_INIT()
      CALL CWORKC_INIT()
      CALL CWORKX_INIT()
      CALL CRANK_INIT()
      CALL CCMATR_INIT()
      CALL STAR_INIT()
      CALL UICORR_INIT()
      CALL CHRCRD_INIT()
      CALL OBSTR_INIT()
      CALL PDMAT_INIT()
      CALL FIRSTS_INIT()
      CALL LOCALVARS_INIT()
C
c
c     Initialize LHS error flag as False
      KLLERR = .False.
c     re-initialize LRUN to zero as it is not reset in LHS_RUN
c     (because of testing on it in LHS_COROUT to return matrices)
      LRUN = 0
c     re-initialize LPOSDEF to zero as it is not reset in LHS_RUN
c     (because of setting IPOSDEF in LHS_COROUT to return value)
      LPOSDEF = 0
c     Initialize correlation variable NCV found in UICORR module
c     NCV is initialized in RDPAR2 in standard LHS
      NCV = 0
c     Initialize variable ICM found in CPARAMS module
c     ICM is initialized in RDPAR2 in standard LHS
      ICM = 0
c     Initialize variable NNames found in InByCall module
c     NNames was local to RDPAR2; it is used by most LHS_xxx routines
      NNames = 0
c
c
c
cc     Note:  Error file now opened in LHS_INIT or LHS_INIT_MEM
cc            rather than LHS driver
c
c      moved setup for LHS error file 99 to here -- Error file is
c      now opened/closed in driver routine rather than subroutine FILEOC
c
c
c        -- Processing to open the LHS Error File
c
c        -- We want to be sure that the error file is created anew for this
c        -- run in order to assure that it has the proper date and time stamp
c        -- information.  Thus, we open it with Status=Unknown and immediately
c        -- close it with Status=Delete.  Then we re-open the file and write
c        -- an error message to it.  We then close it with Status=Keep in
c        -- order to assure that the program buffers are flushed and the file
c        -- will really exist on the disk.  Finally, we re-open the file and
c        -- leave it open until the program executes a normal termination.
c        -- If the program crashes, the error file will remain in existance.
c        -- However, on normal termination, the error file will be deleted.
c        -- Thus, an external program can check for normal termination by
c        -- checking for the existance of this file.
         Open (99, File='LHS.ERR', Status='UNKNOWN', Form='FORMATTED')
         Write (99,*) 'One line into the file just to be sure...'
         Close  (99, Status='DELETE')
         Open (99, File='LHS.ERR', Status='NEW', Form='FORMATTED')
         Write (99,*) 'An error occurred during LHS processing.'
         Write (99,*) 'Consult the message file for additional ',
     1      'information.'
         Close (99, Status='KEEP')
         Open (99, File='LHS.ERR', Status='OLD', Form='FORMATTED')
c
c     Check that LHS_INIT was not also called for this LHS run sequence:
      IF (LINIT == 1) THEN
         Write (*,9019)
         Write (99,9019)
         IError = 1
         KLLERR = .TRUE.
         RETURN
      END IF
c
      IF (LPREP /= 0) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9016)
         WRITE (99,9016)
         Return
      END IF

c From LHS subroutine Read
c     --Set Default Parameter Values (called from READ in standard LHS)
      Call SetDef
      NamOut = 1
C     --default width of output file depends on the value
C     --of ISamW, as set in the PARMS module.  If it is set
C     --narrow there (ISamW=0), then user cannot override it.
C     --Otherwise, set it for multi-column printing here unless
C     --the user overrides it with LHSSCOL in the input.
c     Default ISamW=1 value set in PARMS module
      IF (ISamW == 0) THEN
         I1Col = 1
      Else
         I1Col = 0
      END IF
c     IPtVal can be changed by keyword LHSPVAL
c     if IPtVal = 0, ALL point values MUST be provided by user
      IPtVal = 1
c  End of section from READ.for

c From LHS subroutine LHS.for
c     Open scratch files: 2,3,7,8,9
      CALL FILEOC(1)
c
c
      N = LHSOBS
c     following test from READ routine
      If (N < 1) Then
         Write (*,9006) N
         Write (99,9006) N
         KLLERR = .TRUE.
c        Added IError flag setting not in original READ code
         IError = 1
         RETURN
      End If
      Call ChkDim(1,N,NMax,PNObs,LObs)
cccc      If(KLLERR) Return
      IF (KLLERR) THEN
         IError = 1
         Return
      END IF
cccc  End of section for READ routine

      ISeed = LHSSEED
      ISEEDSV = ISeed
c     following test from READ routine
      If (ISeed <= 0 ) Then
         Write (*,9002) ISeed
         Write (99,9002) ISeed
         KLLERR = .TRUE.
c        Added IError flag setting not in original READ code
         IError = 1
         RETURN
      End If
c           -- Set the random number generator here if necessary
      IRSet=1
cccc  End of section for READ routine


cc    Initialize file labels (Set in Read in Standard LHS)
cc       used in Banner
      CmdLin = ' '
      Mfile = ' '
      Sfile = ' '
      TreeFl = ' '
cc
c     Set indicator to determine LHS_INIT routine was called:
      LINIT = 1
      Return
c
 9016 FORMAT('1',5X,'LHS_INIT or LHS_INIT_MEM must be called ',
     x'before all other LHS Input-By-Call Subroutines')
 9002 Format('1',5X,'The random number generator seed value must ',
     1       'be positive.',/,5X,'The following value was found: ',I12)
 9006 Format('1',5X,'The number of observations requested ',
     1            'is less than one:',I5)
 9019 FORMAT(//,5X,'User must call LHS_INIT_MEM or LHS_INIT, but not',
     x' both for the same LHS run sequence')
      END SUBROUTINE
