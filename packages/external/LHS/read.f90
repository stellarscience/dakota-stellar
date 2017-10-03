C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  29 May 101   12:19 pm
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::Read
      Subroutine Read
cc    READ is called from routine:  LHS                                 sld01
cc    READ calls routines:  SETDEF,GETCL,NEWCRD,RDPAR,CHKDIM,OUTCRD,    sld01
cc                          LJUST,RDPAR2                                sld01
c
c     This routine opens the input file as specified by the user and
c     determines whether it is of the "old" format (as identified by the
c     presence of the NOBS keyword) or the "new" format.  If the "old"
c     format is found, the file is passed off to RDPAR for processing of
c     all input data.  If the "new" format is found, the data processing
c     keywords are read in this routine.  It then prepares the file from
c     which RDPAR2 will read the distribution information, and calls
c     RDPAR2 to perform that input and verification function.
c
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARAMS provides:  LENC,ISamW,IPtVal                               sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  I1Col,ISeed, NREP,N,NamOut,IRSet                sld01
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER PNOBS*5,LOBS*12                                         sld01
      PARAMETER (PNOBS='NOBS ',LOBS='OBSERVATIONS')                     sld01
c
      Character Card*(LENC)
      Logical SamOpn, MsgOpn
c
c     --Set Default Parameter Values
c
ccc***The following information is duplicated in LHS_INIT
      Call SetDef
      If(KLLERR) Return
      NamOut = 1
C     --default width of output file depends on the value
C     --of ISamW, as set in the PARMS module.  If it is set
C     --narrow there (ISamW=0), then user cannot override it.
C     --Otherwise, set it for multi-column printing here unless
C     --the user overrides it with LHSSCOL in the input.
      IF (ISamW == 0) THEN
         I1Col = 1
      Else
         I1Col = 0
      END IF
      IPtVal = 1
ccc***End of information duplicated in LHS_INIT
      TreeFl = ' '
      MFile = 'The screen'
c
c   NOTE:  Keyword file was previously obtained in the READ subroutine.
c   The location was moved to accomodate LHS conversion to DLL for use
c   in the Visual Basic program; putting control of keyword inout in driver
c   routine allows all the other routines to remain the same for both the
c   Visual Basic and stanad-alone versions.
c
c   Following section commented out 11-96:
c     -- Get the keyword file name from the command line.
c     -- If none specified, offer the user the opportunity to input
c     -- one from the keyboard or accept the default KEYWORD.DAT
c
c      Print *, 'Welcome to LHS - The Latin Hypercube Sampling Program'
c      Print *, ' '
c
c      Call GetCL(CmdLin)
c      If(KLLERR) Return
c      If (CmdLin /= ' ') Then
c         Open (5, File=CmdLin, Status="OLD", Err=10)
c         Go To 100
c 10      Print *, 'Error opening file specified on command line.'
c         Print *
c      End If
c 15   Print *, 'Enter the name of the LHS input file to be read,'
c      Print *, 'enter / to exit LHS, or enter . to accept the '
c      Print *, 'default input file name (KEYWORD.DAT):  '
c      Print *
c      Read (*,9001) CmdLin
c      If (CmdLin == '/') Then
c         Print *, 'Program Terminated'
c        Stop
c         KLLERR = .TRUE.
c         RETURN
c      Else If (CmdLin == '.') Then
c         Open (5, File='KEYWORD.DAT', Status="OLD", Err=20)
c         CmdLin = 'KEYWORD.DAT'
c         Go To 100
c      Else
c         Open (5, File=CmdLin, Status="OLD", Err=20)
c         Go To 100
c      End If
c
c 20   Print *, 'Error opening the file ', CmdLin
c      Print *
c      Go To 15
c  End of 11-96 commenting out.
c
c     -- Now ready to start reading the keyword file.
c
c     -- Read the keyword file.  If the keyword file contains the
c     -- keyword NOBS, then rewind the file and read it using RDPAR
c     -- since it is assumed to correspond to the original LHS input
c     -- format.  Otherwise, rewind the file and continue to read in
c     -- this routine.
c
 100  Call NewCrd(Card,5,IEnd)
cc      If(KLLERR) Return   NewCrd has no error conditions              sld01
      If (IEnd == 1) Then
         Write (99,*) 'The LHS input file is empty!'
         KLLERR = .TRUE.
         RETURN
      End If
c
      IFound = 0
      Do While (IEnd == 0)
         If ( Card(1:5) == 'NOBS ' ) Then
            IFound = 1
            Exit
         End If
         Call NewCrd(Card,5,IEnd)
cc         If(KLLERR) Return    NewCrd has no error conditions          sld01
      End Do
c
      Rewind 5
      If (IFound == 1) Then
c        -- Old format found - prepare for RDPAR.  Note that the
c        -- old output format is default for the old input format,
c        -- so set NamOut to zero as a default before calling RDPAR.
         NamOut = 0
         Call RDPAR
         If(KLLERR) Return
         Return
      End If
c
c     -- Now we know that we are reading the new input format.
c     -- Look for the new keywords and set the values associated with them.
c
      Call NewCrd(Card,5,IEnd)
      If(KLLERR) Return
c
      Do While (IEnd == 0)
c
         If ( Card(1:7) == 'LHSOBS ' ) Then
c           -- set the number of observations to be generated
            Read(Card(8:),*) N
ccc***The following information is duplicated in LHS_INIT
            If (N < 1) Then
               Write (*,9006) N
               Write (99,9006) N
               Write (6,9006) N
               KLLERR = .TRUE.
               RETURN
            End If
            Call ChkDim(1,N,NMax,PNObs,LObs)
            If(KLLERR) Return
c
         Else If ( Card(1:8) == 'LHSSEED ' ) Then
c           -- set the value for the random number generator seed
            Read(Card(9:),*) ISeed
ccc   Variable ISeedSv added for second random number generator seed    SLD
            ISEEDSV = ISeed                                             SLD
ccc                                                                     SLD
            If (ISeed <= 0 ) Then
               Write (*,9002) ISeed
               Write (99,9002) ISeed
               Write (6,9002) ISeed
               KLLERR = .TRUE.
               RETURN
            End If
c           -- Set the random number generator here if necessary
            IRSet=1
ccc***End of information duplicated in LHS_INIT
c
         Else If ( Card(1:8) == 'LHSREPS ' ) Then
c           -- set the number of complete samples (repetitions) to be done
            Read(Card(9:),*)NRep
ccc***The following information is duplicated in LHS_OPTIONS
            If (NRep < 1) Then
               Write (*,9007) NRep
               Write (99,9007) NRep
               Write (6,9007) NRep
               KLLERR = .TRUE.
               RETURN
            End If
ccc***End of information duplicated in LHS_OPTIONS
c
ccc***The following information is duplicated in LHS_FILES
         Else If ( Card(1:8) == 'LHSTITL ' ) Then
c           -- set and left-justify the title of the LHS run
            Title(1:97) = Card(9:)
c
         Else If ( Card(1:8) == 'LHSRPTS ' ) Then
c           -- set the LHS reporting (output) options
            Call OutCrd(Card)
            If(KLLERR) Return
ccc***End of information duplicated in LHS_FILES
c
ccc***The following information is duplicated in LHS_OPTIONS
         Else If ( Card(1:8) == 'LHSPVAL ' ) Then
c           -- set the type of calculation to be done when setting the point
c           -- estimate to be put in the new format output file
            Read(Card(9:),*)IPTVAL
ccc***End of information duplicated in LHS_OPTIONS
c
ccc***The following information is duplicated in LHS_FILES
         Else If ( Card(1:8) == 'LHNONAM ' ) Then
c           -- LHS is to print its results in the old format (no names)
            NamOut = 0
c
         Else If ( Card(1:8) == 'LHSSCOL ' ) Then
c           -- LHS is to print its results in a single column
            I1Col = 1
ccc***End of information duplicated in LHS_FILES
c
ccc***The following information is duplicated in LHS_FILES
         Else If ( Card(1:8) == 'LHSWCOL ' ) Then
c           -- LHS is to print its results in a wide format
            ISamW = 2
c           -- Sample output file may already be open with
c           -- the wrong format.  If it's open, close it and
c           -- reopen it with the right format to ensure that
c           -- the user can specify the keywords in any order.
            Inquire (Unit = 1, Opened = SamOpn)
            If (SamOpn) Then
               CLOSE (1)
               Open(1, File=SFile, RECL=32000)
            End If
ccc***End of information duplicated in LHS_FILES
c
ccc***The following information is duplicated in LHS_OPTIONS
         Else If ( Card(1:8) == 'LHSOPTS ' ) Then
c           -- set the LHS processing options
            If ( Index(Card,'RANDOM SAMPLE')  > 0 ) Then
               IRS = 1
               Title(98:) = 'RANDOM SAMPLE'
            End If
            If ( Index(Card,'RANDOM PAIRING') > 0 ) IRP = 1
ccc***End of information duplicated in LHS_OPTIONS
c
ccc***The following information is duplicated in LHS_FILES
         Else If ( Card(1:7) == 'LHSOUT '  ) Then
c           -- set the name of the file to which the sample will be written
            SFile = Card(8:)
            Call LJust(SFile)
            If(KLLERR) Return
            Inquire (Unit = 1, Opened = SamOpn)
            If (SamOpn) Then
               Write (99,*) 'The keyword LHSOUT appears twice in the ',
     1            'keyword input file!'
               KLLERR = .TRUE.
               RETURN
            End If
            IF (ISamW == 2) THEN
               Open(1, File=SFile, RECL=32000)
            Else
               Open(1, File=SFile)
            END IF
c
         Else If ( Card(1:7) == 'LHSMSG '  ) Then
c           -- set the name of the message (program output) file
            MFile = Card(8:)
            Call LJust(MFile)
            If(KLLERR) Return
            Open(6, File=MFile, Form='FORMATTED')
c            Open(6, File=MFile, Form='FORMATTED',
c     1            Carriage Control='FORTRAN')
ccc***End of information duplicated in LHS_FILES
c
         Else If ( Card(1:8) == 'PRETRIN ' ) Then
c           -- set the name of the ETPRE/EVNTRE tree input file from
c           -- which the LHS input is to be read
            TreeFl = Card(9:)
            Call LJust(TreeFl)
            If(KLLERR) Return
c
         Else If ( Card(1:7) == 'ENDKEY '  .OR.
     1             Card(1:8) == 'DATASET:'      ) Then
c           -- stop processing keyword file
            Exit
c
         Else
c           -- all other keywords are ignored
c
         End If
c
c        -- read a new record and repeat the loop
         Call NewCrd(Card,5,IEnd)
         If(KLLERR) Return
      End Do
c
ccc***The following information is duplicated in LHS_FILES
      Inquire (Unit = 1, Opened = SamOpn)
      Inquire (Unit = 6, Opened = MsgOpn)
c
      If (.NOT. SamOpn  .OR.  .NOT. MsgOpn) Then
         Write (99,*) 'The required keywords LHSOUT and LHSMSG were not'
         Write (99,*) 'both present in the keyword input file.'
         KLLERR = .TRUE.
         RETURN
      End If
ccc***End of information duplicated in LHS_FILES
c
c     -- Check to be sure that all mandatory keywords were specified
c     -- is performed with calls to CHKZRO and CHKDIM at the end of
c     -- subroutine RDPAR2.
c
c     -- We are now done reading the keyword file, so close it and
c     -- open the "tree input file" as Unit 5 and read the distributions
c     -- and correlations.
c     -- Modified 5/14/97 to eliminate the close and re-open as this seems to
c     -- give problems on some operating systems.  Changed it to a Rewind only.
c
      If (TreeFl == ' ') Then
         REWIND 5
      Else
         CLOSE (5)
c         Open (5, TreeFl, Status="OLD", Err=999)
         Open (5, File=TreeFl, Status="OLD", Err=999)
      END If
C
      Call RDPAR2
      If(KLLERR) Return
c
      Return
c
c     -- Failure while opening the tree file
c
 999  Write (99,*)
      Write (99,*)  'Error opening the tree input file:', TreeFl
      Write (99,*)
      KLLERR = .TRUE.
      RETURN
c
c Format changed as per Nelson Deane request on May 22, 2001            SLD
c 9001 Format (A50)                                                     SLD
 9001 Format (A)                                                        SLD
c End of May 22, 2001 requested change                                  SLD
 9002 Format('1',5X,'The random number generator seed value must ',
     1       'be positive.',/,5X,'The following value was found: ',I12)
 9006 Format('1',5X,'The number of observations requested ',
     1            'is less than one:',I5)
 9007 Format('1',5X,'The number of repetitions requested ',
     1            'is less than one:',I5)
c
      End
