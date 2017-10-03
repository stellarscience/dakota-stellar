C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   10:34 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_OPTIONS
      SUBROUTINE LHS_OPTIONS(LHSREPS,LHSPVAL,LHSOPTS,IError)
c     This routine provides optional parameters when using the LHS
c     Input-by-Call Version of LHS
c     LHS_OPTIONS does not call any other external routines
c     Integer value input:
c     LHSREPS = LHS keyword for multiple number of samples
c     LHSPVAL = LHS keyword for selection of point value
c     String input:
c     LHSOPTS = allows user control of sampling and pairing methods
c      This string may contain "RANDOM PAIRING","RANDOM SAMPLE", or both

c
      USE KILLFILE
      USE CPARAM
c     CPARAM provides: NRep,IptVal,IRS,IRP
      USE InByCall
cc    InByCall provides subroutine flags: LINT,LPREP
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(LEN=*) :: LHSOPTS
      INTEGER :: LHSREPS,LHSPVAL
c
      CHARACTER(LEN=LENC) :: LCard
      LOGICAL :: MsgOpn
c
cc    IRP is random pairing switch, IRP=1 for random pairing
cc        flag may be set to one in:  READ,RDPAR,LHS_OPTIONS
cc        it is initialized to zero in: LHS,SETDEF,LHS_PREP,LHS_OPTIONS
cc    IRS is random sampling switch:
cc        IRS=0 for no random sampling,
cc        IRS=1 to do random sampling.
cc        IRS is initialized to zero in SETDEF and LHS_OPTIONS
cc        IRS may be set to one in READ, RDPAR or LHS_OPTIONS
c
c     Protect against multiple calls
c        to LHS_OPTIONS by reseting variables:
      IPtVal = 1
      NRep = 1
      IRS = 0
      IRP = 0
cc    Note: NRep,IRS,IRP are "normally" initialized in SetDef
cc          which is called from LHS_INIT, or RDPAR or READ
cc          IPtVal=1 is "normally" initialized in LHS_INIT or READ
c
c
c     Test that LHS_INIT has been called
      IF (LINIT /= 1) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9006)
         WRITE (99,9006)
         INQUIRE (Unit = 4, Opened = MsgOpn)
         IF (MsgOpn) THEN
            WRITE(4,9006)
         END IF
         Return
      END IF
c     Test that LHS_PREP has not been called prematurely
      IF (LPREP /= 0) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9001)
         WRITE (99,9001)
         INQUIRE (Unit = 4, Opened = MsgOpn)
         IF (MsgOpn) THEN
            WRITE(4,9001)
         END IF
         Return
      END IF
c
c     Test that LHSREPS > 0
      NRep = LHSREPS
c     Information in this section from READ.for
      If (NRep < 1) Then
         KLLERR = .TRUE.
         IError = 1
         Write (*,9007) NRep
         Write (99,9007) NRep
         INQUIRE (Unit = 4, Opened = MsgOpn)
         IF (MsgOpn) THEN
            Write (4,9007) NRep
         END IF
         RETURN
      End If
c     End of section from READ.fo

c     Test that LHSPVAL has a value of 0,1,or 2
c       If 0, point values must be supplied for ALL data items
c             in the (x)DIST Subroutines
      IptVal = LHSPVAL
      IF (IptVal < 0 .OR. IptVal > 2) THEN
         KLLERR = .TRUE.
         IError = 1
         Write (*,9008) IptVal
         Write (99,9008) IptVal
         INQUIRE (Unit = 4, Opened = MsgOpn)
         IF (MsgOpn) THEN
            Write (4,9008) IptVal
         END IF
         RETURN
      END IF
c     This section from READ.for
cc    Make sure all LHSOPTS are in upper case for testing
      LCard = LHSOPTS
c     -- Now convert Lcard to upper case
      Maxi = Len(LCard)
      Do i=1, Maxi
         If ( IChar(LCard(i:i)) > 96  .AND.  IChar(LCard(i:i)) < 123 )
     1         LCard(i:i) = Char( IChar(LCard(i:i)) - 32 )
      End Do
c     Determine which options are requested:
c       "RANDOM PAIRING","RANDOM SAMPLE", or both
      IF ( INDEX(LCard,'RANDOM SAMPLE') > 0 ) THEN
         IRS = 1
         Title(98:) = 'RANDOM SAMPLE'
      END IF
      IF ( INDEX(LCard,'RANDOM PAIRING') > 0 ) IRP = 1
c     End of Section from READ.for
      LOPTS = 1
      Return

 9001 FORMAT('1',5X,'LHS_PREP has been called prematurely ',/,5X,
     x'Call LHS_PREP just before call to LHS_RUN')
 9006 FORMAT(//,5x,'LHS_INIT or LHS_INIT_MEM must be called before ',
     x'any other LHS Input-By-Call Subroutines')
 9007 Format('1',5X,'The number of repetitions requested ',
     1            'is less than one:',I5)
 9008 FORMAT('1',5x,'The value of LHSPVAL is < 0 or > 2, LHSPVAL = ',I5
     1)
c
      END SUBROUTINE




