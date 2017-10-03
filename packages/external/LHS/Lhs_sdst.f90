C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   11:06 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_SDIST
      SUBROUTINE LHS_SDIST(NAMVAR,IPTFLAG,PTVAL,DISTYPE,
     x               NUMINTV,NOBSpI,ENDPNTS,IError,IDISTNO,IPVNO)
c     LHS_SDIST inputs data for user-defined n-part histogram
c     LHS_SDIST calls routines: LJUST, CHKSTR, WRTCRD
c
c   Descriptions of call list parameters:
c   Inputs:
c     NAMVAR = name of the variable
c     IPTFLAG = 0, if PTVAL is NOT to be used
c     IPTFLAG = 1, if PTVAL is to be used
c     PTVAL = point value associated with NAMVAR
c     DISTYPE = character string naming distribution type
c     NUMINTV = integer number of intervals "n" in a n-part histogram
c     NOBSpI = integer series of "n" integer values "obs", each "obs(i)"
c              represents the number of observations LHS will draw
c              from the i-th interval
c     ENDPNTS = array (type real) consisting of:
c                 first_point (the distribution minimum)
c                 followed by "n" interval endpoints
c         Note:  The sum of all interval observations must equal the number of
c                observations requested for the LHS run

c   Outputs:
c     IError = error flag, returned = 1 to indicate some error occurred
c     IDISTNO = current distribution number
c     IPVNO = point value index number
c
      USE KILLFILE
      USE PARMS
c     PARMS provides: NamLen,MAXTB,LENC
      USE CPARAM
c     CPARAM provides: IptVal, List, IVarNm, IDIST arrays
      USE InByCall
cc    InByCall provides subroutine flags: LINT,LPREP,LDIST,NNames,IScrh6
      USE DISTNM
c     DISTNM provides:  DIST,IDSST,IDSEND,IDSPAR,LEND,MAXPAR
C     COMMON/STAR/NSUBOB(NINTMX),SUBINT(NINTMX+1),NINT
      USE STAR
c     STAR provides: NINT and arrays NSUBOB, SUBINT
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     Call List Variables:
      INTEGER :: IPTFLAG,NUMINTV,IDISTNO,IPVNO,IError
      INTEGER :: NOBSpI(1)
      DOUBLE PRECISION :: PTVAL,ENDPNTS(1)
      CHARACTER(LEN=*) :: NAMVAR
      CHARACTER(LEN=*) :: DISTYPE
c
c     LHS Internal Variables:
      CHARACTER(LEN=NamLen) :: Name
      CHARACTER(LEN=35) :: NamVal
cc    CHARACTER(LEN=32768) :: LCard, > dimension was needed in RDPAR2
      CHARACTER(LEN=40) :: LCard
      Character Card*(LENC)
c
c
c  check to see if message file is open,
c    if not then open as a scratch file
      IF (IScrh6 == 0) THEN
c        Open scratch file, S6
         OPEN(4, FILE='S4', Form='FORMATTED')
c    MOdified by SFW 3/18/2002
c    1            Carriage Control='FORTRAN')
         IScrh6 = 1
      END IF

c Check to see that INIT has been called and PREP has not
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
c  check to be sure NAMVAR is <= 16 characters and not all blanks
c  terminate trailing blanks; eliminate any leading blanks
      LCard = NAMVAR
      LCard = TRIM(LCard)
      Call LJust(LCard)
      IL = LEN_TRIM(LCard)
      IF (IL == 0) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9040)
         WRITE(99,9040)
         WRITE(4,9040)
         Return
      END IF
      IF (IL > NamLen) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9010) NAMVAR
         WRITE(99,9010) NAMVAR
         WRITE(4,9010) NAMVAR
         Return
      END IF
      NAME = LCard      
c
c     IPtVal can be changed by keyword LHSPVAL -- (in LHS_OPTIONS)
c        default is 1 (i.e. mean value) set in LHS_INIT
c     If IPtVal = 0, ALL point values MUST be provided by user
      IVFnd = 0
      IF (IPtVal == 0) THEN
         If (IPTFLAG == 1) then
            Value = PTVAL
            IVFnd = 1
         Else
c        Error condition as point value must be provided
            KLLERR = .TRUE.
            IError = 1
            WRITE(*,9013) NAMVAR
            WRITE (99,9013) NAMVAR
            WRITE(4,9013) NAMVAR
            Return
         End If
      ELSE
        IF (IPTFLAG == 1) THEN
c        Write warning that point value provided will not be used
            WRITE(*,9014) PTVAL,NAMVAR
            WRITE (99,9014) PTVAL,NAMVAR
            WRITE(4,9014) PTVAL,NAMVAR
            IVFnd = 0
         END IF
      END IF
c
cccc Following section from RDPAR2
cccc these variables appear in section of code:
c       Err = flag used for tracking when editing multiple errors
c       Err is not needed for LHS_DIST routine
c       IVFnd = 0 or 1, Flag for point value present=1, or not=0
c       NNames is initialized in LHS_INIT and used by other
c          LHS_Distribution subroutines through module InByCall
c
c        -- check the name against the existing list of names
c
         IFound = 0
         Do i=1, NNames
            If ( Name == List(i) ) Then
               IFound = 1
               IList = i
               If (IVarNm(i) /= 0) Then
cccc  Following section modified:
c                 -- duplicate definition found
cccc                  Print *, 'Multiple definitions found for ', Name
c****** add prints to message and error files
cccc                  Write(99,*) 'Multiple definitions found for ', Name
cccc                  Write(4,*) 'Multiple definitions found for ', Name
c******
cccc                  Err = .True.
                  IError = 1
                  KLLERR = .TRUE.
                  WRITE(*,9015) Name
                  WRITE(4,9015) Name
                  WRITE(99,9015) Name
                  Return
cccc  end of modification
               End If
               Exit
            End If
         End Do
c        -- If name not found in list, add to the list
         If (IFound == 0) Then
            NNames = NNames + 1
            IList = NNames
            List(IList) = Name
         End If
c
c        -- If no point value found and LHSPVAL 0 is specified, then
c        -- declare an error and stop.
c
c
cccc following section omitted as it was checked on input to this routine
cccc         If ( IPtVal == 0  .AND.  IVFnd == 0 ) Then
cccc            Print *, 'Error: When the keyword LHSPVAL 0 is present, ',
cccc     1               'all point values must be specified.'
cccc            Print *, 'Point value missing for ', Name
c****** add prints to message and error files
cccc            Write(99,*) 'Error: When the keyword LHSPVAL 0 is ',
cccc     1               'present, all point values must be specified.'
cccc            Write(99,*) 'Point value missing for ', Name
cccc            Write(4,*) 'Error: When the keyword LHSPVAL 0 is ',
cccc     1               'present, all point values must be specified.'
cccc            Write(4,*) 'Point value missing for ', Name
c******
cccc           Err = .True.
cccc         End If
         If ( IVFnd /= 0 ) PValue(IList) = Value
c
c
      NamVal = Name
c
c     Determine distribution name
c       terminate trailing blanks; eliminate any leading blanks
      LCard = TRIM(DISTYPE)
      Call LJust(LCard)

c     -- Now convert Lcard to upper case
      Maxi = Len(LCard)
      Do i=1, Maxi
         If ( IChar(LCard(i:i)) > 96  .AND.  IChar(LCard(i:i)) < 123 )
     1         LCard(i:i) = Char( IChar(LCard(i:i)) - 32 )
ccc         If ( IChar(Card(i:i)) == 9 )  Card(i:i) = ' '
ccc         If ( Card(i:i) == ',' )       Card(i:i) = ' '
      End Do
c
ccc   Place data normally read from input into appropriate arrays
      NINT = NUMINTV
      NINTP1 = NINT + 1
      DO K = 1,NINT
        NSUBOB(K) = NOBSpI(K)
        SUBINT(K) = ENDPNTS(K)
      END DO
      SUBINT(NINTP1) = ENDPNTS(NINTP1)
ccc
c
c   This is another section from RDPAR2
c
c        -- Read and check the * distributions
c
        IFound = 0
         Do ID=1, LEND
            If ( IDSPAR(ID) /= -1 ) Cycle
            IDL=IDSEND(ID)-IDSST(ID)+1
            IDL1=IDL+1
            IF ( LCard(1:IDL) .EQ. DIST(IDSST(ID):IDSEND(ID))) THEN
cccc following section modified for LHS_DIST
cc               Read (LCard(IDL1:),*,ERR=9000) NINT                    sld01
cc             Read (LCard(IDL1:),*,ERR=9000) NINT,(NSUBOB(I),I=1,NINT),sld01
cc     1                            (SUBINT(I),I=1,NINT+1)              sld01
cccc               Read (LCard(IDL1:),*,ERR=9000,END=9000) NINT             sld01
cccc               Read (LCard(IDL1:),*,ERR=9000,END=9000) NINT,            sld01
cccc     1               (NSUBOB(I),I=1,NINT),(SUBINT(I),I=1,NINT+1)        sld01

csld determine what should be in CARD (lenth is 80)
c   Card would contain 80 columns with Name(which has been removed at
c    at this point),Dist,NumObs,Obs(i),1,NumObs,Endpoints(i),1,NumObs+1
c
      CARD(1:IDL) = LCard(1:IDL)
      Write(CARD(IDL1:),*,ERR=9090) NINT,(NSUBOB(I),I=1,NINT),
     x               (SUBINT(I),I=1,NINTP1)
c
c
               CALL CHKSTR ( DIST(IDSST(ID):IDSEND(ID)), CARD)
cccc end of section modifications
               If(KLLERR) Return
               CALL WRTCRD ( ID, NamVal)
CC               If(KLLERR) Return    WRTCRD has no error conditions	sld01
               IVarNm(IList) = NV
               IFound = 1
               Exit
            END IF
         End Do
c
c    end of RDPAR2 section
c
c        -- If IFound = 1 then we have successfully completed a Distribution
         IF (IFound /= 1) THEN
               KLLERR = .TRUE.
               IError = 1
               WRITE(*,9012) DISTYPE
               WRITE(4,9012) DISTYPE
               WRITE(99,9012) DISTYPE
               RETURN
         END IF
c     Set current distribution number for return to LHS_SDIST
      IDISTNO = NV
c     Set the point value index number for return to LHS_SDIST
      IPVNO = IList
c     Set flag to indicate that LHS_SDIST has been called
      LDIST = 1
      Return
c
 9001 FORMAT('1',5X,'LHS_PREP has been called prematurely ',/,5X,
     x'Call LHS_PREP just before call to LHS_RUN')
 9006 FORMAT(//,5x,'LHS_INIT or LHS_INIT_MEM must be called before ',
     x'any other LHS Input-By-Call Subroutines')
 9010 FORMAT('1',5X,'Variable Name exceeds 16 characters, NAMVAR = '
     x,A)
 9012 FORMAT('1',5X,'Distribution type not found, distribution name: '
     x,A)
 9013 FORMAT('1',5X,'Error: When the keyword LHSPVAL 0 is present, ',
     x/,10X,'all point values must be specified.',/,10X,
     x'Point value missing for ',A)
 9014 FORMAT(//,5X,'Point Value = ',1PE10.3,' provided for variable ',A
     x,/,8X,' will not be used, option "LHSPVAL" must be set to zero by'
     x,/,8X, ' calling LHS_OPTIONS prior to distribution call')
 9015 FORMAT('1',5X,'Multiple definitions found for ', A)
 9040 FORMAT(//,5X,'Variable Name is all blanks')
 9090 WRITE(9091)  Name
 9091 FORMAT(//,5X,'Error in SDist on interal read of Card, Name=',A)
      KLLERR = .True.
      IError = 1
      Return
c
      END SUBROUTINE
