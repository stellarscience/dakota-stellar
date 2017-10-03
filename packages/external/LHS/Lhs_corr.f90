C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   11:10 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_CORR
      SUBROUTINE LHS_CORR(NAM1,NAM2,CORRVAL,IError)
c     LHS_CORR specfies a correlation matrix entry
c     LHS_CORR calls routines:
c   Descriptions of call list parameters:
c   Inputs:
c     NAM1 = first variable name, character
c     NAM2 = second variable name, character
c     CORRVAL = correlation value parameter, real #  between -1 and 1
c   Outputs:
c     IError = error flag, returned = 1 to indicate some error occurred
c
c
      USE KILLFILE
      USE PARMS
c     PARMS provides: NamLen,MAXTB,LENC
      USE CPARAM
c     CPARAM provides: IptVal, List, IVarNm, IDIST arrays,
C                      N,NV,ICM,ICORR
      USE InByCall
cc    InByCall provides subroutine flags: LINT,LPREP,LDIST,NNames,IScrh6
      USE DISTNM
c     DISTNM provides:  DIST,IDSST,IDSEND,IDSPAR,LEND,MAXPAR
C     COMMON/UICORR/ICVAR(NCVAR),JCVAR(NCVAR),CVAR(NCVAR),NCV
      USE UICORR
c     UICORR provides: NCV and arrays ICVAR,JCVAR,CVAR
C     INCLUDE 'CCMATR.INC'                                              GDW-96
      USE CCMATR                        
cc    CCMATR provides:  NCM, and arrays LCM and CORR                                      sld01
C     INCLUDE 'CWORKC.INC'                                              GDW-96  
      USE CWORKC                        
cc    CWORKC provides:  Q and S arrays                                  sld01
C====================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     Call List Variables:
      INTEGER :: IError
      CHARACTER(LEN=*) :: NAM1,NAM2
      DOUBLE PRECISION :: CORRVAL
c     LHS Internal Variables:
      CHARACTER(LEN=NamLen) :: NAME, Name2
CCCC      CHARACTER(LEN=35) :: NamVal
cc    CHARACTER(LEN=32768) :: LCard, > dimension was needed in RDPAR2
      CHARACTER(LEN=40) :: LCard
CCCC      Character Card*(LENC)
cccc      LOGICAL Err  used in RDPAR2 to track multiple error conditions
c
c
c  check to see if message file is open,
c    if not then open as a scratch file
      IF (IScrh6 == 0) THEN
c        Open scratch file, S6
         OPEN(4, FILE='S4', Form='FORMATTED')
c Modified by SFW on 3/18/2002 
c    1            Carriage Control='FORTRAN')
         IScrh6 = 1
      END IF

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
         WRITE(*,9051)
         WRITE (99,9051)
         WRITE(4,9051)
         Return
      END IF
c
c  check to be sure NAM1 is <= 16 characters and not all blanks
c  terminate trailing blanks; eliminate any leading blanks
      LCard = NAM1
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
         WRITE(*,9010) NAM1
         WRITE(99,9010) NAM1
         WRITE(4,9010) NAM1
         Return
      END IF
      NAME = LCard
c  check to be sure NAM2 is < 16 characters and not all blanks
c  terminate trailing blanks; eliminate any leading blanks
      LCard = NAM2
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
         WRITE(*,9010) NAM2
         WRITE(99,9010) NAM2
         WRITE(4,9010) NAM2
         Return
      END IF
      Name2 = LCard
c
c    Set correlation value from input
      Value = CORRVAL
c    note: Name,Name2,and Value are read from input in Standard LHS
c    following section is from RDPAR2
CCCC         If ( Card(1:10) == 'CORRELATE ' ) Then
c
c           -- Process correlation information
c
c           -- strip off the keyword and left-justify
CCCC            Card = Card(10:)
CCCC            Call LJust(Card)
cc            If(KLLERR) Return   LJust has no error conditions         sld01
c
c           -- first variable name is now left justified - get it
CCCC            ISpc = Index(Card,' ')
CCCC           Name = Card(1:ISpc-1)
CCCC            Card = Card(ISpc:)
CCCC            Call LJust(Card)
cc            If(KLLERR) Return    LJust has no error conditions        sld01
c
c           -- second variable name is now left justified - get it
CCCC           ISpc = Index(Card,' ')
CCCC            Name2 = Card(1:ISpc-1)
CCCC           Card = Card(ISpc:)
CCCC            Call LJust(Card)
cc            If(KLLERR) Return      LJust has no error conditions 
c
c           -- check to see that neither name is blank
            If ( Name == ' '  .OR.  Name2 == ' ' ) Then
               Print *,'Error: two non-blank names must be specified'
               Print *,'when specifying correlation between variables.'
               Print *,'The following names were found: ',Name
               Print *,'                           and: ',Name2
c****** add prints to message and error files
               Write(99,*)'Error: two non-blank names must be specified'
               Write(99,*)'when specifying correlation between '
     1 ,           'variables.'
               Write(99,*)'The following names were found: ',Name
               Write(99,*)'                           and: ',Name2
               Write(4,*)'Error: two non-blank names must be specified'
               Write(4,*)'when specifying correlation between '
     1 ,           'variables.'
               Write(4,*)'The following names were found: ',Name
               Write(4,*)'                           and: ',Name2
c******
CCCC              Err = .True.
            End If
c
c           -- Correlation value is now left justified - get it and check
c           -- that it is in the range -1.0 < Value < 1.0
CCCC            Read (Card,*,Err=151) Value

            If ( Name /= Name2  .AND.  Abs(Value) >= 1.0 ) Then
               Print *, 'Error: The absolute value of all correlation',
     1            ' values must be less than 1.0'
               Print *, 'The following correlation was found:'
               Print *, 'CORRELATE ', Name, Name2, Value
c****** add prints to message and error files
               Write(99,*) 'Error: The absolute value of all ',
     1            'correlation values must be less than 1.0'
               Write(99,*) 'The following correlation was found:'
               Write(99,*) 'CORRELATE ', Name, Name2, Value
               Write(4,*) 'Error: The absolute value of all ',
     1            'correlation values must be less than 1.0'
               Write(4,*) 'The following correlation was found:'
               Write(4,*) 'CORRELATE ', Name, Name2, Value
c******
CCCC              Err = .True.
            End If
c
c           -- Now increment the correlation counter, determine which
c           -- variable numbers correspond to the given names, and
c           -- store this pair's correlation data in ICVAR, JCVAR and CVAR
c
            IFound = 0
            JFound = 0
            Do i=1, NNames
               If ( Name == List(i) ) Then
                  IFound = 1
                  IList = i
                  If (JFound == 1) Exit
               End If
               If ( Name2 == List(i) ) Then
                  JFound = 1
                  JList = i
                  If (IFound == 1) Exit
               End If
            End Do
c
c           -- If name not found in list, add to the list
            If (IFound == 0) Then
               NNames = NNames + 1
               IList = NNames
               List(IList) = Name
            End If
            If (Name == Name2) Then
c              -- both names are the same - correlation must be 1.00
               JList = IList
               If ( Abs(1.0-Value) > 1.0E-06 ) Then
                  Print *, 'Error: If a variable is to be correlated ',
     1              'with itself, the correlation must be 1.0'
                  Print *, 'CORRELATE ', Name, Name2, Value
c****** add prints to message and error files
                  Write(99,*) 'Error: If a variable is to be correlated'
     1,              ' with itself, the correlation must be 1.0'
                  Write(99,*) 'CORRELATE ', Name, Name2, Value
                  Write(4,*) 'Error: If a variable is to be correlated',
     1              ' with itself, the correlation must be 1.0'
                  Write(4,*) 'CORRELATE ', Name, Name2, Value
c******
CCCC                 Err = .True.
               End If
            Else If (JFound == 0) Then
               NNames = NNames + 1
               JList = NNames
               List(JList) = Name2
            End If
c
c           -- Store this pair's correlation data in ICVAR, JCVAR and CVAR.
c           -- These will be checked in CMCRD.
            ICM = 1
            NCV = NCV + 1
            ICVAR(NCV) = IList
            JCVAR(NCV) = JList
            CVAR(NCV) = Value
c
c           -- now done with the correlation record so cycle and repeat the
c           -- whole process with a new record
CCCC            Cycle
c
c           -- Error reading value
CCCC 151        Print *, 'Correlation Value read was not a number'
CCCC            Print *, Name, Name2, Card
c****** add prints to message and error files
CCCC            Write(99,*) 'Correlation Value read was not a number'
CCCC            Write(99,*) Name, Name2, Card
CCCC            Write(4,*) 'Correlation Value read was not a number'
CCCC            Write(4,*) Name, Name2, Card
c******
CCCC           Err = .True.
c
CCCC         End If
c     end of section from RDPAR2
c
      RETURN
c
 9040 FORMAT(//,5X,'Variable Name is all blanks')
 9051 FORMAT('1',5X,'LHS_PREP has been called prematurely ',/,5X,
     x'Call LHS_PREP just before call to LHS_RUN, in CORR')
 9006 FORMAT(//,5x,'LHS_INIT or LHS_INIT_MEM must be called before ',
     x'any other LHS Input-By-Call Subroutines')
 9010 FORMAT('1',5X,'Variable Name exceeds 16 characters, NAMVAR = '
     x,A)
c
      END SUBROUTINE
