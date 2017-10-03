C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   10:31 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_SAMEAS
      SUBROUTINE LHS_SAMEAS(NEWNAM,ORGNAM,IError,IDISTREF,IPVNO)
c     LHS_SAMEAS allows the user to specify an alias for a named distribution
c     LHS_SAMEAS calls routine: LJUST
c
c   Descriptions of call list parameters:
c   Inputs:
c     NEWNAM = new variable name, character
c     ORGNAM = old variable name, character
c   Outputs:
c     IError = error flag, returned = 1 to indicate some error occurred
c     IDISTREF = distribution reference number
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

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     Call List Variables:
      INTEGER :: IDISTREF,IPVNO,IError
      CHARACTER(LEN=*) :: NEWNAM,ORGNAM

c     LHS Internal Variables:
      CHARACTER(LEN=NamLen) :: Name
CCCC      CHARACTER(LEN=35) :: NamVal
cc    CHARACTER(LEN=32768) :: LCard, > dimension was needed in RDPAR2
      CHARACTER(LEN=40) :: LCard
CCCC      Character Card*(LENC)
cccc      LOGICAL Err  used in RDPAR2 to track multiple error conditions

c
c  check to see if message file is open,
c    if not then open as a scratch file
      IF (IScrh6 == 0) THEN
c        Open scratch file, S6
         OPEN(4, FILE='S4', Form='FORMATTED')
c    Modified by SFW on 3/18/2002
c    1            Carriage Control='FORTRAN')
         IScrh6 = 1
      END IF
c
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
c
c  check to be sure NEWNAM is <= 16 characters and not all blanks
c  terminate trailing blanks; eliminate any leading blanks
      LCard = NEWNAM
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
         WRITE(*,9010) NEWNAM
         WRITE(99,9010) NEWNAM
         WRITE(4,9010) NEWNAM
         Return
      END IF
      NAME = LCard
c
c test to be sure new name is not already in the list
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
cccc                  Write(6,*) 'Multiple definitions found for ', Name
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
c set LCard to ORGNAM
c  check to be sure ORGNAM is <= 16 characters and not all blanks
c  terminate trailing blanks; eliminate any leading blanks
      LCard = ORGNAM
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
         WRITE(*,9010) ORGNAM
         WRITE(99,9010) ORGNAM
         WRITE(4,9010) ORGNAM
         Return
      END IF
      NAME = LCard
c take "SAME AS" loop from RDPAR2
c        -- Check for the "Same As" distribution type
csld go through the following code:
CCCC         If ( LCard(1:8) == 'SAME AS ' ) Then
CCCC            Card = LCard(9:)
CCCC            Call LJust(Card)
cc            If(KLLERR) Return     LJust has no error conditions       sld01
cc            Read (Card,*,Err=9000) Name                               sld01
CCCC            Read (Card,*,Err=9000,END=9000) Name                        sld01
c           -- check that there is only one name on the line
CCCC            ISp = Index(Card,' ')
CCCC            Card = Card(ISp:)
CCCC            Call LJust(Card)
cc            If(KLLERR) Return   LJust has no error conditions         sld01
CCCC            If ( Card /= ' ' ) Then
CCCC               Card = LCard
CCCC               Print *, 'Error: too many names found on a Same As ',
CCCC     1                  'distribution definition record.'
CCCC               Print *, List(IList), ' SAME AS ', Card(1:50)
c****** add prints to message and error files
CCCC               Write(99,*) 'Error: too many names found on a Same',
CCCC     1                  ' As distribution definition record.'
CCCC               Write(99,*) List(IList), ' SAME AS ', Card(1:50)
CCCC               Write(6,*) 'Error: too many names found on a Same',
CCCC     1                  ' As distribution definition record.'
CCCC               Write(6,*) List(IList), ' SAME AS ', Card(1:50)
c******
CCCC               Err = .True.
CCCC            End If
c           -- check the name against the existing list of names
            IFound = 0
            Do i=1, NNames
               If ( Name == List(i) ) Then
                  IFound = 1
                  ISame = i
                  Exit
               End If
            End Do
c           -- If name not found in list, add to the list
            If (IFound == 0) Then
               NNames = NNames + 1
               ISame = NNames
               List(ISame) = Name
            End If
c           -- Now set the distribution information for LHS
            IVarNm(IList) = -ISame
CCCC            Cycle
CCCC         End If
c     Set the point value index number for return to LHS_SAMEAS
      IPVNO = IList
c     Set the distribution reference number for return to LHS_SAMEAS
      IDISTREF = ISame
c
      RETURN
c
 9040 FORMAT(//,5X,'Variable Name is all blanks')
 9001 FORMAT('1',5X,'LHS_PREP has been called prematurely ',/,5X,
     x'Call LHS_PREP just before call to LHS_RUN')
 9006 FORMAT(//,5x,'LHS_INIT or LHS_INIT_MEM must be called before ',
     x'any other LHS Input-By-Call Subroutines')
 9010 FORMAT('1',5X,'Variable Name exceeds 16 characters, NAMVAR = '
     x,A)
c
 9015 FORMAT(//,5X, 'Multiple definitions found for ', A16)
c
      END SUBROUTINE
