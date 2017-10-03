C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   11:27 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_PREP
      SUBROUTINE LHS_PREP(IError,NUMNAM,NUMVAR)
c     Subroutine LHS_PREP must be called prior to LHS_RUN
c       It checks that all necessary information has been provided
c       to make an LHS run
c     LHS_PREP returns:
c     IError = 1, to inidicate some error was found
c     NUMNAM - the integer number of variable names including "same as"
c     NUMVAR- the number of variables for which observations
c                are generated (with no duplication for same as)
c

c     LHS_PREP calls routines:  WRTPAR,FileOC,BANNER,ChkZro,ChkDim
c                               CMCRD, PMTRX,POSDEF,CHLSKY
c
      USE KILLFILE
      USE PARMS
cc    PARMS provides:  IPrint,NMax
      USE CPARAM
c     CPARAMS provides:  N,NV,IRSet,ICM,Mfile,Sfile,TreeFl,Cmdlin
      USE InByCall
cc    InByCall provides subroutine flags: NNames,LINIT,LPREP,LFILES,
cc                                        IScrh1,IScrh6,LRUN,LPOSDEF
      USE UICORR
c     UICORR provides: NCV and arrays ICVAR,JCVAR,CVAR
C     INCLUDE 'CCMATR.INC'                                              GDW-96
      USE CCMATR                        
cc    CCMATR provides:  NCM and arrays CORR,LCM                         sld01
C     INCLUDE 'CWORKC.INC'                                              GDW-96  
      USE CWORKC                        
cc    CWORKC provides:  Q and S arrays                                  sld01

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C====================
      CHARACTER*11 TIMRES
C====================
c     next lines added 1-17-96 to convert date/time routine to LF90
      integer :: dt(8)
      character (len=5) :: zone
      character * 8 DATRES
      character * 2 Chr2
C
      CHARACTER CDUM1*5, CDUM2*12
      DATA CDUM1,CDUM2/' ',' '/
cccc  LOGICAL Err  used in RDPAR2 to track multiple error conditions
      LOGICAL Err
C
      LOC1(I,J)=J+(I*I-I)/2
C
      Err = .False.
c
c     see if unit 6 is open, if not open a scratch for message info
      IF (IScrh6 == 0) THEN
c        open scratch file
         OPEN(1, FILE='S1')
         Mfile = ' '
         IScrh6 = 1
      END IF
c
c Check to see that INIT has been called and PREP has not
c     Test that LHS_INIT has been called
      IF (LINIT /= 1) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9026)
         WRITE (99,9026)
         WRITE(4,9026)
         Return
      END IF
c     Test that LHS_PREP has not been called previously
      IF (LPREP /= 0) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9022)
         WRITE (99,9022)
         WRITE(4,9022)
         Return
      END IF
c     Test that at least one call has been made to a LHS_xDIST routine
      IF (LDIST /= 1) THEN
         KLLERR = .TRUE.
         IError = 1
         WRITE(*,9023)
         WRITE (99,9023)
         WRITE(4,9023)
         Return
      END IF
c

cccc Following section is from stand LHS routine:
      Call date_and_time (DATRES,TIMRES,zone,dt)
      Chr2 = TIMRES(8:9)
C
      If (IPrint > 0 ) Write (*,9010) 'START TIME:            ',
     1                                dt(5), dt(6), dt(7), Chr2
C====================
C
c     see if unit 1 is open, if not open a scratch for banner info
      IF (IScrh1 == 0) THEN
c        open scratch file with default format, ISamW = 1
         OPEN(1, FILE='S1')
         IScrh1 = 1
      END IF

      CALL BANNER(1)

c
cccc Following section from RDPAR2
c
c     -- If LHS was not required to run in this case (NNames=0), the
c     -- LHS program should die gracefully here.
      If (NNames == 0) Then
         Print *, 'No information for LHS to process.'
         Print *, 'LHS termination.'
         Write(99,*) 'No information for LHS to process.'
         Write(99,*) 'LHS termination.'
         Write(4,*) 'No information for LHS to process.'
         Write(4,*) 'LHS termination.'
         Call FileOC(0)
         If(KLLERR) Return
         KLLERR = .TRUE.
c        Added IError flag setting not in original RDPAR2 code
         IError = 1
c
         RETURN
      End If
c
c     -- Now check that each named data item found has a valid distribution
c     -- associated with it.  If not, it appeared on a Correlate record
c     -- but either had no definition or was defined in a file or
c     -- post-processor record.  Either way, it is an error and LHS can
c     -- not continue processing.  Also, one "Same As" record can not point
c     -- to another "Same As" record - it must point to a defined
c     -- distribution, so check "Same As" records to assure that the
c     -- IVarNm pointed to is positive.
c
      Do i=1, NNames
         If (IVarNm(i) == 0) Then
            Print *, 'Distribution required, not found for ', List(i)
c****** add prints to message and error files
            Write(99,*) 'Distribution required, not found for '
     1, List(i)
            Write(4,*) 'Distribution required, not found for '
     1, List(i)
c******
            Err = .True.
         End If
         If (IVarNm(i) < 0  .AND.  IVarNm(i) /= -9999999) Then
            j = -IVarNm(i)
            If (IVarNm(j) < 0  .AND.  IVarNm(j) /= -9999999) Then
               Print *, 'One distribution of type Same As points to ',
     1            'another Same As distribution.'
               Print *, 'This is invalid.  Each Same As distribution ',
     1            'may only point to an'
               Print *, 'explicitly defined distribution.'
               Print *, List(i), ' Same As ', List(j)
c****** add prints to message and error files
               Write(99,*) 'One distribution of type Same As points',
     1            ' to another Same As distribution.'
               Write(99,*) 'This is invalid.  Each Same As ',
     1            'distribution may only point to an'
               Write(99,*) 'explicitly defined distribution.'
               Write(99,*) List(i), ' Same As ', List(j)
               Write(4,*) 'One distribution of type Same As points',
     1            ' to another Same As distribution.'
               Write(4,*) 'This is invalid.  Each Same As ',
     1            'distribution may only point to an'
               Write(4,*) 'explicitly defined distribution.'
               Write(4,*) List(i), ' Same As ', List(j)
c******
               Err = .True.
            End If
         End If
      End Do
c
c     -- Now replace the current contents of ICVAR and JCVAR
c     -- They were originally set to point to the List array contents,
c     -- but they need to correspond to the order in which the distributions
c     -- were defined (IVarNm(i)).  This transformation made here to prepare
c     -- the data for the CMCRD routine.
c
cccc  Note: NCV is incremented in LHS_CORR, initialized in LHS_INIT
      If (ICM == 1) Then
         Do i=1, NCV
            ij = ICVar(i)
            ICVar(i) = IVarNm(ij)
            jj = JCVar(i)
            JCVar(i) = IVarNm(jj)
c
            If ( ICVar(i) == -9999999 ) Then
               Print *, 'Error: Distributions of type Constant can ',
     1            'not be correlated to other distributions.'
               Print *, List(ij),' is CONSTANT'
c****** add prints to message and error files
               Write(99,*) 'Error: Distributions of type Constant can',
     1            ' not be correlated to other distributions.'
               Write(99,*) List(ij),' is CONSTANT'
               Write(4,*) 'Error: Distributions of type Constant can',
     1            ' not be correlated to other distributions.'
               Write(4,*) List(ij),' is CONSTANT'
c******
               Err = .True.
            Else If ( ICVar(i) <= 0 ) Then
               Print *, 'Error: Distributions of type Same As can ',
     1            'not be correlated to other distributions.'
               Print *, List(ij),' SAME AS ', List(-ICVar(i))
c****** add prints to message and error files
               Write(99,*) 'Error: Distributions of type Same As can',
     1            ' not be correlated to other distributions.'
               Write(99,*) List(ij),' SAME AS ', List(-ICVar(i))
               Write(4,*) 'Error: Distributions of type Same As can',
     1            ' not be correlated to other distributions.'
               Write(4,*) List(ij),' SAME AS ', List(-ICVar(i))
c******
               Err = .True.
            End If
c
            If ( JCVar(i) == -9999999 ) Then
               Print *, 'Error: Distributions of type Constant can ',
     1            'not be correlated to other distributions.'
               Print *, List(jj),' is CONSTANT'
c****** add prints to message and error files
               Write(99,*) 'Error: Distributions of type Constant can',
     1            ' not be correlated to other distributions.'
               Write(99,*) List(jj),' is CONSTANT'
               Write(4,*) 'Error: Distributions of type Constant can',
     1            ' not be correlated to other distributions.'
               Write(4,*) List(jj),' is CONSTANT'
c******
               Err = .True.
            Else If ( JCVar(i) <= 0 ) Then
               Print *, 'Error: Distributions of type Same As can ',
     1            'not be correlated to other distributions.'
               Print *, List(jj),' SAME AS ', List(-JCVar(i))
c****** add prints to message and error files
               Write(99,*) 'Error: Distributions of type Same As can',
     1            ' not be correlated to other distributions.'
               Write(99,*) List(jj),' SAME AS ', List(-JCVar(i))
               Write(4,*) 'Error: Distributions of type Same As can',
     1            ' not be correlated to other distributions.'
               Write(4,*) List(jj),' SAME AS ', List(-JCVar(i))
c******
               Err = .True.
            End If
c
         End Do
      End If
c
c     -- Stop if an error was encountered during data checking
      If (Err) Then
         Write(99,*) ' Error was encountered during data checking'
         Write(4,*) ' Error was encountered during data checking'
         KLLERR = .TRUE.
c        Added IError flag setting not in original RDPAR2 code
         IError = 1
c
         RETURN
      Endif

c
c
c     -- Now check for a few more errors
c
      Call ChkZro(N,NV,IRSet)
cccc
cccc      If(KLLERR) Return
c        Added IError flag setting not in original RDPAR2 code
      IF (KLLERR) THEN
         IError = 1
         Return
      END IF
cccc
      Call ChkDim(2,NV,NVAR,CDUM1,CDUM2)
cccc      If(KLLERR) Return
c        Added IError flag setting not in original RDPAR2 code
      IF (KLLERR) THEN
         IError = 1
         Return
      END IF
cccc
      If (ICM == 1) then
         Call CMCRD
cccc         If(KLLERR) Return
c        Added IError flag setting not in original RDPAR2 code
         IF (KLLERR) THEN
            IError = 1
            Return
         END IF
cccc
      Endif
c

c     Following section is from Standard LHS.for
c     -- Start to set up the sampling process.
c
      CALL WRTPAR
cc    IF(KLLERR) Return   WRTPAR has no error conditions                sld01
C
C     -- IF THE USER HAS SPECIFIED A CORRELATION STRUCTURE (ICM=1) THEN
C     -- THE CORRELATION MATRIX IS ECHOED AND CHECKED TO MAKE SURE THAT
C     -- IT IS POSITIVE DEFINITE. THE CHOLESKY FACTORIZATION IS COMPUTED
C     -- TO BE USED LATER AS PART OF THE PROCESS FOR INDUCING THE DESIRED
C     -- CORRELTION STRUCTURE
C
      IF(ICM.NE.0)THEN
         CALL PMTRX(NCM,3)
cc    store input correlation data for retrieval by LHS_COROUT
      kkk = (NCM*(NCM+1))/2
      do jjj = 1,kkk
         VCTR1(jjj) = CORR(jjj)
      end do
cc    save row\column location information
      do jjj = 1,NCM
         LCMSav(jjj) = LCM(jjj)
      end do
cc
cc         If(KLLERR) Return  PMTRX has no error conditions             sld01
         IF(IRP.EQ.1)THEN
            IRP=0
            WRITE(4,9003)
         ENDIF
         IF(N.LE.NV)THEN
            ICM=0
            WRITE(4,9001)
            GO TO 10
         ENDIF
         CALL POSDEF(ITER)
         IF (ITER > 1) THEN
           LPOSDEF= 1
         Else
           LPOSDEF= 0
         END IF

cc    store adjusted correlation data for retrieval by LHS_COROUT
      do jjj = 1,kkk
         VCTR1(kkk+jjj) = CORR(jjj)
      end do
cccsld
         If(KLLERR) Return
         IF(ITER.GT.1)THEN
            WRITE(4,9002)
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
C
c     End of section from Standard LHS.for


c     -- Input has been successfully read - let's do the sampling!
cccc  End of RDPAR2 section

c     Set flag to indicate that LHS_PREP has been called
      LPREP = 1
c     Return number of variable names and number of variables to be
c     sampled to the user
      NUMNAM = NNames
      NUMVAR = NV
      Return
c
c
c     -- Process all errors that occur reading parameters from the input file
c
cccc 9000 Write (4,9005) NamVal, LCard(1:90), LCARD(91:180)
cccc      Write (99,9005) NamVal, LCard(1:90), LCARD(91:180)
cccc      KLLERR = .TRUE.
cccc      RETURN
c
c     -- Format Statements
c
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
 9010 FORMAT (5x,A,I2,":",I2,":",I2,".",A)
c
 9022 FORMAT('1',5X,'LHS_PREP should only be called once prior to ',
     x'calling LHS_RUN')
 9023 FORMAT('1',5x,'LHS_DIST, LHS_UDIST, or LHS_SDIST must be called',
     x' prior to calling LHS_PREP')
 9026 FORMAT(//,5x,'LHS_INIT or LHS_INIT_MEM must be called before ',
     x'any other LHS Input-By-Call Subroutines')
      END SUBROUTINE
