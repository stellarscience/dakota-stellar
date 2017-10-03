C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101   11:19 am
      PROGRAM TestRun
c     This program is a test driver for the Input-by-Subroutine Call
c     Version of LHS  (S.L.Daniel April-May 2001)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(LEN=80) :: LHSOPTS,LHSMSG,LHSTITL,LHSOUT,FLOPTS
      CHARACTER(LEN=16) :: NAMVAR,NEWNAM,ORGNAM,NAM1,NAM2
      CHARACTER(LEN=16) :: LSTDNAM(200)
      CHARACTER(LEN=30) :: DISTYPE
      INTEGER :: LHSOBS,LHSSEED,LHSREPS,LHSPVAL,IPTFLAG,NUMPRMS
      INTEGER :: NUMINTV,NOBSpI(20),NUMPTS,INDXNAM(200),NUMVAR
      INTEGER :: NumCor
      INTEGER :: RFLAG
      DOUBLE PRECISION :: PTVAL,APRAMS(4),XVAL(50),YVAL(50),ENDPNTS(20)
      DOUBLE PRECISION :: CORRVAL,PTVALST(200),SMATX(170,250),
     x C1MATX(170,170),C2MATX(170,170),RMATX(170,250)
c     SMATX is sample matrix (50,500) stored as a single vector
c     use statement function:  LOC(I,J)=(J-1)*N+I where N is number
c                              of samples, I is index forsample number,
c                              and J is index for number of variables
c     To print use:     where NV is number of variables
c            Do I=1, N
c               Write(1,*)I,NV,(X(LOC(I,J)),J=1,NV)
c            End Do
      PARAMETER (MAXNAM = 200, MAXVAR = 170, MAXOBS = 250,
     x           MAXNUMC = 50)
      OPEN(15, FILE='TestRun.out')
      OPEN(16, FILE='InpData.txt')
      IError = 0
c
ccc      LHSOBS = 500    read from input
      LHSSEED = 1
c
      LHSREPS = 1
CC
c
      LHSOUT = '1NewIn.Out'
      LHSMSG = '2NewIn.Out'
      LHSTITL = ' Test New Input Routines'
c     FLOPTS available are:  CORR HIST DATA LHNONAM LHSSCOL
c     (as described in manual for keywords:  LHSRPTS LHSNONAM and LHSSCOL)
cc      FLOPTS = ' CORR '
      FLOPTS = ' '
c
      IPTFLAG = 0
      PTVAL = 0
cc use with "fixed" point values:
cc      IPTFLAG = 1
cc and PTVAL's read with input data
c
      READ(16,*) LHSOBS
      READ(16,9085) FLOPTS
 9085 FORMAT(A80)

c     Call to initialize the LHS code with NAMELIST variables included
      LNMAX =  25000
      LMAXNNV = 190000
      LNVAR = 1500
      LNINTMX = 75
      LNCVAR = 1500
      LMAXTB =  6001
      LIPrint =  -1
      LISamW =  2
      Call LHS_INIT_MEM(LHSOBS,LHSSEED,LNMAX,LMAXNNV,LNVAR,
     xLNINTMX,LNCVAR,LMAXTB,LIPrint,LISamW,IError)
      WRITE(15,*) ' returned from LHS_INIT_MEM'

c     Call to initialize the LHS code for Data Input:
c     (User must call either LHS_INIT or LHS_INIT_MEM but not both)
cc      CALL LHS_INIT(LHSOBS,LHSSEED,IError)
cc      WRITE(15,*) ' returned from LHS_INIT'

cc use with options
      LHSPVAL = 1
cc      LHSPVAL = 2
cc for user specified point value:
cc      LHSPVAL = 0
      LHSOPTS = ' '
cc      LHSOPTS = ' RANDOM Pairing RANDOM SAMPLE'
c     Call to set LHS Options:
c     (This call is optional)
      CALL LHS_OPTIONS(LHSREPS,LHSPVAL,LHSOPTS,IError)
cc      WRITE(15,*) ' returned from LHS_OPTIONS'

c
c     Call to provide file information and file options:
c     (This call is optional)
      CALL LHS_FILES(LHSOUT,LHSMSG,LHSTITL,FLOPTS,IError)
cc      WRITE(15,*) ' returned from LHS_FILES'
c
c
      WRITE(15,9088) LHSREPS,LHSPVAL,LHSOPTS,FLOPTS,LHSTITL,
     xLHSOUT,LHSMSG,LHSOBS,LHSSEED,LHSREPS,LHSPVAL
c  Set up Case structure for inputing distribution data:
c
      JSTOP = 0
      DO WHILE (JSTOP == 0)
c
c        read case indiactor, IVAL
         READ(16,*) IVAL
         WRITE(15,*) IVal
         SELECT CASE (IVAL)
c
c        IVAL = -1 for calls to LHS_DIST
         CASE (-1)
           READ(16,*) NAMVAR
           WRITE(15,*) ' NAMVAR = ',NAMVAR
c
cc           READ(16,*) PtVal
cc           WRITE(15,*) ' Point Value = ',PtVal
c
           READ(16,9080) DISTYPE
           WRITE(15,*) ' DISTYPE = ',DISTYPE
           READ(16,*) NUMPRMS
           WRITE(15,*) ' NUMPRMS = ',NUMPRMS
           READ(16,*) (APRAMS(K),K=1,NUMPRMS)
           WRITE(15,*) ' APRAMS = ',(APRAMS(K),K=1,NUMPRMS)
c     Calls to specify distribution information (At least one DIST
c        routine must be called:
c     Call to input distributions with known number of parameter:
      CALL LHS_DIST(NAMVAR,IPTFLAG,PTVAL,DISTYPE,APRAMS,NUMPRMS,
     x IError,IDISTNO,IPVNO)
      WRITE(15,9089)  IDISTNO,IPVNO
 9089 FORMAT(5X,' RETURNED FROM LHS_DIST',//,5X,'IDISTNO = ',I5,/
     X,5X,'IPVNO = ',I5)
c
c        IVAL = -5 for calls to LHS_CONST
         CASE (-5)
           READ(16,*) NAMVAR
           WRITE(15,*) ' NAMVAR = ',NAMVAR
c
cc           READ(16,*) PtVal
cc           WRITE(15,*) ' Point Value = ',PtVal
c
           READ(16,*) PTVAL
           WRITE(15,*) ' PTVAL = ',PTVAL
         CALL LHS_CONST(NAMVAR,PTVAL,IError,IPVNO)
         WRITE(15,*) ' RETURNED FROM LHS_CONST, IPVNO = ',IPVNO
c
c
c        IVAL = -4 for calls to LHS_SAMEAS
          CASE (-4)
            READ(16,*) NEWNAM
            WRITE(15,*) ' NEWNAM = ',NEWNAM
            READ(16,*) ORGNAM
            WRITE(15,*) ' ORGNAM = ',ORGNAM
c     Call to input an alias for a named distribution (i.e."Same As"):
c     (This call is optional)
      CALL LHS_SAMEAS(NEWNAM,ORGNAM,IError,IDISTREF,IPVNO)
      WRITE(15,9084) IDISTREF,IPVNO
 9084 FORMAT(//,' RETURNED FROM LHS_SAMEAS',/,5X,'IDISTREF = ',I5,
     x/,5X,'IPVNO = ',I5)
c
c
c        IVAL = -3 for calls to LHS_SDIST
          CASE (-3)
            READ(16,*) NAMVAR
            WRITE(15,*) ' NAMVAR = ',NAMVAR
c
cc           READ(16,*) PtVal
cc           WRITE(15,*) ' Point Value = ',PtVal
c
            READ(16,9080) DISTYPE
            WRITE(15,*) ' DISTYPE = ',DISTYPE
            READ(16,*) NUMINTV
            WRITE(15,*) ' NUMINTV = ',NUMINTV
            READ(16,*) (NOBSpI(K),K=1,NUMINTV)
            WRITE(15,*) ' NOBSpI(i) = ',(NOBSpI(K),K=1,NUMINTV)
            READ(16,*) (ENDPNTS(K),K=1,NUMINTV+1)
            WRITE(15,*) ' ENDPNTS(i) = ',(ENDPNTS(K),K=1,NUMINTV+1)
c***************
c     Call to input distributions with the "*" number of parameters:
c     (i.e.  UNIFORM*  or LOGUNIFORM*)
      CALL LHS_SDIST(NAMVAR,IPTFLAG,PTVAL,DISTYPE,
     x               NUMINTV,NOBSpI,ENDPNTS,IError,IDISTNO,IPVNO)
 9082 FORMAT(5X,' RETURNED FROM LHS_SDIST',//,5X,'IDISTNO = ',I5,/
     X,5X,'IPVNO = ',I5)
       WRITE(15,9082) IDISTNO,IPVNO
c

c        IVAL = -2 for calls to LHS_UDIST
          CASE (-2)
            READ(16,*) NAMVAR
            WRITE(15,*) ' NAMVAR = ',NAMVAR
c
cc           READ(16,*) PtVal
cc           WRITE(15,*) ' Point Value = ',PtVal
c
            READ(16,9080) DISTYPE
            WRITE(15,*) ' DISTYPE = ',DISTYPE
            READ(16,*) NUMPTS
            WRITE(15,*) ' NUMPTS = ',NUMPTS
            READ(16,*) (XVAL(K),YVAL(K),K=1,NUMPTS)
            WRITE(15,*) ' XVAL,YVAL = ',(XVAL(K),YVAL(K),K=1,NUMPTS)
c***************
c     Call to input distributions with user defined types:
      CALL LHS_UDIST(NAMVAR,IPTFLAG,PTVAL,DISTYPE,
     x               NUMPTS,XVAL,YVAL,IError,IDISTNO,IPVNO)
 9081 FORMAT(5X,' RETURNED FROM LHS_UDIST',//,5X,'IDISTNO = ',I5,/
     X,5X,'IPVNO = ',I5)
       WRITE(15,9081) IDISTNO,IPVNO

c
c     IVAL = -6 correlation of distributions
         CASE (-6)
           READ(16,*) NAM1
           WRITE(15,*) ' NAM1 = ',NAM1
           READ(16,*) NAM2
           WRITE(15,*) ' NAM2 = ',NAM2
           READ(16,*) CORRVAL
           WRITE(15,*) ' CORRVAL = ',CORRVAL
c     Call to specify a correlation matrix entry:
c     (This call is optional)
      CALL LHS_CORR(NAM1,NAM2,CORRVAL,IError)
cc      WRITE(15,*) '  Returned from LHS_CORR '

c     IVAL = -99 to stop distribution input
         CASE (-99)
            JSTOP = 1
            Exit
c
         CASE DEFAULT
c        something is wrong; print and stop
            WRITE(15,*)' Error on Data input, Case index incorrect',Ival
            stop
         END SELECT
c
      END DO
c
c End of Distribution information input

c     Call to input an alias for a named distribution (i.e."Same As"):
c     (This call is optional)
ccc      CALL LHS_SAMEAS(NEWNAM,ORGNAM,IError,NUMVAR)

c     Call to specify a "constant" distribution:
c     (This call is optional)
ccc      CALL LHS_CONST(NAMVAR,PTVAL,IError,NUMVAR)
c
c     Call to specify a correlation matrix entry:
c     (This call is optional)
ccc      CALL LHS_CORR(NAM1,NAM2,CORRVAL,IError)
c
c     Call to check data input prior to running LHS
c     (This call is mandatory)
      CALL LHS_PREP(IError,NUMNAM,NUMVAR)
      WRITE(15,*) ' Have returned from LHS_PREP '

c     LPSwiler addition 9/27/06, to reflect capability to accept and 
c     return rank matrix.  RFLAG value of 0 indicates normal operation.      
      RFLAG = 0
c
c     Call to run LHS sampling routine:
c     (This call is mandatory)

      CALL LHS_RUN(MAXVAR,MAXOBS,MAXNAM,IError,
     x LSTDNAM,INDXNAM,PTVALST,NUMNAM,SMATX,NUMVAR,RMATX,RFLAG)


c     Print the list of variables, their point values and vector index
      WRITE(15,9090)
 9090 FORMAT(/////,9X,'I',2X,'NAME',19X,'PT. VALUE',6X,'INDEX'/)
      do I = 1,NUMNAM
         WRITE (15,9091) I,LSTDNAM(I),PTVALST(I),INDXNAM(I)
 9091 FORMAT (5X,I5,2X,A16,2X,1PE15.4,2X,I8)
      end do
c     Print the sample matrix
ccc      WRITE(15,9092)
ccc 9092 FORMAT (/////,5X,'MATRIX OF OBSERVATIONS:',//)
ccc      DO I = 1,LHSOBS
ccc        WRITE(15,*) I,NUMVAR,(SMATX(J,I),J=1,NUMVAR)
ccc      END DO

      IF (FLOPTS .eq. 'CORR') THEN
c     Call to retrieve correlation matrices:
c     (This call is optional)
      CALL LHS_COROUT(MAXVAR,IERROR,C1MATX,C2MATX,
     xNumCor,NumVar,IPosDef)

      WRITE(15,*) ' IPosDef = ',IPosDef
c
      WRITE(15,*) ' Input Rank Correlation Matrix :'
      DO I = 1,NumVar
         WRITE(15,9077) (C1MATX(I,J),J=1,I)
      END DO
c
      WRITE(15,*) ' Full Matrix Print for C1MATX '
      do I = 1,NumVar
         WRITE(15,9077) (C1MATX(I,J),J=1,NumVar)
      end do
c
c
      WRITE(15,*) ' Raw Data Correlation Matrix :'
      DO I = 1,NumVar
         WRITE(15,9077) (C2MATX(I,J),J=1,I)
      END DO
c
      WRITE(15,*) ' Rank Data Correlation Matrix :'
      DO I = 1,NumVar
         WRITE(15,9077) (C2MATX(I,J),J=I,NumVar)
 9077 FORMAT(/,5X,15F8.4)
      END DO
c
      WRITE(15,*) ' Full Matrix Print for C2MATX '
      do I = 1,NumVar
         WRITE(15,9077) (C2MATX(I,J),J=1,NumVar)
      end do
      END IF

c
c     Call to retrieve last random number seed used by LHS
c     (This call is optional)
      Call LHS_RtvSEED(IError,LastSeed)
c
      WRITE(15,*) ' Last random number seed = ', LastSeed
c

c     Call to terminate LHS calls and deallocate storage used by LHS
c     (Although not mantory -- It should be called to terminate LHS properly)
      CALL LHS_CLOSE(IError)

      STOP
 9080    FORMAT(A30)
 9088 FORMAT(//,5X,'LHS KEYWORDS INPUT:',//,10X,'LHSREPS = ',I5,/,
     x 10X,'LHSPVAL = ',I5,/,10X,'LHSOPTS = ',A,/,10X,'FLOPTS = ',A,
     x /,10X,'LHSTITL = ',A,/10X,'LHSOUT = ',A,/,10X,'LHSMSG = ',A,/,
     x 10X,'LHSOBS = ',I5,/,10X,'LHSSEED = ',I5,/10X,'LHSREPS = ',I5,
     x /,10X,'LHSPVAL = ',I5)
c
      END PROGRAM

       DOUBLE PRECISION FUNCTION RNUMLHS1()
cc     This function and the next (RNUMLHS2) were added as part of a
cc     refactor to get random number generation working correctly on
cc     the Mac. Previously, RNUMLHS1 & 2 were defined in separate .f90 files.
cc     According to comments there, they were called by:
cc                           ETA,BINOM,CUMULC,CUMULD,ENTRPY,EXPON,     sld01
cc                           GAMMA,GAMMAB,GAMMAM,GEOM,HGEOM,IGAUS,	sld01
cc                           IGAUSF,MIX,NBINOM,NORMAL,PARETO,POISON,    sld01
cc                           TRIANG,UNIFRM,WEIBUL                       sld01
cc Note:  With addition of second random number generator routine       SLD
cc        IGAUSF,GAMMAB,GAMMAM now call RNUMLHS2 instead of RNUMLHS1          SLD
cc     
Cc     For standalone LHS, this remains true.  When LHS is compiled as a 
cc     library for use in Dakota, Dakota exports rnumlhs1 and 2 for use by LHS.
cc     (See packages/pecos/src/BoostRNG_Monostate.hpp.) These functions returns
cc     random numbers produced either by the f90 function DEFAULTRNUM1 or mt13997.
cc     JAS
         RNUMLHS1 = DEFAULTRNUM1() ! defined in DefaultRnum1.f90
         RETURN
       END

       DOUBLE PRECISION FUNCTION RNUMLHS2()
cc     See comments for RNUMLHS1()
         RNUMLHS2 = DEFAULTRNUM2()  ! defined in DefaultRnum2.f90
         RETURN
       END
