C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  11 Jul 101    8:13 am
      PROGRAM InpDrv
c     This program is contains generic example calls for the
c     Input-by-Subroutine Call Version of LHS
c     (S.L.Daniel April-May 2001)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(LEN=80) :: LHSOPTS,LHSMSG,LHSTITL,LHSOUT,FLOPTS
      CHARACTER(LEN=16) :: NAMVAR,NEWNAM,ORGNAM,NAM1,NAM2
      CHARACTER(LEN=16) :: LSTDNAM(200)
      CHARACTER(LEN=30) :: DISTYPE
      INTEGER :: LHSOBS,LHSSEED,LHSREPS,LHSPVAL,IPTFLAG,NUMPRMS
      INTEGER :: NUMINTV,NOBSpI(20),NUMPTS,INDXNAM(200),NUMVAR
      INTEGER :: NUMNAM,IDISTNO,IPVNO,IDISTREF,NumCor
      DOUBLE PRECISION :: PTVAL,APRAMS(4),XVAL(50),YVAL(50),ENDPNTS(20)
      DOUBLE PRECISION :: CORRVAL,PTVALST(200),SMATX(170,250),
     x C1MATX(170,170),C2MATX(170,170)
      PARAMETER (MAXNAM = 200, MAXVAR = 170, MAXOBS = 250)

c
c
c  Description of input parameters:
c     LHSOBS - integer number of observations requested (LHS keyword)
c     LHSSEED - integer random number seed (LHS keyword)
c     LHSREPS = Integer value for multiple number of samples (LHS keyword)
c     LHSPVAL = Integer value for selection of point value (LHS keyword)
c             = 0, point value MUST be specified by user for each distribution
c             = 1, LHS calculates the mean of the distribution (Default)
c             = 2, LHS calculates the median of the distribution
c     LHSOPTS = character string that allows user control of sampling
c               and pairing methods; String may contain "RANDOM PAIRING",
c               "RANDOM SAMPLE", or both  (LHS keyword)
c     LHSOUT - character string containing the LHS output file name
c              (this file contains the sample information) (LHS keyword)
c     LHSMSG - character string containing the LHS message file name
c              (LHS keyword)
c     LHSTITL - character string containing the Tile of the LHS run
c              (LHS keyword)
c     FLOPTS - character string containing any requested options
c     FLOPTS available are:  CORR HIST DATA LHNONAM LHSSCOL LHSWCOL
c     (as described in manual for keywords:  LHSRPTS, LHNONAM,
C                                            LHSWCOL and LHSSCOL)
c     NAMVAR - name of the variable
c     IPTFLAG - integer flag associated with point value, PTVAL
c     IPTFLAG = 0, if PTVAL is NOT to be used
c     IPTFLAG = 1, if PTVAL is to be used
c     PTVAL = point value associated with NAMVAR
c     DISTYPE = character string naming distribution type
c     APRAMS =  array of real numbers that are the distribution parameters
c     NUMPRMS = number of parameters in array APRAMS
c     NEWNAM = new variable name, character
c     ORGNAM = old variable name, character
c     NUMINTV = integer number of intervals "n" in a n-part histogram
c     NOBSpI = integer series of "n" integer values "obs", each "obs(i)"
c              represents the number of observations LHS will draw
c              from the i-th interval
c     ENDPNTS = array (type real) consisting of:
c                 first_point (the distribution minimum)
c                 followed by "n" interval endpoints
c         Note:  The sum of all interval observations must equal the number of
c                observations requested for the LHS run
c     NUMPTS = number of points in the user-defined distribution
c     XVAL = real array containing the x coordinates for i=1,NUMPTS points
c     YVAL = real array containing the y coordinates for i=1,NUMPTS points
c     NAM1 = first variable name, character
c     NAM2 = second variable name, character
c     CORRVAL = correlation value parameter, real # between -1 and 1
c
c  Description of output parameters:
c     IError - integer error flag
c        If an error is found, IError = 1 is returned
c     IDISTNO = current distribution number
c     IPVNO = point value index number
c     IDISTREF = distribution reference number
c     NUMNAM - the integer number of variable names including "same as"
c     NUMVAR- the number of variables for which observations
c                are generated (with no duplication for same as)
c     LSTDNAM - a list of distribution names (type character)
c     INDXNAM - integer array containing index number(position)
c               of names in sample data
c     PTVALST - an array of the associated point values (type real)
c     SMATX - the sample data matrix dimension (MAXVAR,MAXOBS)
c     C1MATX - original correlation matrix in lower triangle;
c              adjusted (to positive definite) matrix in upper triangle
c     C2MATX - raw data correlation matrix created by sample data in lower
c              triangle; rank data correlation matrix in upper triangle
c     NUMCOR - the total number of correlations requested
c     IPosDef - indicator for positive definite check; returned with value
c               of one if adjustment to matrix was required, zero otherwise
c  NAMELIST parameters used in LHS_INIT_MEM call list
c  Note: Any NAMELIST parameter can be set to -1 to accept LHS default value
c     LNMAX - Maximum number of observations
c     LMAXNNV  - Maximum number of variables * number of samples
c     LNVAR - Maximum number of variables
c     LNINTMX Maximum number of sub-intervals for any uniform* and
c                loguniform* distributions
c     LNCVAR - Maximum number of correlations (pairs)
c     LMAXTB  - Maximum number of entries that can be in a table
c     LIPrint - LPrint is used to control the amount of printing
C           0 = Nothing to the screen
C           1 = Normal
c     LISamW - LISamW controls the width of the Sample Output File
C           0 = Narrow -- Single Column
C           1 = Normal (compiler default -- 80 characters - I think)
C           2 = Very wide (for input into spreadsheets)
c
      OPEN(15, FILE='TestRun.out')
      OPEN(16, FILE='InpData.txt')
c
      IPTFLAG = 0
      PTVAL = 0.0
      IError = 0
c
      LHSOBS = 250
      LHSSEED = 1
c     Call to initialize the LHS code for Data Input:
c     LHS_INIT or LHS_INIT_MEM must be the first Input-by-Call
c     routine called for each run sequence
      CALL LHS_INIT(LHSOBS,LHSSEED,IError)
c OR (but not both)
      LNMAX =  2000
      LMAXNNV = 100000
      LNVAR = 50
      LNINTMX = -1
      LNCVAR = 50
      LMAXTB =  -1
      LIPrint = -1
      LISamW =  2
      Call LHS_INIT_MEM(LHSOBS,LHSSEED,LNMAX,LMAXNNV,LNVAR,
     xLNINTMX,LNCVAR,LMAXTB,LIPrint,LISamW,IError)
c
      LHSREPS = 1
      LHSPVAL = 1
      LHSOPTS = ' '
c     Call to set LHS Options:
c     (This call is optional)
      CALL LHS_OPTIONS(LHSREPS,LHSPVAL,LHSOPTS,IError)
c
c
      LHSOUT = '1NewIn.Out'
      LHSMSG = '2NewIn.Out'
      LHSTITL = ' Test New Input Routines'
      FLOPTS = ' CORR '
c     Call to provide file information and file options:
c     (This call is optional)
      CALL LHS_FILES(LHSOUT,LHSMSG,LHSTITL,FLOPTS,IError)
c
c
      NAMVAR = 'dist1'
      NUMPRMS = 2
      DISTYPE = 'normal'
      APRAMS(1) = 12.0
      APRAMS(2) = 7.0
c     Calls to specify distribution information (At least one DIST
c        routine must be called:
c     Call to input distributions with known number of parameter:
      CALL LHS_DIST(NAMVAR,IPTFLAG,PTVAL,DISTYPE,APRAMS,NUMPRMS,
     x IError,IDISTNO,IPVNO)
c
c
      NAMVAR =  'PI'
      PTVAL = 3.14159
c     Call to specify a "constant" distribution:
c     (This call is optional)
      CALL LHS_CONST(NAMVAR,PTVAL,IError,IPVNO)
c
c
      NEWNAM = 'DistSame4'
      ORGNAM = 'dist4'
c     Call to input an alias for a named distribution (i.e."Same As"):
c     (This call is optional)
      CALL LHS_SAMEAS(NEWNAM,ORGNAM,IError,IDISTREF,IPVNO)
c
c
      PTVAL = 0.0
      NAMVAR = 'Dist103'
      DISTYPE = 'uniform*'
      NUMINTV = 3
      NOBSpI(1) = 50
      NOBSpI(2) = 50
      NOBSpI(3) = 150
      ENDPNTS(1) = 0
      ENDPNTS(2) = 4
      ENDPNTS(3) = 5
      ENDPNTS(4) = 10
c     Call to input distributions with the "*" number of parameters:
c     (i.e.  UNIFORM*  or LOGUNIFORM*)
      CALL LHS_SDIST(NAMVAR,IPTFLAG,PTVAL,DISTYPE,
     x               NUMINTV,NOBSpI,ENDPNTS,IError,IDISTNO,IPVNO)
c
c
      NAMVAR = 'dist62 '
      DISTYPE ='discrete cumulative'
      NUMPTS = 3
      XVAL(1) = 7.8
      YVAL(1) =  .33
      XVAL(2) = 9.2
      YVAL(2) =  .666
      XVAL(3) = 11.4
      YVAL(3) =  1.0
c     Call to input distributions with user defined types:
c     (i.e. discrete cumulative, continuous linear, continuous frequency,
c           discrete histogram, continuous logarithmic)
      CALL LHS_UDIST(NAMVAR,IPTFLAG,PTVAL,DISTYPE,
     x               NUMPTS,XVAL,YVAL,IError,IDISTNO,IPVNO)
c
c
      NAM1 =  'dist115'
      NAM2 =  'dist117'
      CORRVAL =  0.7
c     Call to specify a correlation matrix entry:
c     (This call is optional)
      CALL LHS_CORR(NAM1,NAM2,CORRVAL,IError)
c
c
c     Call to check data input prior to running LHS
c     (This call is mandatory)
      CALL LHS_PREP(IError,NUMNAM,NUMVAR)
c
c
C  SET ABOVE:   MAXNAM = 200, MAXVAR = 170, MAXOBS = 250)
c     Call to run LHS sampling routine:
c     (This call is mandatory)
      CALL LHS_RUN(MAXVAR,MAXOBS,MAXNAM,IError,
     x LSTDNAM,INDXNAM,PTVALST,NUMNAM,SMATX,NUMVAR)
c
c  Example print:
c     Print the list of variables, their point values and vector index
      WRITE(15,9090)
 9090 FORMAT(/////,9X,'I',2X,'NAME',19X,'PT. VALUE',6X,'INDEX'/)
      do I = 1,NUMNAM
         WRITE (15,9091) I,LSTDNAM(I),PTVALST(I),INDXNAM(I)
 9091 FORMAT (5X,I5,2X,A16,2X,1PE15.4,2X,I8)
      end do
c     Print the sample matrix
      WRITE(15,9092)
 9092 FORMAT (/////,5X,'MATRIX OF OBSERVATIONS:',//)
      DO I = 1,LHSOBS
        WRITE(15,*) I,NUMVAR,(SMATX(J,I),J=1,NUMVAR)
      END DO
c     This routine is not implemented in current Beta Version
c     but is being alpha tested by 6410
c     Call to retrieve correlation matrices:
c     (This call is optional)
      CALL LHS_COROUT(MAXVAR,IError,C1MATX,C2MATX,
     x NumCor,NumVar,IPosDef)
c
c     Call to retrieve last random number seed used by LHS
c     (This call is optional)
      call LHS_RtvSEED(IError,LastSeed)
c
c     Call to terminate LHS calls and deallocate storage used by LHS
c     (Although not mantory -- It should be called to terminate LHS properly)
      CALL LHS_CLOSE(IError)

      STOP
c
 9088 FORMAT('1',5X,'LHS KEYWORDS INPUT:',//,10X,'LHSREPS = ',I5,/,
     x 10X,'LHSPVAL = ',I5,/,10X,'LHSOPTS = ',A,/,10X,'FLOPTS = ',A,
     x /,10X,'LHSTITL = ',A,/10X,'LHSOUT = ',A,/,10X,'LHSMSG = ',A,/,
     x 10X,'LHSOBS = ',I5,/,10X,'LHSSEED = ',I5,/10X,'LHSREPS = ',I5,
     x /,10X,'LHSPVAL = ',I5,/,10X,'LHSOPTS = ',A//)
c
      END PROGRAM
