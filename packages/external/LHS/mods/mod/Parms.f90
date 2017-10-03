C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  10 Jul 101    8:35 am
C**********************************************************
C     -- THIS FILE IS AN INCLUDE BLOCK THAT CONTAINS ALL OF
C     -- THE COMMON PARAMETERS THAT ARE USED IN LHS.
C
c      PARAMETER (NMAX=20001, MAXNNV=200000)
c      PARAMETER (NVAR=401, NINTMX=401, NCVAR=801)
c      PARAMETER (LENT=125, LENC=80, NAMLEN=16)
c      PARAMETER (MAXTB=5001)

c above dimensions changed 12-2-96                 gdw
C      PARAMETER (NMAX=24576, MAXNNV=200000)
C      PARAMETER (NVAR=1024, NINTMX=64, NCVAR=1024)
C      PARAMETER (LENT=125, LENC=80, NAMLEN=16)
C      PARAMETER (MAXTB=5001)
c
c   description of parameters:
c       NMAX   - Maximum number of observations
c       MAXNNV - Maximum number of variables * number of samples       
c       NVAR   - Maximum number of variables
c       NINTMX - Maximum number of sub-intervals for any uniform* and
c                loguniform* distributions
c       NCVAR  - Maximum number of correlations (pairs)
c       LENT   - Maximum length of title card
c       LENC   - Maximum length of a single (one card) of input
c       NAMLEN - Maximum length of a name
c       MAXTB  - Maximum number of entries that can be in a table 
C
C     -- END OF PARAMETER DECLARATION
C**********************************************************
C
      MODULE PARMS
C
C       These are the elements of the old include file
C
c        dimensions changed 12-2-96                 gdw
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NMAX
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::MAXNNV
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NVAR
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NINTMX
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::NCVAR
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::MAXTB
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::IPrint
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::ISamW
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER :: NMAX=24576, MAXNNV=200000
        INTEGER :: NVAR=1024, NINTMX=64, NCVAR=1024
        INTEGER, PARAMETER :: LENT=125, LENC=256, NAMLEN=16
        INTEGER :: MAXTB=5001
C
C        IPrint is used to control the amount of printing
C           0 = Nothing to the screen
C           1 = Normal
C           2 = Debug printout (someday)
        INTEGER :: IPrint = 1
C
C        ISamW controls the width of the Sample Output File
C           0 = Narrow -- Single Column
C           1 = Normal (compiler default -- 80 characters - I think)
C           2 = Very wide (for input into spreadsheets)
        INTEGER :: ISamW = 1
C
C       Now here is the initialization routine
      CONTAINS
C
ccc      SUBROUTINE PRAMS_INIT()                                        SLD
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PRAMS_INIT
      SUBROUTINE PRAMS_INIT(UseFile)                                    SLD
C
C       Define these variables to be part of a namelist, then
C       provide default values for them.
C
        NAMELIST /LHS/ NMAX, MAXNNV, NVAR, NINTMX, NCVAR,
     1                 MAXTB, IPrint, ISamW
c                                                                       SLD
ccc        LOGICAL :: YESNO                                             SLD
        LOGICAL :: YESNO,UseFile,LInpMem                                SLD
        OPTIONAL UseFile                                                SLD
c                                                                       SLD
C
        NMAX=24576
        MAXNNV=200000
cc        NVAR=512                                                      sld01
        NVAR=1024
        NINTMX=64
        NCVAR=1024
        MAXTB=5001
        IPrint = 1
        ISamW = 1
c                                                                       SLD
c     UseFile variable added to allow by-pass of SIPRA.INI for the      SLD
c     Input-By-Call version of LHS -- LHS_INIT_MEM allows user to set   SLD
c     variables in the NAMELIST input from SIPRA.INI                    SLD
      IF (PRESENT(UseFile)) THEN                                        SLD
         LInpMem = UseFile                                              SLD
      ELSE                                                              SLD
         LInpMem = .True.                                               SLD
      END IF                                                            SLD
c                                                                       SLD
      IF (LInpMem) THEN                                                 SLD
C
C       Now open the initialization file, read the namelist variables,
C       and close the namelist file to get any changes to these default
C       values.
C
        INQUIRE (FILE="SIPRA.INI", EXIST=YESNO)
C
        IF ( YESNO ) THEN
C
           OPEN (19, FILE="SIPRA.INI", ERR=200, ACTION="READ")
           READ (19, NML=LHS, ERR=100, END=100)
           CLOSE (19)
C
        ELSE
C
           GOTO 200
C
        END IF
      END IF                                                            SLD
C
        RETURN
C
C  = = = = = = = = =  ERROR HANDLING SECTION  = = = = = = = = = = = = =
C
C       An error condition occurred while reading the file.  Close it,
C       write an error message, and go on using the defaults.
C
 100    WRITE (*,*)
        WRITE (*,*) "**** Error reading file SIPRA.INI. ****"
        WRITE (*,*) "****   Default dimensions used.    ****"
        WRITE (*,*)
C
        CLOSE (19)
C
        RETURN
C
C       An error condition occurred while opening the file.
C       Write an error message, and go on using the defaults.
C
 200    WRITE (*,*)
        WRITE (*,*) "**** Error opening file SIPRA.INI. ****"
        WRITE (*,*) "****   Default dimensions used.    ****"
        WRITE (*,*)
C
        RETURN
C
      END SUBROUTINE
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::PRAMS_CLOSE
      SUBROUTINE PRAMS_CLOSE()
C
C       Nothing to deallocate.
C
        RETURN
C
      END SUBROUTINE
C
      END MODULE
