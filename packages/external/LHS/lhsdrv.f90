C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Jun 101   10:23 am
       Program LHSDRV
c      New driver for LHS Program developed to accomodate conversion
c      to DLLs for Visual Basic Interface with stand-alone version of code
c
C      INCLUDE 'KILLFILE.INC'                                           GDW-96  
      USE KILLFILE                      
C      INCLUDE 'PARMS.INC'                                              GDW-96  
      USE PARMS                         
C      INCLUDE 'CPARAM.INC'                                             GDW-96  
      USE CPARAM
C
cc    added because of storage used in COROUT for VCTR1,VCTR2 used by
cc    LHS_XXXX input by call routines
      USE InByCall                                                      SLD
cc
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
cc
        CALL InByCall_INIT()                                            SLD
cc
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
       KLLERR = .False.
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

C
c
c  NOTE:  Keyword file was previously obtained in the READ subroutine.
c  The location was moved to accomodate LHS conversion to DLL for use
c  in the Visual Basic program; putting control of keyword input in driver
c  routine allows all the other routines to remain the same for both the
c  Visual Basic and stand-alone versions.
c
c
c     -- Get the keyword file name from the command line.
c     -- If none specified, offer the user the opportunity to input
c     -- one from the keyboard or accept the default KEYWORD.DAT
c
      If (IPrint > 0 ) Then
        Print *, 'Welcome to LHS - The Latin Hypercube Sampling Program'
        Print *, ' '
      End If
c------------------------------------------------------------------------
c     NOTE: We're no longer using GetCL to get the command line arguments
c     because it is a Lahey Fortran for Windows call which is _NOT_
c     portable; we're using the standard F77 subroutine getarg.
c
c     2/17/2004 - slbrow
c------------------------------------------------------------------------
C$$$  Call GetCL(CmdLin)      
      call getarg(1,CmdLin)

      If (CmdLin /= ' ') Then
          LenFil = LEN_TRIM(CmdLin)
          print *, '***',LenFil,'***'
         Open (5, File=CmdLin(1:LenFil), Status="OLD", Err=10)
         Go To 100
 10      Print *, 'Error opening file specified on command line.'
         Print *
      End If
 15   If (IPrint == 0 ) Then
        Print *, 'Welcome to LHS - The Latin Hypercube Sampling Program'
        Print *, ' '
      End If
      Print *, 'Enter the name of the LHS input file to be read,'
      Print *, 'enter / to exit LHS, or enter . to accept the '
      Print *, 'default input file name (KEYWORD.DAT):  '
      Print *
      Read (*,9001) CmdLin
 9001 FORMAT (A128)
      If (CmdLin == '/') Then
         Print *, 'Program Terminated'
         Stop
      Else If (CmdLin == '.') Then
         Open (5, File='KEYWORD.DAT', Status="OLD", Err=20)
         CmdLin = 'KEYWORD.DAT'
         Go To 100
      Else
         LenFil = LEN_TRIM(CmdLin)
         Open (5, File=CmdLin(1:LenFil), Status="OLD", Err=20)
         Go To 100
      End If
c
 20   Print *, 'Error opening the file ', CmdLin
      Print *
      Go To 15
c
c     -- Now ready to start reading the keyword file.
 100  Continue
      Call LHS
C
      If (KLLERR) then
c     An Error has occurred
         Write (99,*) 'Error was detected during LHS run'
         Close (99, Status = 'Keep')
       Else
c      No error has occurred, close and delete error file
C        -- Since this is a normal program termination, close the LHS
C        -- Error file with Status=Delete so that other processors can
C        -- know that the program terminated successfully.
         Close (99, Status='DELETE')
       Endif
C
       Close (6)
C
        CALL LOCALVARS_CLOSE()
        CALL DISTNM_CLOSE()
        CALL CPARAM_CLOSE()
cc
        CALL InByCall_CLOSE()                                           SLD
cc
        CALL CSAMP_CLOSE()
        CALL CWORKC_CLOSE()
        CALL CWORKX_CLOSE()
        CALL CRANK_CLOSE()
        CALL CCMATR_CLOSE()
        CALL STAR_CLOSE()
        CALL UICORR_CLOSE()
        CALL CHRCRD_CLOSE()
        CALL OBSTR_CLOSE()
        CALL PDMAT_CLOSE()
        CALL FIRSTS_CLOSE()
C
        CALL PRAMS_CLOSE()
C
      If (IPrint > 0 ) Then
        PRINT *
        PRINT *, "***** LHS COMPLETE *****"
        PRINT *
      End If
C
       Stop
C
       end

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
         DOUBLE PRECISION DEFAULTRNUM1
         RNUMLHS1 = DEFAULTRNUM1() ! defined in DefaultRnum1.f90
         RETURN
       END

       DOUBLE PRECISION FUNCTION RNUMLHS2()
cc     See comments for RNUMLHS1()
         DOUBLE PRECISION DEFAULTRNUM2
         RNUMLHS2 = DEFAULTRNUM2()  ! defined in DefaultRnum2.f90
         RETURN
       END




