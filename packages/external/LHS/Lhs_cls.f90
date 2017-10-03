C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  22 Jun 101    8:46 am
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::LHS_CLOSE
      SUBROUTINE LHS_CLOSE(IError)
c     LHS_CLOSE terminates LHS calls and deallocates storage
c
      USE KILLFILE
      USE PARMS                         
cc    PARMS provides:  ISamW,IPtVal
      USE CPARAM                        
cc    CPARAM provides: N,ISeed,I1Col,NamOut
      USE InByCall
cc    InByCall provides subroutine flags: LINT,LPREP,LRUN,LFILES,LDIST
C
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
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      If (KLLERR) then
         IError = 1
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
C      close scratch files:
       call fileoc(0)
c      close message file
       IF (IScrh6 == 1) THEN
          Close (4, Status='DELETE')
       Else
          CLOSE (4)
       END IF

C
        CALL LOCALVARS_CLOSE()
        CALL DISTNM_CLOSE()
        CALL CPARAM_CLOSE()
        CALL InByCall_CLOSE()
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
      RETURN


      END SUBROUTINE
