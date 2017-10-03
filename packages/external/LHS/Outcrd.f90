C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    9:02 am
C****************************************************************
C SUBROUTINE OUTCRD PROCESSES THE OUTPUT PARAMETER OPTIONS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::OUTCRD
      SUBROUTINE OUTCRD(CARD)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
cc    only 2001 sld changes were comments                               sld01
cc    OUTCRD is called from:  RDPAR,READ                                sld01
cc    OUTCRD does not call any other external routines                  sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARMS provides:  LENC                                             sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  IDATA,IHIST,ICORR                               sld01

C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER CARD*(LENC),PDATA*4,PHIST*4,PCORR*4,BLANK
      PARAMETER (PDATA='DATA',PHIST='HIST',PCORR='CORR',BLANK=' ')
C
      IC=7
   10 CONTINUE
      IC=IC+1
      IF(IC.GT.LENC)GO TO 20
      IF(CARD(IC:IC).EQ.BLANK)GO TO 10
      IE=IC+3
      IF(CARD(IC:IE).EQ.PDATA)THEN
        IDATA=1
        IC=IE+1
        GO TO 10
      ELSE IF(CARD(IC:IE).EQ.PHIST)THEN
        IHIST=1
        IC=IE+1
        GO TO 10
      ELSE IF(CARD(IC:IE).EQ.PCORR)THEN
        ICORR=1
        IC=IE+1
        GO TO 10
      ELSE
        WRITE(4,9001)CARD
        WRITE(99,9001)CARD
        KLLERR = .TRUE.
        RETURN
      ENDIF
   20 CONTINUE
      RETURN
 9001 FORMAT('1',5X,'THE FOLLOWING OUTPUT OPTION CARD REQUESTED ',
     1       'AN UNDEFINED OUTPUT OPTION',/,6X,'PLEASE CHECK THE ',
     2       'USER MANUAL FOR THE CORRECT OUTPUT OPTION CARD ',
     3       'SYNTAX',//,3X,'***',A,'***')
      END
