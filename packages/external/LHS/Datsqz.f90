C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    9:15 am
C****************************************************************
C SUBROUTINE DATSQZ PROCESSES PARAMETER CARDS WHICH REQUIRE
C CONVERTING CHARACTER DATA TO INTEGER DATA
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::DATSQZ
      SUBROUTINE DATSQZ(CARD,CRDTYP,TMPCRD)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE
cc    only 2001 sld changes were comments                               sld01
cc    DATSQZ is called from:  RDPAR                                     sld01
cc    DATSQZ does not call any other external routines                  sld01                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LENTT=11)
      CHARACTER CARD*(*),CRDTYP*(*),TMPCRD*(LENTT),BLANK,MINUS
      PARAMETER (BLANK=' ',MINUS='-')
      IZERO=ICHAR('0')
      ININE=ICHAR('9')
      LENC=LEN(CARD)
      IC=0
C
C     -- SEARCH FOR BEGINNING OF NON-BLANK CHARACTER STRING
   10 CONTINUE
      IC=IC+1
      IF(IC.GT.LENC)THEN
        WRITE(4,9001)CRDTYP
        WRITE(99,9001)CRDTYP
        KLLERR = .TRUE.
        RETURN
      ENDIF
      IF(CARD(IC:IC).EQ.BLANK)GO TO 10
      IBEG=IC
      IF(CARD(IC:IC).EQ.MINUS)GO TO 20
      ITEST=ICHAR(CARD(IC:IC))
      IF(ITEST.LT.IZERO.OR.ITEST.GT.ININE)THEN
        WRITE(4,9002)CRDTYP,CARD(IC:IC)
        WRITE(99,9002)CRDTYP,CARD(IC:IC)
        KLLERR = .TRUE.
        RETURN
      ENDIF
C
C     -- SEARCH FOR ENDING OF NON-BLANK CHARACTER STRING
   20 CONTINUE
      IC=IC+1
      IF(IC.GT.LENC)GO TO 30
      IF(CARD(IC:IC).EQ.BLANK)GO TO 30
      ITEST=ICHAR(CARD(IC:IC))
      IF(ITEST.LT.IZERO.OR.ITEST.GT.ININE)THEN
        WRITE(4,9002)CRDTYP,CARD(IC:IC)
        WRITE(99,9002)CRDTYP,CARD(IC:IC)
        KLLERR = .TRUE.
        RETURN
      ENDIF
      GO TO 20
C
C     -- MOVE NON-BLANK CHARACTER STRING INTO TMPCRD RIGHT-JUSTIFIED
   30 CONTINUE
      IEND=IC-1
      ILEN=IEND-IBEG+1
      IF(ILEN.GT.LENTT)THEN
        WRITE(4,9003)CRDTYP,ILEN,LENTT
        WRITE(99,9003)CRDTYP,ILEN,LENTT
        KLLERR = .TRUE.
        RETURN
      ENDIF
      TMPCRD=BLANK
      IT=LENTT-ILEN
      DO 40 I=IBEG,IEND
        IT=IT+1
        TMPCRD(IT:IT)=CARD(I:I)
   40 CONTINUE
      RETURN
 9001 FORMAT('1',5X,'THE PARAMETER CARD ',A,'CONTAINS NO DATA')
 9002 FORMAT('1',5X,'THE PARAMETER CARD ',A,'CONTAINS THE ',
     1       'NON-NUMERIC CHARACTER ',A)
 9003 FORMAT('1',5X,'THE DATA ON PARAMETER CARD ',A,'CONTAINS ',I2,
     1       ' DIGITS',/,6X,'THE MAXIMUM NUMBER OF DIGITS ALLOWED ',
     2       'IS ',I2)
      END
