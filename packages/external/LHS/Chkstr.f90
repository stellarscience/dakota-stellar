C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    8:54 am
C****************************************************************
C SUBROUTINE CHKSTR CHECKS PARAMETERS OF THE UNIFORM* AND THE
C LOGUNIFORM* DISTRIBUTIONS FOR CONSISTENCY
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CHKSTR
      SUBROUTINE CHKSTR(PAR,CARD)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
cc    only 2001 sld changes were comments                               sld01
cc    CHKSTR is called from:  RDPAR,RDPAR2                              sld01
cc    CHKSTR does not call any other external routines                  sld01
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARMS provides:  LENC,NINTMX                                      sld01
C
C     These statements removed to make modules work - GDW-96
c     COMMON/STAR/NSUBOB(NINTMX),SUBINT(NINTMX+1),NINT
      USE STAR
cc    STAR provides:  NINT,SUBINT,NSUBOB                                sld01
c     COMMON/CHRCRD/CRDSTR
      USE CHRCRD
cc    CHRCRD provides:  CRDSTR                                          sld01
c     COMMON/OBSTR/NSTR,NOBSTR(NVAR)
      USE OBSTR
cc    OBSTR provides:  NSTR,NOBSTR                                      sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER PAR*(*),CARD*(LENC),PLOG*3
c     CHARACTER*(LENC) CRDSTR(NVAR)
      PARAMETER (PLOG='LOG')
C
      IF(NINT.EQ.0)THEN
        WRITE(4,9001)PAR
        WRITE(99,9001)PAR
        KLLERR = .TRUE.
        RETURN
C
      ELSE IF(NINT.GT.NINTMX)THEN
        WRITE(4,9002)PAR,NINT,NINTMX
        WRITE(99,9002)PAR,NINT,NINTMX
        KLLERR = .TRUE.
        RETURN
C
      ELSE
        WRITE(8)NINT
C
      ENDIF
C
      NSTR=NSTR+1
      DO 100 I=1,NINT
C
        IF(PAR(1:3).EQ.PLOG.AND.SUBINT(I).LE.0.0)THEN
          WRITE(4,9003)PAR,I,SUBINT(I)
          WRITE(99,9003)PAR,I,SUBINT(I)
          KLLERR = .TRUE.
          RETURN
C
        ELSE IF(SUBINT(I).GE.SUBINT(I+1))THEN
          WRITE(4,9004)PAR,I,SUBINT(I),SUBINT(I+1)
          WRITE(99,9004)PAR,I,SUBINT(I),SUBINT(I+1)
          KLLERR = .TRUE.
          RETURN
C
        ELSE IF(NSUBOB(I).LT.0)THEN
          WRITE(4,9005)PAR,I
          WRITE(99,9005)PAR,I
          KLLERR = .TRUE.
          RETURN
C
        ELSE
          NOBSTR(NSTR)=NOBSTR(NSTR)+NSUBOB(I)
          WRITE(8)NSUBOB(I),SUBINT(I),SUBINT(I+1)
C
        ENDIF
C
  100 CONTINUE
C
      CRDSTR(NSTR)=CARD
C
      RETURN
C
 9001 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE NUMBER OF ',
     1       'SUBINTERVALS IS ZERO')
 9002 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE NUMBER OF ',
     1       'SUBINTERVALS REQUESTED ',I3,/,6X,'IS GREATER THAN THE ',
     2       'MAXIMUM NUMBER OF SUBINTERVALS CURRENTLY PERMITTED ',I3,
     3       /,6X,'PLEASE CONSULT THE USER MANUAL FOR INSTRUCTIONS ',
     4       'ON HOW TO ALLOW MORE SUBINTERVALS')
 9003 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE SUBINTERVAL ',
     1       'LIMIT FOR SUBINTERVAL ',I3,/,6X,'IS LESS THAN OR ',
     2       'EQUAL TO ZERO ',G20.10)
 9004 FORMAT('1',5X,'ON THE ',A,'DISTRIBUTION FOR SUBINTERVAL ',
     1       I3,' THE LOWER LIMIT ',G20.10,/,6X,'IS GREATER ',
     2       'THAN OR EQUAL TO THE UPPER LIMIT ',G20.10)
 9005 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION SUBINTERVAL ',I3,
     1       ' REQUESTED A NEGATIVE NUMBER OF OBSERVATIONS')
C
      END
