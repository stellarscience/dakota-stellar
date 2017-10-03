C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    9:05 am
C****************************************************************
C SUBROUTINE CHKZRO CHECKS TO MAKE SURE THAT THE MINIMUM
C REQUIREMENTS FOR A SAMPLE HAVE BEEN MET
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::CHKZRO
      SUBROUTINE CHKZRO(N,NV,IRSET)
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
      USE KILLFILE                      
cc    only 2001 sld changes were comments                               sld01
cc    CHKZRO is called from:  RDPAR,RDPAR2                              sld01
cc    CHKZRO does not call any other external routines                  sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(N.EQ.0)THEN
        WRITE(4,9001)
        WRITE(99,9001)
        KLLERR = .TRUE.
        RETURN
c      -- it is now legal to have zero variables if all information
c      -- is constants, so remove this check and make modifications
c      -- in the main program to bypass sampling.
c      ELSE IF(NV.EQ.0)THEN
c        WRITE(4,9002)
c        STOP 'CHKZRO'
      ELSE IF(IRSET.EQ.0)THEN
        WRITE(4,9003)
        WRITE(99,9003)
        KLLERR = .TRUE.
        RETURN
      ENDIF
      RETURN
 9001 FORMAT('1',5X,'THE NUMBER OF OBSERVATIONS HAS NOT BEEN ',
     1       'SPECIFIED')
 9002 FORMAT('1',5X,'NO VARIABLES HAVE BEEN SPECIFIED')
 9003 FORMAT('1',5X,'A RANDOM SEED HAS NOT BEEN SPECIFIED')
      END
