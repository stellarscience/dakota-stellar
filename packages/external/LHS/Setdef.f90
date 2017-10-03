C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD  27 Mar 101    9:12 am
C****************************************************************
C SUBROUTINE SETDEF SETS THE DEFAULT VALUES OF THE PARAMETERS
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::SETDEF
      SUBROUTINE SETDEF
C
cc    only 2001 sld changes were comments                               sld01
cc    SETDEF is called from:  RDPAR,READ                                sld01
cc    SETDEF does not call any other external routines                  sld01
C
c      INCLUDE 'KILLFILE.INC'  -- include not needed 10-96
C     INCLUDE 'PARMS.INC'                                               GDW-96  
      USE PARMS                         
cc    PARMS provides:  NVAR,LENT                                        sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,NV,IRS,NREP,IDATA,IHIST,ICOOR,ICM,IRP,IV1,    sld01
cc                      TITLE,IDST,List,IVarNam,PValue			sld01
C
C     These statements removed to make modules work - GDW-96
C     COMMON/OBSTR/NSTR,NOBSTR(NVAR)
      USE OBSTR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DO 10 I=1, LENT
         TITLE(I:I)=' '
   10 CONTINUE
      N=0
      NV=0
      IRS=0
      ICM=0
      NREP=1
      IRP=0
      IV1=0
      IDATA=0
      IHIST=0
      ICORR=0
      NSTR=0
      DO I=1,NVAR
         IDIST(I)=0
         NOBSTR(I)=0
c        -- List is the array that contains all of the names encountered
         List(i) = ' '
c        -- IVarNm(i) holds the distribution number that goes with List(i)
c        -- (the order in which the dist was defined, not the dist type)
         IVarNm(i) = 0
c        -- PValue(i) holds the optional point values read from the input
c        -- for each distribution name
         PValue(i) = 0.0
      End Do
C
      RETURN
      END
