C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  JAS  09 Oct 2013   4:23 pm
C****************************************************************
C FUNCTION DEFAULTRNUM2 IS USED TO CREATE A RANDOM NUMBER.
C
      DOUBLE PRECISION FUNCTION DEFAULTRNUM2()
c     DEFAULTRNUM2 is a second random number generating routine added to
c       provide random numbers used in the Acceptance-Rejection scheme
c       used to generate the GAMMA and INVERSE GAUSSIAN distributions
cc    DEFAULTRNUM2 calls routines:  RMCNP2,RMCNPI2
cc    DEFAULTRNUM2 is called from:  RNUMLHS2                                JAS-13
C
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE    -- not needed (RMCNP2 & RMCNPI2 set no errors)  sld01
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  JSEED                                           sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION RMCNP2
      EXTERNAL RMCNP2
C
      CALL RMCNPI2(JSEED)
cc      If(KLLERR) Return  -- RMCNPI2 nor RMCNP2 do not set errors      sld01
      DEFAULTRNUM2 = RMCNP2(JSEED)
C
      RETURN
      END
