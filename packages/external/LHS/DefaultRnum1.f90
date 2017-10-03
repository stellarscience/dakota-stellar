C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  JAS  09 Oct 2013    3:10 pm
C****************************************************************
C FUNCTION DEFAULTRNUM1 IS USED TO CREATE A RANDOM NUMBER.
C
      DOUBLE PRECISION FUNCTION DEFAULTRNUM1()
cc    DEFAULTRNUM1 calls routines:  RMCNP,RMCNPI                               sld01
cc    DEFAULTRNUM1 is called from RNUMLHS1 in standalone LHS.    JAS-13
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE      -- not needed (RMCNP & RMCNPI set no errors)  sld01
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  ISEED                                           sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION RMCNP
      EXTERNAL RMCNP
C
      CALL RMCNPI(ISEED)
cc      If(KLLERR) Return  -- RMCNPI nor RMCNP do not set errors        sld01
      DEFAULTRNUM1 = RMCNP(ISEED)
C
      RETURN
      END
