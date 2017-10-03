C***********************************************************************
C LHS (Latin Hypercube Sampling) UNIX Library/Standalone. 
C Copyright (c) 2004, Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
C
C This software is distributed under the GNU Lesser General Public License.
C For more information, see the README file in the LHS directory. 
C***********************************************************************
C     Last change:  SLD   5 Apr 101    7:30 am
C****************************************************************
C SUBROUTINE HISTO GENERATES HISTOGRAMS OF THE VARIABLES
C
!LHS_EXPORT_DEC ATTRIBUTES DLLEXPORT::HISTO
      SUBROUTINE HISTO
cc    HISTO is called from routine:  HSTOUT                             sld01
cc    HISTO calls routine:  SIFT                                        sld01
C     INCLUDE 'KILLFILE.INC'                                            GDW-96  
cc      USE KILLFILE            not needed				sld01
C
C     -- THIS IS A MODIFICATION OF DISTAT -- IT MAKES THE HISTOGRAM
C     -- AND PRINTS THE MEAN AND THE MEDIAN
C           N = NO. OF OBSERVATIONS
C           XV = ARRAY NAME
C
C     INCLUDE 'PARMS.INC'                                               GDW-96  
cc      USE PARMS                                                       sld01
C     INCLUDE 'CPARAM.INC'                                              GDW-96  
      USE CPARAM                        
cc    CPARAM provides:  N,IRS                                           sld01
C     INCLUDE 'CRANK.INC'                                               GDW-96  
      USE CRANK                         
cc    CRANK provides:  XV array                                         sld01
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IFLG=20)
      CHARACTER*1 KX
      DATA KX / 'X' /
C
      IF (N-1) 10,20,30
   10 WRITE(4,260)
      RETURN
   20 WRITE(4,270)XV(1)
      RETURN
   30 CONTINUE
      CALL SIFT (XV,N)
cc      If(KLLERR) Return     SIFT has no error conditions              sld01
      SUM=0.0
      SUMSQ=0.0
      DO 40 I=1,N
      SUM=SUM+XV(I)
      SUMSQ=SUMSQ+XV(I)*XV(I)
   40 CONTINUE
      XMEAN=SUM/FLOAT(N)
      XVAR=(SUMSQ-FLOAT(N)*XMEAN*XMEAN)/FLOAT(N)
      IF(IRS.NE.0)XVAR=XVAR*FLOAT(N)/FLOAT(N-1)
      N2=N/2
      LOC=(N+1)/2-N2
      FMED=(XV(N2+LOC)+XV((N+1-LOC)/2+1))/2.0
      R=XV(N)-XV(1)
      IF (R.NE.0.0) GO TO 50
      WRITE(4,280)
      GO TO 250
   50 CONTINUE
      M=1
      SIZE=R/IFLG
      POWER=LOG10(SIZE)
      IF (POWER.GE.0.) GO TO 60
      POWER=AINT(POWER)
      GO TO 70
   60 POWER=AINT(POWER)+1.
   70 SIZE=SIZE/10.**POWER
      CON=.01
   80 IF (SIZE.LE.(CON+.005)) GO TO 90
      CON=CON+.01
      GO TO 80
   90 SIZE=CON*10.**POWER
      TEMP=XV(1)/SIZE
cc    3-Way test on floating ok as it is really a <= condition          sld01
      IF (TEMP) 100,100,110
  100 TEMP=AINT(TEMP-1.0)
      GO TO 120
  110 TEMP=AINT(TEMP)
  120 IF (XV(1)-TEMP*SIZE.GT.0.0) GO TO 130
      TEMP=TEMP-0.5
  130 CELL=TEMP*SIZE+SIZE
      POINT=TEMP*SIZE+SIZE*0.5
      I=0
      KSUM=0
      WRITE(4,290)
      K=-1
      M=1
  140 I=I+1
      K=K+1
  150 IF (I-N) 160,160,190
  160 IF (XV(I)-CELL) 140,140,170
  170 IF (K) 180,180,200
  180 WRITE(4,310)POINT,K
      GO TO 240
  190 M=2
  200 IF (K-90) 210,210,220
  210 KK=K
      GO TO 230
  220 KK=90
  230 WRITE(4,310)POINT,K,(KX,J=1,KK)
      KSUM=KSUM+K
      IF (M.GT.1) GO TO 250
      K=0
  240 CELL=CELL+SIZE
      POINT=POINT+SIZE
      GO TO 150
  250 WRITE(4,320)KSUM
      WRITE(4,300)XV(1),XV(N),R,XMEAN,FMED,XVAR
      RETURN
C
  260 FORMAT(' N is Zero',//)
  270 FORMAT(' One Obs. ',1PE17.8,//)
  280 FORMAT(' No Histogram - Range =0',/)
  290 FORMAT(/,5X,'Midpoint',10X,'Freq.',/)
  300 FORMAT(//,6X,'Min',12X,'Max',11X,'Range',11X,'Mean',10X,
     1       'Median',8X,'Variance',//,1X,6(1PE15.7),/)
  310 FORMAT(1X,1PE15.7,5X,0P,I5,4X,90A1)
  320 FORMAT('0',20X,I5)
      END
