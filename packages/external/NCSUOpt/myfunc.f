C+-----------------------------------------------------------------------+
C| Program       : main.f (subfile myfunc.f)                             |
C| Last modified : 04-12-2001                                            |
C| Written by    : Joerg Gablonsky                                       |
C| Testfunctions from the literature to test global optimization methods.|
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| Wrapper for the test functions. In the last position of the array     |
C| iidata the problem id is stored. According to this number, this       |
C| routine calls the different subroutines containing the different      |
C| test functions.                                                       |
C+-----------------------------------------------------------------------+
      SUBROUTINE myfunc(n,x,f,flag, 
     +                  iidata, iisize, ddata, idsize, cdata, icsize)

      IMPLICIT NONE
      INTEGER n,flag
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f
      INTEGER problem
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)
      EXTERNAL teconst
      EXTERNAL telinear
      EXTERNAL tequadr
      EXTERNAL jonestest
      EXTERNAL branin
      EXTERNAL shekel
      EXTERNAL hartman
      EXTERNAL goldprice
      EXTERNAL sixhump
      EXTERNAL Shubert

      problem = iidata(iisize)
      IF (problem .eq. 0) then
         CALL teconst(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF (problem .eq. 1) then 
         CALL telinear(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
            ELSE IF (problem .eq. 2) then 
         CALL tequadr(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF (problem .eq. 3) then 
         CALL jonestest(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF (problem .eq. 4) then 
         CALL branin(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF ((problem .ge. 5) .and. (problem .LE. 7)) then 
         CALL shekel(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF ((problem .ge. 8) .and. (problem .LE. 9)) then 
         CALL hartman(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF (problem .eq. 10) then 
         CALL goldprice(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF (problem .eq. 11) then 
         CALL sixhump(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      ELSE IF (problem .eq. 12) then 
         CALL shubert(x, n, flag, f, iidata, 
     +        iisize, ddata, idsize, cdata, icsize)
      END IF
      END

      SUBROUTINE teconst(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)

      IMPLICIT NONE
      INTEGER n,flag
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 100
      flag = 0
      END

      SUBROUTINE telinear(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)

      IMPLICIT NONE
      INTEGER n,flag,i
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 0.D0
      f = 2.D0*x(1)
      DO 100, i = 2,n
         f = f + (x(i))*3.D0
100   CONTINUE
      flag = 0
      END

      SUBROUTINE tequadr(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)

      IMPLICIT NONE
      INTEGER n,flag,i
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 10
      DO 100, i = 1,n
         f = f + (x(i)-5.3)*(x(i)-5.3)
100   CONTINUE
      flag = 0
      END

      SUBROUTINE jonestest(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)

      IMPLICIT NONE
      INTEGER n,flag
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f, help1, help2, pi
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 0.D0
      help1 = 0.D0
      pi = acos(0.0D0)
      pi = pi + pi
      help1 = -sin(4*pi*x(1)) + 2*sin(2*pi*x(2))*sin(2*pi*x(2))
      IF (help1 .le. 0) THEN
         flag = 0
         help1 = x(1)*x(1)
         help2 = x(2)*x(2)
         f = (4.D0 - 2.1D0*help2 + (help1*help1)/3.D0)*help1
         f = f + x(1)*x(2) + (-4.D0 + 4.D0*help2)*help2
      ELSE
         flag = 1
      END IF
      END


      SUBROUTINE goldprice(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n,flag
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f, help1, help2, pi
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 0.D0
      help1 = 0.D0
      pi = acos(0.0D0)
      pi = pi + pi
      help1 = x(1)+x(2)+1.D0
      f = 19.D0-14.D0*(x(1)+x(2))+3*(x(1)*x(1)+x(2)*x(2))+6*x(1)*x(2)
      help1 = help1*help1*f+1
      help2 = 2.D0*x(1)-3.D0*x(2)
      f = 18-32*x(1)+12*x(1)*x(1) +48*x(2)-36*x(1)*x(2)+27*x(2)*x(2)
      help2 = help2*help2*f + 30.D0
      f = help1*help2
      flag = 0
      END


      SUBROUTINE branin(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n,flag
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f, help1, pi
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 0.D0
      help1 = 0.D0
      pi = acos(0.0D0)
      pi = pi + pi
      help1 = x(2) - 1.275D0*x(1)*x(1)/(pi*pi)
      help1 = help1 + 5*x(1)/pi - 6.D0
      f = help1*help1 + 10.D0 + (1 - 1/(8.D0*pi))*10.D0*cos(x(1))
      flag = 0
      END

      SUBROUTINE shekel(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n,flag,i,j,m
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f, help1, pi
      DOUBLE PRECISION a(10,4), c(10)
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      m = iidata(1)
      f = 0.D0
      help1 = 0.D0
      pi = acos(0.0D0)
      pi = pi + pi
      a(1,1) = 4.0D0
      a(1,2) = 4.0D0
      a(1,3) = 4.0D0
      a(1,4) = 4.0D0

      a(2,1) = 1.0D0
      a(2,2) = 1.0D0
      a(2,3) = 1.0D0
      a(2,4) = 1.0D0

      a(3,1) = 8.0D0
      a(3,2) = 8.0D0
      a(3,3) = 8.0D0
      a(3,4) = 8.0D0

      a(4,1) = 6.0D0
      a(4,2) = 6.0D0
      a(4,3) = 6.0D0
      a(4,4) = 6.0D0

      a(5,1) = 3.0D0
      a(5,2) = 7.0D0
      a(5,3) = 3.0D0
      a(5,4) = 7.0D0

      a(6,1) = 2.0D0
      a(6,2) = 9.0D0
      a(6,3) = 2.0D0
      a(6,4) = 9.0D0

      a(7,1) = 5.0D0
      a(7,2) = 5.0D0
      a(7,3) = 3.0D0
      a(7,4) = 3.0D0

      a(8,1) = 8.0D0
      a(8,2) = 1.0D0
      a(8,3) = 8.0D0
      a(8,4) = 1.0D0

      a(9,1) = 6.0D0
      a(9,2) = 2.0D0
      a(9,3) = 6.0D0
      a(9,4) = 2.0D0

      a(10,1) = 7.0D0
      a(10,2) = 3.6D0
      a(10,3) = 7.0D0
      a(10,4) = 3.6D0

      c(1)  = .1D0
      c(2)  = .2D0
      c(3)  = .2D0
      c(4)  = .4D0
      c(5)  = .4D0
      c(6)  = .6D0
      c(7)  = .3D0
      c(8)  = .7D0
      c(9)  = .5D0
      c(10) = .5D0
      f = 0.D0

      DO 10, i=1,m
        help1 = 0.D0
        DO 20, j = 1,n
          help1 = help1 + (x(j) - a(i,j))*(x(j) - a(i,j))
20      CONTINUE
        help1 = help1 + c(i)
        f = f - 1/help1
10    CONTINUE

      flag = 0
      END

      SUBROUTINE hartman(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n,flag,i,j,m
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f, help1, pi
      DOUBLE PRECISION a(4,6), c(4), p(4,6)
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)


      m = iidata(1)
      f = 0.D0
      help1 = 0.D0
      pi = acos(0.0D0)
      pi = pi + pi
      IF (n .eq. 3) THEN
        a(1,1) =  3.0D0
        a(1,2) = 10.0D0
        a(1,3) = 30.0D0

        a(2,1) =   .1D0
        a(2,2) = 10.0D0
        a(2,3) = 35.0D0

        a(3,1) =  3.0D0
        a(3,2) = 10.0D0
        a(3,3) = 30.0D0

        a(4,1) =  0.1D0
        a(4,2) = 10.0D0
        a(4,3) = 35.0D0

        c(1)  = 1.0D0
        c(2)  = 1.2D0
        c(3)  = 3.0D0
        c(4)  = 3.2D0

        p(1,1) =  0.3689D0
        p(1,2) =  0.1170D0
        p(1,3) =  0.2673D0

        p(2,1) =  0.4699D0
        p(2,2) =  0.4387D0
        p(2,3) =  0.7470D0

        p(3,1) =  0.1091D0
        p(3,2) =  0.8732D0
        p(3,3) =  0.5547D0

        p(4,1) =  0.03815D0
        p(4,2) =  0.5743D0
        p(4,3) =  0.8828D0
      ELSE
        a(1,1) = 10.0D0
        a(1,2) =  3.0D0
        a(1,3) = 17.0D0
        a(1,4) =  3.5D0
        a(1,5) =  1.7D0
        a(1,6) =  8.0D0

        a(2,1) =  0.05D0
        a(2,2) = 10.0D0
        a(2,3) = 17.0D0
        a(2,4) =  0.1D0
        a(2,5) =  8.0D0
        a(2,6) = 14.0D0

        a(3,1) =  3.0D0
        a(3,2) =  3.5D0
        a(3,3) =  1.7D0
        a(3,4) = 10.0D0
        a(3,5) = 17.0D0
        a(3,6) =  8.0D0

        a(4,1) = 17.0D0
        a(4,2) =  8.0D0
        a(4,3) =  0.05D0
        a(4,4) = 10.0D0
        a(4,5) =  0.1D0
        a(4,6) = 14.0D0

        c(1)  = 1.0D0
        c(2)  = 1.2D0
        c(3)  = 3.0D0
        c(4)  = 3.2D0

        p(1,1) = .1312D0
        p(1,2) = .1696D0
        p(1,3) = .5569D0
        p(1,4) = .0124D0
        p(1,5) = .8283D0
        p(1,6) = .5886D0

        p(2,1) = .2329D0
        p(2,2) = .4135D0
        p(2,3) = .8307D0
        p(2,4) = .3736D0
        p(2,5) = .1004D0
        p(2,6) = .9991D0

        p(3,1) = .2348D0
        p(3,2) = .1451D0
        p(3,3) = .3522D0
        p(3,4) = .2883D0
        p(3,5) = .3047D0
        p(3,6) = .6650D0

        p(4,1) = .4047D0
        p(4,2) = .8828D0
        p(4,3) = .8732D0
        p(4,4) = .5743D0
        p(4,5) = .1091D0
        p(4,6) = .0381D0

      END IF
      f = 0.D0

      DO 10, i=1,m
        help1 = 0.D0
        DO 20, j = 1,n
          help1 = help1 + a(i,j)*(x(j) - p(i,j))*(x(j) - p(i,j))
20      CONTINUE
        help1 = c(i)*dexp(-help1)
        f = f - help1
10    CONTINUE

      flag = 0
      END


      SUBROUTINE sixhump(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n,flag
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f, help1, help2, pi
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 0.D0
      help1 = 0.D0
      pi = acos(0.0D0)
      pi = pi + pi
      help1 = 4.D0 - 2.1D0*x(1)*x(1)+x(1)*x(1)*x(1)*x(1)/3.D0
      help2 = -4.D0 + 4.D0*x(2)*x(2)
      f = help1*x(1)*x(1) + x(1)*x(2) + help2*x(2)*x(2)
      flag = 0
      END

      SUBROUTINE Shubert(x,n,flag,f,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n,flag,i
      DOUBLE PRECISION x(n)
      DOUBLE PRECISION f, help1, help2
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)

      f = 0.D0
      help1 = 0.D0
      help2 = 0.D0
      DO 10, i = 1,5
        help1 = help1+real(i)*dcos((real(i)+1.D0)*x(1) + real(i))
        help2 = help2+real(i)*dcos((real(i)+1.D0)*x(2) + real(i))
10    CONTINUE
      f = help1*help2
      flag = 0
      END
