#ifndef __CBLAS_H__
#define __CBLAS_H__

#ifdef HAVE_CONFIG_H
  #include "OPT++_config.h"
  /* BLAS TYPE DEFINITIONS  */
  #define daxpy F77_FUNC(daxpy, DAXPY)
  #define dswap F77_FUNC(dswap, DSWAP)
  #define ddot F77_FUNC(ddot, DDOT)
  #define dscal F77_FUNC(dscal, DSCAL)
  #define dnrm2 F77_FUNC(dnrm2, DNRM2)
  #define dcopy F77_FUNC(dcopy, DCOPY)
#else
  #include "cblas_config.h"
  /* BLAS TYPE DEFINITIONS  */
  #define daxpy OPTPP_GLOBAL(daxpy, DAXPY)
  #define dswap OPTPP_GLOBAL(dswap, DSWAP)
  #define ddot OPTPP_GLOBAL(ddot, DDOT)
  #define dscal OPTPP_GLOBAL(dscal, DSCAL)
  #define dnrm2 OPTPP_GLOBAL(dnrm2, DNRM2)
  #define dcopy OPTPP_GLOBAL(dcopy, DCOPY)
#endif

#ifdef __cplusplus
extern "C" {
#endif


typedef int  Integer;

typedef struct { 
            float real;
            float imag; 
} Complex;

typedef struct {
	    double real;
            double imag; 
} Zomplex;

typedef enum { NoTranspose,
               Transpose,
               ConjugateTranspose } MatrixTranspose;

typedef enum { UpperTriangle,
               LowerTriangle } MatrixTriangle;

typedef enum { UnitTriangular,
               NotUnitTriangular } MatrixUnitTriangular;

typedef enum { LeftSide,
               RightSide } OperationSide;


/********                                                             ********
 ********                    LEVEL 1 BLAS                             ********
 ********                                                             ********/
/*                             Swap two vectors                              *
 *                                x <-> y                                    */

extern void dswap( Integer *n,  double *x, Integer *incx,  double *y, 
		   Integer *incy );

/*                              Scale a vector                               *
 *                               x <- alpha*x                                */

extern void dscal( Integer *n,  double *alpha,  double *x, Integer *incx );

/*                         Copy one vector to another                        *
 *                                  y <- x                                   */

extern void dcopy( Integer *n,  double *x, Integer *incx,  double *y, 
		   Integer *incy );

/*                 Scale a vector then add to another vector                 *
 *                             y <- alpha*x + y                              */

extern void daxpy( Integer *n,  double *alpha,  double *x, Integer *incx, 
		   double *y, Integer *incy );

/*                         Dot product of two vectors                        * 
 *                                 dot <- xTy                                */

extern double  ddot( Integer *n, double *x, Integer *incx,
		     double *y, Integer *incy );

/*                            2-Norm of a vector                             *
 *                              nrm2 <- ||x||2                               */

extern double  dnrm2( Integer *n,  double *x, Integer *incx );

#ifdef __cplusplus
}
#endif


#endif /* !__CBLAS_H__ */

