/* dv2axy.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dv2axy_(p, w, a, x, y)
integer *p;
doublereal *w, *a, *x, *y;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET W = A*X + Y  --  W, X, Y = P-VECTORS, A = SCALAR  *** */



    /* Parameter adjustments */
    --y;
    --x;
    --w;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	w[i__] = *a * x[i__] + y[i__];
    }
    return 0;
} /* dv2axy_ */

