/* dv7scl.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dv7scl_(n, x, a, y)
integer *n;
doublereal *x, *a, *y;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET X(I) = A*Y(I), I = 1(1)N  *** */



    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = *a * y[i__];
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DV7SCL FOLLOWS  *** */
} /* dv7scl_ */

