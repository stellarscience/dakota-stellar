/* dv7cpy.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dv7cpy_(p, y, x)
integer *p;
doublereal *y, *x;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET Y = X, WHERE X AND Y ARE P-VECTORS  *** */



    /* Parameter adjustments */
    --x;
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = x[i__];
    }
    return 0;
} /* dv7cpy_ */

