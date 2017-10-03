/* i7pnvr.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int i7pnvr_(n, x, y)
integer *n, *x, *y;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;


/*  ***  SET PERMUTATION VECTOR X TO INVERSE OF Y  *** */


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = y[i__];
	x[j] = i__;
/* L10: */
    }

/* L999: */
    return 0;
/*  ***  LAST LINE OF I7PNVR FOLLOWS  *** */
} /* i7pnvr_ */

