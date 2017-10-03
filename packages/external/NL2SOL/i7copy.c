/* i7copy.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int i7copy_(p, y, x)
integer *p, *y, *x;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET Y = X, WHERE X AND Y ARE INTEGER P-VECTORS  *** */



    /* Parameter adjustments */
    --x;
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = x[i__];
    }
/* L999: */
    return 0;
} /* i7copy_ */

