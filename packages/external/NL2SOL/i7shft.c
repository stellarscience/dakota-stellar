/* i7shft.f -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int i7shft_(integer *n, integer *k, integer *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, t, k1, ii, nm1;


/*  ***  SHIFT X(K),...,X(N) LEFT CIRCULARLY ONE POSITION IF K .GT. 0. */
/*  ***  SHIFT X(-K),...,X(N) RIGHT CIRCULARLY ONE POSITION IF K .LT. 0. */



    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*k < 0) {
	goto L20;
    }
    if (*k >= *n) {
	goto L999;
    }
    nm1 = *n - 1;
    t = x[*k];
    i__1 = nm1;
    for (i__ = *k; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = x[i__ + 1];
    }
    x[*n] = t;
    goto L999;

L20:
    k1 = -(*k);
    if (k1 >= *n) {
	goto L999;
    }
    t = x[*n];
    nm1 = *n - k1;
    i__1 = nm1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n - ii;
	x[i__ + 1] = x[i__];
/* L30: */
    }
    x[k1] = t;
L999:
    return 0;
/*  ***  LAST LINE OF I7SHFT FOLLOWS  *** */
} /* i7shft_ */

