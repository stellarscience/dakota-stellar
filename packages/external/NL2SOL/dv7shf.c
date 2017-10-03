/* dv7shf.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dv7shf_(n, k, x)
integer *n, *k;
doublereal *x;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t;
    static integer nm1;


/*  ***  SHIFT X(K),...,X(N) LEFT CIRCULARLY ONE POSITION  *** */



    /* Parameter adjustments */
    --x;

    /* Function Body */
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
L999:
    return 0;
} /* dv7shf_ */

