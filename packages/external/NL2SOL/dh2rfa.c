/* dh2rfa.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dh2rfa_(n, a, b, x, y, z__)
integer *n;
doublereal *a, *b, *x, *y, *z__;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t;


/*  ***  APPLY 2X2 HOUSEHOLDER REFLECTION DETERMINED BY X, Y, Z TO */
/*  ***  N-VECTORS A, B  *** */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = a[i__] * *x + b[i__] * *y;
	a[i__] += t;
	b[i__] += t * *z__;
/* L10: */
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DH2RFA FOLLOWS  *** */
} /* dh2rfa_ */

