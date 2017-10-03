/* dd7mlp.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dd7mlp_(n, x, y, z__, k)
integer *n;
doublereal *x, *y, *z__;
integer *k;
{
    /* Initialized data */

    static doublereal one = 1.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal t;


/* ***  SET X = DIAG(Y)**K * Z */
/* ***  FOR X, Z = LOWER TRIANG. MATRICES STORED COMPACTLY BY ROW */
/* ***  K = 1 OR -1. */

/* /6S */
/*     DOUBLE PRECISION X(1), Y(N), Z(1) */
/* /7S */
/* / */
    /* Parameter adjustments */
    --y;
    --x;
    --z__;

    /* Function Body */

    l = 1;
    if (*k >= 0) {
	goto L30;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = one / y[i__];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    x[l] = t * z__[l];
	    ++l;
/* L10: */
	}
/* L20: */
    }
    goto L999;

L30:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = y[i__];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    x[l] = t * z__[l];
	    ++l;
/* L40: */
	}
/* L50: */
    }
L999:
    return 0;
/*  ***  LAST CARD OF DD7MLP FOLLOWS  *** */
} /* dd7mlp_ */

