/* dh2rfg.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

doublereal dh2rfg_(a, b, x, y, z__)
doublereal *a, *b, *x, *y, *z__;
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal c__, t, a1, b1;


/*  ***  DETERMINE X, Y, Z SO  I + (1,Z)**T * (X,Y)  IS A 2X2 */
/*  ***  HOUSEHOLDER REFLECTION SENDING (A,B)**T INTO (C,0)**T, */
/*  ***  WHERE  C = -SIGN(A)*SQRT(A**2 + B**2)  IS THE VALUE DH2RFG */
/*  ***  RETURNS. */


/* /+ */
/* / */

/*  ***  BODY  *** */

    if (*b != zero) {
	goto L10;
    }
    *x = zero;
    *y = zero;
    *z__ = zero;
    ret_val = *a;
    goto L999;
L10:
    t = abs(*a) + abs(*b);
    a1 = *a / t;
    b1 = *b / t;
/* Computing 2nd power */
    d__1 = a1;
/* Computing 2nd power */
    d__2 = b1;
    c__ = sqrt(d__1 * d__1 + d__2 * d__2);
    if (a1 > zero) {
	c__ = -c__;
    }
    a1 -= c__;
    *z__ = b1 / a1;
    *x = a1 / c__;
    *y = b1 / c__;
    ret_val = t * c__;
L999:
    return ret_val;
/*  ***  LAST LINE OF DH2RFG FOLLOWS  *** */
} /* dh2rfg_ */

