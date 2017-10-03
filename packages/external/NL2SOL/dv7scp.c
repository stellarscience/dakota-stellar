/* dv7scp.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dv7scp_(p, y, s)
integer *p;
doublereal *y, *s;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*  ***  SET P-VECTOR Y TO SCALAR S  *** */



    /* Parameter adjustments */
    --y;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y[i__] = *s;
    }
    return 0;
} /* dv7scp_ */

