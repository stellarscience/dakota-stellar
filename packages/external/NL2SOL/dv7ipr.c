/* dv7ipr.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dv7ipr_(n, ip, x)
integer *n, *ip;
doublereal *x;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;


/*     PERMUTE X SO THAT X.OUTPUT(I) = X.INPUT(IP(I)). */
/*     IP IS UNCHANGED ON OUTPUT. */


    /* Parameter adjustments */
    --x;
    --ip;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ip[i__];
	if (j == i__) {
	    goto L30;
	}
	if (j > 0) {
	    goto L10;
	}
	ip[i__] = -j;
	goto L30;
L10:
	t = x[i__];
	k = i__;
L20:
	x[k] = x[j];
	k = j;
	j = ip[k];
	ip[k] = -j;
	if (j > i__) {
	    goto L20;
	}
	x[k] = t;
L30:
	;
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DV7IPR FOLLOWS  *** */
} /* dv7ipr_ */

