/* ds7lvm.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int ds7lvm_(p, y, s, x)
integer *p;
doublereal *y, *s, *x;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    extern doublereal dd7tpr_();
    static doublereal xi;
    static integer im1;


/*  ***  SET  Y = S * X,  S = P X P SYMMETRIC MATRIX.  *** */
/*  ***  LOWER TRIANGLE OF  S  STORED ROWWISE.         *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION S(P*(P+1)/2) */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  NO INTRINSIC FUNCTIONS  *** */

/*  ***  EXTERNAL FUNCTION  *** */


/* -----------------------------------------------------------------------
 */

    /* Parameter adjustments */
    --x;
    --y;
    --s;

    /* Function Body */
    j = 1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = dd7tpr_(&i__, &s[j], &x[1]);
	j += i__;
/* L10: */
    }

    if (*p <= 1) {
	goto L999;
    }
    j = 1;
    i__1 = *p;
    for (i__ = 2; i__ <= i__1; ++i__) {
	xi = x[i__];
	im1 = i__ - 1;
	++j;
	i__2 = im1;
	for (k = 1; k <= i__2; ++k) {
	    y[k] += s[j] * xi;
	    ++j;
/* L30: */
	}
/* L40: */
    }

L999:
    return 0;
/*  ***  LAST CARD OF DS7LVM FOLLOWS  *** */
} /* ds7lvm_ */

