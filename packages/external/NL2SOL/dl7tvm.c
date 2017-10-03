/* dl7tvm.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dl7tvm_(n, x, l, y)
integer *n;
doublereal *x, *l, *y;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, i0, ij;
    static doublereal yi;


/*  ***  COMPUTE  X = (L**T)*Y, WHERE  L  IS AN  N X N  LOWER */
/*  ***  TRIANGULAR MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY */
/*  ***  OCCUPY THE SAME STORAGE.  *** */

/*     DIMENSION L(N*(N+1)/2) */
/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --y;
    --x;
    --l;

    /* Function Body */
    i0 = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yi = y[i__];
	x[i__] = 0.;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ij = i0 + j;
	    x[j] += yi * l[ij];
/* L10: */
	}
	i0 += i__;
/* L20: */
    }
/* L999: */
    return 0;
/*  ***  LAST CARD OF DL7TVM FOLLOWS  *** */
} /* dl7tvm_ */

