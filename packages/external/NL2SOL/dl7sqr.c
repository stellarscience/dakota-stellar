/* dl7sqr.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dl7sqr_(n, a, l)
integer *n;
doublereal *a, *l;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer i0, j0, ii, ij, ik, jj, jk, ip1, np1;


/*  ***  COMPUTE  A = LOWER TRIANGLE OF  L*(L**T),  WITH BOTH */
/*  ***  L  AND  A  STORED COMPACTLY BY ROWS.  (BOTH MAY OCCUPY THE */
/*  ***  SAME STORAGE. */

/*  ***  PARAMETERS  *** */

/*     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2) */

/*  ***  LOCAL VARIABLES  *** */


    /* Parameter adjustments */
    --l;
    --a;

    /* Function Body */
    np1 = *n + 1;
    i0 = *n * (*n + 1) / 2;
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = np1 - ii;
	ip1 = i__ + 1;
	i0 -= i__;
	j0 = i__ * (i__ + 1) / 2;
	i__2 = i__;
	for (jj = 1; jj <= i__2; ++jj) {
	    j = ip1 - jj;
	    j0 -= j;
	    t = 0.;
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		ik = i0 + k;
		jk = j0 + k;
		t += l[ik] * l[jk];
/* L10: */
	    }
	    ij = i0 + j;
	    a[ij] = t;
/* L20: */
	}
/* L30: */
    }
/* L999: */
    return 0;
} /* dl7sqr_ */

