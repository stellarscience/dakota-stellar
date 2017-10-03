/* dl7tsq.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dl7tsq_(n, a, l)
integer *n;
doublereal *a, *l;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, m, i1, ii;
    static doublereal lj, lii;
    static integer iim1;


/*  ***  SET A TO LOWER TRIANGLE OF (L**T) * L  *** */

/*  ***  L = N X N LOWER TRIANG. MATRIX STORED ROWWISE.  *** */
/*  ***  A IS ALSO STORED ROWWISE AND MAY SHARE STORAGE WITH L.  *** */

/*     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2) */


    /* Parameter adjustments */
    --l;
    --a;

    /* Function Body */
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ii + 1;
	ii += i__;
	m = 1;
	if (i__ == 1) {
	    goto L30;
	}
	iim1 = ii - 1;
	i__2 = iim1;
	for (j = i1; j <= i__2; ++j) {
	    lj = l[j];
	    i__3 = j;
	    for (k = i1; k <= i__3; ++k) {
		a[m] += lj * l[k];
		++m;
/* L10: */
	    }
/* L20: */
	}
L30:
	lii = l[ii];
	i__2 = ii;
	for (j = i1; j <= i__2; ++j) {
/* L40: */
	    a[j] = lii * l[j];
	}
/* L50: */
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF DL7TSQ FOLLOWS  *** */
} /* dl7tsq_ */

