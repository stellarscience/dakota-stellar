/* dl7srt.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dl7srt_(n1, n, l, a, irc)
integer *n1, *n;
doublereal *l, *a;
integer *irc;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer i0, j0, ij, ik, jk;
    static doublereal td;
    static integer im1, jm1;


/*  ***  COMPUTE ROWS N1 THROUGH N OF THE CHOLESKY FACTOR  L  OF */
/*  ***  A = L*(L**T),  WHERE  L  AND THE LOWER TRIANGLE OF  A  ARE BOTH 
*/
/*  ***  STORED COMPACTLY BY ROWS (AND MAY OCCUPY THE SAME STORAGE). */
/*  ***  IRC = 0 MEANS ALL WENT WELL.  IRC = J MEANS THE LEADING */
/*  ***  PRINCIPAL  J X J  SUBMATRIX OF  A  IS NOT POSITIVE DEFINITE -- */
/*  ***  AND  L(J*(J+1)/2)  CONTAINS THE (NONPOS.) REDUCED J-TH DIAGONAL. 
*/

/*  ***  PARAMETERS  *** */

/*     DIMENSION L(N*(N+1)/2), A(N*(N+1)/2) */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

/*  ***  BODY  *** */

    /* Parameter adjustments */
    --a;
    --l;

    /* Function Body */
    i0 = *n1 * (*n1 - 1) / 2;
    i__1 = *n;
    for (i__ = *n1; i__ <= i__1; ++i__) {
	td = 0.;
	if (i__ == 1) {
	    goto L40;
	}
	j0 = 0;
	im1 = i__ - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    t = 0.;
	    if (j == 1) {
		goto L20;
	    }
	    jm1 = j - 1;
	    i__3 = jm1;
	    for (k = 1; k <= i__3; ++k) {
		ik = i0 + k;
		jk = j0 + k;
		t += l[ik] * l[jk];
/* L10: */
	    }
L20:
	    ij = i0 + j;
	    j0 += j;
	    t = (a[ij] - t) / l[j0];
	    l[ij] = t;
	    td += t * t;
/* L30: */
	}
L40:
	i0 += i__;
	t = a[i0] - td;
	if (t <= 0.) {
	    goto L60;
	}
	l[i0] = sqrt(t);
/* L50: */
    }

    *irc = 0;
    goto L999;

L60:
    l[i0] = t;
    *irc = i__;

L999:
    return 0;

/*  ***  LAST CARD OF DL7SRT  *** */
} /* dl7srt_ */

