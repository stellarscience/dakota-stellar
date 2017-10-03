/* do7prd.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int do7prd_(l, ls, p, s, w, y, z__)
integer *l, *ls, *p;
doublereal *s, *w, *y, *z__;
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer y_dim1, y_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal wk, yi;


/*  ***  FOR I = 1..L, SET S = S + W(I)*Y(.,I)*(Z(.,I)**T), I.E., */
/*  ***        ADD W(I) TIMES THE OUTER PRODUCT OF Y(.,I) AND Z(.,I). */

/*     DIMENSION S(P*(P+1)/2) */

    /* Parameter adjustments */
    --w;
    --s;
    z_dim1 = *p;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    y_dim1 = *p;
    y_offset = y_dim1 + 1;
    y -= y_offset;

    /* Function Body */

    i__1 = *l;
    for (k = 1; k <= i__1; ++k) {
	wk = w[k];
	if (wk == zero) {
	    goto L30;
	}
	m = 1;
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    yi = wk * y[i__ + k * y_dim1];
	    i__3 = i__;
	    for (j = 1; j <= i__3; ++j) {
		s[m] += yi * z__[j + k * z_dim1];
		++m;
/* L10: */
	    }
/* L20: */
	}
L30:
	;
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF DO7PRD FOLLOWS  *** */
} /* do7prd_ */

