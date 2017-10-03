/* dr7tvm.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dr7tvm_(n, p, y, d__, u, x)
integer *n, *p;
doublereal *y, *d__, *u, *x;
{
    /* System generated locals */
    integer u_dim1, u_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal t;
    extern doublereal dd7tpr_();
    static integer ii, pl, pp1;


/*  ***  SET Y TO R*X, WHERE R IS THE UPPER TRIANGULAR MATRIX WHOSE */
/*  ***  DIAGONAL IS IN D AND WHOSE STRICT UPPER TRIANGLE IS IN U. */

/*  ***  X AND Y MAY SHARE STORAGE. */



/*  ***  LOCAL VARIABLES  *** */


/*  ***  BODY  *** */

    /* Parameter adjustments */
    --x;
    u_dim1 = *n;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --d__;
    --y;

    /* Function Body */
/* Computing MIN */
    pl = min(*n,*p);
    pp1 = pl + 1;
    i__1 = pl;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = pp1 - ii;
	t = x[i__] * d__[i__];
	if (i__ > 1) {
	    i__2 = i__ - 1;
	    t += dd7tpr_(&i__2, &u[i__ * u_dim1 + 1], &x[1]);
	}
	y[i__] = t;
/* L10: */
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DR7TVM FOLLOWS  *** */
} /* dr7tvm_ */

