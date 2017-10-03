/* drldst.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

doublereal drldst_(p, d__, x, x0)
integer *p;
doublereal *d__, *x, *x0;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static doublereal emax, xmax;
    static integer i__;
    static doublereal t;


/*  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  *** */
/*  ***  NL2SOL VERSION 2.2  *** */


/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */

/*  ***  BODY  *** */

    /* Parameter adjustments */
    --x0;
    --x;
    --d__;

    /* Function Body */
    emax = 0.;
    xmax = 0.;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = (d__1 = d__[i__] * (x[i__] - x0[i__]), abs(d__1));
	if (emax < t) {
	    emax = t;
	}
	t = d__[i__] * ((d__1 = x[i__], abs(d__1)) + (d__2 = x0[i__], abs(
		d__2)));
	if (xmax < t) {
	    xmax = t;
	}
/* L10: */
    }
    ret_val = 0.;
    if (xmax > 0.) {
	ret_val = emax / xmax;
    }
/* L999: */
    return ret_val;
/*  ***  LAST CARD OF DRLDST FOLLOWS  *** */
} /* drldst_ */

