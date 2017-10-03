/* dc7vfn.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dc7vfn_(iv, l, lh, liv, lv, n, p, v)
integer *iv;
doublereal *l;
integer *lh, *liv, *lv, *n, *p;
doublereal *v;
{
    /* Initialized data */

    static doublereal half = .5;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int dv7scl_(), dl7nvr_(), dl7tsq_();
    static integer cov;


/*  ***  FINISH COVARIANCE COMPUTATION FOR  DRN2G,  DRNSG  *** */



/*  ***  LOCAL VARIABLES  *** */


/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/* /6 */
/*     DATA CNVCOD/55/, COVMAT/26/, F/10/, FDH/74/, H/56/, MODE/35/, */
/*    1     RDREQ/57/, REGD/67/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --l;
    --iv;
    --v;

    /* Function Body */

/*  ***  BODY  *** */

    iv[1] = iv[55];
    i__ = iv[35] - *p;
    iv[35] = 0;
    iv[55] = 0;
    if (iv[74] <= 0) {
	goto L999;
    }
/* Computing 2nd power */
    i__1 = i__ - 2;
    if (i__1 * i__1 == 1) {
	iv[67] = 1;
    }
    if (iv[57] % 2 != 1) {
	goto L999;
    }

/*     ***  FINISH COMPUTING COVARIANCE MATRIX = INVERSE OF F.D. HESSIAN. 
*/

    cov = abs(iv[56]);
    iv[74] = 0;

    if (iv[26] != 0) {
	goto L999;
    }
    if (i__ >= 2) {
	goto L10;
    }
    dl7nvr_(p, &v[cov], &l[1]);
    dl7tsq_(p, &v[cov], &v[cov]);

L10:
/* Computing MAX */
    i__1 = 1, i__2 = *n - *p;
    d__1 = v[10] / (half * (real) max(i__1,i__2));
    dv7scl_(lh, &v[cov], &d__1, &v[cov]);
    iv[26] = cov;

L999:
    return 0;
/*  ***  LAST LINE OF DC7VFN FOLLOWS  *** */
} /* dc7vfn_ */

