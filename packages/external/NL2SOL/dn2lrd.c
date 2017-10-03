/* dn2lrd.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b4 = -1.;
static integer c__1 = 1;

/* Subroutine */ int dn2lrd_(dr, iv, l, lh, liv, lv, nd, nn, p, r__, rd, v)
doublereal *dr;
integer *iv;
doublereal *l;
integer *lh, *liv, *lv, *nd, *nn, *p;
doublereal *r__, *rd, *v;
{
    /* Initialized data */

    static doublereal onev[1] = { 1. };

    /* System generated locals */
    integer dr_dim1, dr_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer step1;
    static doublereal a;
    static integer i__, j, m;
    static doublereal s, t;
    extern /* Subroutine */ int do7prd_(), dl7ivm_();
    extern doublereal dd7tpr_();
    extern /* Subroutine */ int dv7scp_(), dl7itv_();
    static doublereal ff;
    static integer cov;


/*  ***  COMPUTE REGRESSION DIAGNOSTIC AND DEFAULT COVARIANCE MATRIX FOR 
*/
/*        DRN2G  *** */

/*  ***  PARAMETERS  *** */


/*  ***  CODED BY DAVID M. GAY (WINTER 1982, FALL 1983)  *** */

/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  CONSTANTS  *** */


/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */

/*  ***  IV AND V SUBSCRIPTS  *** */

/* /6 */
/*     DATA F/10/, H/56/, MODE/35/, RDREQ/57/, STEP/40/ */
/* /7 */
/* / */
/* /6 */
/*     DATA NEGONE/-1.D+0/, ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --l;
    --iv;
    --v;
    --rd;
    --r__;
    dr_dim1 = *nd;
    dr_offset = dr_dim1 + 1;
    dr -= dr_offset;

    /* Function Body */

/* ++++++++++++++++++++++++++++++++  BODY  +++++++++++++++++++++++++++++++
 */

    step1 = iv[40];
    i__ = iv[57];
    if (i__ <= 0) {
	goto L999;
    }
    if (i__ % 4 < 2) {
	goto L30;
    }
    ff = 1.;
    if (v[10] != 0.) {
	ff = 1. / sqrt((abs(v[10])));
    }
    dv7scp_(nn, &rd[1], &c_b4);
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = r__[i__];
	a = d__1 * d__1;
	m = step1;
	i__2 = *p;
	for (j = 1; j <= i__2; ++j) {
	    v[m] = dr[i__ + j * dr_dim1];
	    ++m;
/* L10: */
	}
	dl7ivm_(p, &v[step1], &l[1], &v[step1]);
	s = dd7tpr_(p, &v[step1], &v[step1]);
	t = 1. - s;
	if (t <= 0.) {
	    goto L20;
	}
	a = a * s / t;
	rd[i__] = sqrt(a) * ff;
L20:
	;
    }

L30:
    if (iv[35] - *p < 2) {
	goto L999;
    }

/*  ***  COMPUTE DEFAULT COVARIANCE MATRIX  *** */

    cov = abs(iv[56]);
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = step1;
	i__2 = *p;
	for (j = 1; j <= i__2; ++j) {
	    v[m] = dr[i__ + j * dr_dim1];
	    ++m;
/* L40: */
	}
	dl7ivm_(p, &v[step1], &l[1], &v[step1]);
	dl7itv_(p, &v[step1], &l[1], &v[step1]);
	do7prd_(&c__1, lh, p, &v[cov], onev, &v[step1], &v[step1]);
/* L50: */
    }

L999:
    return 0;
/*  ***  LAST LINE OF DN2LRD FOLLOWS  *** */
} /* dn2lrd_ */

