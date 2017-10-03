/* dg7qsb.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c_n1 = -1;

/* Subroutine */ int dg7qsb_(b, d__, dihdi, g, ipiv, ipiv1, ipiv2, ka, l, lv, 
	p, p0, pc, step, td, tg, v, w, x, x0)
doublereal *b, *d__, *dihdi, *g;
integer *ipiv, *ipiv1, *ipiv2, *ka;
doublereal *l;
integer *lv, *p, *p0, *pc;
doublereal *step, *td, *tg, *v, *w, *x, *x0;
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer step_dim1, step_offset;

    /* Local variables */
    static doublereal nred, pred;
    static integer k, kinit, p1;
    extern /* Subroutine */ int ds7bqn_();
    extern doublereal dd7tpr_();
    extern /* Subroutine */ int dv7scp_(), ds7ipr_(), dg7qts_(), dv7ipr_(), 
	    dv7cpy_();
    static integer kb, p10;
    extern /* Subroutine */ int dv7vmp_();
    static integer ns;
    static doublereal ds0, rad;


/*  ***  COMPUTE HEURISTIC BOUNDED NEWTON STEP  *** */

/*     DIMENSION DIHDI(P*(P+1)/2), L(P*(P+1)/2) */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  V SUBSCRIPTS  *** */


/* /6 */
/*     DATA DST0/3/, DSTNRM/2/, GTSTEP/4/, NREDUC/6/, PREDUC/7/, */
/*    1     RADIUS/8/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --dihdi;
    --l;
    --v;
    --x0;
    --x;
    --w;
    --tg;
    --td;
    step_dim1 = *p;
    step_offset = step_dim1 + 1;
    step -= step_offset;
    --ipiv2;
    --ipiv1;
    --ipiv;
    --g;
    --d__;
    b -= 3;

    /* Function Body */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
 */

    p1 = *pc;
    if (*ka < 0) {
	goto L10;
    }
    nred = v[6];
    ds0 = v[3];
    goto L20;
L10:
    *p0 = 0;
    *ka = -1;

L20:
    kinit = -1;
    if (*p0 == p1) {
	kinit = *ka;
    }
    dv7cpy_(p, &x[1], &x0[1]);
    pred = zero;
    rad = v[8];
    kb = -1;
    v[2] = zero;
    if (p1 > 0) {
	goto L30;
    }
    nred = zero;
    ds0 = zero;
    dv7scp_(p, &step[step_offset], &zero);
    goto L60;

L30:
    dv7cpy_(p, &td[1], &d__[1]);
    dv7ipr_(p, &ipiv[1], &td[1]);
    dv7vmp_(p, &tg[1], &g[1], &d__[1], &c_n1);
    dv7ipr_(p, &ipiv[1], &tg[1]);
L40:
    k = kinit;
    kinit = -1;
    v[8] = rad - v[2];
    dg7qts_(&td[1], &tg[1], &dihdi[1], &k, &l[1], &p1, &step[step_offset], &v[
	    1], &w[1]);
    *p0 = p1;
    if (*ka >= 0) {
	goto L50;
    }
    nred = v[6];
    ds0 = v[3];

L50:
    *ka = k;
    v[8] = rad;
    p10 = p1;
    ds7bqn_(&b[3], &d__[1], &step[(step_dim1 << 1) + 1], &ipiv[1], &ipiv1[1], 
	    &ipiv2[1], &kb, &l[1], lv, &ns, p, &p1, &step[step_offset], &td[1]
	    , &tg[1], &v[1], &w[1], &x[1], &x0[1]);
    if (ns > 0) {
	ds7ipr_(&p10, &ipiv1[1], &dihdi[1]);
    }
    pred += v[7];
    if (ns != 0) {
	*p0 = 0;
    }
    if (kb <= 0) {
	goto L40;
    }

L60:
    v[3] = ds0;
    v[6] = nred;
    v[7] = pred;
    v[4] = dd7tpr_(p, &g[1], &step[step_offset]);

/* L999: */
    return 0;
/*  ***  LAST LINE OF DG7QSB FOLLOWS  *** */
} /* dg7qsb_ */

