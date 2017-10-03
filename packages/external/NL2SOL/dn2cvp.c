/* dn2cvp.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <stdio.h>
#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dn2cvp_(integer *iv, integer *liv, integer *lv, integer *
	p, doublereal *v)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer i1, ii, pu, cov1;


/*  ***  PRINT COVARIANCE MATRIX FOR  DRN2G  *** */


/*  ***  LOCAL VARIABLES  *** */


/*     ***  IV SUBSCRIPTS  *** */


/* /6 */
/*     DATA COVMAT/26/, COVPRT/14/, COVREQ/15/, NEEDHD/36/, NFCOV/52/, */
/*    1     NGCOV/53/, PRUNIT/21/, RCOND/53/, REGD/67/, STATPR/23/ */
/* /7 */
/* / */
/*  ***  BODY  *** */

    /* Parameter adjustments */
    --iv;
    --v;

    /* Function Body */
    if (iv[1] > 8) {
	goto L999;
    }
    pu = iv[21];
    if (pu == 0) {
	goto L999;
    }
    if (iv[23] == 0) {
	goto L30;
    }
    if (iv[52] > 0)
	printf("\n%3d EXTRA FUNC. EVALS FOR COVARIANCE AND DIAGNOSTICS.\n", iv[52]);
    if (iv[53] > 0)
	printf("%3d EXTRA GRAD. EVALS FOR COVARIANCE AND DIAGNOSTICS.\n", iv[53]);

L30:
    if (iv[14] <= 0) {
	goto L999;
    }
    cov1 = iv[26];
    if (iv[67] <= 0 && cov1 <= 0) {
	goto L70;
    }
    iv[36] = 1;
/* Computing 2nd power */
    d__1 = v[53];
    t = d__1 * d__1;
    if (abs(iv[15]) > 2) {
	goto L50;
    }

	printf("\nRECIPROCAL CONDITION OF F.D. HESSIAN = AT MOST %# -9.2g\n", t);
    goto L70;

L50:
	printf("\nRECIPROCAL CONDITION OF (J**T)*J = AT LEAST %# -9.2g\n", t);

L70:
    if (iv[14] % 2 == 0) {
	goto L999;
    }
    iv[36] = 1;
    if (cov1 < 0) {
	goto L80;
    } else if (cov1 == 0) {
	goto L110;
    } else {
	goto L130;
    }
L80:
    if (-1 == cov1)
	printf("\n++++++ INDEFINITE COVARIANCE MATRIX ++++++\n");
    if (-2 == cov1)
	printf("\n++++++ OVERSIZE STEPS IN COMPUTING COVARIANCE +++++\n");
    goto L999;

L110:
	printf("\n++++++ COVARIANCE MATRIX NOT COMPUTED ++++++\n");
    goto L999;

L130:
    i__ = abs(iv[15]);
    if (i__ <= 1)
	printf("\nCOVARIANCE = SCALE * H**-1 * (J**T * J) * H**-1\nWHERE H = F.D. HESSIAN\n\n");
    if (i__ == 2)
	printf("\nCOVARIANCE = H**-1, WHERE H = FINITE-DIFFERENCE HESSIAN\n\n");
    if (i__ > 2)
	printf("\nCOVARIANCE = SCALE * (J**T * J)**-1\n\n");
    ii = cov1 - 1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i1 = ii + 1;
	ii += i__;
	printf("ROW %2d    %# -12.3g", i__, v[i1]);
	while(++j <= ii)
		printf((j - i1) % 5 ? " %# -11.3g" : "\n          %# -12.3g", v[j]);
	printf("\n");
/* L170: */
    }

L999:
    return 0;
/*  ***  LAST CARD OF DN2CVP FOLLOWS  *** */
} /* dn2cvp_ */

