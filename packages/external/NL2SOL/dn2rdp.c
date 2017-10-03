/* dn2rdp.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <stdio.h>
#include "f2c.h"

/* Subroutine */ int dn2rdp_(integer *iv, integer *liv, integer *lv, integer *
	n, doublereal *rd, doublereal *v)
{

    /* System generated locals */
    integer rd_dim1, i__1;

    /* Local variables */
    static integer pu;

/*  ***  PRINT REGRESSION DIAGNOSTICS FOR MLPSL AND NL2S1 *** */


/*     ***  NOTE -- V IS PASSED FOR POSSIBLE USE BY REVISED VERSIONS OF */
/*     ***  THIS ROUTINE. */


/*  ***  IV AND V SUBSCRIPTS  *** */


/* /6 */
/*     DATA COVPRT/14/, F/10/, NEEDHD/36/, PRUNIT/21/, REGD/67/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --v;
    rd_dim1 = *n;
    --rd;

    /* Function Body */
    pu = iv[21];
    if (pu == 0) {
	goto L999;
    }
    if (iv[14] < 2) {
	goto L999;
    }
    if (iv[67] <= 0) {
	goto L999;
    }
    iv[36] = 1;
    if (v[10] != 0.) {
	goto L10;
    } else {
	goto L30;
    }
L10:
	printf("\nREGRESSION DIAGNOSTIC = SQRT( G(I)**T * H(I)**-1 * G(I) / ABS(F) )...");
	for(i__1 = 1; i__1 <= rd_dim1; i__1++)
		printf(i__1 % 6 == 1 ? "\n  %# -11.3g" : " %# -11.3g", rd[i__1]);
	printf("\n");
    goto L999;
L30:
	printf("\nREGRESSION DIAGNOSTIC = SQRT( G(I)**T * H(I)**-1 * G(I) )...");
	for(i__1 = 1; i__1 <= rd_dim1; i__1++)
		printf(i__1 % 6 == 1 ? "\n  %# -11.3g" : " %# -11.3g", rd[i__1]);
	printf("\n");

L999:
    return 0;
/*  ***  LAST LINE OF DN2RDP FOLLOWS  *** */
} /* dn2rdp_ */

