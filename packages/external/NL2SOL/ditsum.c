/* ditsum.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <stdio.h>
#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int ditsum_(doublereal *d__, doublereal *g, integer *iv,
	integer *liv, integer *lv, integer *p, doublereal *v, doublereal *x)
{
    /* Initialized data */

    static char model1[4*6] = "    " "    " "    " "    " "  G " "  S ";
    static char model2[4*6] = " G  " " S  " "G-S " "S-G " "-S-G" "-G-S";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, m, nf, ng, ol, pu, iv1, alg;
    static doublereal oldf, reldf, nreldf, preldf;


/*  ***  PRINT ITERATION SUMMARY FOR ***SOL (VERSION 2.3)  *** */

/*  ***  PARAMETER DECLARATIONS  *** */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */

/* /6S */
/*     REAL MODEL1(6), MODEL2(6) */
/* /7S */
/* / */

/*  ***  NO EXTERNAL FUNCTIONS OR SUBROUTINES  *** */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA ALGSAV/51/, NEEDHD/36/, NFCALL/6/, NFCOV/52/, NGCALL/30/, */
/*    1     NGCOV/53/, NITER/31/, OUTLEV/19/, PRNTIT/39/, PRUNIT/21/, */
/*    2     SOLPRT/22/, STATPR/23/, SUSED/64/, X0PRT/24/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DSTNRM/2/, F/10/, F0/13/, FDIF/11/, NREDUC/6/, PREDUC/7/, */
/*    1     RELDX/17/, STPPAR/5/ */
/* /7 */
/* / */

/* /6 */
/*     DATA ZERO/0.D+0/ */
/* /7 */
/* / */
/* /6S */
/*     DATA MODEL1(1)/4H    /, MODEL1(2)/4H    /, MODEL1(3)/4H    /, */
/*    1     MODEL1(4)/4H    /, MODEL1(5)/4H  G /, MODEL1(6)/4H  S /, */
/*    2     MODEL2(1)/4H G  /, MODEL2(2)/4H S  /, MODEL2(3)/4HG-S /, */
/*    3     MODEL2(4)/4HS-G /, MODEL2(5)/4H-S-G/, MODEL2(6)/4H-G-S/ */
/* /7S */
    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --g;
    --d__;

    /* Function Body */
/* / */

/* -------------------------------  BODY  -------------------------------- */

    pu = iv[21];
    if (pu == 0) {
	goto L999;
    }
    iv1 = iv[1];
    if (iv1 > 62) {
	iv1 += -51;
    }
    ol = iv[19];
    alg = (iv[51] - 1) % 2 + 1;
    if (iv1 < 2 || iv1 > 15) {
	goto L370;
    }
    if (iv1 >= 12) {
	goto L120;
    }
    if (iv1 == 2 && iv[31] == 0) {
	goto L390;
    }
    if (ol == 0) {
	goto L120;
    }
    if (iv1 >= 10 && iv[39] == 0) {
	goto L120;
    }
    if (iv1 > 2) {
	goto L10;
    }
    ++iv[39];
    if (iv[39] < abs(ol)) {
	goto L999;
    }
L10:
    nf = iv[6] - abs(iv[52]);
    iv[39] = 0;
    reldf = 0.;
    preldf = 0.;
/* Computing MAX */
    d__1 = abs(v[13]), d__2 = abs(v[10]);
    oldf = max(d__1,d__2);
    if (oldf <= 0.) {
	goto L20;
    }
    reldf = v[11] / oldf;
    preldf = v[7] / oldf;
L20:
    if (ol > 0) {
	goto L60;
    }

/*        ***  PRINT SHORT SUMMARY LINE  *** */

    if (iv[36] == 1 && alg == 1)
	printf("\n  IT  NF    F         RELDF   PRELDF   RELDX   MODEL  STPPAR\n");
    if (iv[36] == 1 && alg == 2)
	printf("\n   IT   NF     F           RELDF    PRELDF    RELDX   STPPAR\n");
    iv[36] = 0;
    if (alg == 2) {
	goto L50;
    }
    m = iv[64];
	printf("%5d %4d %# -10.3g %# -8.2g %# -8.2g %# -7.1g%.3s%.4s %# -7.1g\n",
		iv[31], nf, v[10], reldf, preldf, v[17],
		model1 + (m - 1 << 2), model2 + (m - 1 << 2), v[5]);
    goto L120;

L50:
	printf("%5d %4d %# -10.3g %# -9.2g %# -9.2g %# -8.1g %# -8.1g\n",
		iv[31], nf, v[10], reldf, preldf, v[17], v[5]);
    goto L120;

/*     ***  PRINT LONG SUMMARY LINE  *** */

L60:
    if (iv[36] == 1 && alg == 1)
	printf("\n   IT   NF      F       RELDF   PRELDF   RELDX   MODEL  STPPAR  D*STEP  NPRELDF\n");
    if (iv[36] == 1 && alg == 2)
	printf("\n   IT   NF       F        RELDF    PRELDF    RELDX   STPPAR   D*STEP   NPRELDF\n");
    iv[36] = 0;
    nreldf = 0.;
    if (oldf > 0.) {
	nreldf = v[6] / oldf;
    }
    if (alg == 2) {
	goto L90;
    }
    m = iv[64];
	printf("%5d %4d %# -10.3g %# -8.2g %# -8.2g %# -7.1g%.3s%.4s %# -7.1g %# -7.1g %# -8.2g\n",
		iv[31], nf, v[10], reldf, preldf, v[17],
		model1 + (m - 1 << 2), model2 + (m - 1 << 2), v[5],
		v[2], nreldf);
    goto L120;

L90:
	printf("%5d %4d %# -10.3g %# -9.2g %# -9.2g %# -8.1g %# -8.1g %# -8.1g %# -9.2g\n",
		iv[31], nf, v[10], reldf, preldf, v[17], v[5], v[2], nreldf);

L120:
    if (iv1 <= 2) {
	goto L999;
    }
    i__ = iv[23];
    if (i__ == -1) {
	goto L460;
    }
    if (i__ + iv1 < 0) {
	goto L460;
    }
    switch (iv1) {
	case 1:  goto L999;
	case 2:  goto L999;
	case 3:  goto L130;
	case 4:  goto L150;
	case 5:  goto L170;
	case 6:  goto L190;
	case 7:  goto L210;
	case 8:  goto L230;
	case 9:  goto L250;
	case 10:  goto L270;
	case 11:  goto L290;
	case 12:  goto L310;
	case 13:  goto L330;
	case 14:  goto L350;
	case 15:  goto L500;
    }

L130:
	printf("\n***** X-CONVERGENCE *****\n");
    goto L430;

L150:
	printf("\n***** RELATIVE FUNCTION CONVERGENCE *****\n");
    goto L430;

L170:
	printf("\n***** X- AND RELATIVE FUNCTION CONVERGENCE *****\n");
    goto L430;

L190:
	printf("\n***** ABSOLUTE FUNCTION CONVERGENCE *****\n");
    goto L430;

L210:
	printf("\n***** SINGULAR CONVERGENCE *****\n");
    goto L430;

L230:
	printf("\n***** FALSE CONVERGENCE *****\n");
    goto L430;

L250:
	printf("\n***** FUNCTION EVALUATION LIMIT *****\n");
    goto L430;

L270:
	printf("\n***** ITERATION LIMIT *****\n");
    goto L430;

L290:
	printf("\n***** STOPX *****\n");
    goto L430;

L310:
	printf("\n***** INITIAL F(X) CANNOT BE COMPUTED *****\n");

    goto L390;

L330:
	printf("\n***** BAD PARAMETERS TO ASSESS *****\n");
    goto L999;

L350:
	printf("\n***** GRADIENT COULD NOT BE COMPUTED *****\n");
    if (iv[31] > 0) {
	goto L460;
    }
    goto L390;

L370:
	printf("\n***** IV(1) = %4d *****\n", iv[1]);
    goto L999;

/*  ***  INITIAL CALL ON DITSUM  *** */

L390:
    if (iv[24] != 0) {
	printf("\n    I     INITIAL X(I)       D(I)\n\n");
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__)
		printf("%5d    %# -13.6g      %# .3g\n", i__, x[i__], d__[i__]);
    }
/*     *** THE FOLLOWING ARE TO AVOID UNDEFINED VARIABLES WHEN THE */
/*     *** FUNCTION EVALUATION LIMIT IS 1... */
    v[2] = 0.;
    v[11] = 0.;
    v[6] = 0.;
    v[7] = 0.;
    v[17] = 0.;
    if (iv1 >= 12) {
	goto L999;
    }
    iv[36] = 0;
    iv[39] = 0;
    if (ol == 0) {
	goto L999;
    }
    if (ol < 0 && alg == 1)
	printf("\n   IT   NF     F       RELDF    PRELDF   RELDX   MODEL  STPPAR\n");
    if (ol < 0 && alg == 2)
	printf("\n   IT   NF     F          RELDF    PRELDF    RELDX   STPPAR\n");
    if (ol > 0 && alg == 1)
	printf("\n   IT   NF     F       RELDF    PRELDF   RELDX   MODEL  STPPAR  D*STEP  NPRELDF\n");
    if (ol > 0 && alg == 2)
	printf("\n   IT   NF     F          RELDF    PRELDF    RELDX   STPPAR   D*STEP   NPRELDF\n");
    if (alg == 1)
	printf("\n    0 %4d %# -9.3g\n", iv[6], v[10]);
    if (alg == 2)
	printf("\n    0 %4d %# -10.3g\n", iv[6], v[10]);
    goto L999;

/*  ***  PRINT VARIOUS INFORMATION REQUESTED ON SOLUTION  *** */

L430:
    iv[36] = 1;
    if (iv[23] <= 0) {
	goto L460;
    }
/* Computing MAX */
    d__1 = abs(v[13]), d__2 = abs(v[10]);
    oldf = max(d__1,d__2);
    preldf = 0.;
    nreldf = 0.;
    if (oldf <= 0.) {
	goto L440;
    }
    preldf = v[7] / oldf;
    nreldf = v[6] / oldf;
L440:
    nf = iv[6] - iv[52];
    ng = iv[30] - iv[53];
	printf("\nFUNCTION    %# -13.6g   RELDX       %# .3g\nFUNC. EVALS %6d          GRAD. EVALS %6d\nPRELDF      %# -10.3g      NPRELDF     %# .3g\n",
		v[10], v[17], nf, ng, preldf, nreldf);

L460:
    if (iv[22] == 0) {
	goto L999;
    }
    iv[36] = 1;
    if (iv[51] > 2) {
	goto L999;
    }
	printf("\n    I      FINAL X(I)        D(I)          G(I)\n\n");
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__)
	printf("%5d      %# -13.6g    %# -10.3g  %# .3g\n", i__, x[i__], d__[i__], g[i__]);
    goto L999;

L500:
	printf("INCONSISTENT DIMENSIONS\n");
L999:
    return 0;
/*  ***  LAST CARD OF DITSUM FOLLOWS  *** */
} /* ditsum_ */

