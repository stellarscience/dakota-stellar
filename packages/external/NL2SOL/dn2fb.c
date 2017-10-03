/* dn2fb.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dn2fb_(integer *n, integer *p, doublereal *x, doublereal 
	*b, S_fp calcr, integer *iv, integer *liv, integer *lv, doublereal *v,
	 integer *ui, doublereal *ur, U_fp uf)
{
    /* Initialized data */

    static doublereal hlim = .1;
    static doublereal negpt5 = -.5;
    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal h__;
    static integer i__, k;
    static doublereal t;
    static integer d1;
    static doublereal h0;
    static integer n1, n2, r1, dk, nf, ng, rn;
    static doublereal xk;
    static integer j1k, dr1, rd1, iv1;
    static doublereal xk1;
    extern /* Subroutine */ int drn2gb_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dv7scp_(integer *, doublereal *, 
	    doublereal *), divset_(integer *, integer *, integer *, integer *,
	     doublereal *);


/*  ***  MINIMIZE A NONLINEAR SUM OF SQUARES USING RESIDUAL VALUES ONLY.. */
/*  ***  THIS AMOUNTS TO   DN2G WITHOUT THE SUBROUTINE PARAMETER CALCJ. */

/*  ***  PARAMETERS  *** */

/* /6 */
/*     INTEGER IV(LIV), UI(1) */
/*     DOUBLE PRECISION X(P), B(2,P), V(LV), UR(1) */
/* /7 */
/* / */

/* -----------------------------  DISCUSSION  ---------------------------- */

/*        THIS AMOUNTS TO SUBROUTINE NL2SNO (REF. 1) MODIFIED TO HANDLE */
/*     SIMPLE BOUNDS ON THE VARIABLES... */
/*           B(1,I) .LE. X(I) .LE. B(2,I), I = 1(1)P. */
/*        THE PARAMETERS FOR  DN2FB ARE THE SAME AS THOSE FOR  DN2GB */
/*     (WHICH SEE), EXCEPT THAT CALCJ IS OMITTED.  INSTEAD OF CALLING */
/*     CALCJ TO OBTAIN THE JACOBIAN MATRIX OF R AT X,  DN2FB COMPUTES */
/*     AN APPROXIMATION TO IT BY FINITE (FORWARD) DIFFERENCES -- SEE */
/*     V(DLTFDJ) BELOW.   DN2FB DOES NOT COMPUTE A COVARIANCE MATRIX. */
/*        THE NUMBER OF EXTRA CALLS ON CALCR USED IN COMPUTING THE JACO- */
/*     BIAN APPROXIMATION ARE NOT INCLUDED IN THE FUNCTION EVALUATION */
/*     COUNT IV(NFCALL), BUT ARE RECORDED IN IV(NGCALL) INSTEAD. */

/* V(DLTFDJ)... V(43) HELPS CHOOSE THE STEP SIZE USED WHEN COMPUTING THE */
/*             FINITE-DIFFERENCE JACOBIAN MATRIX.  FOR DIFFERENCES IN- */
/*             VOLVING X(I), THE STEP SIZE FIRST TRIED IS */
/*                       V(DLTFDJ) * MAX(ABS(X(I)), 1/D(I)), */
/*             WHERE D IS THE CURRENT SCALE VECTOR (SEE REF. 1).  (IF */
/*             THIS STEP IS TOO BIG, I.E., IF CALCR SETS NF TO 0, THEN */
/*             SMALLER STEPS ARE TRIED UNTIL THE STEP SIZE IS SHRUNK BE- */
/*             LOW 1000 * MACHEP, WHERE MACHEP IS THE UNIT ROUNDOFF. */
/*             DEFAULT = MACHEP**0.5. */

/*  ***  REFERENCE  *** */

/* 1.  DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1981), AN ADAPTIVE */
/*             NONLINEAR LEAST-SQUARES ALGORITHM, ACM TRANS. MATH. */
/*             SOFTWARE, VOL. 7, NO. 3. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */

/* +++++++++++++++++++++++++++  DECLARATIONS  +++++++++++++++++++++++++++ */

/*  ***  EXTERNAL SUBROUTINES  *** */


/* DIVSET.... PROVIDES DEFAULT IV AND V INPUT COMPONENTS. */
/* DRN2GB... CARRIES OUT OPTIMIZATION ITERATIONS. */
/* DN2RDP... PRINTS REGRESSION DIAGNOSTICS. */
/* DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  IV AND V COMPONENTS  *** */

/* /6 */
/*     DATA COVREQ/15/, D/27/, DINIT/38/, DLTFDJ/43/, J/70/, MODE/35/, */
/*    1     NEXTV/47/, NFCALL/6/, NFGCAL/7/, NGCALL/30/, NGCOV/53/, */
/*    2     R/61/, REGD0/82/, TOOBIG/2/, VNEED/4/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    b -= 3;
    --x;
    --iv;
    --v;
    --ui;
    --ur;

    /* Function Body */

/* ---------------------------------  BODY  ------------------------------ */

    if (iv[1] == 0) {
	divset_(&c__1, &iv[1], liv, lv, &v[1]);
    }
    iv[15] = 0;
    iv1 = iv[1];
    if (iv1 == 14) {
	goto L10;
    }
    if (iv1 > 2 && iv1 < 12) {
	goto L10;
    }
    if (iv1 == 12) {
	iv[1] = 13;
    }
    if (iv[1] == 13) {
	iv[4] = iv[4] + *p + *n * (*p + 2);
    }
    drn2gb_(&b[3], &x[1], &v[1], &iv[1], liv, lv, n, n, &n1, &n2, p, &v[1], &
	    v[1], &v[1], &x[1]);
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  STORAGE ALLOCATION  *** */

    iv[27] = iv[47];
    iv[61] = iv[27] + *p;
    iv[82] = iv[61] + *n;
    iv[70] = iv[82] + *n;
    iv[47] = iv[70] + *n * *p;
    if (iv1 == 13) {
	goto L999;
    }

L10:
    d1 = iv[27];
    dr1 = iv[70];
    r1 = iv[61];
    rn = r1 + *n - 1;
    rd1 = iv[82];

L20:
    drn2gb_(&b[3], &v[d1], &v[dr1], &iv[1], liv, lv, n, n, &n1, &n2, p, &v[r1]
	    , &v[rd1], &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L30;
    } else if (i__1 == 0) {
	goto L50;
    } else {
	goto L999;
    }

/*  ***  NEW FUNCTION VALUE (R VALUE) NEEDED  *** */

L30:
    nf = iv[6];
    (*calcr)(n, p, &x[1], &nf, &v[r1], &ui[1], &ur[1], (U_fp)uf);
    if (nf > 0) {
	goto L40;
    }
    iv[2] = 1;
    goto L20;
L40:
    if (iv[1] > 0) {
	goto L20;
    }

/*  ***  COMPUTE FINITE-DIFFERENCE APPROXIMATION TO DR = GRAD. OF R  *** */

/*     *** INITIALIZE D IF NECESSARY *** */

L50:
    if (iv[35] < 0 && v[38] == zero) {
	dv7scp_(p, &v[d1], &one);
    }

    j1k = dr1;
    dk = d1;
    ng = iv[30] - 1;
    if (iv[1] == -1) {
	--iv[53];
    }
    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {
	if (b[(k << 1) + 1] >= b[(k << 1) + 2]) {
	    goto L110;
	}
	xk = x[k];
/* Computing MAX */
	d__1 = abs(xk), d__2 = one / v[dk];
	h__ = v[43] * max(d__1,d__2);
	h0 = h__;
	++dk;
	t = negpt5;
	xk1 = xk + h__;
	if (xk - h__ >= b[(k << 1) + 1]) {
	    goto L60;
	}
	t = -t;
	if (xk1 > b[(k << 1) + 2]) {
	    goto L80;
	}
L60:
	if (xk1 <= b[(k << 1) + 2]) {
	    goto L70;
	}
	t = -t;
	h__ = -h__;
	xk1 = xk + h__;
	if (xk1 < b[(k << 1) + 1]) {
	    goto L80;
	}
L70:
	x[k] = xk1;
	nf = iv[7];
	(*calcr)(n, p, &x[1], &nf, &v[j1k], &ui[1], &ur[1], (U_fp)uf);
	++ng;
	if (nf > 0) {
	    goto L90;
	}
	h__ = t * h__;
	xk1 = xk + h__;
	if ((d__1 = h__ / h0, abs(d__1)) >= hlim) {
	    goto L70;
	}
L80:
	iv[2] = 1;
	iv[30] = ng;
	goto L20;
L90:
	x[k] = xk;
	iv[30] = ng;
	i__2 = rn;
	for (i__ = r1; i__ <= i__2; ++i__) {
	    v[j1k] = (v[j1k] - v[i__]) / h__;
	    ++j1k;
/* L100: */
	}
	goto L120;
/*        *** SUPPLY A ZERO DERIVATIVE FOR CONSTANT COMPONENTS... */
L110:
	dv7scp_(n, &v[j1k], &zero);
	j1k += *n;
L120:
	;
    }
    goto L20;

L999:
    return 0;

/*  ***  LAST CARD OF  DN2FB FOLLOWS  *** */
} /* dn2fb_ */

