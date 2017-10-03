/* dparck.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <stdio.h>
#include "f2c.h"

 static void
s_copy(char *a, char *b, ftnlen la, ftnlen lb)	/* from libf2c */
{
	char *aend, *bend;

	aend = a + la;

	if(la <= lb)
		if (a <= b || a >= b + la)
			while(a < aend)
				*a++ = *b++;
		else
			for(b += la; a < aend; )
				*--aend = *--b;

	else {
		bend = b + lb;
		if (a <= b || a >= bend)
			while(b < bend)
				*a++ = *b++;
		else {
			a += lb;
			while(b < bend)
				*--a = *--bend;
			a += lb;
			}
		while(a < aend)
			*a++ = ' ';
		}
	}

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__5 = 5;

/* Subroutine */ int dparck_(integer *alg, doublereal *d__, integer *iv,
	integer *liv, integer *lv, integer *n, doublereal *v)
{
    /* Initialized data */

    static doublereal big = 0.;
    static char cngd[4*3] = "---C" "HANG" "ED V";
    static char dflt[4*3] = "NOND" "EFAU" "LT V";
    static integer ijmp = 33;
    static integer jlim[4] = { 0,24,0,24 };
    static integer ndflt[4] = { 32,25,32,25 };
    static integer miniv[4] = { 82,59,103,103 };
    static doublereal machep = -1.;
    static doublereal tiny = 1.;
    static doublereal zero = 0.;
    static char vn[4*2*34] = "EPSL" "ON.." "PHMN" "FC.." "PHMX" "FC.." "DECF"
	    "AC.." "INCF" "AC.." "RDFC" "MN.." "RDFC" "MX.." "TUNE" "R1.."
	    "TUNE" "R2.." "TUNE" "R3.." "TUNE" "R4.." "TUNE" "R5.." "AFCT"
	    "OL.." "RFCT" "OL.." "XCTO" "L..." "XFTO" "L..." "LMAX" "0..."
	    "LMAX" "S..." "SCTO" "L..." "DINI" "T..." "DTIN" "IT.." "D0IN"
	    "IT.." "DFAC" "...." "DLTF" "DC.." "DLTF" "DJ.." "DELT" "A0.."
	    "FUZZ" "...." "RLIM" "IT.." "COSM" "IN.." "HUBE" "RC.." "RSPT"
	    "OL.." "SIGM" "IN.." "ETA0" "...." "BIAS" "....";
    static doublereal vm[34] = { .001,-.99,.001,.01,1.2,.01,1.2,0.,0.,.001,
	    -1.,0.0,0.,0.0,0.,0.,0.0,0.0,0.,-10.,0.,0.,0.,0.0,0.0,0.0,1.01,
	    1e10,0.0,0.,0.,0.,0.0,0. };
    static doublereal vx[34] = { .9,-.001,10.,.8,100.,.8,100.,.5,.5,1.,1.,0.0,
	    0.0,.1,1.,1.,0.0,0.0,1.,0.0,0.0,0.0,1.,1.,1.,1.,1e10,0.0,1.,0.0,
	    1.,1.,1.,1. };
    static char varnm[1*2] = "P" "P";
    static char sh[1*2] = "S" "H";


    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m, ii;
    static doublereal vk;
    static integer pu, iv1, alg1, miv1, miv2;
    static char which[4*3];
    extern doublereal dr7mdc_(integer *);
    extern /* Subroutine */ int dv7dfl_(integer *, integer *, doublereal *),
	    dv7cpy_(integer *, doublereal *, doublereal *);
    static integer parsv1, ndfalt;
    extern /* Subroutine */ int divset_(integer *, integer *, integer *,
	    integer *, doublereal *);


/*  ***  CHECK ***SOL (VERSION 2.3) PARAMETERS, PRINT CHANGED VALUES  *** */

/*  ***  ALG = 1 FOR REGRESSION, ALG = 2 FOR GENERAL UNCONSTRAINED OPT. */


/* DIVSET  -- SUPPLIES DEFAULT VALUES TO BOTH IV AND V. */
/* DR7MDC -- RETURNS MACHINE-DEPENDENT CONSTANTS. */
/* DV7CPY  -- COPIES ONE VECTOR TO ANOTHER. */
/* DV7DFL  -- SUPPLIES DEFAULT PARAMETER VALUES TO V ALONE. */

/*  ***  LOCAL VARIABLES  *** */

/* /6S */
/*     INTEGER VARNM(2), SH(2) */
/*     REAL CNGD(3), DFLT(3), VN(2,34), WHICH(3) */
/* /7S */
/* / */

/*  ***  IV AND V SUBSCRIPTS  *** */



/* /6 */
/*     DATA ALGSAV/51/, DINIT/38/, DTYPE/16/, DTYPE0/54/, EPSLON/19/, */
/*    1     INITS/25/, IVNEED/3/, LASTIV/44/, LASTV/45/, LMAT/42/, */
/*    2     NEXTIV/46/, NEXTV/47/, NVDFLT/50/, OLDN/38/, PARPRT/20/, */
/*    3     PARSAV/49/, PERM/58/, PRUNIT/21/, VNEED/4/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --iv;
    --v;
    --d__;

    /* Function Body */
/* /6S */
/*     DATA VN(1,1),VN(2,1)/4HEPSL,4HON../ */
/*     DATA VN(1,2),VN(2,2)/4HPHMN,4HFC../ */
/*     DATA VN(1,3),VN(2,3)/4HPHMX,4HFC../ */
/*     DATA VN(1,4),VN(2,4)/4HDECF,4HAC../ */
/*     DATA VN(1,5),VN(2,5)/4HINCF,4HAC../ */
/*     DATA VN(1,6),VN(2,6)/4HRDFC,4HMN../ */
/*     DATA VN(1,7),VN(2,7)/4HRDFC,4HMX../ */
/*     DATA VN(1,8),VN(2,8)/4HTUNE,4HR1../ */
/*     DATA VN(1,9),VN(2,9)/4HTUNE,4HR2../ */
/*     DATA VN(1,10),VN(2,10)/4HTUNE,4HR3../ */
/*     DATA VN(1,11),VN(2,11)/4HTUNE,4HR4../ */
/*     DATA VN(1,12),VN(2,12)/4HTUNE,4HR5../ */
/*     DATA VN(1,13),VN(2,13)/4HAFCT,4HOL../ */
/*     DATA VN(1,14),VN(2,14)/4HRFCT,4HOL../ */
/*     DATA VN(1,15),VN(2,15)/4HXCTO,4HL.../ */
/*     DATA VN(1,16),VN(2,16)/4HXFTO,4HL.../ */
/*     DATA VN(1,17),VN(2,17)/4HLMAX,4H0.../ */
/*     DATA VN(1,18),VN(2,18)/4HLMAX,4HS.../ */
/*     DATA VN(1,19),VN(2,19)/4HSCTO,4HL.../ */
/*     DATA VN(1,20),VN(2,20)/4HDINI,4HT.../ */
/*     DATA VN(1,21),VN(2,21)/4HDTIN,4HIT../ */
/*     DATA VN(1,22),VN(2,22)/4HD0IN,4HIT../ */
/*     DATA VN(1,23),VN(2,23)/4HDFAC,4H..../ */
/*     DATA VN(1,24),VN(2,24)/4HDLTF,4HDC../ */
/*     DATA VN(1,25),VN(2,25)/4HDLTF,4HDJ../ */
/*     DATA VN(1,26),VN(2,26)/4HDELT,4HA0../ */
/*     DATA VN(1,27),VN(2,27)/4HFUZZ,4H..../ */
/*     DATA VN(1,28),VN(2,28)/4HRLIM,4HIT../ */
/*     DATA VN(1,29),VN(2,29)/4HCOSM,4HIN../ */
/*     DATA VN(1,30),VN(2,30)/4HHUBE,4HRC../ */
/*     DATA VN(1,31),VN(2,31)/4HRSPT,4HOL../ */
/*     DATA VN(1,32),VN(2,32)/4HSIGM,4HIN../ */
/*     DATA VN(1,33),VN(2,33)/4HETA0,4H..../ */
/*     DATA VN(1,34),VN(2,34)/4HBIAS,4H..../ */
/* /7S */
/* / */


/* /6S */
/*     DATA VARNM(1)/1HP/, VARNM(2)/1HP/, SH(1)/1HS/, SH(2)/1HH/ */
/*     DATA CNGD(1),CNGD(2),CNGD(3)/4H---C,4HHANG,4HED V/, */
/*    1     DFLT(1),DFLT(2),DFLT(3)/4HNOND,4HEFAU,4HLT V/ */
/* /7S */
/* / */

/* ...............................  BODY  ................................ */

    pu = 0;
    if (21 <= *liv) {
	pu = iv[21];
    }
    if (51 > *liv) {
	goto L20;
    }
    if (*alg == iv[51]) {
	goto L20;
    }
    if (pu != 0)
	printf("\nTHE FIRST PARAMETER TO DIVSET SHOULD BE %d RATHER THAN %d\n",
		*alg, iv[51]);
    iv[1] = 67;
    goto L999;
L20:
    if (*alg < 1 || *alg > 4) {
	goto L340;
    }
    miv1 = miniv[*alg - 1];
    if (iv[1] == 15) {
	goto L360;
    }
    alg1 = (*alg - 1) % 2 + 1;
    if (iv[1] == 0) {
	divset_(alg, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 != 13 && iv1 != 12) {
	goto L30;
    }
    if (58 <= *liv) {
/* Computing MAX */
	i__1 = miv1, i__2 = iv[58] - 1;
	miv1 = max(i__1,i__2);
    }
    if (3 <= *liv) {
	miv2 = miv1 + max(iv[3],0);
    }
    if (44 <= *liv) {
	iv[44] = miv2;
    }
    if (*liv < miv1) {
	goto L300;
    }
    iv[3] = 0;
    iv[45] = max(iv[4],0) + iv[42] - 1;
    iv[4] = 0;
    if (*liv < miv2) {
	goto L300;
    }
    if (*lv < iv[45]) {
	goto L320;
    }
L30:
    if (iv1 < 12 || iv1 > 14) {
	goto L60;
    }
    if (*n >= 1) {
	goto L50;
    }
    iv[1] = 81;
    if (pu == 0) {
	goto L999;
    }
	printf("\n/// BAD%.1s = %d\n", varnm + (alg1 - 1), *n);
    goto L999;
L50:
    if (iv1 != 14) {
	iv[46] = iv[58];
    }
    if (iv1 != 14) {
	iv[47] = iv[42];
    }
    if (iv1 == 13) {
	goto L999;
    }
    k = iv[49] - 19;
    i__1 = *lv - k;
    dv7dfl_(&alg1, &i__1, &v[k + 1]);
    iv[54] = 2 - alg1;
    iv[38] = *n;
    s_copy(which, dflt, (ftnlen)4, (ftnlen)4);
    s_copy(which + 4, dflt + 4, (ftnlen)4, (ftnlen)4);
    s_copy(which + 8, dflt + 8, (ftnlen)4, (ftnlen)4);
    goto L110;
L60:
    if (*n == iv[38]) {
	goto L80;
    }
    iv[1] = 17;
    if (pu == 0) {
	goto L999;
    }
	printf("\n/// %.1s CHANGED FROM %d TO %d\n", varnm + (alg1 - 1), iv[38], *n);
    goto L999;

L80:
    if (iv1 <= 11 && iv1 >= 1) {
	goto L100;
    }
    iv[1] = 80;
    if (pu != 0)
	printf("\n///  IV(1) = %d SHOULD BE BETWEEN 0 and 14.\n", iv1);
    goto L999;

L100:
    s_copy(which, cngd, (ftnlen)4, (ftnlen)4);
    s_copy(which + 4, cngd + 4, (ftnlen)4, (ftnlen)4);
    s_copy(which + 8, cngd + 8, (ftnlen)4, (ftnlen)4);

L110:
    if (iv1 == 14) {
	iv1 = 12;
    }
    if (big > tiny) {
	goto L120;
    }
    tiny = dr7mdc_(&c__1);
    machep = dr7mdc_(&c__3);
    big = dr7mdc_(&c__6);
    vm[11] = machep;
    vx[11] = big;
    vx[12] = big;
    vm[13] = machep;
    vm[16] = tiny;
    vx[16] = big;
    vm[17] = tiny;
    vx[17] = big;
    vx[19] = big;
    vx[20] = big;
    vx[21] = big;
    vm[23] = machep;
    vm[24] = machep;
    vm[25] = machep;
    vx[27] = dr7mdc_(&c__5);
    vm[28] = machep;
    vx[29] = big;
    vm[32] = machep;
L120:
    m = 0;
    i__ = 1;
    j = jlim[alg1 - 1];
    k = 19;
    ndfalt = ndflt[alg1 - 1];
    i__1 = ndfalt;
    for (l = 1; l <= i__1; ++l) {
	vk = v[k];
	if (vk >= vm[i__ - 1] && vk <= vx[i__ - 1]) {
	    goto L140;
	}
	m = k;
	if (pu != 0)
		printf("\n///  %.4s%.4s.. V(%d) = %#.3g SHOULD BE BETWEEN %#.3g AND %#.3g\n",
			vn + ((i__ << 1) - 2 << 2), vn + ((i__ << 1) - 1 << 2),
			k, vk, vm[i__ - 1], vx[i__ - 1]);
L140:
	++k;
	++i__;
	if (i__ == j) {
	    i__ = ijmp;
	}
/* L150: */
    }

    if (iv[50] == ndfalt) {
	goto L170;
    }
    iv[1] = 51;
    if (pu == 0) {
	goto L999;
    }
	printf("\nIV(NVDFLT) = %d RATHER THAN %d\n", iv[50], ndfalt);
    goto L999;
L170:
    if ((iv[16] > 0 || v[38] > zero) && iv1 == 12) {
	goto L200;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] > zero) {
	    goto L190;
	}
	m = 18;
	if (pu != 0)
		printf("\n///  D(%d) = %#.3g SHOULD BE POSITIVE\n", i__, d__[i__]);
L190:
	;
    }
L200:
    if (m == 0) {
	goto L210;
    }
    iv[1] = m;
    goto L999;

L210:
    if (pu == 0 || iv[20] == 0) {
	goto L999;
    }
    if (iv1 != 12 || iv[25] == alg1 - 1) {
	goto L230;
    }
    m = 1;
	printf("\nNONDEFAULT VALUES....\nINIT%.1s..... IV(25) = %d\n", sh + (alg1 - 1), iv[25]);
L230:
    if (iv[16] == iv[54]) {
	goto L250;
    }
    if (m == 0)
	printf("\n%.12sALUES...\n\n", which);
    m = 1;
	printf("DTYPE..... IV(16) = %d\n", iv[16]);
L250:
    i__ = 1;
    j = jlim[alg1 - 1];
    k = 19;
    l = iv[49];
    ndfalt = ndflt[alg1 - 1];
    i__1 = ndfalt;
    for (ii = 1; ii <= i__1; ++ii) {
	if (v[k] == v[l]) {
	    goto L280;
	}
	if (m == 0)
		printf("\n%.12sALUES...\n\n", which);
	m = 1;
	printf("%.8s.. V(%d) = %#.7g\n", vn + ((i__ << 1) - 2 << 2), k, v[k]);
L280:
	++k;
	++l;
	++i__;
	if (i__ == j) {
	    i__ = ijmp;
	}
/* L290: */
    }

    iv[54] = iv[16];
    parsv1 = iv[49];
    dv7cpy_(&iv[50], &v[parsv1], &v[19]);
    goto L999;

L300:
    iv[1] = 15;
    if (pu == 0) {
	goto L999;
    }
	printf("\n/// LIV = %d MUST BE AT LEAST %d\n", *liv, miv2);
    if (*liv < miv1) {
	goto L999;
    }
    if (*lv < iv[45]) {
	goto L320;
    }
    goto L999;

L320:
    iv[1] = 16;
    if (pu != 0)
	printf("\n/// LV = %d MUST BE AT LEAST %d\n", *lv, iv[45]);
    goto L999;

L340:
    iv[1] = 67;
    if (pu != 0)
	printf("\n/// ALG = %d MUST BE 1, 2, 3, or 4\n", *alg);
    goto L999;
L360:
    if (pu != 0)
	printf("\n/// LIV = %d MUST BE AT LEAST %d TO COMPUTE TRUE MIN. LIV AND MIN. LV\n",
		*liv, miv1);
    if (44 <= *liv) {
	iv[44] = miv1;
    }
    if (45 <= *liv) {
	iv[45] = 0;
    }

L999:
    return 0;
/*  ***  LAST LINE OF DPARCK FOLLOWS  *** */
} /* dparck_ */

