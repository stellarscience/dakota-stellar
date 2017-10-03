/* dr7mdc.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__4 = 4;

doublereal dr7mdc_(k)
integer *k;
{
    /* Initialized data */

    static doublereal big = 0.;
    static doublereal eta = 0.;
    static doublereal machep = 0.;
    static doublereal zero = 0.;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern doublereal d1mach_();


/*  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  *** */


/*  ***  THE CONSTANT RETURNED DEPENDS ON K... */

/*  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS. */
/*  ***        K = 2... SQUARE ROOT OF ETA. */
/*  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH */
/*  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1. */
/*  ***        K = 4... SQUARE ROOT OF MACHEP. */
/*  ***        K = 5... SQUARE ROOT OF BIG (SEE K = 6). */
/*  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS. */

/* /+ */
/* / */

    if (big > zero) {
	goto L1;
    }
    big = d1mach_(&c__2);
    eta = d1mach_(&c__1);
    machep = d1mach_(&c__4);
L1:

/* -------------------------------  BODY  --------------------------------
 */

    switch ((int)*k) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	case 5:  goto L50;
	case 6:  goto L60;
    }

L10:
    ret_val = eta;
    goto L999;

L20:
    ret_val = sqrt(eta * 256.) / 16.;
    goto L999;

L30:
    ret_val = machep;
    goto L999;

L40:
    ret_val = sqrt(machep);
    goto L999;

L50:
    ret_val = sqrt(big / 256.) * 16.;
    goto L999;

L60:
    ret_val = big;

L999:
    return ret_val;
/*  ***  LAST CARD OF DR7MDC FOLLOWS  *** */
} /* dr7mdc_ */

