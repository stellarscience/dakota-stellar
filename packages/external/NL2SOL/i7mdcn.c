/* i7mdcn.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

integer i7mdcn_(k)
integer *k;
{
    /* Initialized data */

    static integer mdperm[3] = { 2,4,1 };

    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer i1mach_();



/*  ***  RETURN INTEGER MACHINE-DEPENDENT CONSTANTS  *** */

/*     ***  K = 1 MEANS RETURN STANDARD OUTPUT UNIT NUMBER.   *** */
/*     ***  K = 2 MEANS RETURN ALTERNATE OUTPUT UNIT NUMBER.  *** */
/*     ***  K = 3 MEANS RETURN  INPUT UNIT NUMBER.            *** */
/*          (NOTE -- K = 2, 3 ARE USED ONLY BY TEST PROGRAMS.) */

/*  +++  PORT VERSION FOLLOWS... */
    ret_val = i1mach_(&mdperm[*k - 1]);
/*  +++  END OF PORT VERSION  +++ */

/*  +++  NON-PORT VERSION FOLLOWS... */
/*     INTEGER MDCON(3) */
/*     DATA MDCON(1)/6/, MDCON(2)/8/, MDCON(3)/5/ */
/*     I7MDCN = MDCON(K) */
/*  +++  END OF NON-PORT VERSION  +++ */

/* L999: */
    return ret_val;
/*  ***  LAST CARD OF I7MDCN FOLLOWS  *** */
} /* i7mdcn_ */

