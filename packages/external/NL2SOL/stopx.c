/* stopx.f -- translated by f2c (version 19970211).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

logical stopx_(idummy)
integer *idummy;
{
    /* System generated locals */
    logical ret_val;

/*     *****PARAMETERS... */

/*     .................................................................. 
*/

/*     *****PURPOSE... */
/*     THIS FUNCTION MAY SERVE AS THE STOPX (ASYNCHRONOUS INTERRUPTION) */
/*     FUNCTION FOR THE NL2SOL (NONLINEAR LEAST-SQUARES) PACKAGE AT */
/*     THOSE INSTALLATIONS WHICH DO NOT WISH TO IMPLEMENT A */
/*     DYNAMIC STOPX. */

/*     *****ALGORITHM NOTES... */
/*     AT INSTALLATIONS WHERE THE NL2SOL SYSTEM IS USED */
/*     INTERACTIVELY, THIS DUMMY STOPX SHOULD BE REPLACED BY A */
/*     FUNCTION THAT RETURNS .TRUE. IF AND ONLY IF THE INTERRUPT */
/*     (BREAK) KEY HAS BEEN PRESSED SINCE THE LAST CALL ON STOPX. */

/*     .................................................................. 
*/

    ret_val = FALSE_;
    return ret_val;
} /* stopx_ */

