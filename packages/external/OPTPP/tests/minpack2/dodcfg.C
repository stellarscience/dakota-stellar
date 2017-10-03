// Optimal Design with Composite Materials from MINPACK-2 collection


#ifdef ANSI_HEADERS
#include <cstring>
#include <cmath>
#else
#include <string.h>
#include <math.h>
#endif

#include "OptNewton.h"
#include "NLF.h"
#include "minpack2.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

void dodcps(double t, double mu1, double mu2, double t1, double t2, 
	   double *result, int option, double lambda);

void dodcfg(int nx, int ny, ColumnVector& x, double &f, ColumnVector& fgrad, 
	    MINPACK_TASK task, double lambda)

{
    /* System generated locals */
    int i_1, i_2, i_3, i_4;
    double d_1, d_2;

    /* Local variables */
    double area, dpsi, dvdx, dvdy, temp, hxhy;
    int i, j, k;
    double v;
    int feval, geval;
    double gradv, dpsip, t1, t2, vb, hx, hy, vl, vr, vt;

    /* Function Body */
/*     ********** */

/*     Subroutine dodcfg */

/*     This subroutine computes the function and gradient of the */
/*     optimal design with composite materials problem. */

/*     The subroutine statement is */

/*       subroutine dodcfg(nx,ny,x,f,fgrad,task,lambda) */

/*     where */

/*       nx is an int variable. */
/*         On entry nx is the number of grid points in the first */
/*            coordinate direction. */
/*         On exit nx is unchanged. */

/*       ny is an int variable. */
/*         On entry ny is the number of grid points in the second */
/*            coordinate direction. */
/*         On exit ny is unchanged. */

/*       x is a double precision array of dimension nx*ny. */
/*         On entry x specifies the vector x if task = 'F', 'G', or 'FG'. 
*/
/*            Otherwise x need not be specified. */
/*         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise 
*/
/*            x is set according to task. */

/*       f is a double precision variable. */
/*         On entry f need not be specified. */
/*         On exit f is set to the function evaluated at x if task = 'F' 
*/
/*            or 'FG'. */

/*       fgrad is a double precision array of dimension nx*ny. */
/*         On entry fgrad need not be specified. */
/*         On exit fgrad contains the gradient evaluated at x if */
/*            task = 'G' or 'FG'. */

/*       task is a character variable. */
/*         On entry task specifies the action of the subroutine: */

/*            task               action */
/*            ----               ------ */
/*             'F'     Evaluate the function at x. */
/*             'G'     Evaluate the gradient at x. */
/*             'FG'    Evaluate the function and the gradient at x. */
/*             'XS'    Set x to the standard starting point xs. */

/*         On exit task is unchanged. */

/*       lambda is a double precision variable. */
/*         On entry lambda is the Lagrange multiplier. */
/*         On exit lambda is unchanged. */

/*     Subprograms called */

/*       MINPACK-supplied   ...   dodcps */

/*     MINPACK-2 Project. November 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick. */

/*     ********** */
/*     Initialization. */
    hx = 1. / (double) (nx + 1);
    hy = 1. / (double) (ny + 1);
    hxhy = hx * hy;
    area = hxhy * .5;
/*     Compute the break points. */
    t1 = sqrt(lambda * 2. * 1. / 2.);
    t2 = sqrt(lambda * 2. * 2. / 1.);
/*     Compute the standard starting point if task = 'XS'. */
    if (task == XStandard) {
	i_1 = ny;
	for (j = 1; j <= i_1; ++j) {
/* Computing MAX */
	    i_2 = j, i_3 = ny - j + 1;
	    temp = (double) min(i_3,i_2) * hy;
	    i_2 = nx;
	    for (i = 1; i <= i_2; ++i) {
		k = nx * (j - 1) + i;
/* Computing MAX */
/* Computing MAX */
		i_3 = i, i_4 = nx - i + 1;
		d_2 = (double) min(i_4,i_3) * hx;
/* Computing 2nd power */
		d_1 = min(temp,d_2);
		x(k) = -(d_1 * d_1);
/* L10: */
	    }
/* L20: */
	}
	return;
    }
    if (task == Function || task == FuncGrad) {
	feval = true;
    } else {
	feval = false;
    }
    if (task == Gradient || task == FuncGrad) {
	geval = true;
    } else {
	geval = false;
    }
/*     Evaluate the function if task = 'F', the gradient if task = 'G', */

/*     or both if task = 'FG'. */
    if (feval) {
	f = 0.;
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = 0.;
/* L30: */
	}
    }
/*     Computation of the function and the gradient over the lower */
/*     triangular elements. */
    i_1 = ny;
    for (j = 0; j <= i_1; ++j) {
	i_2 = nx;
	for (i = 0; i <= i_2; ++i) {
	    k = nx * (j - 1) + i;
	    v = 0.;
	    vr = 0.;
	    vt = 0.;
	    if (j >= 1 && i >= 1) {
		v = x(k);
	    }
	    if (i < nx && j > 0) {
		vr = x(k+1);
	    }
	    if (i > 0 && j < ny) {
		vt = x(k+nx);
	    }
	    dvdx = (vr - v) / hx;
	    dvdy = (vt - v) / hy;
/* Computing 2nd power */
	    d_1 = dvdx;
/* Computing 2nd power */
	    d_2 = dvdy;
	    gradv = d_1 * d_1 + d_2 * d_2;
	    if (feval) {
		dodcps(gradv, 1.0, 2.0, t1, t2, &dpsi, 0, 
			lambda);
		f += dpsi;
	    }
	    if (geval) {
		dodcps(gradv, 1.0, 2.0, t1, t2, &dpsip, 1, 
			lambda);
		if (i >= 1 && j >= 1) {
		    fgrad(k) -= (dvdx / hx + dvdy / hy) * 2. * dpsip;
		}
		if (i < nx && j > 0) {
		    fgrad(k+1) += dvdx / hx * 2. * dpsip;
		}
		if (i > 0 && j < ny) {
		    fgrad(k+nx) += dvdy / hy * 2. * dpsip;
		}
	    }
/* L40: */
	}
/* L50: */
    }
/*     Computation of the function and the gradient over the upper */
/*     triangular elements. */
    i_1 = ny + 1;
    for (j = 1; j <= i_1; ++j) {
	i_2 = nx + 1;
	for (i = 1; i <= i_2; ++i) {
	    k = nx * (j - 1) + i;
	    vb = 0.;
	    vl = 0.;
	    v = 0.;
	    if (i <= nx && j > 1) {
		vb = x(k-nx);
	    }
	    if (i > 1 && j <= ny) {
		vl = x(k-1);
	    }
	    if (i <= nx && j <= ny) {
		v = x(k);
	    }
	    dvdx = (v - vl) / hx;
	    dvdy = (v - vb) / hy;
/* Computing 2nd power */
	    d_1 = dvdx;
/* Computing 2nd power */
	    d_2 = dvdy;
	    gradv = d_1 * d_1 + d_2 * d_2;
	    if (feval) {
		dodcps(gradv, 1.0, 2.0, t1, t2, &dpsi, 0, 
			lambda);
		f += dpsi;
	    }
	    if (geval) {
		dodcps(gradv, 1.0, 2.0, t1, t2, &dpsip, 1, 
			lambda);
		if (i <= nx && j > 1) {
		    fgrad(k-nx) -= dvdy / hy * 2. * dpsip;
		}
		if (i > 1 && j <= ny) {
		    fgrad(k-1) -= dvdx / hx * 2. * dpsip;
		}
		if (i <= nx && j <= ny) {
		    fgrad(k) += (dvdx / hx + dvdy / hy) * 2. * dpsip;
		}
	    }
/* L60: */
	}
/* L70: */
    }
/*     Scale the function. */
    if (feval) {
	f = area * f;
    }
/*     Integrate v over the domain. */
    if (feval) {
	temp = 0.;
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    temp += x(k);
/* L80: */
	}
	f += hxhy * temp;
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = area * fgrad(k) + hxhy;
/* L90: */
	}
    }
  }
void dodcps(double t, double mu1, double mu2, double t1, double t2, 
	    double *result, int option, double lambda)
{
    /* Local variables */
    double sqrtt;

/*     ********** */

/*     This subroutine computes the function psi(t) and the scaled */
/*     functions psi'(t)/t and psi''(t)/t for the optimal design */
/*     with composite materials problem. */

/*     The subroutine statement is */

/*       subroutine dodcps(t,mu1,mu2,t1,t2,result,option,lambda) */

/*     where */

/*       t is a double precision variable. */
/*         On entry t is the variable t */
/*         On exit t is unchanged */

/*       mu1 is a double precision variable. */
/*         On entry mu1 is the reciprocal shear modulus of material 1. */
/*         On exit mu1 is unchanged. */

/*       mu2 is a double precision variable. */
/*         On entry mu2 is the reciprocal shear modulus of material 2. */
/*         On exit mu2 is unchanged. */

/*       t1 is a double precision variable. */
/*         On entry t1 is the first breakpoint. */
/*         On exit t1 is unchanged. */

/*       t2 is a double precision variable. */
/*         On entry t2 is the second breakpoint. */
/*         On exit t2 is unchanged. */

/*       result is a double precision variable. */
/*         On entry result need not be specified. */
/*         On exit result is set according to task. */

/*       option is an int variable. */
/*         On entry option specifies the action of the subroutine: */

/*            if option = 0 then evaluate the function psi(t). */
/*            if option = 1 then evaluate the scaled function psi
'(t)/t. */
/*            if option = 2 then evaluate the scaled function psi''(t)/t. 
*/

/*        On option task is unchanged. */

/*       lambda is a double precision variable */
/*         On entry lambda is the Lagrange multiplier. */
/*         On exit lambda is unchanged. */

/*     MINPACK-2 Project. November 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick. */

/*     ********** */
    sqrtt = sqrt(t);
    if (option == 0) {
	if (sqrtt <= t1) {
	    *result = mu2 * .5 * t;
	} else if (sqrtt > t1 && sqrtt < t2) {
	    *result = mu2 * t1 * sqrtt - lambda * mu1;
	} else if (sqrtt >= t2) {
	    *result = mu1 * .5 * t + lambda * (mu2 - mu1);
	}
    } else if (option == 1) {
	if (sqrtt <= t1) {
	    *result = mu2 * .5;
	} else if (sqrtt > t1 && sqrtt < t2) {
	    *result = mu2 * .5 * t1 / sqrtt;
	} else if (sqrtt >= t2) {
	    *result = mu1 * .5;
	}
    } else if (option == 2) {
	if (sqrtt <= t1) {
	    *result = 0.;
	} else if (sqrtt > t1 && sqrtt < t2) {
	    *result = mu2 * -.25 * t1 / (sqrtt * t);
	} else if (sqrtt >= t2) {
	    *result = 0.;
	}
    }
}

