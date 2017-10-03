// Steady-State Combustion from MINPACK-2 collection


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

void dsscfg(int nx, int ny, ColumnVector& x, double &f, ColumnVector& fgrad, 
	    MINPACK_TASK task, double lambda)
{
/*
 *     Subroutine dsscfg 
 *
 *     This subroutine computes the function and gradient of the 
 *     steady state combustion problem. 
 *     The subroutine statement is 
 *       subroutine dsscfg(nx,ny,x,f,fgrad,task,lambda) 
 *     where 
 *       nx is an int variable. 
 *         On entry nx is the number of grid points in the first 
 *            coordinate direction. 
 *         On exit nx is unchanged. 
 *       ny is an int variable. 
 *         On entry ny is the number of grid points in the second 
 *            coordinate direction. 
 *         On exit ny is unchanged. 
 *       x is a double precision array of dimension nx*ny. 
 *         On entry x specifies the vector x if task = 'F', 'G', or 'FG'. 
 *            Otherwise x need not be specified. 
 *         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise 
 *            x is set according to task. 
 *       f is a double precision variable. 
 *         On entry f need not be specified. 
 *         On exit f is set to the function evaluated at x if task = 'F' 
 *            or 'FG'. 
 *       fgrad is a double precision array of dimension nx*ny. 
 *         On entry fgrad need not be specified. 
 *         On exit fgrad contains the gradient evaluated at x if 
 *            task = 'G' or 'FG'. 
 *       task is a character variable. 
 *         On entry task specifies the action of the subroutine: 
 *            task               action 
 *            ----               ------ 
 *            'F'      Evaluate the function at x. 
 *            'G'      Evaluate the gradient at x. 
 *            'FG'     Evaluate the function and the gradient at x. 
 *            'XS'     Set x to the standard starting point xs. 
 *            'XL'     Set x to the lower bound xl. 
 *            'XU'     Set x to the upper bound xu. 
 *         On exit task is unchanged. 
 *       lambda is a double precision variable. 
 *         On entry lambda is a nonnegative Frank-Kamenetski parameter. 
 *         On exit lambda is unchanged. 
 *     MINPACK-2 Project. November 1993. 
 *     Argonne National Laboratory and University of Minnesota. 
 *     Brett M. Averick. 
 */
    /* System generated locals */
    int i_1, i_2, i_3, i_4;
    double d_1, d_2;

    /* Local variables */
    double area, fexp, dvdx, dvdy, temp, expv, temp1;
    int i, j, k;
    double v;
    int feval, geval;
    double fquad, expvb, expvl, expvr, expvt, vb, hx, hy, vl, vr, 
	    vt;

/***************/

    hx = 1. / (double) (nx + 1);
    hy = 1. / (double) (ny + 1);
    area = hx * .5 * hy;
/*     Compute the standard starting point if task = 'XS'. */
    if (task == XStandard) {
	temp1 = lambda / (lambda + 1.);
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
		d_1 = (double) min(i_4,i_3) * hx;
		x(k) = temp1 * sqrt((min(temp,d_1)));
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
/*     Compute the function if task = 'F', the gradient if task = 'G', or */
/*     both if task = 'FG'. */
    if (feval) {
	fquad = 0.;
	fexp = 0.;
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = 0.;
/* L30: */
	}
    }
/*     Computation of the function and the gradient over the lower */
/*     triangular elements.  The trapezoidal rule is used to estimate */
/*     the integral of the exponential term. */
    i_1 = ny;
    for (j = 0; j <= i_1; ++j) {
	i_2 = nx;
	for (i = 0; i <= i_2; ++i) {
	    k = nx * (j - 1) + i;
	    v = 0.;
	    vr = 0.;
	    vt = 0.;
	    if (i != 0 && j != 0) {
		v = x(k);
	    }
	    if (i != nx && j != 0) {
		vr = x(k + 1);
	    }
	    if (i != 0 && j != ny) {
		vt = x(k + nx);
	    }
	    dvdx = (vr - v) / hx;
	    dvdy = (vt - v) / hy;
	    expv = exp(v);
	    expvr = exp(vr);
	    expvt = exp(vt);
	    if (feval) {
/* Computing 2nd power */
		d_1 = dvdx;
/* Computing 2nd power */
		d_2 = dvdy;
		fquad = fquad + d_1 * d_1 + d_2 * d_2;
		fexp -= lambda * (expv + expvr + expvt) / 3.;
	    }
	    if (geval) {
		if (i != 0 && j != 0) {
		    fgrad(k) = fgrad(k) - dvdx / hx - dvdy / hy - lambda * 
			    expv / 3.;
		}
		if (i != nx && j != 0) {
		    fgrad(k + 1) = fgrad(k + 1) + dvdx / hx - lambda * expvr 
			    / 3.;
		}
		if (i != 0 && j != ny) {
		    fgrad(k + nx) = fgrad(k + nx) + dvdy / hy - lambda * 
			    expvt / 3.;
		}
	    }
/* L40: */
	}
/* L50: */
    }
/*     Computation of the function and the gradient over the upper */
/*     triangular elements.  The trapezoidal rule is used to estimate */
/*     the integral of the exponential term. */
    i_1 = ny + 1;
    for (j = 1; j <= i_1; ++j) {
	i_2 = nx + 1;
	for (i = 1; i <= i_2; ++i) {
	    k = nx * (j - 1) + i;
	    vb = 0.;
	    vl = 0.;
	    v = 0.;
	    if (i != nx + 1 && j != 1) {
		vb = x(k - nx);
	    }
	    if (i != 1 && j != ny + 1) {
		vl = x(k - 1);
	    }
	    if (i != nx + 1 && j != ny + 1) {
		v = x(k);
	    }
	    dvdx = (v - vl) / hx;
	    dvdy = (v - vb) / hy;
	    expvb = exp(vb);
	    expvl = exp(vl);
	    expv = exp(v);
	    if (feval) {
/* Computing 2nd power */
		d_1 = dvdx;
/* Computing 2nd power */
		d_2 = dvdy;
		fquad = fquad + d_1 * d_1 + d_2 * d_2;
		fexp -= lambda * (expvb + expvl + expv) / 3.;
	    }
	    if (geval) {
		if (i != nx + 1 && j != 1) {
		    fgrad(k - nx) = fgrad(k - nx) - dvdy / hy - lambda * 
			    expvb / 3.;
		}
		if (i != 1 && j != ny + 1) {
		    fgrad(k - 1) = fgrad(k - 1) - dvdx / hx - lambda * expvl 
			    / 3.;
		}
		if (i != nx + 1 && j != ny + 1) {
		    fgrad(k) = fgrad(k) + dvdx / hx + dvdy / hy - lambda * 
			    expv / 3.;
		}
	    }
/* L60: */
	}
/* L70: */
    }
/*     Scale the result. */
    if (feval) {
	f = area * (fquad * .5 + fexp);
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = area * fgrad(k);
/* L80: */
	}
    }
}

