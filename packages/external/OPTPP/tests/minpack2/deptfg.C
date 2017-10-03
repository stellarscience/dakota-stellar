// Elastic-Plastic Torsion function from MINPACK-2 collection

#ifdef ANSI_HEADERS
#include <cmath>
#include <cstring>
#else
#include <math.h>
#include <string.h>
#endif

#include "OptNewton.h"
#include "NLF.h"
#include "minpack2.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

double d_sign_(double *a, double *b)
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}

void deptfg(int nx, int ny, ColumnVector& x, double &f, ColumnVector& fgrad, 
	    MINPACK_TASK task, double c)
{
/*
 *     Subroutine deptfg 
 *
 *     This subroutine computes the function and gradient of the 
 *     elastic-plastic torsion problem. 
 *     The subroutine statement is 
 *       subroutine deptfg(nx,ny,x,f,fgrad,task,c) 
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
 *             'F'     Evaluate the function at x. 
 *             'G'     Evaluate the gradient at x. 
 *             'FG'    Evaluate the function and the gradient at x. 
 *             'XS'    Set x to the standard starting point xs. 
 *             'XL'    Set x to the lower bound xl. 
 *             'XU'    Set x to the upper bound xu. 
 *         On exit task is unchanged. 
 *       c is a double precision variable. 
 *         On entry c is the angle of twist per unit length. 
 *         On exit c is unchanged. 
 *     MINPACK-2 Project. November 1993. 
 *     Argonne National Laboratory and University of Minnesota. 
 *     Brett M. Averick and Jorge J. More'. */

    /* System generated locals */
    int i_1, i_2, i_3, i_4;
    double d_1, d_2;

    /* Local variables */
    double area, flin, dvdx, dvdy, temp, cdiv3, temp1;
    int i, j, k;
    double v;
    int feval, geval;
    double fquad, vb, hx, hy, vl, vr, vt;

/***************/

    hx = 1. / (double) (nx + 1);
    hy = 1. / (double) (ny + 1);
    area = hx * .5 * hy;
    cdiv3 = c / 3.;
/*     Compute a lower bound for x if task = 'XL' or an upper bound if */
/*     task = 'XU'. */
    if (task == XLower || task == XUpper) {
	if (task == XLower) {
	    temp1 = -1.;
	}
	if (task == XUpper) {
	    temp1 = 1.;
	}
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
		d_1 = min(temp,d_2);
		x(k) = d_sign_(&d_1, &temp1);
/* L10: */
	    }
/* L20: */
	}
	return;
    }
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
		d_1 = (double) min(i_4,i_3) * hx;
		x(k) = min(temp,d_1);
/* L30: */
	    }
/* L40: */
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
/*     or both if taskn = 'FG'. */

    if (feval) {
	fquad = 0.;
	flin = 0.;
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = 0.;
/* L50: */
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
	    if (i >= 1 && j >= 1) {
		v = x(k);
	    }
	    if (i < nx && j > 0) {
		vr = x(k + 1);
	    }
	    if (i > 0 && j < ny) {
		vt = x(k + nx);
	    }
	    dvdx = (vr - v) / hx;
	    dvdy = (vt - v) / hy;
	    if (feval) {
/* Computing 2nd power */
		d_1 = dvdx;
/* Computing 2nd power */
		d_2 = dvdy;
		fquad = fquad + d_1 * d_1 + d_2 * d_2;
		flin -= cdiv3 * (v + vr + vt);
	    }
	    if (geval) {
		if (i != 0 && j != 0) {
		    fgrad(k) = fgrad(k) - dvdx / hx - dvdy / hy - cdiv3;
		}
		if (i != nx && j != 0) {
		    fgrad(k + 1) = fgrad(k + 1) + dvdx / hx - cdiv3;
		}
		if (i != 0 && j != ny) {
		    fgrad(k + nx) = fgrad(k + nx) + dvdy / hy - cdiv3;
		}
	    }
/* L60: */
	}
/* L70: */
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
		vb = x(k - nx);
	    }
	    if (i > 1 && j <= ny) {
		vl = x(k - 1);
	    }
	    if (i <= nx && j <= ny) {
		v = x(k);
	    }
	    dvdx = (v - vl) / hx;
	    dvdy = (v - vb) / hy;
	    if (feval) {
/* Computing 2nd power */
		d_1 = dvdx;
/* Computing 2nd power */
		d_2 = dvdy;
		fquad = fquad + d_1 * d_1 + d_2 * d_2;
		flin -= cdiv3 * (vb + vl + v);
	    }
	    if (geval) {
		if (i != nx + 1 && j != 1) {
		    fgrad(k - nx) = fgrad(k - nx) - dvdy / hy - cdiv3;
		}
		if (i != 1 && j != ny + 1) {
		    fgrad(k - 1) = fgrad(k - 1) - dvdx / hx - cdiv3;
		}
		if (i != nx + 1 && j != ny + 1) {
		    fgrad(k) = fgrad(k) + dvdx / hx + dvdy / hy - cdiv3;
		}
	    }
/* L80: */
	}
/* L90: */
    }
/*     Scale the result. */
    if (feval) {
	f = area * (fquad * .5 + flin);
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = area * fgrad(k);
/* L100: */
	}
    }
}

