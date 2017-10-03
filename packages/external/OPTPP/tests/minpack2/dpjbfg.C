// Journal Bearing function from MINPACK-2 collection


#ifdef ANSI_HEADERS
#include <cstring>
#include <cmath>
#else
#include <string.h>
#include <math.h>
#endif

#include "NLF.h"
#include "minpack2.h"

using NEWMAT::ColumnVector;
using namespace OPTPP;

void dpjbfg(int nx, int ny, ColumnVector& x, double &f, ColumnVector& fgrad, 
	    MINPACK_TASK task, double ecc, double b)
{
/*
 *     Subroutine dpjbfg 
 *     This subroutine computes the function and gradient of the 
 *     pressure distribution in a journal bearing problem. 
 *
 *     The subroutine statement is 
 *
 *       subroutine dpjbfg(nx,ny,x,f,fgrad,task,ecc,b) 
 *
 *     where 
 *
 *       nx is an int variable. 
 *         On entry nx is the number of grid points in the first 
 *            coordinate direction. 
 *         On exit nx is unchanged. 
 *
 *       ny is an int variable. 
 *         On entry ny is the number of grid points in the second 
 *            coordinate direction. 
 *         On exit ny is unchanged. 
 *
 *       x is a double precision array of dimension nxny. 
 *         On entry x specifies the vector x if task = 'F', 'G', or 'FG'. 
 *
 *            Otherwise x need not be specified. 
 *         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise 
 *
 *            x is set according to task. 
 *
 *       f is a double precision variable. 
 *         On entry f need not be specified. 
 *         On exit f is set to the function evaluated at x if task = 'F' 
 *
 *            or 'FG'. 
 *
 *       fgrad is a double precision array of dimension nxny. 
 *         On entry fgrad need not be specified. 
 *         On exit fgrad contains the gradient evaluated at x if 
 *            task = 'G' or 'FG'. 
 *
 *       task is a character variable. 
 *         On entry task specifies the action of the subroutine: 
 *
 *            task               action 
 *            ----               ------ 
 *             'F'     Evaluate the function at x. 
 *             'G'     Evaluate the gradient at x. 
 *             'FG'    Evaluate the function and the gradient at x. 
 *             'XS'    Set x to the standard starting point xs. 
 *             'XL'    Set x to the lower bound xl. 
 *
 *         On exit task is unchanged. 
 *
 *       ecc is a double precision variable 
 *         On entry ecc is the eccentricity in (0,1). 
 *         On exit ecc is unchanged 
 *
 *       b is a double precision variable 
 *         On entry b defines the domain as D = (0,2*pi) X (0,2*b). 
 *         On exit b is unchanged. 
 *
 *     MINPACK-2 Project. November 1993. 
 *     Argonne National Laboratory and University of Minnesota. 
 *     Brett M. Averick and Jorge J. More'. 
 */

    /* System generated locals */
    int i_1, i_2;
    double d_1, d_2, d_3, d_4, d_5, d_6, d_7;

    /* Local variables */
    double flin, dvdx, dvdy, temp, hxhy;
    int i, j, k;
    double v;
    double fquad, ehxhy, trule, vb, pi, hx, hy, vl, xi, vr, vt;
    bool    feval, geval;

/*     Initialization. */
    pi = atan(1.) * 4.;
    hx = pi * 2. / (double) (nx + 1);
    hy = b  * 2. / (double) (ny + 1);
    hxhy = hx * hy;
    ehxhy = ecc * hxhy;
/*     Compute the lower bound xl for x if task = 'XL'. */
    if (task == XLower) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    x(k) = 0.;
	}
	return;
    }
/*     Compute the standard starting point if task = 'XS'. */
    if (task == XStandard) {
	i_1 = nx;
	for (i = 1; i <= i_1; ++i) {
/* Computing MAX */
	    d_1 = sin((double) i * hx);
	    temp = max(0.,d_1);
	    i_2 = ny;
	    for (j = 1; j <= i_2; ++j) {
		k = nx * (j - 1) + i;
		x(k) = temp;
	    }
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
	flin = 0.;
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = 0.;
	}
    }
/*     Computation of the quadratic part of the function and */
/*     corresponding components of the gradient over the */
/*     lower triangular elements. */
    i_1 = nx;
    for (i = 0; i <= i_1; ++i) {
	xi = (double) i * hx;
	d_1 = xi + hx;
/* Computing 3rd power */
	d_2 = ecc * cos(xi) + 1, d_3 = d_2;
/* Computing 3rd power */
	d_4 = ecc * cos(d_1) + 1, d_5 = d_4;
/* Computing 3rd power */
	d_6 = ecc * cos(xi) + 1, d_7 = d_6;
	trule = hxhy * (d_3 * (d_2 * d_2) + d_5 * (d_4 * d_4) + d_7 * (d_6 * 
		d_6)) / 6.;
	i_2 = ny;
	for (j = 0; j <= i_2; ++j) {
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
	    if (feval) {
/* Computing 2nd power */
		d_1 = dvdx;
/* Computing 2nd power */
		d_2 = dvdy;
		fquad += trule * (d_1 * d_1 + d_2 * d_2);
	    }
	    if (geval) {
		if (i != 0 && j != 0) {
		    fgrad(k) = fgrad(k) - trule * (dvdx / hx + dvdy / hy);
		}
		if (i != nx && j != 0) {
		    fgrad(k + 1) = fgrad(k+1) + trule * dvdx / hx;
		}
		if (i != 0 && j != ny) {
		    fgrad(k + nx) = fgrad(k+nx) + trule * dvdy / hy;
		}
	    }
	}
    }
/*     Computation of the quadratic part of the function and */
/*     corresponding components of the gradient over the upper */
/*     triangular elements. */
    i_1 = nx + 1;
    for (i = 1; i <= i_1; ++i) {
	xi = (double) i * hx;
	d_1 = xi - hx;
/* Computing 3rd power */
	d_2 = ecc * cos(xi) + 1, d_3 = d_2;
/* Computing 3rd power */
	d_4 = ecc * cos(d_1) + 1, d_5 = d_4;
/* Computing 3rd power */
	d_6 = ecc * cos(xi) + 1, d_7 = d_6;
	trule = hxhy * (d_3 * (d_2 * d_2) + d_5 * (d_4 * d_4) + d_7 * (d_6 * 
		d_6)) / 6.;
	i_2 = ny + 1;
	for (j = 1; j <= i_2; ++j) {
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
	    if (feval) {
/* Computing 2nd power */
		d_1 = dvdx;
/* Computing 2nd power */
		d_2 = dvdy;
		fquad += trule * (d_1 * d_1 + d_2 * d_2);
	    }
	    if (geval) {
		if (i <= nx && j > 1) {
		    fgrad(k - nx) = fgrad(k - nx) - trule * dvdy / hy;
		}
		if (i > 1 && j <= ny) {
		    fgrad(k - 1) = fgrad(k - 1) - trule * dvdx / hx;
		}
		if (i <= nx && j <= ny) {
		    fgrad(k) = fgrad(k) +trule * (dvdx / hx + dvdy / hy);
		}
	    }
	}
    }
/*     Computation of the linear part of the function and */
/*     corresponding components of the gradient. */
    i_1 = nx;
    for (i = 1; i <= i_1; ++i) {
	temp = sin((double) i * hx);
	if (feval) {
	    i_2 = ny;
	    for (j = 1; j <= i_2; ++j) {
		k = nx * (j - 1) + i;
		flin += temp * x(k);
	    }
	}
	if (geval) {
	    i_2 = ny;
	    for (j = 1; j <= i_2; ++j) {
		k = nx * (j - 1) + i;
		fgrad(k) = fgrad(k) - ehxhy * temp;
	    }
	}
    }
/*     Finish off the function. */
    if (feval) {
	f = fquad * .5 - ehxhy * flin;
    }
}

