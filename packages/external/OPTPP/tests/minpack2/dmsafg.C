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

void dmsafg(int nx, int ny, ColumnVector& x, double &f, 
	    ColumnVector& fgrad, 
	    MINPACK_TASK task, ColumnVector& bottom, ColumnVector& top, 
            ColumnVector& left, ColumnVector& right)

{
    /* System generated locals */
    int i_1, i_2;
    double d_1, d_2;

    /* Local variables */
    double area, dvdx, dvdy;
    int i, j, k;
    double betai, v;
    int feval, geval;
    double xline, yline, fl, vb, fu, alphaj, hx, hy, vl, vr, vt;


/*     Subroutine dmsafg 

 *     This subroutine computes the function and gradient of the 
 *     minimal surface area problem. 

 *     The subroutine statement is 

 *       subroutine dmsafg(nx,ny,x,f,fgrad,task,bottom,top,left,right) 

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

 *         On exit task is unchanged. 

 *       bottom is a double precision array of dimension nx + 2. 
 *         On entry bottom must contain boundary data beginning 
 *            with the lower left corner of the domain. 
 *         On exit bottom is unchanged. 

 *       top is a double precision array of dimension nx + 2. 
 *         On entry top must contain boundary data beginning with 
 *            the upper left corner of the domain. 
 *         On exit top is unchanged. 

 *       left is a double precision array of dimension ny + 2. 
 *         On entry left must contain boundary data beginning with 
 *            the lower left corner of the domain. 
 *         On exit left is unchanged. 

 *       right is a double precision array of dimension ny + 2. 
 *         On entry right must contain boundary data beginning with 
 *            the lower right corner of the domain. 
 *         On exit right is unchanged. 

 *     MINPACK-2 Project. November 1993. 
 *     Argonne National Laboratory and University of Minnesota. 
 *     Brett M. Averick. */

/*     ********** */
/*     Initialize. */
    hx = 1. / (double) (nx + 1);
    hy = 1. / (double) (ny + 1);
    area = hx * .5 * hy;
/*     Compute the standard starting point if task = 'XS'. */
    if (task == XStandard) {
	i_1 = ny;
	for (j = 1; j <= i_1; ++j) {
	    alphaj = (double) j * hy;
	    i_2 = nx;
	    for (i = 1; i <= i_2; ++i) {
		k = nx * (j - 1) + i;
		betai = (double) i * hx;
		yline = alphaj * top(i + 1) + (1. - alphaj) * bottom(i + 1);
		xline = betai * right(j + 1) + (1. - betai) * left(j + 1);
		x(k) = (yline + xline) / 2.;
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
/*     Computation of the function and gradient over the lower */
/*     triangular elements. */
    i_1 = ny;
    for (j = 0; j <= i_1; ++j) {
	i_2 = nx;
	for (i = 0; i <= i_2; ++i) {
	    k = nx * (j - 1) + i;
	    if (i >= 1 && j >= 1) {
		v = x(k);
	    } else {
		if (j == 0) {
		    v = bottom(i + 1);
		}
		if (i == 0) {
		    v = left(j + 1);
		}
	    }
	    if (i < nx && j > 0) {
		vr = x(k + 1);
	    } else {
		if (i == nx) {
		    vr = right(j + 1);
		}
		if (j == 0) {
		    vr = bottom(i + 2);
		}
	    }
	    if (i > 0 && j < ny) {
		vt = x(k + nx);
	    } else {
		if (i == 0) {
		    vt = left(j + 2);
		}
		if (j == ny) {
		    vt = top(i + 1);
		}
	    }
	    dvdx = (vr - v) / hx;
	    dvdy = (vt - v) / hy;
/* Computing 2nd power */
	    d_1 = dvdx;
/* Computing 2nd power */
	    d_2 = dvdy;
	    fl = sqrt(d_1 * d_1 + 1. + d_2 * d_2);
	    if (feval) {
		f += fl;
	    }
	    if (geval) {
		if (i >= 1 && j >= 1) {
		    fgrad(k) -= (dvdx / hx + dvdy / hy) / fl;
		}
		if (i < nx && j > 0) {
		    fgrad(k + 1) += dvdx / hx / fl;
		}
		if (i > 0 && j < ny) {
		    fgrad(k + nx) += dvdy / hy / fl;
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
	    if (i <= nx && j > 1) {
		vb = x(k - nx);
	    } else {
		if (j == 1) {
		    vb = bottom(i + 1);
		}
		if (i == nx + 1) {
		    vb = right(j);
		}
	    }
	    if (i > 1 && j <= ny) {
		vl = x(k - 1);
	    } else {
		if (j == ny + 1) {
		    vl = top(i);
		}
		if (i == 1) {
		    vl = left(j + 1);
		}
	    }
	    if (i <= nx && j <= ny) {
		v = x(k);
	    } else {
		if (i == nx + 1) {
		    v = right(j + 1);
		}
		if (j == ny + 1) {
		    v = top(i + 1);
		}
	    }
	    dvdx = (v - vl) / hx;
	    dvdy = (v - vb) / hy;
/* Computing 2nd power */
	    d_1 = dvdx;
/* Computing 2nd power */
	    d_2 = dvdy;
	    fu = sqrt(d_1 * d_1 + 1. + d_2 * d_2);
	    if (feval) {
		f += fu;
	    }
	    if (geval) {
		if (i <= nx && j > 1) {
		    fgrad(k - nx) -= dvdy / hy / fu;
		}
		if (i > 1 && j <= ny) {
		    fgrad(k - 1) -= dvdx / hx / fu;
		}
		if (i <= nx && j <= ny) {
		    fgrad(k) += (dvdx / hx + dvdy / hy) / fu;
		}
	    }
/* L60: */
	}
/* L70: */
    }
/*     Scale the function and the gradient. */
    if (feval) {
	f = area * f;
    }
    if (geval) {
	i_1 = nx * ny;
	for (k = 1; k <= i_1; ++k) {
	    fgrad(k) = area * fgrad(k);
/* L80: */
	}
    }
}
void dmsabc(int nx, int ny, ColumnVector& bottom, ColumnVector& top, 
            ColumnVector& left, ColumnVector& right)
{
    /* System generated locals */
    int i_1;
    double d_1, d_2, d_3;

    /* Local variables */
    ColumnVector njac(4)	/* was [2][2] */;
    int i, j, k;
    ColumnVector u(2), nf(2);
    int limit;
    double fnorm, hx, hy, xt, yt, det;

    /* Function Body */
/*     ********** */

/*     Subroutine dmsabc */

/*     This subroutine computes Enneper's boundary conditions for the */
/*     minimal surface area problem on the unit square centered at the */
/*     origin. */

/*     The subroutine statement is */

/*       subroutine dmsabc(nx,ny,hx,hy,bottom,top,left,right) */

/*     where */

/*       nx is an int variable. */
/*         On entry nx is the number of grid points in the first */
/*            coordinate direction. */
/*         On exit nx is unchanged. */

/*       ny is an int variable. */
/*         On entry ny is the number of grid points in the second */
/*            coordinate direction. */
/*         On exit ny is unchanged. */

/*       bottom is a double precision array of dimension nx + 2. */
/*         On entry bottom need not be specified. */
/*         On exit bottom contains boundary values for the bottom */
/*            boundary of the domain. */

/*       top is a double precision array of dimension nx + 2. */
/*         On entry top need not be specified. */
/*         On exit top contains boundary values for the top boundary of */

/*            the domain. */
/*       left is a double precision array of dimension ny + 2. */
/*         On entry left need not be specified. */
/*         On exit left contains boundary values for the left boundary */
/*            of the domain. */

/*       right is a double precision array of dimension ny + 2. */
/*         On entry right need not be specified. */
/*         On exit right contains boundary values for the right boundary 
*/
/*            of the domain. */

/*     MINPACK-2 Project. November 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick. */

/*     ********** */
/*     Compute Enneper's boundary conditions: bottom, top, left, then */
/*     right.  Enneper's boundary values are obtained by defining */
/*     bv(x,y) = u**2 - v**2 where u and v are the unique solutions of */
/*     x = u + u*(v**2) - (u**3)/3, y = -v - (u**2)*v + (v**3)/3. */

    hx = 1. / (double) (nx + 1);
    hy = 1. / (double) (ny + 1);
    for (j = 1; j <= 4; ++j) {
	if (j == 1) {
	    yt = -.5;
	    xt = -.5;
	    limit = nx + 2;
	} else if (j == 2) {
	    yt = .5;
	    xt = -.5;
	    limit = nx + 2;
	} else if (j == 3) {
	    yt = -.5;
	    xt = -.5;
	    limit = ny + 2;
	} else if (j == 4) {
	    yt = -.5;
	    xt = .5;
	    limit = ny + 2;
	}
/*        Use Newton's method to solve xt = u + u*(v**2) - (u**3)/3, */
/*        yt = -v - (u**2)*v + (v**3)/3. */
	i_1 = limit;
	for (i = 1; i <= i_1; ++i) {
	    u(1) = xt;
	    u(2) = -yt;
	    for (k = 1; k <= 5; ++k) {
              /* Computing 2nd power */
		d_1 = u(1);
                /* Computing 3rd power */
		d_2 = u(1), d_3 = d_2;
		nf(1) = u(1) + u(1) * (d_1 * d_1) - d_3 * (d_2 * d_2) / 3. - 
			xt;
                /* Computing 2nd power */
		d_1 = u(1);
                /* Computing 3rd power */
		d_2 = u(2), d_3 = d_2;
		nf(2) = -u(2) - d_1 * d_1 * u(2) + d_3 * (d_2 * d_2) / 3. - 
			yt;
		fnorm = sqrt(nf(1) * nf(1) + nf(2) * nf(2));
		if (fnorm <= 1e-10) {
		    break;
		}
                /* Computing 2nd power */
		d_1 = u(2);
                /* Computing 2nd power */
		d_2 = u(1);
		njac(1) = d_1 * d_1 + 1. - d_2 * d_2;
		njac(3) = u(1) * 2. * u(2);
		njac(2) = u(1) * -2. * u(2);
                /* Computing 2nd power */
		d_1 = u(1);
                /* Computing 2nd power */
		d_2 = u(2);
		njac(4) = -1. - d_1 * d_1 + d_2 * d_2;
		det = njac(1) * njac(4) - njac(3) * njac(2);
		u(1) -= (njac(4) * nf(1) - njac(3) * nf(2)) / det;
		u(2) -= (njac(1) * nf(2) - njac(2) * nf(1)) / det;
	    }

            /* Converged */
	    if (j == 1) {
		bottom(i) = u(1) * u(1) - u(2) * u(2);
		xt += hx;
	    } else if (j == 2) {
		top(i) = u(1) * u(1) - u(2) * u(2);
		xt += hx;
	    } else if (j == 3) {
		left(i) = u(1) * u(1) - u(2) * u(2);
		yt += hy;
	    } else if (j == 4) {
		right(i) = u(1) * u(1) - u(2) * u(2);
		yt += hy;
	    }
	}
    }
}



