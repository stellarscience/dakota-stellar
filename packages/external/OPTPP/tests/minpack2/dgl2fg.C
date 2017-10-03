// Ginzburg-Landau Superconductivity function from MINPACK-2 collection

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

#define true  1
#define false 0

void dgl2fg(int nx, int ny, ColumnVector& x, double &f, ColumnVector& fgrad, 
	    Task task, int vornum)
{
    /* Local variables */
    int i, j, k, itemp;
    void dgl2fc();
    int ctr;

/*     Subroutine dgl2fg 
 *     This subroutine computes the function and gradient of the 
 *     Ginzburg-Landau (2-dimensional) superconductivity problem. 
 *     The subroutine statement is 
 *       subroutine dgl2fg(nx,ny,x,f,fgrad,task,w,vornum) 
 *     where 
 *       nx is an integer variable. 
 *         On entry nx is the number of grid points in the first 
 *            coordinate direction. 
 *         On exit nx is unchanged. 
 *       ny is an integer variable. 
 *         On entry ny is the number of grid points in the second 
 *            coordinate direction. 
 *         On exit ny is unchanged. 
 *       x is a double precision array of dimension 4*nxny. 
 *         On entry x specifies the vector x if task = 'F', 'G', or 'FG'. 
 *            Otherwise x need not be specified. 
 *         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise 
 *            x is set according to task. 
 *       f is a double precision variable. 
 *         On entry f need not be specified. 
 *         On exit f is set to the function evaluated at x if task = 'F' 
 *            or 'FG'. 
 *       fgrad is a double precision array of dimension 4*nxny. 
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
 *       w is a double precision work array of dimension 4*(nx+1)(ny+1). 
 *       vornum is an integer variable. 
 *         On entry vornum specifies the number of vortices. 
 *         On exit vornum is unchanged. 
 *     Subprograms called 
 *       MINPACK-supplied   ...   dgl2fc 
 *     MINPACK-2 Project. November 1993. 
 *     Argonne National Laboratory and University of Minnesota. 
 *     Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright. */

/*     ********** */

    itemp = (nx + 1) * (ny + 1);

/*     Pack work array. */
    ctr = 1;
    for (j = 1; j <= ny; ++j) {
      for (i = 1; i <= nx; ++i) {
        k = (j-1)*nx + i;
        w(ctr) = x(k);
        w(itemp+ctr) = x(nxny+k);
        w(2*itemp+ctr) = x(2*nxny+k);
        w(3*itemp+ctr) = x(3*nxny+k);
        ++ctr;
      }
      w(ctr)         = 0.;
      w(itemp+ctr)   = 0.;
      w(2*itemp+ctr) = 0.;
      w(3*itemp+ctr) = 0.;
      ++ctr;
    }

    dgl2fc(nx, ny, &w(1), &w(itemp + 1), &w(2*itemp + 1), &w(3*itemp + 1), 
           f, &fgrad(1), &fgrad(nx * ny + 1), &fgrad(2*nx*ny + 1), 
           &fgrad(3*nx*ny + 1), task, vornum);

/*     Unpack work array */

    ctr = 1;
    for (j = 1; j <= ny; ++j) {
	for (i = 1; i <= nx; ++i) {
	    k = (j - 1) * nx + i;
            x(k) = w(ctr);
            x(nx*ny+k) = w(itemp+ctr);
            x(2*nx*ny+k) = w(2*itemp+ctr);
            x(3*nx*ny+k) = w(3*itemp+ctr);
	    ++ctr;
	}
	++ctr;
    }
}

void dgl2fc(int nx, int ny, double *x, double *y, double *vpotx, 
	    double *vpoty, double *f, double *gradx, double *grady, 
	    double *gradax, double *graday, Task task, int *vornum)
{
    /* System generated locals */
    int x_dim1, x_offset, y_dim1, y_offset, vpotx_dim1, vpotx_offset, 
	    vpoty_dim1, vpoty_offset, gradx_dim1, gradx_offset, grady_dim1, 
	    grady_offset, gradax_dim1, gradax_offset, graday_dim1, 
	    graday_offset, i_1, i_2;
    double d_1, d_2, d_3;

    /* Builtin functions */
    double sqrt(), atan();
    double sin(), cos();

    /* Local variables */
    double sfac, bave, fkin;
    int i, j;
    double fcond, delsq, x1, x2, sqrtv, fkinx1, fkinx2, fkiny1, 
	    fkiny2, ffield, pi, hx, hy, tkappa, xy, arg, sqn, xpt, ypt, cfac;

    /* Parameter adjustments */
    x_dim1 = nx + 1;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    y_dim1 = nx + 1;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    vpotx_dim1 = nx + 1;
    vpotx_offset = vpotx_dim1 + 1;
    vpotx -= vpotx_offset;
    vpoty_dim1 = nx + 1;
    vpoty_offset = vpoty_dim1 + 1;
    vpoty -= vpoty_offset;
    gradx_dim1 = nx;
    gradx_offset = gradx_dim1 + 1;
    gradx -= gradx_offset;
    grady_dim1 = nx;
    grady_offset = grady_dim1 + 1;
    grady -= grady_offset;
    gradax_dim1 = nx;
    gradax_offset = gradax_dim1 + 1;
    gradax -= gradax_offset;
    graday_dim1 = nx;
    graday_offset = graday_dim1 + 1;
    graday -= graday_offset;

/*     Subroutine dgl2fc 
 *     This subroutine computes the function and gradient of the 
 *     Ginzburg-Landau (2-dimensional) superconductivity problem. 
 *     The subroutine statement is 
 *       subroutine dgl2fc(nx,ny,x,y,vpotx,vpoty,f, 
 *    +                  gradx,grady,gradax,graday,task,vornum) 
 *     where 
 *       nx is an integer variable. 
 *         On entry nx is the number of grid points in the first 
 *            coordinate direction. 
 *         On exit nx is unchanged. 
 *       ny is an integer variable. 
 *         On entry ny is the number of grid points in the second 
 *            coordinate direction. 
 *         On exit ny is unchanged. 
 *       x is a double precision array of dimension nxny. 
 *         On entry x specifies the real part of the order parameter 
 *            if task = 'F', 'G', or 'FG'. 
 *            Otherwise x need not be specified. 
 *         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise 
 *            x is set according to task. 
 *       y is a double precision array of dimension nxny. 
 *         On entry y specifies the imaginary part of the order parameter 
 *            if task = 'F', 'G', or 'FG'. 
 *            Otherwise y need not be specified. 
 *         On exit y is unchanged if task = 'F', 'G', or 'FG'. Otherwise 
 *            y is set according to task. 
 *       vpotx is a double precision array of dimension nxny. 
 *         On entry vpotx specifies the x component of the vector 
 *            potential if task = 'F', 'G', or 'FG'. 
 *            Otherwise vpotx need not be specified. 
 *         On exit vpotx is unchanged if task = 'F', 'G', or 'FG'. 
 *            Otherwise vpotx is set according to task. 
 *       vpoty is a double precision array of dimension nxny. 
 *         On entry vpoty specifies the y component of the vector 
 *            potential if task = 'F', 'G', or 'FG'. 
 *            Otherwise vpoty need not be specified. 
 *         On exit vpoty is unchanged if task = 'F', 'G', or 'FG'. 
 *            Otherwise vpoty is set according to task. 
 *       f is a double precision variable. 
 *         On entry f need not be specified. 
 *         On exit f is set to the function evaluated at x if task = 'F' 
 *            or 'FG'. 
 *       gradx is a double precision array of dimension nxny. 
 *         On entry gradx need not be specified. 
 *         On exit gradx contains the gradient with respect to x 
 *            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'. 
 *       grady is a double precision array of dimension nxny. 
 *         On entry grady need not be specified. 
 *         On exit grady contains the gradient with respect to y 
 *            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'. 
 *            if task = 'G' or 'FG'. 
 *       gradax is a double precision array of dimension nxny. 
 *         On entry gradax need not be specified. 
 *         On exit gradax contains the gradient with respect to vpotx 
 *            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'. 
 *            if task = 'G' or 'FG'. 
 *       graday is a double precision array of dimension nxny. 
 *         On entry graday need not be specified. 
 *         On exit graday contains the gradient with respect to vpoty 
 *            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'. 
 *            if task = 'G' or 'FG'. 
 *       task is a character variable. 
 *         On entry task specifies the action of the subroutine: 
 *            task               action 
 *            ----               ------ 
 *             'F'     Evaluate the function at (x,y,vpotx,vpoty). 
 *             'G'     Evaluate the gradient at (x,y,vpotx,vpoty). 
 *             'FG'    Evaluate the function and the gradient at 
 *                         (x,y,vpotx,vpoty). 
 *             'XS'    Set (x,y,vpotx,vpoty) to the standard starting 
 *                         point xs. 
 *         On exit task is unchanged. 
 *       vornum is an integer variable. 
 *         On entry vornum is the number of vortices. 
 *         On exit vornum is unchanged. 
 *     MINPACK-2 Project. November 1993. 
 *     Argonne National Laboratory and University of Minnesota. 
 *     Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright. */


/*     Initialize. */
    tkappa = 5.;
    hx = sqrt(*vornum / 2.) * 3. / (double) (nx);
    hy = sqrt(*vornum / 2.) * 3. * sqrt(3.) / (double) (ny);
    sqn = (double) (nx * ny);
    pi = atan(1.) * 4.;
    bave = pi * 2. * *vornum * tkappa / (sqn * hx * hy);
    sqrtv = sqrt((double) (*vornum)) * pi;
    if (task == XStandard) {
/*        Initial Order Parameter. */
	i_1 = ny + 1;
	for (j = 1; j <= i_1; ++j) {
	    ypt = ((double) j - 1.) * hy;
	    i_2 = nx + 1;
	    for (i = 1; i <= i_2; ++i) {
		xpt = ((double) i - 1.) * hx;
/* Computing 2nd power */
		d_1 = sin(sqrtv * xpt / 6.) * sin(sqrtv * ypt / (sqrt(3.) * 
			2. * 3.));
		x[i + j * x_dim1] = 1. - d_1 * d_1;
		y[i + j * y_dim1] = 0.;
/* L10: */
	    }
/* L20: */
	}
/*        Initial Vector Potential. */
	i_1 = ny + 1;
	for (j = 1; j <= i_1; ++j) {
	    i_2 = nx + 1;
	    for (i = 1; i <= i_2; ++i) {
		xpt = ((double) i - 1.) * hx;
		vpotx[i + j * vpotx_dim1] = 0.;
		vpoty[i + j * vpoty_dim1] = bave * xpt / tkappa;
/* L30: */
	    }
/* L40: */
	}
	return 0;
    }
/*     Enforce vortex constraint and boundary conditions. */
/*     Right face for order parameter and vector potential. */
    i_1 = ny + 1;
    for (j = 1; j <= i_1; ++j) {
	arg = pi * 2. * *vornum * ((double) j - 1.) / (double) (ny);
	x[nx + 1 + j * x_dim1] = x[j * x_dim1 + 1] * cos(arg) - y[j * y_dim1 
		+ 1] * sin(arg);
	y[nx + 1 + j * y_dim1] = x[j * x_dim1 + 1] * sin(arg) + y[j * y_dim1 
		+ 1] * cos(arg);
	vpotx[nx + 1 + j * vpotx_dim1] = vpotx[j * vpotx_dim1 + 1];
	vpoty[nx + 1 + j * vpoty_dim1] = vpoty[j * vpoty_dim1 + 1] + pi * 2. 
		* *vornum / ((double) (ny) * hy);
/* L50: */
    }
/*     Top face for order parameter and vector potential. */
    i_1 = nx + 1;
    for (i = 1; i <= i_1; ++i) {
	x[i + (ny + 1) * x_dim1] = x[i + x_dim1];
	y[i + (ny + 1) * y_dim1] = y[i + y_dim1];
	vpotx[i + (ny + 1) * vpotx_dim1] = vpotx[i + vpotx_dim1];
	vpoty[i + (ny + 1) * vpoty_dim1] = vpoty[i + vpoty_dim1];
/* L60: */
    }
    if (task == Function || task == FuncGrad) {
/*        Compute the Condensation Energy Density */
	fcond = 0.;
	i_1 = nx;
	for (i = 1; i <= i_1; ++i) {
	    i_2 = ny;
	    for (j = 1; j <= i_2; ++j) {
/* Computing 2nd power */
		d_1 = x[i + j * x_dim1];
/* Computing 2nd power */
		d_2 = y[i + j * y_dim1];
		delsq = d_1 * d_1 + d_2 * d_2;
/* Computing 2nd power */
		d_1 = delsq;
		fcond = fcond - delsq + d_1 * d_1 / 2.;
/* L70: */
	    }
/* L80: */
	}
	fcond /= sqn;
/*        Compute the Kinetic Energy Density. */
	fkin = 0.;
	i_1 = nx;
	for (i = 1; i <= i_1; ++i) {
	    i_2 = ny;
	    for (j = 1; j <= i_2; ++j) {
		x1 = x[i + 1 + j * x_dim1] - x[i + j * x_dim1] * cos(hx * 
			vpotx[i + j * vpotx_dim1]) + y[i + j * y_dim1] * sin(
			hx * vpotx[i + j * vpotx_dim1]);
		x2 = y[i + 1 + j * y_dim1] - y[i + j * y_dim1] * cos(hx * 
			vpotx[i + j * vpotx_dim1]) - x[i + j * x_dim1] * sin(
			hx * vpotx[i + j * vpotx_dim1]);
/* Computing 2nd power */
		d_1 = x1;
/* Computing 2nd power */
		d_2 = x2;
/* Computing 2nd power */
		d_3 = hx;
		fkin += (d_1 * d_1 + d_2 * d_2) / (d_3 * d_3);
		x1 = x[i + (j + 1) * x_dim1] - x[i + j * x_dim1] * cos(hy * 
			vpoty[i + j * vpoty_dim1]) + y[i + j * y_dim1] * sin(
			hy * vpoty[i + j * vpoty_dim1]);
		x2 = y[i + (j + 1) * y_dim1] - y[i + j * y_dim1] * cos(hy * 
			vpoty[i + j * vpoty_dim1]) - x[i + j * x_dim1] * sin(
			hy * vpoty[i + j * vpoty_dim1]);
/* Computing 2nd power */
		d_1 = x1;
/* Computing 2nd power */
		d_2 = x2;
/* Computing 2nd power */
		d_3 = hy;
		fkin += (d_1 * d_1 + d_2 * d_2) / (d_3 * d_3);
/* L90: */
	    }
/* L100: */
	}
	fkin /= sqn;
/*        Compute the Magnetic Field Energy Density. */
	ffield = 0.;
	i_1 = nx;
	for (i = 1; i <= i_1; ++i) {
	    i_2 = ny;
	    for (j = 1; j <= i_2; ++j) {
		xy = (vpoty[i + 1 + j * vpoty_dim1] - vpoty[i + j * 
			vpoty_dim1]) / hx - (vpotx[i + (j + 1) * vpotx_dim1] 
			- vpotx[i + j * vpotx_dim1]) / hy;
/* Computing 2nd power */
		d_1 = xy;
		ffield += d_1 * d_1;
/* L110: */
	    }
/* L120: */
	}
/* Computing 2nd power */
	d_1 = tkappa;
	ffield = ffield * (d_1 * d_1) / sqn;
	*f = fcond + fkin + ffield;
    }
    if (task == Gradient || task == FuncGrad) {
	i_1 = ny;
	for (j = 1; j <= i_1; ++j) {
	    i_2 = nx;
	    for (i = 1; i <= i_2; ++i) {
/* Computing 2nd power */
		d_1 = x[i + j * x_dim1];
/* Computing 2nd power */
		d_2 = y[i + j * y_dim1];
		gradx[i + j * gradx_dim1] = x[i + j * x_dim1] * (d_1 * d_1 - 
			1. + d_2 * d_2);
		gradx[i + j * gradx_dim1] = gradx[i + j * gradx_dim1] * 2. / 
			sqn;
/* Computing 2nd power */
		d_1 = x[i + j * x_dim1];
/* Computing 2nd power */
		d_2 = y[i + j * y_dim1];
		grady[i + j * grady_dim1] = y[i + j * y_dim1] * (d_1 * d_1 - 
			1. + d_2 * d_2);
		grady[i + j * grady_dim1] = grady[i + j * grady_dim1] * 2. / 
			sqn;
		gradax[i + j * gradax_dim1] = 0.;
		graday[i + j * graday_dim1] = 0.;
/* L130: */
	    }
/* L140: */
	}
/*        Kinetic Energy Part, Interior Points */
	i_1 = nx;
	for (i = 2; i <= i_1; ++i) {
	    i_2 = ny;
	    for (j = 2; j <= i_2; ++j) {
		fkinx1 = 2. / (hx * hx * sqn) * (x[i + 1 + j * x_dim1] - x[i 
			+ j * x_dim1] * cos(hx * vpotx[i + j * vpotx_dim1]) + 
			y[i + j * y_dim1] * sin(hx * vpotx[i + j * vpotx_dim1]
			));
		fkinx2 = 2. / (hx * hx * sqn) * (y[i + 1 + j * y_dim1] - y[i 
			+ j * y_dim1] * cos(hx * vpotx[i + j * vpotx_dim1]) - 
			x[i + j * x_dim1] * sin(hx * vpotx[i + j * vpotx_dim1]
			));
		fkiny1 = 2. / (hy * hy * sqn) * (x[i + (j + 1) * x_dim1] - x[
			i + j * x_dim1] * cos(hy * vpoty[i + j * vpoty_dim1]) 
			+ y[i + j * y_dim1] * sin(hy * vpoty[i + j * 
			vpoty_dim1]));
		fkiny2 = 2. / (hy * hy * sqn) * (y[i + (j + 1) * y_dim1] - y[
			i + j * y_dim1] * cos(hy * vpoty[i + j * vpoty_dim1]) 
			- x[i + j * x_dim1] * sin(hy * vpoty[i + j * 
			vpoty_dim1]));
		ffield = (vpotx[i + j * vpotx_dim1] - vpotx[i + (j + 1) * 
			vpotx_dim1]) / hy + (vpoty[i + 1 + j * vpoty_dim1] - 
			vpoty[i + j * vpoty_dim1]) / hx;
/* Computing 2nd power */
		d_1 = tkappa;
		ffield = d_1 * d_1 * 2. / sqn * ffield;
		gradx[i + j * gradx_dim1] = gradx[i + j * gradx_dim1] - cos(
			hx * vpotx[i + j * vpotx_dim1]) * fkinx1 - sin(hx * 
			vpotx[i + j * vpotx_dim1]) * fkinx2 - cos(hy * vpoty[
			i + j * vpoty_dim1]) * fkiny1 - sin(hy * vpoty[i + j *
			 vpoty_dim1]) * fkiny2;
		grady[i + j * grady_dim1] = grady[i + j * grady_dim1] + sin(
			hx * vpotx[i + j * vpotx_dim1]) * fkinx1 - cos(hx * 
			vpotx[i + j * vpotx_dim1]) * fkinx2 + sin(hy * vpoty[
			i + j * vpoty_dim1]) * fkiny1 - cos(hy * vpoty[i + j *
			 vpoty_dim1]) * fkiny2;
		gradax[i + j * gradax_dim1] = gradax[i + j * gradax_dim1] + 
			ffield / hy + fkinx1 * (hx * x[i + j * x_dim1] * sin(
			hx * vpotx[i + j * vpotx_dim1]) + hx * y[i + j * 
			y_dim1] * cos(hx * vpotx[i + j * vpotx_dim1])) + 
			fkinx2 * (hx * y[i + j * y_dim1] * sin(hx * vpotx[i + 
			j * vpotx_dim1]) - hx * x[i + j * x_dim1] * cos(hx * 
			vpotx[i + j * vpotx_dim1]));
		graday[i + j * graday_dim1] = graday[i + j * graday_dim1] - 
			ffield / hx + fkiny1 * (hy * x[i + j * x_dim1] * sin(
			hy * vpoty[i + j * vpoty_dim1]) + hy * y[i + j * 
			y_dim1] * cos(hy * vpoty[i + j * vpoty_dim1])) + 
			fkiny2 * (hy * y[i + j * y_dim1] * sin(hy * vpoty[i + 
			j * vpoty_dim1]) - hy * x[i + j * x_dim1] * cos(hy * 
			vpoty[i + j * vpoty_dim1]));
		fkinx1 = 2. / (hx * hx * sqn) * (x[i + j * x_dim1] - x[i - 1 
			+ j * x_dim1] * cos(hx * vpotx[i - 1 + j * vpotx_dim1]
			) + y[i - 1 + j * y_dim1] * sin(hx * vpotx[i - 1 + j *
			 vpotx_dim1]));
		fkinx2 = 2. / (hx * hx * sqn) * (y[i + j * y_dim1] - y[i - 1 
			+ j * y_dim1] * cos(hx * vpotx[i - 1 + j * vpotx_dim1]
			) - x[i - 1 + j * x_dim1] * sin(hx * vpotx[i - 1 + j *
			 vpotx_dim1]));
		fkiny1 = 2. / (hy * hy * sqn) * (x[i + j * x_dim1] - x[i + (j 
			- 1) * x_dim1] * cos(hy * vpoty[i + (j - 1) * 
			vpoty_dim1]) + y[i + (j - 1) * y_dim1] * sin(hy * 
			vpoty[i + (j - 1) * vpoty_dim1]));
		fkiny2 = 2. / (hy * hy * sqn) * (y[i + j * y_dim1] - y[i + (j 
			- 1) * y_dim1] * cos(hy * vpoty[i + (j - 1) * 
			vpoty_dim1]) - x[i + (j - 1) * x_dim1] * sin(hy * 
			vpoty[i + (j - 1) * vpoty_dim1]));
		gradx[i + j * gradx_dim1] = gradx[i + j * gradx_dim1] + 
			fkinx1 + fkiny1;
		grady[i + j * grady_dim1] = grady[i + j * grady_dim1] + 
			fkinx2 + fkiny2;
		ffield = (vpotx[i + (j - 1) * vpotx_dim1] - vpotx[i + j * 
			vpotx_dim1]) / hy + (vpoty[i + 1 + (j - 1) * 
			vpoty_dim1] - vpoty[i + (j - 1) * vpoty_dim1]) / hx;
/* Computing 2nd power */
		d_1 = tkappa;
		ffield = d_1 * d_1 * 2. / sqn * ffield;
		gradax[i + j * gradax_dim1] -= ffield / hy;
		ffield = (vpotx[i - 1 + j * vpotx_dim1] - vpotx[i - 1 + (j + 
			1) * vpotx_dim1]) / hy + (vpoty[i + j * vpoty_dim1] - 
			vpoty[i - 1 + j * vpoty_dim1]) / hx;
/* Computing 2nd power */
		d_1 = tkappa;
		ffield = d_1 * d_1 * 2. / sqn * ffield;
		graday[i + j * graday_dim1] += ffield / hx;
/* L150: */
	    }
/* L160: */
	}
/*        Kinetic Energy Part, Boundary Points. */
/*        Bottom J = 1 */
	i_1 = nx;
	for (i = 2; i <= i_1; ++i) {
	    fkinx1 = 2. / (hx * hx * sqn) * (x[i + 1 + x_dim1] - x[i + x_dim1]
		     * cos(hx * vpotx[i + vpotx_dim1]) + y[i + y_dim1] * sin(
		    hx * vpotx[i + vpotx_dim1]));
	    fkinx2 = 2. / (hx * hx * sqn) * (y[i + 1 + y_dim1] - y[i + y_dim1]
		     * cos(hx * vpotx[i + vpotx_dim1]) - x[i + x_dim1] * sin(
		    hx * vpotx[i + vpotx_dim1]));
	    fkiny1 = 2. / (hy * hy * sqn) * (x[i + (x_dim1 << 1)] - x[i + 
		    x_dim1] * cos(hy * vpoty[i + vpoty_dim1]) + y[i + y_dim1] 
		    * sin(hy * vpoty[i + vpoty_dim1]));
	    fkiny2 = 2. / (hy * hy * sqn) * (y[i + (y_dim1 << 1)] - y[i + 
		    y_dim1] * cos(hy * vpoty[i + vpoty_dim1]) - x[i + x_dim1] 
		    * sin(hy * vpoty[i + vpoty_dim1]));
	    ffield = (vpotx[i + vpotx_dim1] - vpotx[i + (vpotx_dim1 << 1)]) / 
		    hy + (vpoty[i + 1 + vpoty_dim1] - vpoty[i + vpoty_dim1]) /
		     hx;
/* Computing 2nd power */
	    d_1 = tkappa;
	    ffield = d_1 * d_1 * 2. / sqn * ffield;
	    gradx[i + gradx_dim1] = gradx[i + gradx_dim1] - cos(hx * vpotx[i 
		    + vpotx_dim1]) * fkinx1 - sin(hx * vpotx[i + vpotx_dim1]) 
		    * fkinx2 - cos(hy * vpoty[i + vpoty_dim1]) * fkiny1 - sin(
		    hy * vpoty[i + vpoty_dim1]) * fkiny2;
	    grady[i + grady_dim1] = grady[i + grady_dim1] + sin(hx * vpotx[i 
		    + vpotx_dim1]) * fkinx1 - cos(hx * vpotx[i + vpotx_dim1]) 
		    * fkinx2 + sin(hy * vpoty[i + vpoty_dim1]) * fkiny1 - cos(
		    hy * vpoty[i + vpoty_dim1]) * fkiny2;
	    gradax[i + gradax_dim1] = gradax[i + gradax_dim1] + ffield / hy + 
		    fkinx1 * (hx * x[i + x_dim1] * sin(hx * vpotx[i + 
		    vpotx_dim1]) + hx * y[i + y_dim1] * cos(hx * vpotx[i + 
		    vpotx_dim1])) + fkinx2 * (hx * y[i + y_dim1] * sin(hx * 
		    vpotx[i + vpotx_dim1]) - hx * x[i + x_dim1] * cos(hx * 
		    vpotx[i + vpotx_dim1]));
	    graday[i + graday_dim1] = graday[i + graday_dim1] - ffield / hx + 
		    fkiny1 * (hy * x[i + x_dim1] * sin(hy * vpoty[i + 
		    vpoty_dim1]) + hy * y[i + y_dim1] * cos(hy * vpoty[i + 
		    vpoty_dim1])) + fkiny2 * (hy * y[i + y_dim1] * sin(hy * 
		    vpoty[i + vpoty_dim1]) - hy * x[i + x_dim1] * cos(hy * 
		    vpoty[i + vpoty_dim1]));
	    fkinx1 = 2. / (hx * hx * sqn) * (x[i + x_dim1] - x[i - 1 + x_dim1]
		     * cos(hx * vpotx[i - 1 + vpotx_dim1]) + y[i - 1 + y_dim1]
		     * sin(hx * vpotx[i - 1 + vpotx_dim1]));
	    fkinx2 = 2. / (hx * hx * sqn) * (y[i + y_dim1] - y[i - 1 + y_dim1]
		     * cos(hx * vpotx[i - 1 + vpotx_dim1]) - x[i - 1 + x_dim1]
		     * sin(hx * vpotx[i - 1 + vpotx_dim1]));
	    fkiny1 = 2. / (hy * hy * sqn) * (x[i + (ny + 1) * x_dim1] - x[i 
		    + ny * x_dim1] * cos(hy * vpoty[i + ny * vpoty_dim1]) + 
		    y[i + ny * y_dim1] * sin(hy * vpoty[i + ny * vpoty_dim1]
		    ));
	    fkiny2 = 2. / (hy * hy * sqn) * (y[i + (ny + 1) * y_dim1] - y[i 
		    + ny * y_dim1] * cos(hy * vpoty[i + ny * vpoty_dim1]) - 
		    x[i + ny * x_dim1] * sin(hy * vpoty[i + ny * vpoty_dim1]
		    ));
	    gradx[i + gradx_dim1] = gradx[i + gradx_dim1] + fkinx1 + fkiny1;
	    grady[i + grady_dim1] = grady[i + grady_dim1] + fkinx2 + fkiny2;
	    ffield = (vpotx[i + ny * vpotx_dim1] - vpotx[i + (ny + 1) * 
		    vpotx_dim1]) / hy + (vpoty[i + 1 + ny * vpoty_dim1] - 
		    vpoty[i + ny * vpoty_dim1]) / hx;
/* Computing 2nd power */
	    d_1 = tkappa;
	    ffield = d_1 * d_1 * 2. / sqn * ffield;
	    gradax[i + gradax_dim1] -= ffield / hy;
	    ffield = (vpotx[i - 1 + vpotx_dim1] - vpotx[i - 1 + (vpotx_dim1 <<
		     1)]) / hy + (vpoty[i + vpoty_dim1] - vpoty[i - 1 + 
		    vpoty_dim1]) / hx;
/* Computing 2nd power */
	    d_1 = tkappa;
	    ffield = d_1 * d_1 * 2. / sqn * ffield;
	    graday[i + graday_dim1] += ffield / hx;
/* L170: */
	}
/*        Left I = 1. */
	i_1 = ny;
	for (j = 2; j <= i_1; ++j) {
	    fkinx1 = 2. / (hx * hx * sqn) * (x[j * x_dim1 + 2] - x[j * x_dim1 
		    + 1] * cos(hx * vpotx[j * vpotx_dim1 + 1]) + y[j * y_dim1 
		    + 1] * sin(hx * vpotx[j * vpotx_dim1 + 1]));
	    fkinx2 = 2. / (hx * hx * sqn) * (y[j * y_dim1 + 2] - y[j * y_dim1 
		    + 1] * cos(hx * vpotx[j * vpotx_dim1 + 1]) - x[j * x_dim1 
		    + 1] * sin(hx * vpotx[j * vpotx_dim1 + 1]));
	    fkiny1 = 2. / (hy * hy * sqn) * (x[(j + 1) * x_dim1 + 1] - x[j * 
		    x_dim1 + 1] * cos(hy * vpoty[j * vpoty_dim1 + 1]) + y[j * 
		    y_dim1 + 1] * sin(hy * vpoty[j * vpoty_dim1 + 1]));
	    fkiny2 = 2. / (hy * hy * sqn) * (y[(j + 1) * y_dim1 + 1] - y[j * 
		    y_dim1 + 1] * cos(hy * vpoty[j * vpoty_dim1 + 1]) - x[j * 
		    x_dim1 + 1] * sin(hy * vpoty[j * vpoty_dim1 + 1]));
	    ffield = (vpotx[j * vpotx_dim1 + 1] - vpotx[(j + 1) * vpotx_dim1 
		    + 1]) / hy + (vpoty[j * vpoty_dim1 + 2] - vpoty[j * 
		    vpoty_dim1 + 1]) / hx;
/* Computing 2nd power */
	    d_1 = tkappa;
	    ffield = d_1 * d_1 * 2. / sqn * ffield;
	    gradx[j * gradx_dim1 + 1] = gradx[j * gradx_dim1 + 1] - cos(hx * 
		    vpotx[j * vpotx_dim1 + 1]) * fkinx1 - sin(hx * vpotx[j * 
		    vpotx_dim1 + 1]) * fkinx2 - cos(hy * vpoty[j * vpoty_dim1 
		    + 1]) * fkiny1 - sin(hy * vpoty[j * vpoty_dim1 + 1]) * 
		    fkiny2;
	    grady[j * grady_dim1 + 1] = grady[j * grady_dim1 + 1] + sin(hx * 
		    vpotx[j * vpotx_dim1 + 1]) * fkinx1 - cos(hx * vpotx[j * 
		    vpotx_dim1 + 1]) * fkinx2 + sin(hy * vpoty[j * vpoty_dim1 
		    + 1]) * fkiny1 - cos(hy * vpoty[j * vpoty_dim1 + 1]) * 
		    fkiny2;
	    gradax[j * gradax_dim1 + 1] = gradax[j * gradax_dim1 + 1] + 
		    ffield / hy + fkinx1 * (hx * x[j * x_dim1 + 1] * sin(hx * 
		    vpotx[j * vpotx_dim1 + 1]) + hx * y[j * y_dim1 + 1] * cos(
		    hx * vpotx[j * vpotx_dim1 + 1])) + fkinx2 * (hx * y[j * 
		    y_dim1 + 1] * sin(hx * vpotx[j * vpotx_dim1 + 1]) - hx * 
		    x[j * x_dim1 + 1] * cos(hx * vpotx[j * vpotx_dim1 + 1]));
	    graday[j * graday_dim1 + 1] = graday[j * graday_dim1 + 1] - 
		    ffield / hx + fkiny1 * (hy * x[j * x_dim1 + 1] * sin(hy * 
		    vpoty[j * vpoty_dim1 + 1]) + hy * y[j * y_dim1 + 1] * cos(
		    hy * vpoty[j * vpoty_dim1 + 1])) + fkiny2 * (hy * y[j * 
		    y_dim1 + 1] * sin(hy * vpoty[j * vpoty_dim1 + 1]) - hy * 
		    x[j * x_dim1 + 1] * cos(hy * vpoty[j * vpoty_dim1 + 1]));
	    fkinx1 = 2. / (hx * hx * sqn) * (x[nx + 1 + j * x_dim1] - x[nx 
		    + j * x_dim1] * cos(hx * vpotx[nx + j * vpotx_dim1]) + y[
		    nx + j * y_dim1] * sin(hx * vpotx[nx + j * vpotx_dim1]))
		    ;
	    fkinx2 = 2. / (hx * hx * sqn) * (y[nx + 1 + j * y_dim1] - y[nx 
		    + j * y_dim1] * cos(hx * vpotx[nx + j * vpotx_dim1]) - x[
		    nx + j * x_dim1] * sin(hx * vpotx[nx + j * vpotx_dim1]))
		    ;
	    fkiny1 = 2. / (hy * hy * sqn) * (x[j * x_dim1 + 1] - x[(j - 1) * 
		    x_dim1 + 1] * cos(hy * vpoty[(j - 1) * vpoty_dim1 + 1]) + 
		    y[(j - 1) * y_dim1 + 1] * sin(hy * vpoty[(j - 1) * 
		    vpoty_dim1 + 1]));
	    fkiny2 = 2. / (hy * hy * sqn) * (y[j * y_dim1 + 1] - y[(j - 1) * 
		    y_dim1 + 1] * cos(hy * vpoty[(j - 1) * vpoty_dim1 + 1]) - 
		    x[(j - 1) * x_dim1 + 1] * sin(hy * vpoty[(j - 1) * 
		    vpoty_dim1 + 1]));
	    sfac = sin(pi * 2. * *vornum * (j - 1.) / (double) (ny));
	    cfac = cos(pi * 2. * *vornum * (j - 1.) / (double) (ny));
	    gradx[j * gradx_dim1 + 1] = gradx[j * gradx_dim1 + 1] + cfac * 
		    fkinx1 + sfac * fkinx2 + fkiny1;
	    grady[j * grady_dim1 + 1] = grady[j * grady_dim1 + 1] - sfac * 
		    fkinx1 + cfac * fkinx2 + fkiny2;
	    ffield = (vpotx[(j - 1) * vpotx_dim1 + 1] - vpotx[j * vpotx_dim1 
		    + 1]) / hy + (vpoty[(j - 1) * vpoty_dim1 + 2] - vpoty[(j 
		    - 1) * vpoty_dim1 + 1]) / hx;
/* Computing 2nd power */
	    d_1 = tkappa;
	    ffield = d_1 * d_1 * 2. / sqn * ffield;
	    gradax[j * gradax_dim1 + 1] -= ffield / hy;
	    ffield = (vpotx[nx + j * vpotx_dim1] - vpotx[nx + (j + 1) * 
		    vpotx_dim1]) / hy + (vpoty[nx + 1 + j * vpoty_dim1] - 
		    vpoty[nx + j * vpoty_dim1]) / hx;
/* Computing 2nd power */
	    d_1 = tkappa;
	    ffield = d_1 * d_1 * 2. / sqn * ffield;
	    graday[j * graday_dim1 + 1] += ffield / hx;
/* L180: */
	}
/*        Kinetic Energy Part, at origin (only needed in zero field). 
*/
	fkinx1 = 2. / (hx * hx * sqn) * (x[x_dim1 + 2] - x[x_dim1 + 1] * cos(
		hx * vpotx[vpotx_dim1 + 1]) + y[y_dim1 + 1] * sin(hx * vpotx[
		vpotx_dim1 + 1]));
	fkinx2 = 2. / (hx * hx * sqn) * (y[y_dim1 + 2] - y[y_dim1 + 1] * cos(
		hx * vpotx[vpotx_dim1 + 1]) - x[x_dim1 + 1] * sin(hx * vpotx[
		vpotx_dim1 + 1]));
	fkiny1 = 2. / (hy * hy * sqn) * (x[(x_dim1 << 1) + 1] - x[x_dim1 + 1] 
		* cos(hy * vpoty[vpoty_dim1 + 1]) + y[y_dim1 + 1] * sin(hy * 
		vpoty[vpoty_dim1 + 1]));
	fkiny2 = 2. / (hy * hy * sqn) * (y[(y_dim1 << 1) + 1] - y[y_dim1 + 1] 
		* cos(hy * vpoty[vpoty_dim1 + 1]) - x[x_dim1 + 1] * sin(hy * 
		vpoty[vpoty_dim1 + 1]));
	ffield = (vpotx[vpotx_dim1 + 1] - vpotx[(vpotx_dim1 << 1) + 1]) / hy 
		+ (vpoty[vpoty_dim1 + 2] - vpoty[vpoty_dim1 + 1]) / hx;
/* Computing 2nd power */
	d_1 = tkappa;
	ffield = d_1 * d_1 * 2. / sqn * ffield;
	gradx[gradx_dim1 + 1] = gradx[gradx_dim1 + 1] - cos(hx * vpotx[
		vpotx_dim1 + 1]) * fkinx1 - sin(hx * vpotx[vpotx_dim1 + 1]) * 
		fkinx2 - cos(hy * vpoty[vpoty_dim1 + 1]) * fkiny1 - sin(hy * 
		vpoty[vpoty_dim1 + 1]) * fkiny2;
	grady[grady_dim1 + 1] = grady[grady_dim1 + 1] + sin(hx * vpotx[
		vpotx_dim1 + 1]) * fkinx1 - cos(hx * vpotx[vpotx_dim1 + 1]) * 
		fkinx2 + sin(hy * vpoty[vpoty_dim1 + 1]) * fkiny1 - cos(hy * 
		vpoty[vpoty_dim1 + 1]) * fkiny2;
	gradax[gradax_dim1 + 1] = gradax[gradax_dim1 + 1] + ffield / hy + 
		fkinx1 * (hx * x[x_dim1 + 1] * sin(hx * vpotx[vpotx_dim1 + 1])
		 + hx * y[y_dim1 + 1] * cos(hx * vpotx[vpotx_dim1 + 1])) + 
		fkinx2 * (hx * y[y_dim1 + 1] * sin(hx * vpotx[vpotx_dim1 + 1])
		 - hx * x[x_dim1 + 1] * cos(hx * vpotx[vpotx_dim1 + 1]));
	graday[graday_dim1 + 1] = graday[graday_dim1 + 1] - ffield / hx + 
		fkiny1 * (hy * x[x_dim1 + 1] * sin(hy * vpoty[vpoty_dim1 + 1])
		 + hy * y[y_dim1 + 1] * cos(hy * vpoty[vpoty_dim1 + 1])) + 
		fkiny2 * (hy * y[y_dim1 + 1] * sin(hy * vpoty[vpoty_dim1 + 1])
		 - hy * x[x_dim1 + 1] * cos(hy * vpoty[vpoty_dim1 + 1]));
	fkinx1 = 2. / (hx * hx * sqn) * (x[nx + 1 + x_dim1] - x[nx + x_dim1]
		 * cos(hx * vpotx[nx + vpotx_dim1]) + y[nx + y_dim1] * sin(
		hx * vpotx[nx + vpotx_dim1]));
	fkinx2 = 2. / (hx * hx * sqn) * (y[nx + 1 + y_dim1] - y[nx + y_dim1]
		 * cos(hx * vpotx[nx + vpotx_dim1]) - x[nx + x_dim1] * sin(
		hx * vpotx[nx + vpotx_dim1]));
	fkiny1 = 2. / (hy * hy * sqn) * (x[(ny + 1) * x_dim1 + 1] - x[ny * 
		x_dim1 + 1] * cos(hy * vpoty[ny * vpoty_dim1 + 1]) + y[ny * 
		y_dim1 + 1] * sin(hy * vpoty[ny * vpoty_dim1 + 1]));
	fkiny2 = 2. / (hy * hy * sqn) * (y[(ny + 1) * y_dim1 + 1] - y[ny * 
		y_dim1 + 1] * cos(hy * vpoty[ny * vpoty_dim1 + 1]) - x[ny * 
		x_dim1 + 1] * sin(hy * vpoty[ny * vpoty_dim1 + 1]));
	gradx[gradx_dim1 + 1] = gradx[gradx_dim1 + 1] + fkinx1 + fkiny1;
	grady[grady_dim1 + 1] = grady[grady_dim1 + 1] + fkinx2 + fkiny2;
	ffield = (vpotx[ny * vpotx_dim1 + 1] - vpotx[(ny + 1) * vpotx_dim1 
		+ 1]) / hy + (vpoty[ny * vpoty_dim1 + 2] - vpoty[ny * 
		vpoty_dim1 + 1]) / hx;
/* Computing 2nd power */
	d_1 = tkappa;
	ffield = d_1 * d_1 * 2. / sqn * ffield;
	gradax[gradax_dim1 + 1] -= ffield / hy;
	ffield = (vpotx[nx + vpotx_dim1] - vpotx[nx + (vpotx_dim1 << 1)]) / 
		hy + (vpoty[nx + 1 + vpoty_dim1] - vpoty[nx + vpoty_dim1]) /
		 hx;
/* Computing 2nd power */
	d_1 = tkappa;
	ffield = d_1 * d_1 * 2. / sqn * ffield;
	graday[graday_dim1 + 1] += ffield / hx;
    }
} /* dgl2fc_ */

