/*********************************************************************
*   File:  lbfgs.c
*
*   The source code in this file minimizes an unconstrained objective
*   using a limited memory BFGS approximation for the Hessian.
*   The algorithm is globalized by trust regions, and each trust
*   region subproblem is solved approximately by Powell's dogleg
*   method.  It must be possible to obtain function and gradient
*   information for the objective.
*
*   The routines in this file use the following external routines:
*      potential              - energy.c
*      force_to_grad0         - coorsw_fix.c
*      daxpy                  - blas1.c
*      dcopy                  - blas1.c
*      ddot                   - blas1.c
*      dnrm2                  - blas1.c
*      dscal                  - blas1.c
*
*      updtbh_                - etr/baseline/etr_subs/lbfgs_cmp.f
*      multbv_                - etr/baseline/etr_subs/lbfgs_cmp.f
*      multhv_                - etr/baseline/etr_subs/lbfgs_cmp.f
*
*   Code development:
*      12 Oct 94 - Originated by T. Plantenga, Sandia National Labs.
**********************************************************************/

#include "cblas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <malloc.h>

typedef double real;

#define T_LBFGS  5            /*-- Number of L-BFGS vectors to store */
#define MMIN     1.0e-300
#define MCHEPS   1.1e-16
#define DEBUG


/*-- External function declarations. */
extern void daxpy(int, real, real*, int, real*, int);
extern void dcopy (int, real*, int, real*, int);
extern real ddot(int, real*, int, real*, int);
extern real dnrm2(int, real*, int);
extern void dscal (int, real, real*, int);

void Dogleg (int, real*, real, real, int, int, real*, real*,
             real*, real*, real*);


/*-- External function declaration for FORTRAN subroutines. */
extern void updtbh_ (int*, int*, real*, real*, int*, real*, real*);
extern void multbv_ (int*, int*, real*, real*, int*, real*, real*);
extern void multhv_ (int*, int*, real*, real*, int*, real*, real*);

/*-- Internal function declarations. */
void Update_lbfgs (int, int, real*, real*, real*, int*, real*, real*);

extern int iupdate;

char *status_file = {"lbfgs.status"};
FILE *bfgs_fp;


/********************************************************************/
void lbfgs
   (int          n,              /* I  num unknowns = 3 * num atoms */
    real         stop_tol,       /* I  tol for ||g||/sqrt(n)        */
    int          itmax,          /* I  max num iterations allowed   */
    int          itmax_line,
    int*         iter,           /* IO iters required to find min   */
    real*        fret,           /*  O minimum value                */
    int          iprint,
    int          last_call)      /*    not used                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*   Called by exec_minimization in "options.c".
*   This routine is modelled on mm_nlcg in "nlcg.c", by JC Meza.
*   It finds the minimum of the unconstrained molecular potential
*   energy using trust region methods and a limited memory BFGS
*   approximation of the Hessian.  Each trust region subproblem is
*   solved by Powell's dogleg method.
*
*   The L-BFGS structures are currently FORTRAN subroutines that
*   require memory allocation of CMPS and CMPY for work space.
*   These variables must be passed to the FORTRAN subroutines and
*   not altered anywhere else.
*
*   The outline of the algorithm is as follows:
*      Allocate memory, initialize parameters
*      Compute the gradient g_vec
*      LOOP
*        IF ||g_vec|| < tol THEN RETURN
*        Update L-BFGS matrices B and H
*        Compute the dogleg step d_vec
*        Compute the predicted reduction of the quadratic model
*        Compute the actual potential energy reduction
*        IF ared/pred > eta
*          THEN x_vec = x_vec + d_vec
*               TR_size >= TR_size
*               Compute the gradient g_vec at the new point
*          ELSE TR_size < ||d_vec||
*      CONTINUE
*********************************************************************/
{
  real         *x_vec, *new_x_vec, *g_vec, *d_vec, *y_vec;
  int          i;
  int          iter_num, NUPDT;
  real         obj_val, new_obj_val;
  real         dd1, dd2;
  real         eta, pred, ared;
  real         TR_size;
  real         gnorm, xnorm, dnorm;
  real         *CMPS, *CMPY, *tmp_vec;


/*-- Open the status file for saving output. */
  bfgs_fp = fopen (status_file, "w");
  if (bfgs_fp == NULL) {
    fprintf (stderr, "lbfgs: cannot open %s\n", status_file);
    printf ("*** lbfgs: cannot open %s\n", status_file);
    exit (1);
  }

/*--------------------------------------------------------------------
 *   Allocate memory and set up files.
 *-------------------------------------------------------------------*/

  fprintf (stderr, "Allocate space in lbfgs n=%d\n",n);
  x_vec     = (real *) calloc ((n+1) , sizeof(real));
  new_x_vec = (real *) calloc ((n+1) , sizeof(real));
  g_vec     = (real *) calloc ((n+1) , sizeof(real));
  d_vec     = (real *) calloc ((n+1) , sizeof(real));
  y_vec     = (real *) calloc ((n+1) , sizeof(real));
  CMPS      = (real *) calloc ((n*T_LBFGS) , sizeof(real));
  CMPY      = (real *) calloc ((n*T_LBFGS) , sizeof(real));
  tmp_vec   = (real *) calloc ((n+1) , sizeof(real));

/*--------------------------------------------------------------------
 *   Evaluate the objective and its gradient at the start point.
 *-------------------------------------------------------------------*/

/*  obj_val = potential (str, coor, force);
  force_to_grad0 (str, force, g_vec);
 */

/*--------------------------------------------------------------------
 *   Initialize the trust region algorithm parameters.
 *-------------------------------------------------------------------*/

  eta = 0.3;                   /*-- step acceptance threshold */

  TR_size = 1.0;

  iter_num = 0;                /*-- number inner iterations   */
  NUPDT    = 1;                /*-- number L-BFGS B updates   */

  gnorm = dnrm2 (n, g_vec, 1);

  fprintf(bfgs_fp, "\n\t\t Steepest descent with Trust Regions\n");
  fprintf(bfgs_fp, "     Iter     f(x)       ||grad||/n      Delta      ||step||       ared         pred\n");
  fprintf(bfgs_fp,"    %5d %12.4e %12.4e %12.4e\n",
          iter_num, obj_val, gnorm/sqrt(n), TR_size);

/*--------------------------------------------------------------------
 *   Begin the main loop.
 *-------------------------------------------------------------------*/

  while (iter_num < itmax) {

    iter_num++;

    dd1 = g_vec[0];
    for (i=1; i<n; i++)
      if (g_vec[i] > dd1)  dd1 = g_vec[i];
    if (dd1 <= stop_tol) {
/*-- Exit main loop IF ||g||_2 sufficiently small. */
/*    if (gnorm <= (stop_tol * sqrt(n))) {*/
      xnorm = dnrm2 (n, x_vec, 1);
      fprintf (stderr,
               " step_size: %5d   energy: %12.4f  |Grad|: %10.4f |X|: %10.4f\n",
               iter_num, obj_val, gnorm/sqrt((real)(n)),
               xnorm/sqrt((real)(n)));
      break;
    }

/*-- Solve the dogleg subproblem. */
    Dogleg (n, g_vec, gnorm, TR_size, NUPDT, T_LBFGS, CMPS, CMPY,
            tmp_vec, d_vec, &dnorm);

/*-- Compute pred = - g'd - 1/2 d'Bd. */

    dd1 = ddot (n, g_vec, 1, d_vec, 1);
    i = T_LBFGS;
    multbv_ (&n, &i, d_vec, tmp_vec, &NUPDT, CMPS, CMPY);
    dd2 = ddot (n, d_vec, 1, tmp_vec, 1);

    pred = -dd1 - (0.5 * dd2);

    if (pred < sqrt(MMIN)) {
      printf ("*** Predicted reduction not positive.  <lbfgs>\n");
      printf ("    pred = %15.6e\n", pred);
      exit (1);
    }

/*-- Compute ared by evaluating the objective at the trial point. */

    dcopy (n, x_vec, 1, new_x_vec, 1);
    daxpy (n, 1.0, d_vec, 1, new_x_vec, 1);

/*  call evalf
    p_to_coor0 (str, new_x_vec, coor);
    new_obj_val = potential (str, coor, force);
*/
    ared = obj_val - new_obj_val;

/*-- Decide whether to take the step.
 *-- If yes, then increase TR_size, compute the gradient, and
 *-- update the L-BFGS approximations.
 *-- If no, then decrease TR_size to a fraction of ||d||_2. */

    if ((ared / pred) >= eta) {

/*-- Increase the trust region size. */
      if (ared / pred >= 0.9) {
        dd1 = 10.0 * dnorm;
        if (dd1 > TR_size) TR_size = dd1;
      } else {
        dd1 = 2.0 * dnorm;
        if (dd1 > TR_size) TR_size = dd1;
      }
      if (TR_size > 1.0e3)  TR_size = 1.0e3;

/*-- Get the gradient from the previously calculated force vector.
 *-- Set y_vec = new gradient - old gradient. */
      dcopy (n, g_vec, 1, y_vec, 1);
      dscal (n, -1.0, y_vec, 1);
/* get grad
      force_to_grad0 (str, force, g_vec);
 */
      daxpy (n, 1.0, g_vec, 1, y_vec, 1);

/*-- Update the L-BFGS and inverse L-BFGS approximations. */
      Update_lbfgs (n, T_LBFGS, d_vec, y_vec, tmp_vec, &NUPDT, CMPS, CMPY);

      dcopy (n, new_x_vec, 1, x_vec, 1);

      gnorm = dnrm2 (n, g_vec, 1);
      obj_val = new_obj_val;

/*      p_to_coor0 (str, x_vec, coor); */
      fprintf(bfgs_fp,"    %5d %12.4e %12.4e %12.4e %12.4e %12.3e %12.3e\n",
              iter_num, obj_val, gnorm/sqrt(n), TR_size, dnorm, ared, pred);

    } else {

/*-- Decrease the trust region size by linear interpolation. */
      dd1 = (1.0 - eta) / (1.0 - (ared / pred));
      if (dd1 < 0.1)  dd1 = 0.1;
      if (dd1 > 0.5)  dd1 = 0.5;
      TR_size = dd1 * dnorm;

      fprintf(bfgs_fp,"rej %5d %12.4e %12.4e %12.4e %12.4e %12.3e %12.3e\n",
              iter_num, obj_val, gnorm/sqrt(n), TR_size, dnorm, ared, pred);

      if (TR_size < (100.0 * MCHEPS)) {
        printf ("*** Trust region too small to continue.  <lbfgs>\n");
        fprintf (bfgs_fp, "*** Trust region too small to continue.\n");
        fflush (bfgs_fp);
        exit (1);
      }
    }

  }  /***  end of while(iter_num < itmax)  ***/

/*--------------------------------------------------------------------
 *   Clean up and exit.
 *-------------------------------------------------------------------*/

  *iter = iter_num;
  *fret = obj_val;

  free (x_vec);
  free (new_x_vec);
  free (g_vec);
  free (d_vec);

  return;
}

/********************************************************************/
void Update_lbfgs
   (int   n,               /* I  num unknowns              */
    int   t,               /* I  num vectors to store      */
    real  *s_vec,          /* I  x_new - x_old             */
    real  *y_vec,          /* I  g_new - g_old             */
    real  *Bs_vec,         /* IO scratch space             */
    int   *NUPDT,          /* IO num updates made to B/H   */
    real  *CMPS,           /* IO storage for t s_vecs      */
    real  *CMPY)           /* IO storage for t y_vecs      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*   This routine updates the limited memory BFGS approximations for
*   the Hessian B and its inverse H.  The update uses s_vec and
*   y_vec, whose product should be negative to preserve the positive
*   definiteness of the BFGS approximations.  If s'y is not a descent
*   direction, then y is "damped", an idea due to Powell.
*   Then y_vec becomes:
*        y_damped = alpha*y + (1-alpha)*Bs
*
*                                      0.8*s'Bs
*                       where alpha = ----------
*                                     s'Bs - s'y
*
*   If 0 <= s'y <= 1.0E-8 y'y then the update is skipped to prevent
*   division by a small number.
*
*   The L-BFGS structures are currently FORTRAN subroutines that
*   require the variables NUPDT, CMPS, and CMPY.  These are used as
*   work space storage and should not be altered between FORTRAN calls.
*********************************************************************/
{
  real  s_dot_y, sBs, y_dot_y, alpha;


/*-- Damp the estimate if necessary. */

  s_dot_y = ddot (n, s_vec, 1, y_vec, 1);
  multbv_ (&n, &t, s_vec, Bs_vec, NUPDT, CMPS, CMPY);
  sBs = ddot (n, s_vec, 1, Bs_vec, 1);

  if (s_dot_y < (0.2 * sBs)) {
    fprintf (bfgs_fp, "--- damping L-BFGS update\n");
    alpha = 0.8 * sBs / (sBs - s_dot_y);
    dscal (n, alpha, y_vec, 1);
    daxpy (n, (1.0 - alpha), Bs_vec, 1, y_vec, 1);
    s_dot_y = ddot (n, s_vec, 1, y_vec, 1);
  }

/*-- Decide whether to skip the update. */

  y_dot_y = ddot (n, y_vec, 1, y_vec, 1);
  if ((s_dot_y >= 0.0) && (s_dot_y <= (sqrt(MCHEPS) * y_dot_y))) {
    fprintf (bfgs_fp, "--- skipping L-BFGS update\n");
    return;
  }

/*-- Make the updates. */
  updtbh_ (&n, &t, s_vec, y_vec, NUPDT, CMPS, CMPY);

  return;
}

