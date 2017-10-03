#include "NLP.h"
#include "dims.h"
#include "cblas.h"
#include "pds.h"

using NEWMAT::ColumnVector;

using namespace OPTPP;

extern double alpha;
static int fcn_count = 0;
extern "C" {
extern void get_resid_twafer(int, double *, int *, double *);
}
/* Example file to demonstrate the calling sequence to a 
 * simple NLF1 function
 */
void init_twaf (int ndim, ColumnVector& x)
{
  if (ndim != 4)
  {
    exit (1);
  }
  x(1) = 200.0;
  x(2) = 800.0;
  x(3) = 500.0;
  x(4) = 200.0;
}
void twaf(int n, const ColumnVector& x, double& fx, int& result)
{
  double x1[MAX_NUM_OPT_VARS], resid[MAX_NUM_RESID_POINTS];
  int i, num_nodes;
  
  if (n != 4) return;

  /*  x1 = (double *) calloc(n, DBL_SIZE);
  if (x1 == NULL)
    cerr << "can't allocate x1 for function call"; */

  for (i = 0; i < n; i++)
    x1[i] = x(i+1);

  get_resid_twafer(n, x1, &num_nodes, resid);
  fx = dnrm2(num_nodes, resid, 1);
}
