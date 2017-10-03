#ifndef tstfcn_h
#define tstfcn_h

/*
 *
 * Declarations for TWAFER
 */

/* Initializer for TWAFER */

void init_twaf(int n, ColumnVector& x);

/* TWAFER with no analytic derivative */

void twaf(int n, const ColumnVector& x, double& fx, int& result);

/*
 *
 * Declarations for rosenbrock's functions
 */

/* Initializer for Rosenbrock */

void init_rosen(int n, ColumnVector& x);

/* Rosenbrock with no analytic derivative */

void rosen0(int n, const ColumnVector& x, double& fx, int& result);

/* Rosenbrock with analytic derivative */
void rosen(int mode, int n, const ColumnVector& x, double& fx, ColumnVector& g,
	   int& result);

/* Rosenbrock with analytic derivative and Hessian*/

void rosen2(int mode, int n, const ColumnVector& x, double& fx, 
	    ColumnVector& g, SymmetricMatrix& H, int& result);

/* Scaled version of Rosenbrock with analytic derivative */
void srosen(int mode, int n, const ColumnVector& x, double& fx, 
	    ColumnVector& g, int& result);

/* Scaled version of Rosenbrock with analytic derivative and Hessian*/

void srosen2(int mode, int n, const ColumnVector& x, double& fx, 
	     ColumnVector& g, SymmetricMatrix& H, int& result);

void   init_illum(int n, ColumnVector& x);

void   illum(int, int, const ColumnVector&, double&, ColumnVector&, int&);

void   illum2(int mode, int n, const ColumnVector& x, double& fx, 
	      ColumnVector& g, SymmetricMatrix& H, int& result);

void init_erosen (int ndim, ColumnVector& x);
void erosen(int ndim, const ColumnVector& x, double& fx, int& result);
void erosen1(int mode, int ndim, const ColumnVector& x, double& fx, 
	     ColumnVector&, int& result);

void init_penalty1 (int ndim, ColumnVector& x);
void penalty1(int ndim, const ColumnVector& x, double& fx, int& result);

void init_penalty2 (int ndim, ColumnVector& x);
void penalty2(int ndim, const ColumnVector& x, double& fx, int& result);

void init_epowell (int ndim, ColumnVector& x);
void epowell(int ndim, const ColumnVector& x, double& fx, int& result);

void init_trig (int ndim, ColumnVector& x);
void trig(int ndim, const ColumnVector& x, double& fx, int& result);

void init_gen_brown (int ndim, ColumnVector& x);
void gen_brown(int ndim, const ColumnVector& x, double& fx, int& result);

void init_vardim (int ndim, ColumnVector& x);
void vardim(int ndim, const ColumnVector& x, double& fx, int& result);

void init_chain_singular (int ndim, ColumnVector& x);
void chain_singular(int ndim, const ColumnVector& x, double& fx, int& result);

void init_gen_wood (int ndim, ColumnVector& x);
void gen_wood(int ndim, const ColumnVector& x, double& fx, int& result);

void init_chain_wood (int ndim, ColumnVector& x);
void chain_wood(int ndim, const ColumnVector& x, double& fx, int& result);

void init_broyden (int ndim, ColumnVector& x);
void broyden_tridiag(int ndim, const ColumnVector& x, double& fx, int& result);
void broyden1a(int ndim, const ColumnVector& x, double& fx, int& result);
void broyden1b(int ndim, const ColumnVector& x, double& fx, int& result);
void broyden2a(int ndim, const ColumnVector& x, double& fx, int& result);
void broyden2b(int ndim, const ColumnVector& x, double& fx, int& result);
void tointbroy(int ndim, const ColumnVector& x, double& fx, int& result);

void init_gen_cragg_levy (int ndim, ColumnVector& x);
void gen_cragg_levy(int ndim, const ColumnVector& x, double& fx, int& result);

void init_toint_trig (int ndim, ColumnVector& x);
void toint_trig(int ndim, const ColumnVector& x, double& fx, int& result);

void init_chebyquad (int ndim, ColumnVector& x);
void chebyquad(int ndim, const ColumnVector& x, double& fx, int& result);

void init_nelder (int ndim, ColumnVector& x);
void nelder(int ndim, const ColumnVector& x, double& fx, int& result);


#endif
