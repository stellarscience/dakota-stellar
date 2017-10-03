#ifndef tstfcn_h
#define tstfcn_h

/*
 *
 * Declarations for TWAFER
 */

/* Initializer for TWAFER */

void init_twaf(int n, Teuchos::SerialDenseVector<int,double>& x);

/* TWAFER with no analytic derivative */

void twaf(int n, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

/*
 *
 * Declarations for rosenbrock's functions
 */

/* Initializer for Rosenbrock */

void init_rosen(int n, Teuchos::SerialDenseVector<int,double>& x);

/* Rosenbrock with no analytic derivative */

void rosen0(int n, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

/* Rosenbrock least squares formulation with no analytic derivative */

void rosen0_least_squares(int n, const Teuchos::SerialDenseVector<int,double>& x, Teuchos::SerialDenseVector<int,double>& fx, 
	int& result);

/* Rosenbrock with analytic derivative */

void rosen(int mode, int n, const Teuchos::SerialDenseVector<int,double>& x, double& fx, Teuchos::SerialDenseVector<int,double>& g,
	   int& result);

/* Rosenbrock least squares formulation with analytic derivative*/

void rosen_least_squares(int mode, int n, const Teuchos::SerialDenseVector<int,double>& x, 
	Teuchos::SerialDenseVector<int,double>& fx, Teuchos::SerialDenseMatrix<int,double>& gx, int& result);

/* Rosenbrock with analytic derivative and Hessian*/

void rosen2(int mode, int n, const Teuchos::SerialDenseVector<int,double>& x, double& fx, 
	    Teuchos::SerialDenseVector<int,double>& g, Teuchos::SerialSymDenseMatrix<int,double>& H, int& result);

/* Scaled version of Rosenbrock with analytic derivative */
void srosen(int mode, int n, const Teuchos::SerialDenseVector<int,double>& x, double& fx, 
	    Teuchos::SerialDenseVector<int,double>& g, int& result);

/* Scaled version of Rosenbrock with analytic derivative and Hessian*/

void srosen2(int mode, int n, const Teuchos::SerialDenseVector<int,double>& x, double& fx, 
	     Teuchos::SerialDenseVector<int,double>& g, Teuchos::SerialSymDenseMatrix<int,double>& H, int& result);

void   init_illum(int n, Teuchos::SerialDenseVector<int,double>& x);

void   illum(int, int, const Teuchos::SerialDenseVector<int,double>&, double&, Teuchos::SerialDenseVector<int,double>&, int&);

void   illum2(int mode, int n, const Teuchos::SerialDenseVector<int,double>& x, double& fx, 
	      Teuchos::SerialDenseVector<int,double>& g, Teuchos::SerialSymDenseMatrix<int,double>& H, int& result);

void init_erosen (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void erosen(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);
void erosen1(int mode, int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, 
	     Teuchos::SerialDenseVector<int,double>&, int& result);

void init_penalty1 (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void penalty1(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_penalty2 (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void penalty2(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_epowell (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void epowell(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_trig (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void trig(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_gen_brown (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void gen_brown(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_vardim (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void vardim(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_chain_singular (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void chain_singular(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_gen_wood (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void gen_wood(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_chain_wood (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void chain_wood(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_broyden (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void broyden_tridiag(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);
void broyden1a(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);
void broyden1b(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);
void broyden2a(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);
void broyden2b(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);
void tointbroy(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_gen_cragg_levy (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void gen_cragg_levy(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_toint_trig (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void toint_trig(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_chebyquad (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void chebyquad(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);

void init_nelder (int ndim, Teuchos::SerialDenseVector<int,double>& x);
void nelder(int ndim, const Teuchos::SerialDenseVector<int,double>& x, double& fx, int& result);


#endif
