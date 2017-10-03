#ifndef minpack2_h
#define minpack2_h
/*------------------------------------------------------------------------
// OPT++ 1.4
// Copyright (C) 1994:
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//
// Header file for MINPACK-2 Test problems
//
// Problems:
//   elastic        Elastic Plastic Torsion
//   journal        Pressure Distribution in a Journal Bearing
//   minsurf        Minimal Surfaces
//   optdesign      Optimal Design with Composite Materials
//   inhomo_super   Inhomogeneous Superconductors
//   lennard_jones  Lennard-Jones Cluster
//   steady_comb    Steady-State Combustion
//   homo_super     Homogeneous Superconductors
//------------------------------------------------------------------------
*/

typedef enum { Function, Gradient, FuncGrad, XStandard, XLower, XUpper} MINPACK_TASK;


/* Elastic Plastic Torsion Problem */

void init_elastic(int n, NEWMAT::ColumnVector& x);
void elastic(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, int& result);
void deptfg(int nx, int ny,  NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, MINPACK_TASK task, double c);

/* Pressure Distribution in Journal Bearing Problem */

void journal(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, int& result);
void init_journal(int n, NEWMAT::ColumnVector& x);
void dpjbfg(int nx, int ny,  NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, MINPACK_TASK task, double ecc, double b);

/* Minimal Surfaces Problem */

void minsurf(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, int& result);
void init_minsurf(int n, NEWMAT::ColumnVector& x);
void dmsafg(int nx, int ny,  NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, MINPACK_TASK task, NEWMAT::ColumnVector& bottom, 
    NEWMAT::ColumnVector& top, NEWMAT::ColumnVector& left, 
    NEWMAT::ColumnVector& right);
void dmsabc(int nx, int ny, NEWMAT::ColumnVector& bottom, 
    NEWMAT::ColumnVector& top, NEWMAT::ColumnVector& left, 
    NEWMAT::ColumnVector& right); 
/* Optimal Design with Composite Materials Problem */

void optdesign(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, int& result);
void init_optdesign(int n, NEWMAT::ColumnVector& x);
void dodcfg(int nx, int ny,  NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, MINPACK_TASK task, double lambda);

/* Steady-State Combustion Problem */

void steady_comb(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, int& result);
void init_steady_comb(int n, NEWMAT::ColumnVector& x);
void dsscfg(int nx, int ny,  NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, MINPACK_TASK task, double lambda);

/* Ginzburg-Landau Superconductivity Problem */

void superconduc(int mode, int n, const NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, int& result);
void init_superconduc(int n, NEWMAT::ColumnVector& x);
void dgl2fg(int nx, int ny,  NEWMAT::ColumnVector& x, double& fx, 
    NEWMAT::ColumnVector& g, MINPACK_TASK task, NEWMAT::ColumnVector& w, 
    int vornum);

#endif
