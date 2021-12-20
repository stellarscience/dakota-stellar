/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <ctype.h>
#include <string>
#include <cmath>
#include <iostream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_SerialSpdDenseSolver.hpp>

#include <teuchos_data_types.hpp>
#include <LinearSystemSolver.hpp>
#include <Function.hpp>
#include <RandomVariables.hpp>

using namespace Surrogates;

namespace {

  const size_t NUM_PARAMS = 1;
  const size_t NUM_RESP   = 1;
  const size_t NUM_OBS = 20;

  //----------------------------------------------------------------------------
  // Helper functions

  void
  fill_data(RealMatrix& vars, RealMatrix& resps)
  {
    vars.shapeUninitialized(NUM_OBS, NUM_PARAMS);
    resps.shapeUninitialized(NUM_OBS, NUM_RESP);

    // Variables
    vars(0 ,0) = 0.4356122655; //vars(0 ,1) =  1.402203228;
    vars(1 ,0) =  0.846735555; //vars(1 ,1) = 0.7146972414;
    vars(2 ,0) = 0.6370206252; //vars(2 ,1) = 0.4323235434;
    vars(3 ,0) = 0.9820371062; //vars(3 ,1) =  1.201303441;
    vars(4 ,0) =  1.394309654; //vars(4 ,1) = 0.5527851693;
    vars(5 ,0) = 0.9700363914; //vars(5 ,1) =  1.026058852;
    vars(6 ,0) =  1.102804294; //vars(6 ,1) = 0.6996376666;
    vars(7 ,0) =    0.5567185; //vars(7 ,1) = 0.8682234433;
    vars(8 ,0) = 0.6044407372; //vars(8 ,1) =  1.297886238;
    vars(9 ,0) = 0.2400064561; //vars(9 ,1) =  1.275490984;
    vars(10,0) = 0.5301719731; //vars(10,1) =  1.231657289;
    vars(11,0) =  1.149139043; //vars(11,1) =  1.496428965;
    vars(12,0) = 0.7649387327; //vars(12,1) =  1.068050129;
    vars(13,0) = 0.5894615341; //vars(13,1) =  1.080760183;
    vars(14,0) = 0.2916886869; //vars(14,1) =  1.258397242;
    vars(15,0) = 0.9515128152; //vars(15,1) = 0.6185770598;
    vars(16,0) =  1.213305147; //vars(16,1) = 0.6065278133;
    vars(17,0) =  1.291005599; //vars(17,1) = 0.9494996162;
    vars(18,0) =  1.259308813; //vars(18,1) = 0.7482905864;
    vars(19,0) = 0.4734308531; //vars(19,1) = 0.9235273973;

    // Responses
    resps(0 ,0) = 1219.127632; //resps(0 ,1) = -0.5113435681; resps(0 ,2) =  1.748367759;
    resps(1 ,0) = 1524.007177; //resps(1 ,1) =  0.3596124794; resps(1 ,2) =  0.08742436939;
    resps(2 ,0) = 3395.121208; //resps(2 ,1) =  0.1896335052; resps(2 ,2) =  -0.1316066664;
    resps(3 ,0) = 2674.001642; //resps(3 ,1) =  0.3637451577; resps(3 ,2) =  0.9521114034;
    resps(4 ,0) = 4115.064175; //resps(4 ,1) =   1.667706825; resps(4 ,2) =  -0.3915833834;
    resps(5 ,0) = 418.0000013; //resps(5 ,1) =  0.4279411749; resps(5 ,2) =  0.5677785712;
    resps(6 ,0) = 1636.008251; //resps(6 ,1) =  0.8663584775; resps(6 ,2) = -0.06190928242;
    resps(7 ,0) = 8435.038913; //resps(7 ,1) = -0.1241762334; resps(7 ,2) =   0.4754526975;
    resps(8 ,0) = 15379.03236; //resps(8 ,1) = -0.2835945142; resps(8 ,2) =    1.382288318;
    resps(9 ,0) = 7364.339371; //resps(9 ,1) = -0.5801423932; resps(9 ,2) =    1.506874023;
    resps(10,0) = 10756.05161; //resps(10,1) = -0.3347463234; resps(10,2) =    1.251893691;
    resps(11,0) = 13635.06123; //resps(11,1) =  0.5723060564; resps(11,2) =    1.664730127;
    resps(12,0) = 3283.003074; //resps(12,1) = 0.05110620039; resps(12,2) =   0.7582617115;
    resps(13,0) = 4179.028449; //resps(13,1) = -0.1929151911; resps(13,2) =   0.8733118053;
    resps(14,0) = 6755.256166; //resps(14,1) = -0.5441163311; resps(14,2) =    1.437719276;
    resps(15,0) = 2020.021171; //resps(15,1) =  0.5960881076; resps(15,2) = -0.09311882868;
    resps(16,0) =  2994.02604; //resps(16,1) =   1.168845472; resps(16,2) =   -0.238776585;
    resps(17,0) = 3779.007178; //resps(17,1) =   1.191945649; resps(17,2) =   0.2560467216;
    resps(18,0) = 5491.008536; //resps(18,1) =   1.211713392; resps(18,2) = -0.06971560457;
    resps(19,0) = 6659.076916; //resps(19,1) = -0.2376269259; resps(19,2) =   0.6161874269;
  }

  //----------------------------------------------------------------------------

  void center_training_data(const RealMatrix& trainPoints, RealMatrix& ctrPoints, RealVector& trainMeans, RealVector& trainStdvs)
  {
    size_t num_v = trainPoints.numCols();
    size_t num_obs = trainPoints.numRows();

    trainMeans.sizeUninitialized(num_v);
    trainStdvs.sizeUninitialized(num_v);

    ctrPoints = trainPoints;

    // get means of each column for X input
    for( size_t i=0; i<num_v; i++ )
    {
      Real sum = 0.;
      for( size_t j=0; j<num_obs; j++ )
        sum += ctrPoints(j,i);
      trainMeans(i) = sum/(double(num_obs));
    }

    // subtract means from each value in each column
    for( size_t i=0; i<num_v; i++ )
    {
      trainStdvs(i) = 0.0;
      for( size_t j=0; j<num_obs; j++ )
      {
        ctrPoints(j,i) -= trainMeans(i);
        trainStdvs(i) += std::pow(ctrPoints(j,i),2);
      }
      trainStdvs(i) = std::sqrt(trainStdvs(i)/double(num_obs-1));
    }

    // divide by standard deviation
    for( size_t i=0; i<num_v; i++ )
      for( size_t j=0; j<num_obs; j++ )
        ctrPoints(j,i) /= trainStdvs(i);

    for( size_t j=0; j<num_obs; j++ )
    {
      for( size_t i=0; i<num_v; i++ )
        std::cout <<" input  " << ctrPoints(j,i) << '\n';
      std::cout <<" output " << trainPoints(j,0) << '\n';
    }
  }

  //----------------------------------------------------------------------------

  void compute_trend_func(const RealMatrix& ctrPoints, RealMatrix& trendFunc, unsigned short order)
  {
    size_t num_v = ctrPoints.numCols();
    size_t num_obs = ctrPoints.numRows();

    switch (order)
    {
      case 0:
        trendFunc.shapeUninitialized(num_obs,1);
        break;
      case 1:
        trendFunc.shapeUninitialized(num_obs,num_v+1);
        break;
      case 2:
        trendFunc.shapeUninitialized(num_obs,2*num_v+1);
        break;
    }

    // all orders require ones in first column
    for( size_t i=0; i<num_obs; i++ )
      trendFunc(i,0) = 1.;

    if( order > 0)
    {
      for (size_t j=0; j<num_v; j++) {
        for (size_t i=0; i<num_obs; i++) {
          trendFunc(i,j+1) = ctrPoints(i,j);
          if (order == 2)
            trendFunc(i,num_v+j+1) = ctrPoints(i,j)*ctrPoints(i,j);
        }
      }
    }

  }

  //----------------------------------------------------------------------------

  void optimize_theta_global(const RealMatrix& ctrPoints, RealMatrix& trendFunc, unsigned short order)
  {
//    // Need to set this up by bringing in NCSU - RWH
//    GPinstance = this;
//    Iterator nll_optimizer; // empty envelope
//
//    // bounds for non-log transformation - ie, no exp(theta)
//    //RealVector theta_lbnds(num_v,1.e-5), theta_ubnds(num_v,150.);
//    // bounds for log transformation of correlation parameters
//    size_t num_v = sharedDataRep->numVars;
//    RealVector theta_lbnds(num_v, false), theta_ubnds(num_v, false);
//    theta_lbnds = -9.; theta_ubnds = 5.;
//#ifdef HAVE_NCSU
//    // NCSU DIRECT optimize of Negative Log Likelihood
//    // Uses default convergence tolerance settings in NCSUOptimizer wrapper!
//    int max_iterations = 1000, max_fn_evals = 10000;
//    nll_optimizer.assign_rep(
//        new NCSUOptimizer(theta_lbnds,theta_ubnds,max_iterations,max_fn_evals,
//          negloglikNCSU),false);
//    nll_optimizer.run(); // no pl_iter needed for this optimization
//    const Variables& vars_star = nll_optimizer.variables_results();
//    const Response&  resp_star = nll_optimizer.response_results();
//    copy_data(vars_star.continuous_variables(), thetaParams);
//
//#ifdef DEBUG
//    Cout << "Optimal NLL = " << resp_star.function_value(0) << '\n';
//#endif //DEBUG
//#endif //HAVE_NCSU
  }


  //----------------------------------------------------------------------------

  void compute_covariance_matrix(const RealMatrix& ctrPoints, const RealVector& thetaParams, RealSymMatrix& cov_matrix)
  {
    // Note, this gets only the lower triangle of covmatrix, which is all
    // that is needed by the Teuchos SPD solvers.  However, if this matrix
    // is to be used in full matrix multiplication, or with a non-SPD
    // solver, then the lower part should be copied to the upper part.

    size_t num_v = ctrPoints.numCols();
    size_t num_obs = ctrPoints.numRows();

    cov_matrix.shapeUninitialized(num_obs);

    RealVector expThetaParms(num_v);

    for (size_t i=0; i<num_v; i++)
      expThetaParms[i]=std::exp(thetaParams[i]);

    Real sume, delta;
    for (size_t j=0; j<num_obs; j++) {
      for (size_t k=j; k<num_obs; k++) {
        sume = 0.0;
        for (size_t i=0; i<num_v; i++) {
          Real pt_diff = ctrPoints(j,i) - ctrPoints(k,i);
          sume += expThetaParms[i]*pt_diff*pt_diff;
        }
        cov_matrix(k,j) = std::exp(-1.0*sume);
      }
    }

    std::cout << "covariance matrix" << '\n';
    for (size_t j=0; j<num_obs; j++){
      for (size_t k=0; k<num_obs; k++)
        std::cout << cov_matrix(j,k) << " ";
      std::cout << '\n';
    }
    std::cout << '\n';
  }

  //----------------------------------------------------------------------------

  int compute_cholesky_factor(const RealSymMatrix& covMatrix, const RealMatrix& ctrPoints,
                              const RealVector& thetaParams, RealSymMatrix& factoredMatrix)
  {
    // Gets the Cholesky factorization of the covariance matrix.  If this
    // is called by prediction routines, then they are basically "stuck"
    // with the current covariance matrix, so if it turns out to be
    // singular to wp, then we iteratively condition it with very small
    // values until it is no longer singular.  However, if this routine is
    // called by the likelihood estimation procedures, we want them to
    // "know" about singular covariance as opposed to "working around it"
    // by conditioning.  Thus, this routine returns an integer of 0 if the
    // covariance matrix is positive definite, and returns 1 if it is
    // singular to wp.

    // Preserve our original matrix
    factoredMatrix = covMatrix;

    int ok;
    size_t num_v = ctrPoints.numCols();
    size_t num_obs = ctrPoints.numRows();
    Real nugget = 1.0e-15;
    Teuchos::SerialSpdDenseSolver<int, Real> cov_solver;
    cov_solver.setMatrix( rcp(&factoredMatrix, false) );
    cov_solver.factorWithEquilibration(true);
    ok = cov_solver.factor();
    if (ok > 0)
    {
      do {
        compute_covariance_matrix(ctrPoints, thetaParams, factoredMatrix);
        for (size_t i=0; i<num_obs; i++) factoredMatrix(i,i) += nugget;
        cov_solver.setMatrix( rcp(&factoredMatrix, false) );
        cov_solver.factorWithEquilibration(true);
        ok = cov_solver.factor();
        nugget *= 3.0;
      } while (ok > 0);
      //#ifdef DEBUG_FULL
      std::cout << "COV matrix corrected with nugget: " << nugget/3. << std::endl;
      //#endif
      std::cout << "factored covariance matrix" << '\n';
      for (size_t j=0; j<num_obs; j++){
        for (size_t k=0; k<num_obs; k++)
          std::cout << factoredMatrix(j,k) << " ";
        std::cout << '\n';
      }
      std::cout << '\n';
      return 1;
    }
    else {
      std::cout << "factored covariance matrix" << '\n';
      for (size_t j=0; j<num_obs; j++){
        for (size_t k=0; k<num_obs; k++)
          std::cout << factoredMatrix(j,k) << " ";
        std::cout << '\n';
      }
      std::cout << '\n';
      return 0;
    }
  }

  //----------------------------------------------------------------------------

  void compute_beta_coefficients(const RealSymMatrix& matrix, const RealMatrix& trainPoints,
                                 const RealMatrix& trendFunc, RealMatrix& betaCoeffs)
  {
    // The choice between generalized least squares (GLS) and ordinary
    // least squares (OLS) is determined simply by which matrix gets
    // passed in, e.g. covariance --> GLS, identity --> OLS bool ordinary

    size_t num_v = trainPoints.numCols();
    size_t num_obs = trainPoints.numRows();
    size_t trend_dim = trendFunc.numCols(); // 1 + trendOrder * num_v;

    // These are copies of data to protect the incoming objects - needed? (RWH)
    RealSymMatrix nonconst_matrix = matrix;
    RealMatrix nonconst_points = trainPoints;
    RealMatrix nonconst_trend = trendFunc;

    Teuchos::SerialSpdDenseSolver<int, Real> solver;
    solver.setMatrix( rcp(&nonconst_matrix, false) );

    RealMatrix Rinv_Y(num_obs, 1, false);
    solver.setVectors( rcp(&Rinv_Y, false), rcp(&nonconst_points, false) );
    solver.solve();

    RealMatrix FT_Rinv_Y(trend_dim, 1, false);
    FT_Rinv_Y.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, trendFunc, Rinv_Y, 0.0);

    RealMatrix Rinv_F(num_obs, trend_dim, false);
    solver.setVectors( rcp(&Rinv_F, false), rcp(&nonconst_trend, false) );
    solver.solve();

    RealMatrix FT_Rinv_F(trend_dim, trend_dim, false);
    FT_Rinv_F.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., trendFunc, Rinv_F, 0.);
    // LPS TO DO:
    //RealSymMatrix FT_Rinv_F(trend_dim, false);
    //Teuchos::symMatTripleProduct(Teuchos::TRANS, 1., trendFunc,
    //                             Rinv, FT_Rinv_F);
    //RealSpdSolver Temp_slvr;

    RealMatrix temphold5(trend_dim, 1, false);
    Teuchos::SerialDenseSolver<int, Real> temp_slvr;
    temp_slvr.setMatrix( rcp(&FT_Rinv_F, false) );
    temp_slvr.setVectors( rcp(&temphold5, false), rcp(&FT_Rinv_Y, false) );
    temp_slvr.factorWithEquilibration(true);
    temp_slvr.factor();
    temp_slvr.solve();

    betaCoeffs.shapeUninitialized(trend_dim, 1);
    for (size_t i=0; i<trend_dim; i++)
    {
      betaCoeffs(i,0) = temphold5(i,0);
      std::cout << "Beta coefficient " << i << " : " << betaCoeffs(i,0) << std::endl;
    }

    if (betaCoeffs(0,0) != betaCoeffs(0,0))
      throw std::runtime_error("Nan for beta at exit of get_beta in GPApproximation\n");
  }

  //----------------------------------------------------------------------------

  Real get_process_variance(const RealSymMatrix& matrix, const RealMatrix& trendFunc,
                            const RealMatrix& trainPoints, const RealMatrix& betaCoeffs)
  {
    size_t num_obs = matrix.numRows();

    RealMatrix YFb(num_obs, 1, false);
    RealMatrix temphold3(1, 1, false);

    RealMatrix Rinv_YFb;
    Rinv_YFb.shape(num_obs, 1);
    YFb.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1., trendFunc, betaCoeffs, 0.);
    YFb.scale(-1);
    YFb += trainPoints;

    Teuchos::SerialSpdDenseSolver<int, Real> cov_solver;
    RealSymMatrix nonconst_matrix = matrix;
    cov_solver.setMatrix( rcp(&nonconst_matrix, false) );
    cov_solver.setVectors( rcp(&Rinv_YFb, false), rcp(&YFb, false) );
    cov_solver.solve();

    temphold3.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1., YFb, Rinv_YFb, 0.);

    return temphold3(0,0)/double(num_obs);
  }



//----------------------------------------------------------------------------
// Unit tests

/**
\brief Generate a gaussian process interpolant using ...
to solve a ...
*/
TEUCHOS_UNIT_TEST(models, gp_approx1)
{
  RealMatrix vars;
  RealMatrix resps;
  fill_data(vars, resps);

  std::cout << vars << std::endl;

  RealMatrix ctrPoints;
  RealVector trainMeans;
  RealVector trainStdvs;
  center_training_data(vars, ctrPoints, trainMeans, trainStdvs);

  /// Theta is the vector of covariance parameters for the GP.
  /// We determine the values of theta by optimization
  /// Currently, the covariance function is
  /// theta[0]*exp(-0.5*sume)+delta*pow(sige,2).  sume is the
  /// sum squared of weighted distances; it involves a sum of
  /// theta[1](Xi(1)-Xj(1))^2 + theta[2](Xi(2)-Xj(2))^2 + ...
  /// where Xi(1) is the first dimension value of multi-dimensional
  /// variable Xi.  delta*pow(sige,2) is a jitter term used to
  /// improve matrix computations. delta is zero for the covariance
  /// between different points and 1 for the covariance between the
  /// same point.  sige is the underlying process error.
  RealVector thetaParams(vars.numCols());
  //thetaParams.putScalar(1.0);
  // Use theta parameters taken from A previous Dakota study using the same test data as above
  thetaParams(0) = 4.9982218158e+00;
  RealSymMatrix cov_matrix;
  compute_covariance_matrix(ctrPoints, thetaParams, cov_matrix);

  RealSymMatrix factored_matrix;
  compute_cholesky_factor(cov_matrix, ctrPoints, thetaParams, factored_matrix);

  RealMatrix trendFunc;
  compute_trend_func(ctrPoints, trendFunc, 2 /* order */);

  // Note that we are currently forced to create GPs sequentially for each variable
  RealMatrix beta_coeffs;
  RealMatrix train_vals(vars.numRows(), 1);
  Teuchos::setCol( Teuchos::getCol( Teuchos::Copy, vars, 1), 0, train_vals);
  compute_beta_coefficients(cov_matrix, train_vals, trendFunc, beta_coeffs);

  std::cout << "Process variance = " <<
    get_process_variance(cov_matrix, trendFunc, train_vals, beta_coeffs)
                                     << std::endl;

  /*
    Steps to implement - RWH

    optimize_theta_global();
  */

  TEST_ASSERT( true );
}

}
