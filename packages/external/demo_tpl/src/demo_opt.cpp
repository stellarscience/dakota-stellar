/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       Demo_Opt TPL
//- Description: A demo optimization TPL used as a pattern and for testing Dakota adapters
//- Owner:       Russell Hooper
//- Checked by:  ...

#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <limits>
#include <random>

// Demo TPL headers
#include "demo_opt.hpp"

// -----------------------------------------------------------------

Demo_Opt::Demo_Opt() :
  nln_eq_fn_callback_(NULL),
  nln_ineq_fn_callback_(NULL)
{
}

// -----------------------------------------------------------------

bool
Demo_Opt::set_solver_options(const std::string & filename, bool verbose)
{
  if( verbose )
    std::cout << "Setting Demo_Opt solver options using file \""<<filename<<"\"" << std::endl;

  options_file_ = filename;

  return true;
}

// -----------------------------------------------------------------

bool
Demo_Opt::initialize(bool verbose)
{
  if( verbose )
  {
    std::cout << "Doing Demo_Opt::initialize." << std::endl;
    std::cout << "Registered parameters :\n";
    for( auto ip : int_params_ )
      std::cout << ip.first << " = " << ip.second << std::endl;
    for( auto dp : dbl_params_ )
      std::cout << dp.first << " = " << dp.second << std::endl;
  }

  return true;
}

// -----------------------------------------------------------------

void
Demo_Opt::set_problem_data(const std::vector<double> & init ,
                           const std::vector<double> & lower,
                           const std::vector<double> & upper )
{
  assert( init.size() == lower.size() );
  assert( init.size() == upper.size() );

  init_vals_  = init;
  lower_bnds_ = lower;
  upper_bnds_ = upper;
}

// -----------------------------------------------------------------

bool
Demo_Opt::execute(bool verbose)
{
  if( verbose )
    std::cout << "Doing Demo_Opt::execute." << std::endl;

  int max_evals = 100;
  if( int_params_.count("Maximum Evaluations") > 0 )
    max_evals = int_params_["Maximum Evaluations"];

  double fn_tol = 1.e-4;
  if( dbl_params_.count("Function Tolerance") > 0 )
    fn_tol = dbl_params_["Function Tolerance"];

  int num_params = (int)init_vals_.size();
  best_x_.clear();
  best_f_ = std::numeric_limits<double>::max();

  //assert( dbl_params_.count("Objective Target") > 0 );
  //double target = int_params_["Objective Target"];
  double target = 0.0;

  std::default_random_engine generator;
  std::vector< std::uniform_real_distribution<double> > distributions;
  for( size_t i=0; i<init_vals_.size(); ++i )
    distributions.push_back(std::uniform_real_distribution<double>(lower_bnds_[i],upper_bnds_[i]));

  // Crude "optimization" based on random sampling over parameter space
  std::vector<double> x(num_params);

  double fn;
  std::vector<double> nln_eqs;
  if( nln_eq_fn_callback_ )
    nln_eqs.resize( nln_eq_fn_callback_->get_num_nln_eq() );
  std::vector<double> nln_ineqs;
  if( nln_ineq_fn_callback_ )
    nln_ineqs.resize( nln_ineq_fn_callback_->get_num_nln_ineq() );

  int i = 0;
  while( i<=max_evals && best_f_>fn_tol )
  {
    if( i == 0 )
      x = init_vals_;
    else
      for( int np=0; np<num_params; ++np )
        x[np] = distributions[np](generator);
    // Get objective fn
    fn = obj_fn_callback_->compute_obj(x, false);

    // Get nonlinear equality constraints values if applicable
    if( nln_eq_fn_callback_ )
    {
      nln_eq_fn_callback_->compute_nln_eq(nln_eqs, x, false);
      for( auto eqval : nln_eqs )
        fn += fabs(eqval);
    }

    // Get nonlinear inequality constraints values if applicable
    if( nln_ineq_fn_callback_ )
    {
      nln_ineq_fn_callback_->compute_nln_ineq(nln_ineqs, x, false);
      for( auto eqval : nln_ineqs )
        fn += fabs(eqval);
    }

    if( fabs(fn-target) < best_f_ )
    {
      best_x_         = x;
      best_f_         = fabs(fn-target);
      best_nln_eqs_   = nln_eqs;
      best_nln_ineqs_ = nln_ineqs;
    }
    ++i;
  }

  if( verbose )
    std::cout << "Found best_f_ = " << best_f_ << std::endl;

  return true;
}

// -----------------------------------------------------------------
