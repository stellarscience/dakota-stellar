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

#include <teuchos_data_types.hpp>
#include <CppFunction.hpp>
#include <BoundedVariables.hpp>
#include <AffineVariableTransformation.hpp>
#include <RegressionBuilder.hpp>

#include <DirectANN.hpp>

using namespace Pecos;
using namespace Pecos::util;
using namespace Pecos::surrogates;

namespace {

//--------------------------------------
// Helper functions

/**
\brief Evaluate a quadratic polynomial with no interaction terms, i.e.
*  \[q_i = (i+1) + \sum_j (i+1)*(2*x_j-1)^2+(i+2)*(2*x_j-1)\]
*/
void additive_quadratic_function(const Real* sample, Real* func_vals, const OptionsList &opts){
   const int num_vars = opts.get<int>("num_vars");
   const int num_qoi =  opts.get<int>("num_qoi");

   for (int i=0; i<num_qoi; i++){
     func_vals[i] = (Real)(i+1);
     for (int j=0; j<num_vars; j++){
       Real x = (2*sample[j]-1);
       func_vals[i] += (Real)(i+1)*x*x+(Real)(i+2)*x;
     }
   }
}
/**
\brief Evaluate the "textbook" function of two real variables, i.e.
*  \[q = (x_1 - 1)^4 + (x_2 - 2)^4 \]
*/

void textbook_function(const Real* sample, Real* func_vals, const OptionsList &opts){
   Real x1 = sample[0];
   Real x2 = sample[1];
   func_vals[0] = pow(x1 - 1., 4.) + pow(x2 - 1., 4.);
}

// Compute the mean pointwise absolute error
double mean_abs_error(RealMatrix &pred, RealMatrix &truth) {
  double K = pred.numRows();
  if (K != truth.numRows())
    throw(std::runtime_error("Mismatch between predicted and truth vector dimensions"));
  double eval = 0.0;
  for (int i = 0; i < K; i++) {
    eval += std::abs(pred(i,0) - truth(i,0));
  }
  eval /= K;
  return eval;
}

// Compute the relative mean pointwise absolute error
double relative_mean_abs_error(RealMatrix &pred, RealMatrix &truth) {
  double K = pred.numRows();
  if (K != truth.numRows())
    throw(std::runtime_error("Mismatch between predicted and truth vector dimensions"));
  double eval = 0.0;
  double tmp;
  for (int i = 0; i < K; i++) {
    eval += std::abs((pred(i,0) - truth(i,0))/truth(i,0));
  }
  eval /= K;
  return eval;
}

int test_directANN(Real atol){
  // Data generation parameters
  int num_vars = 2;
  int num_samples = 100;
  int training_data_seed = 3105;
  int test_data_seed = 1000;

  // direct ANN parameters
  int num_qoi = 1; // must be 1 for directANN
  double range_param = 2.0;
  int directANN_seed = 8105;
  int nodes = num_samples - 1;

  // Define the function variables
  //RealVector ranges;
  //define_homogeneous_ranges(num_vars, 0., 1., ranges);
  RealVector ranges(4);
  ranges[0] = 0.;
  ranges[1] = 1.;
  ranges[2] = -1.;
  ranges[3] = 1.;
  std::shared_ptr<Variables> variables(new BoundedVariables);
  std::dynamic_pointer_cast<BoundedVariables>(variables)->set_ranges(ranges);

  // Define the variable transformation
  std::shared_ptr<VariableTransformation>
    var_transform(new AffineVariableTransformation);
  var_transform->set_variables(variables);

  // Define the function to approximate
  CppFunction model;
  OptionsList model_opts;
  model_opts.set("num_vars", num_vars);
  model_opts.set("num_qoi", num_qoi);
  model.set_options(model_opts);
  model.set_function(&additive_quadratic_function);
  //model.set_function(&textbook_function);

  // Generate training and test data
  RealMatrix training_samples, training_vals;
  generate_uniform_samples(num_vars, num_samples, training_data_seed, *var_transform,
                           training_samples);
  model.value(training_samples, training_vals);

  RealMatrix test_samples, test_vals;
  generate_uniform_samples(num_vars, num_samples, test_data_seed, *var_transform,
                           test_samples);
  model.value(test_samples, test_vals);

  // Build the Direct ANN approximation for 3 different scalers
  // and compare to "gold" values
  double mean_normalization_gold = 1.06246e-5;
  double min_max_normalization_gold = 3.26542e-5;
  double standardization_gold = 3.04479e-4;

  RealMatrix pred(num_samples, 1);
  double error_eval;

  std::string scaler_type = "mean normalization";
  DirectANN directANN_mean_nm(nodes, range_param, training_samples, training_vals,
                         scaler_type, directANN_seed);
  directANN_mean_nm.value(test_samples, pred);
  error_eval = relative_mean_abs_error(pred, test_vals);
  if (!(std::abs(error_eval - mean_normalization_gold < atol))) {
    std::cout << "1\n";
    return 1;
  }

  scaler_type = "min-max normalization";
  DirectANN directANN_mm_nm(nodes, range_param, training_samples, training_vals,
                            scaler_type, directANN_seed);
  directANN_mm_nm.value(test_samples, pred);
  error_eval = relative_mean_abs_error(pred, test_vals);
  if (!(std::abs(error_eval - min_max_normalization_gold < atol))) {
    std::cout << "2\n";
    return 2;
  }

  scaler_type = "standardization";
  DirectANN directANN_sd(nodes, range_param, training_samples, training_vals,
                         scaler_type, directANN_seed);
  directANN_sd.value(test_samples, pred);
  error_eval = relative_mean_abs_error(pred, test_vals);
  if (!(std::abs(error_eval - standardization_gold < atol))) {
    std::cout << "3\n";
    return 3;
  }
  else {
    return 0;
  }
}

//----------------------------------------------------------------
// Unit tests

/**
\brief Create a direct ANN from sample data and test performance
       for 3 data scaler types
*/
TEUCHOS_UNIT_TEST(dry_run, test_directANN)
{
  TEST_ASSERT(!test_directANN(1e-10));
}

}
