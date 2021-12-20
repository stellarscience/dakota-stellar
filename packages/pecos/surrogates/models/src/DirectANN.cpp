/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <DirectANN.hpp>
#include <BoundedVariables.hpp>
#include <AffineVariableTransformation.hpp>
#include <RegressionBuilder.hpp>
#include "linear_solvers.hpp"
#include "OptionsList.hpp"

namespace Pecos {
namespace surrogates {
 
DirectANNBasisSet::DirectANNBasisSet(const RealMatrix &weights_in)
  : weights(weights_in) {
}

Real DirectANNBasisSet::nodeSum(int index, const RealVector &x) const
{
  if (index >= weights.numRows())
    throw(std::runtime_error("basis set index is out of bounds"));
  if (x.length() + 1 !=  weights.numCols())
    throw(std::runtime_error("size of input to hidden layer weights matrix is incorrect"));
  Real sum = 0.0;
  for (int i = 0; i < x.length(); i++) {
    sum += weights(index,i)*x[i];
  }
  sum += weights(index, x.length()); // bias term
  return sum;
}

Real DirectANNBasisSet::eval(int index, const RealVector &x) const
{
  return tanh(nodeSum(index, x));
}

/*
 * no longer valid because output node has no tanh(.) ...
Real DirectANNBasisSet::deriv(int index, const RealVector &x, const IntVector &vars) const
{
  if (vars.length() != 1) 
    throw(std::runtime_error("vars vector too long"));
  if (vars[0] >= x.length())
    throw(std::runtime_error("vars[0] entry too large"));
  Real sum = nodeSum(index, x);
  Real tanhsum = tanh(sum);
  return (1.0 - tanhsum*tanhsum)*weights(index, vars[0]);
}
*/

// default constructor
DirectANN::DirectANN(){}
// default destructor
DirectANN::~DirectANN(){}

// main constructor
DirectANN::DirectANN(const int num_nodes, const Real range_param, 
                     const RealMatrix &samples, const RealMatrix &response,
                     const std::string scaler_type, const int seed) {

  const int num_vars = samples.numRows();
  const int num_samples = samples.numCols();
  nodes = std::min(num_nodes, maxNodes); // max of 100 allowed 
  if (num_samples < num_nodes + 1) {
    throw(std::runtime_error("DirectANN needs at least"
          " as many samples are the number of nodes + 1"));
  }

  // Scale the data
  if (scaler_type == "mean normalization")
    ds = std::make_shared<NormalizationScaler>(samples, true);
  else if (scaler_type == "min-max normalization")
    ds = std::make_shared<NormalizationScaler>(samples, false);
  else if (scaler_type == "standardization")
    ds = std::make_shared<StandardizationScaler>(samples);
  else
    throw(std::string("Invalid scaler type"));

  RealMatrix scaled_samples = ds->getScaledFeatures();

  // Randomly generate weights for the first layer
  RealVector ranges;
  define_homogeneous_ranges(nodes, -range_param/2.0, range_param/2.0, ranges);
  std::shared_ptr<Variables> variables(new BoundedVariables);
  std::dynamic_pointer_cast<BoundedVariables>(variables)->set_ranges(ranges);
  std::shared_ptr<VariableTransformation>
    var_transform(new AffineVariableTransformation);
  var_transform->set_variables(variables);

  RealMatrix random_weights;
  generate_uniform_samples(nodes, num_vars + 1, seed, *var_transform,
                             random_weights);
  bs = std::make_shared<DirectANNBasisSet>(random_weights);

  // construct the linear system
  RealMatrix A(num_samples, nodes + 1);
  RealVector b(num_samples);

  for (int s = 0; s < num_samples; s++) {
    for (int n = 0; n < nodes; n++) {
      RealVector sample = Teuchos::getCol<int,Real> (Teuchos::View, 
                                                     scaled_samples, s);
      A(s,n) = bs->eval(n, sample);
    }
    A(s, nodes) = 1.0; // bias term
    b[s] = response(s,0);
  }

  // set solver options and solve
  util::OptionsList params;
  params.set("regression_type",util::ORTHOG_MATCH_PURSUIT);
  std::shared_ptr<util::LinearSystemSolver> solver = regression_solver_factory(params);

  solver->solve(A, b, params);
  solver->get_final_solutions(coeffs);
}

// constructor for existing (i.e. calibrated) ANN -- stills needs
// scaler info
/*
DirectANN::DirectANN(const RealMatrix& bs_matrix, const RealVector& coeffs_in)
{
  bs = std::make_shared<DirectANNBasisSet>(bs_matrix);
  setCoefficients(coeffs_in);
}
*/

// constructor for known 1st layer weights, unknown output layer weights
// initially used to compare output from pecos ANN ... no longer suitable
// for this in current form.
// 
/*
DirectANN::DirectANN(int num_nodes, const RealMatrix& bs_matrix,
                     RealMatrix &samples, RealMatrix& response)
{
  bs = std::make_shared<DirectANNBasisSet>(bs_matrix);

  const int num_samples = samples.numCols();
  nodes = std::min(num_nodes, maxNodes); // max of 100 allowed 
  if (num_samples < num_nodes + 1) {
    throw(std::runtime_error("DirectANN needs at least"
          " as many samples are the number of nodes + 1"));
  }

  // Scale the data
  const Real norm_factor = 0.8;
  const int num_vars = samples.numRows();
  ds = std::make_shared<NormalizationScaler>(samples, false);
  //ds = std::make_shared<StandardizationScaler>(samples);
  RealMatrix scaled_samples = ds->getScaledFeatures();

  // construct the linear system
  RealMatrix A(num_samples, nodes + 1);
  RealVector b(num_samples);

  for (int s = 0; s < num_samples; s++) {
    for (int n = 0; n < nodes; n++) {
      RealVector sample = Teuchos::getCol<int,Real> (Teuchos::View, 
                                                     scaled_samples, s);

      A(s,n) = bs->eval(n, sample);
    }
    A(s, nodes) = 1.0; // bias term
    b(s) = response(s,0);
  }

  // set solver options and solve
  util::OptionsList params;
  params.set("regression_type",util::ORTHOG_MATCH_PURSUIT);
  //params.set("regression_type",util::QR_LEAST_SQ_REGRESSION);
  std::shared_ptr<util::LinearSystemSolver> solver = regression_solver_factory(params);

  solver->solve(A, b, params);
  solver->get_final_solutions(coeffs);
}
*/


Real DirectANN::evaluate(const RealVector &x) const
{
  RealVector u;
  if (ds->has_scaling) {
    u = ds->scaleFeatures(x);
  }
  else {
    u = x;
  }
  const int num_features = u.length();
  const int num_coeffs = coeffs.length();
  if (num_features + 1 != bs->weights.numCols()) {
    throw(std::runtime_error("DirectANN input to hidden layer weights do"
          " not match the number of features"));
  }
  if (nodes != bs->weights.numRows()) {
    throw(std::runtime_error("DirectANN input to hidden layer weights do"
          " not match the number of nodes"));
  }
  Real sum = 0.0;
  for (int i = 0; i < nodes; i++) {
    sum += coeffs[i]*bs->eval(i, u);
  }
  sum += coeffs[num_coeffs - 1]; // bias term
  return sum;
}

void DirectANN::value(const RealMatrix &samples, RealMatrix &approx_values) {
  const int num_features = samples.numRows();
  const int num_samples = samples.numCols();
  if (num_samples != approx_values.numRows()) {
    throw(std::runtime_error("Direct ANN value inputs are not consistent."
          " Number of samples and approximation sizes do not match"));
  }
  RealVector x(num_features);
  // can't use Teuchos::getCol because it doesn't compile with const samples
  for (int i = 0; i < num_samples; i++) {
    for (int j = 0; j < num_features; j++) {
      x(j) = samples(j,i);
    }
    approx_values(i,0) = this->evaluate(x);
  }
}

void DirectANN::gradient(const RealMatrix &samples, int qoi, RealMatrix &gradients) {
  throw(std::string("This Function type does not provide gradients"));
}

void DirectANN::jacobian(const RealVector &sample, RealMatrix &jacobian) {
  throw(std::string("This Function type does not provide a jacobian"));
}

void DirectANN::hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians) {
  throw(std::string("This Function type does not provide hessians"));
}

}  // namespace surrogates
}  // namespace Pecos
