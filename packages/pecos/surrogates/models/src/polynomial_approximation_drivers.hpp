/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_POLYNOMIAL_APPROXIMATION_DRIVERS_HPP
#define PECOS_SURROGATES_POLYNOMIAL_APPROXIMATION_DRIVERS_HPP

#include "PolyApproximation.hpp"
#include "RegressionBuilder.hpp"
#include "linear_algebra.hpp"
#include "math_tools.hpp"
// TODO: to be removed
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>

namespace Pecos {
namespace surrogates {

  enum SolverTypes{QR_LSQ,SVD_LSQ};
  

void regression_solve(const RealMatrix &matrix, const RealMatrix &rhs,
		      const util::OptionsList &opts,
		      RealMatrix &coeffs, util::OptionsList &status){
  switch (opts.get<int>("solver_type")){
    case QR_LSQ:
      //resize_if_needed(coeffs, matrix.numRows(), rhs.numRows());
      util::qr_solve(matrix,rhs,coeffs,Teuchos::NO_TRANS);
      break;
    default:
      util::error("Incorrect regression solver specified.");
      break;
    }
}


void regression_driver(Function &target_function, PolyApproximation &approx, util::OptionsList &opts ){

  // Generate samples to build approximation
  int num_samples = opts.get<int>("num_samples");
  std::string sample_type = opts.get<std::string>("sample_type");

  // Create mc sampler and pass in sample type. For now hack support for
  // uniform mc sampler
  RealMatrix samples;
  int seed = 1337;
  std::shared_ptr<VariableTransformation> var_transform=
    approx.get_variable_transformation();
  generate_uniform_samples(approx.num_vars(), num_samples, seed,
			   *var_transform, samples);

  // Evaluate the function at the build samples
  RealMatrix values;
  target_function.value(samples, values);

  // Generate matrix of the linear sytem to be used in regression solve
  RealMatrix basis_matrix;
  approx.generate_basis_matrix(samples,basis_matrix);

  // Use regressiong to find coefficients
  util::OptionsList regression_status;
  RealMatrix coeffs;
  regression_solve(basis_matrix, values, opts, coeffs,
                   regression_status);

  // Set the approximation coefficients
  approx.set_coefficients(coeffs);

  // figure out how to return class using swig argout
  // Then consider returning a status OptionsList
}

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
