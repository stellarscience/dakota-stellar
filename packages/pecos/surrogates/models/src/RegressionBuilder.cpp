/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "RegressionBuilder.hpp"
#include "PolyApproximation.hpp"

namespace Pecos {
namespace surrogates {

void get_canonical_uniform_samples(int M, int N, int seed, RealMatrix &samples){ 
  const Real range_min = -1.; 
  const Real range_max = 1.; 
  typedef boost::uniform_real<> NumberDistribution; 
  typedef boost::mt19937 RNG; 
  typedef boost::variate_generator<RNG&,NumberDistribution> Generator; 
 
  NumberDistribution distribution( range_min, range_max ); 
  RNG generator; 
  Generator numberGenerator( generator, distribution ); 
  generator.seed( seed );
  
  samples.reshape( M, N );
  for ( int i = 0; i < M; i++ ){
    for ( int j = 0; j < N; j++ )
      samples(i,j) = numberGenerator();
  }
}
  
void generate_uniform_samples(int num_vars, int num_samples, int seed, const VariableTransformation &var_transform, RealMatrix &samples ){
  
  RealMatrix canonical_samples;
  // generate samples on canonical domain of random variable, i.e. [-1,1]
  // for uniform vars
  get_canonical_uniform_samples(num_vars,num_samples,seed,canonical_samples);
  
  // use var transform to map to user variables domain
  var_transform.map_samples_to_user_space(canonical_samples,samples);
}

void solve_regression(const RealMatrix &samples,
		      const RealMatrix &values,
		      util::OptionsList &opts,
		      Approximation &approx){

  PolyApproximation& poly_approx = 
    dynamic_cast<PolyApproximation&>(approx);
    
  // Generate matrix of the linear system to be used in
  // regression solve
  RealMatrix basis_matrix;
  poly_approx.generate_basis_matrix(samples,basis_matrix);
    
  // Solve regression problem to get coefficients
  std::shared_ptr<util::LinearSystemSolver> solver = regression_solver_factory(opts);
  RealMatrix coeffs;
  solver->solve(basis_matrix, values, opts);

  solver->get_final_solutions(coeffs);
  // Set the approximation coefficients
  poly_approx.set_coefficients(coeffs);
}

void RegressionBuilder::
build(util::OptionsList &opts, Approximation &approx, util::OptionsList &result){
  // Generate samples to build approximation
  int num_samples = opts.get<int>("num_samples");
  std::string sample_type = opts.get<std::string>("sample_type");

  // Create mc sampler and pass in sample type.
  // For now hack support for uniform mc sampler
  RealMatrix samples;
  int seed = 1337;
  std::shared_ptr<VariableTransformation> var_transform =
    approx.get_variable_transformation();
  generate_uniform_samples(approx.num_vars(), num_samples, seed,
                           *var_transform, samples);

  // Evaluate the function at the build samples
  RealMatrix values;
  targetFunction_->value(samples, values);

  //\todo consider having opts have multiple parameterLists
  //each associated with a particular aspect of build
  //e.g. opts = (sample_opts, regression_opts)
  solve_regression(samples, values, opts, approx);

  result.set("training_samples",samples);
  result.set("training_values",values);
}

}  // namespace surrogates
}  // namespace Pecos
