/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <AleatoryVariableTransformation.hpp>

namespace Pecos {
namespace surrogates {

AleatoryVariableTransformation::AleatoryVariableTransformation(){}

AleatoryVariableTransformation::~AleatoryVariableTransformation(){}

void AleatoryVariableTransformation::
set_variables(const std::shared_ptr<Surrogates::Variables> vars){
  Surrogates::VariableTransformation::set_variables(vars);
  aleatoryVars_ = Teuchos::rcp_dynamic_cast<Surrogates::AleatoryVariables>(vars_,true);
  if (aleatoryVars_.is_null())
    throw(std::runtime_error("vars is not an object of type AleatoryVariables"));
}

void AleatoryVariableTransformation::
map_samples_to_user_space(const RealMatrix &samples,
			  RealMatrix &transformed_samples) const {
  int num_vars = aleatoryVars_->num_vars();
  if ( samples.numRows() != aleatoryVars_->num_vars() )
    throw( std::runtime_error("Samples have incorrect number of random variables") );

  int num_samples = samples.numCols();
  transformed_samples.shapeUninitialized(num_vars, num_samples);
  for ( int j=0; j<num_samples; j++ ){
    for ( int i=0; i<num_vars; i++ ){
      transformed_samples(j,i) =
	2.*(samples(i,j)-aleatoryVars_->lb(i))/(aleatoryVars_->ub(i)-aleatoryVars_->lb(i))-1.;
    }
  }
}

void AleatoryVariableTransformation::
map_samples_from_user_space(const RealMatrix &samples, RealMatrix &transformed_samples) const {
  int num_vars = aleatoryVars_->num_vars();
  if ( samples.numRows() != aleatoryVars_->num_vars() )
    throw( std::runtime_error("Samples have incorrect number of random variables") );

  int num_samples = samples.numCols();
  transformed_samples.shapeUninitialized(num_vars, num_samples);
  for ( int j=0; j<num_samples; j++ ){
    for ( int i=0; i<num_vars; i++ ){
      transformed_samples(j,i) =
	(samples(i,j)+1)/2.*(aleatoryVars_->ub(i)-aleatoryVars_->lb(i))+aleatoryVars_->lb(i);
    }
  }
}

void AleatoryVariableTransformation::
map_derivatives_to_user_space(const RealVector &derivatives,
			      int dim, RealVector &transformed_derivatives) const{
  transformed_derivatives.shapeUninitialized( derivatives.numRows(),
					      derivatives.numCols() );
  transformed_derivatives.assign( derivatives );
  Real scaling_factor = 2./(aleatoryVars_->ub(dim)-aleatoryVars_->lb(dim));
  transformed_derivatives *= scaling_factor;
}

void AleatoryVariableTransformation::
map_derivatives_from_user_space(const RealVector &derivatives,
				int dim, RealVector &transformed_derivatives) const{
  transformed_derivatives.shapeUninitialized( derivatives.numRows(),
					      derivatives.numCols() );
  transformed_derivatives.assign( derivatives );
  Real scaling_factor = (aleatoryVars_->ub(dim)-aleatoryVars_->lb(dim)) / 2.;
  transformed_derivatives *= scaling_factor;
}

}  // namespace surrogates
}  // namespace Pecos
