/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <AffineVariableTransformation.hpp>

namespace Pecos {
namespace surrogates {

AffineVariableTransformation::AffineVariableTransformation(){}

AffineVariableTransformation::~AffineVariableTransformation(){}

void AffineVariableTransformation::
set_variables(const std::shared_ptr<Variables> &vars){
  VariableTransformation::set_variables(vars);
  //boundedVars_ = std::dynamic_pointer_cast<BoundedVariables>(vars_);
  boundedVars_ = std::dynamic_pointer_cast<BoundedVariables>(vars_);
  if (!boundedVars_)
    throw(std::runtime_error("vars is not an object of BoundedVariables"));
}

void AffineVariableTransformation::
map_samples_from_user_space(const RealMatrix &samples,
			  RealMatrix &transformed_samples) const {
  int num_vars = boundedVars_->num_vars();
  
  if ( samples.numRows() != boundedVars_->num_vars() )
    throw( std::runtime_error("Samples have incorrect number of random variables") );

  int num_samples = samples.numCols();
  transformed_samples.shapeUninitialized(num_vars, num_samples);

  for ( int j=0; j<num_samples; j++ ){
    for ( int i=0; i<num_vars; i++ ){
      transformed_samples(i,j) =
	2.*(samples(i,j)-boundedVars_->lb(i))/(boundedVars_->ub(i)-boundedVars_->lb(i))-1.;
    }
  }
}

void AffineVariableTransformation::
map_samples_to_user_space(const RealMatrix &samples, RealMatrix &transformed_samples) const {
  int num_vars = boundedVars_->num_vars();
  if ( samples.numRows() != boundedVars_->num_vars() )
    throw( std::runtime_error("Samples have incorrect number of random variables") );

  int num_samples = samples.numCols();
  transformed_samples.shapeUninitialized(num_vars, num_samples);
  for ( int j=0; j<num_samples; j++ ){
    for ( int i=0; i<num_vars; i++ ){
      transformed_samples(i,j) =
	(samples(i,j)+1)/2.*(boundedVars_->ub(i)-boundedVars_->lb(i))+boundedVars_->lb(i);
    }
  }
}

void AffineVariableTransformation::
map_derivatives_to_user_space(const RealVector &derivatives,
			      int dim, RealVector &transformed_derivatives) const{
  transformed_derivatives.shapeUninitialized( derivatives.numRows(),
					      derivatives.numCols() );
  transformed_derivatives.assign( derivatives );
  Real scaling_factor = 2./(boundedVars_->ub(dim)-boundedVars_->lb(dim));
  transformed_derivatives *= scaling_factor;
}

void AffineVariableTransformation::
map_derivatives_from_user_space(const RealVector &derivatives,
				int dim, RealVector &transformed_derivatives) const{
  transformed_derivatives.shapeUninitialized( derivatives.numRows(),
					      derivatives.numCols() );
  transformed_derivatives.assign( derivatives );
  Real scaling_factor = (boundedVars_->ub(dim)-boundedVars_->lb(dim)) / 2.;
  transformed_derivatives *= scaling_factor;
}

}  // namespace surrogates
}  // namespace Pecos
