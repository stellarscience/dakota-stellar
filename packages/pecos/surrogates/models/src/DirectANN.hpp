/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_DIRECT_ANN_HPP
#define PECOS_SURROGATES_DIRECT_ANN_HPP

#include <Approximation.hpp>
#include <DataScaler.hpp>
#include "teuchos_data_types.hpp"

namespace Pecos {
namespace surrogates {


/// The basis set is the first layer in the neural network.  It
/// applies random weights and bias to the (previously scaled) inputs
/// and a tanh activation function. gamma(x) = tanh( A0 x + theta0 )

class DirectANNBasisSet
{
public:
  /// random weights and bias [ A0 | theta0 ] for input layer
  RealMatrix weights;

  /// default constructor 
  DirectANNBasisSet(){}

  /// standard constructor accepting weights [ A0 | theta0 ] to apply
  /// to the input layer
  DirectANNBasisSet(const RealMatrix& weights_in);

  /// evaluate the index-th basis function at the point x
  Real eval(int index, const RealVector &x) const;

  /// evaluate the index-th basis function derivative at the point x
  //Real deriv(int index, const RealVector& x, const IntVector& vars) const;

  /// compute the contribution due to the index-th basis function at the point x
  Real nodeSum(int index, const RealVector &x) const;

};
  
class DirectANN: public Approximation {
protected:

  /// basis set that maps the input layer to the hidden layer
  std::shared_ptr<DirectANNBasisSet> bs;

  /// coefficients and bias for hidden layer to the output node
  RealVector coeffs;

  /// number of user-requested hidden layer nodes
  const int maxNodes = 100;

  /// number of nodes in hidden layer
  int nodes;

  /// evaluate the model at the point x
  Real evaluate(const RealVector &x) const;

  /// Coefficients setter
  void setCoefficients(const RealVector &c) {coeffs = c;}

public:

  DirectANN();

  ~DirectANN();

  /// Main constructor
  DirectANN(const int num_nodes, const Real range_param, 
            const RealMatrix &samples, const RealMatrix &response,
            const std::string scaler_type = "mean_normalization", 
            const int seed = 129);


  /// Constructor from existing model. 
  //  not done ... needs scaler info.
  // DirectANN(const RealMatrix& bs_matrix, const RealVector& coeffs_in);

  /// Constructor from existing weights (unknown coeffs)
  //  Needs test coverage ... originally used to compare to Surfpack
  //  but the this ANN implementation is now distinct.
  //  DirectANN(int num_nodes, const RealMatrix &bs_matrix,
  //          RealMatrix &samples, RealMatrix &response);

  /// scaler for the input data
  std::shared_ptr<DataScaler> ds;

  /// Coefficients getter
  RealVector getCoefficients() {return coeffs;}


  /** \copydoc Function::value() */
  virtual void value(const RealMatrix &samples, RealMatrix &approx_values);

  /** \copydoc Function::gradient() */
  virtual void gradient(const RealMatrix &samples, int qoi, RealMatrix &result_0);

  /** \copydoc Function::jacobian() */
  virtual void jacobian(const RealVector &sample, RealMatrix &result_0);

  /** \copydoc Function::hessian() */
  virtual void hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians);

}; // class DirectANN


}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
