/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_FUNCTION_HPP
#define PECOS_SURROGATES_FUNCTION_HPP

#include <string>
#include <vector>
#include <OptionsList.hpp>
#include <Teuchos_SerialDenseHelpers.hpp>
#include <teuchos_data_types.hpp>
#include <math_tools.hpp>

namespace Pecos {
namespace surrogates {

/**
\class Function

\brief Abstract base class representation of a function \f$f(x)\f$ for a
multivariate variable x.

Functions derived from this class can either be \f$C_0, C_1\f$ or \f$C_2\f$
 continuous. Regardlesd of the function regularity any derived class must
implement value(), \f$C_1\f$ functions must implement gradient() and jacobian()
and \f$C_2\f$ functions must implement hessian().

To avoid complexity of the function hierarchy we do not create seperate classes
for  \f$C_0, C_1\f$ and \f$C_2\f$ functions but rather only implement the
functions necessitated by the function regularity. All other functions must
raise an exception if called (these exceptions are implemented in this base
class)
*/
class Function{
protected:
  /// Model specific options
  util::OptionsList opts_;

  /// The number of QoI (outputs) of the vector-valued function
  int numQOI_;

public:
  /// Default constructtor
  Function();

  /// Destructor
  virtual ~Function();

  /**\brief evaluate the vector-valued function \f$f(x)\f$ at a set of
     samples x.

   \param[in] samples (num_vars x num_samples) matrix
       The coordindates of the samples x.

   \param[out] result_0 (num_samples x num_qoi) matrix
       The values of the vectored function at each of the samples.
       Thsi argument is named values_out so that swig knows this is
       an output variable.
   */
  virtual void value(const RealMatrix &samples, RealMatrix &result_0) = 0;

  /**\brief evaluate the gradient of the ith QoI \f$\nabla f_i(x)\f$ of a vector
     valued function \f$f(x)\f$ at a set of samples x.

   \param[in] samples (num_vars x num_samples) matrix
       The coordindates of the samples x.

   \param[in] qoi
       The index of the QoI for which gradients will be computed.

   \param[out] result (num_samples x num_qoi) matrix
       The gradients of a single function output (QoI) at each of the samples.
   */
  virtual void gradient(const RealMatrix &samples, int qoi, RealMatrix &result_0_0);

  /**\brief evaluate the Jacobian \f$\nabla f(x)\f$ of a vector
     valued function \f$f(x)\f$ at a single samples x.

   \param[in] sample (num_vars x 1) vector
       The coordindates of the sample x.

   \param[out] result_0 (num_qoi x num_vars) matrix
       The jacobian of the vector-valued function at the single sample x.
   */
  virtual void jacobian(const RealVector &sample, RealMatrix &result_0);


  /**\brief evaluate the Hessian of the ith QoI \f$\Delta f_i(x)\f$ of a vector
     valued function \f$f(x)\f$ at a set of samples x.

   \param[in] samples (num_vars x num_samples) matrix
       The coordindates of the samples x.

   \param[in] qoi
       The index of the QoI for which gradients will be computed.


   \param[out] hessians (num_samples) list of (num_vars x num_vars) matrices
       The hessians of a single function output (QoI) at each of the samples.
   */
  virtual void hessian(const RealMatrix &samples, int qoi, RealMatrixList &hessians);

  /**\brief set options specific to the model
   */
  virtual void set_options(const util::OptionsList &opts);

  /**\brief get the model specific options
   */
  virtual void get_options(util::OptionsList &opts);

}; // class Function

// Note see http://www.cplusplus.com/doc/tutorial/typecasting/
// for information on type casting and typeid

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
