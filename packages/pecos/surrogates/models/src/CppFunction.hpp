/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_CPP_FUNCTION_HPP
#define PECOS_SURROGATES_CPP_FUNCTION_HPP

#include <Function.hpp>

namespace Pecos {
namespace surrogates {

/**
\class CppFunction

\brief \f$C_0\f$ function that evaluates wraps a function pointer.

This class is very useful for unit testing.
*/
class CppFunction : public Function {
public:

  CppFunction();

  ~CppFunction();

  void set_function(
     void (*function)(const Real* sample, Real* func_vals,
		      const util::OptionsList &opts));

  /** \copydoc Function::value() */
  void value(const RealMatrix &samples, RealMatrix &values_out);

  /** \brief evaluate the function \f$f(x)\f$ at single sample x

  \param[in] sample (num_vars x 1) vector
       The coordindates of the sample x.

   \param[out] result (num_qoi x 1) vector
       The values of the vectored function at the sample x.
   */
  void value(const RealVector &sample, RealVector &result);

protected:
  /// The function \f$f(x)\f$
  void (*targetFunction_)(const Real*, Real*,
		       const util::OptionsList &opts);

}; // class CppFunction

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
