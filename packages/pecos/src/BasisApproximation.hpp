/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef BASIS_APPROXIMATION_HPP
#define BASIS_APPROXIMATION_HPP

#include "pecos_data_types.hpp"


namespace Pecos {

class SharedBasisApproxData;
class SurrogateData;


/// Base class for multivariate basis approximations used for
/// projection of random variables through time or space

/** The base class for basis approximations defined from Fourier
    functions, eigenfunctions, or polynomial functions. */

class BasisApproximation
{
public:

  //
  //- Heading: Constructors, destructor, operator=
  //

  /// default constructor
  BasisApproximation();
  /// standard constructor for envelope
  BasisApproximation(const SharedBasisApproxData& shared_data);
  /// copy constructor
  BasisApproximation(const BasisApproximation& basis_approx);

  /// destructor
  virtual ~BasisApproximation();

  /// assignment operator
  BasisApproximation operator=(const BasisApproximation& basis_approx);

  //
  //- Heading: Virtual functions
  //

  /// retrieve the approximate function value for a given parameter vector
  virtual Real value(const RealVector& x);
  /// retrieve the approximate function gradient for a given parameter vector
  virtual const RealVector& gradient(const RealVector& x);
  /// retrieve the approximate function Hessian for a given parameter vector
  virtual const RealSymMatrix& hessian(const RealVector& x);

  /// set PolynomialApproximation::surrData
  virtual void surrogate_data(const SurrogateData& data);
  /// get PolynomialApproximation::surrData
  virtual const SurrogateData& surrogate_data() const;

  /// return the minimum number of samples (unknowns) required to
  /// build the derived class approximation type in numVars dimensions
  virtual int min_coefficients() const;
  /// calculate the approximation coefficients using the SurrogateData
  virtual void compute_coefficients();
  /// recalculate the approximation coefficients following SurrogateData update
  virtual void increment_coefficients();
  /// restore the approximation coefficients to the state preceding the last
  /// increment
  virtual void decrement_coefficients();
  /// restore the approximation coefficients to a previously incremented state
  /// as identified by the current data increment
  virtual void push_coefficients();
  /// finalize the coefficients by applying all previously evaluated increments
  virtual void finalize_coefficients();

  /// store the current coefficients for later combination
  virtual void store_coefficients(size_t index = _NPOS);
  /// restore a previously stored coefficient state
  virtual void restore_coefficients(size_t index = _NPOS);
  /// swap the current coefficients with a previously stored set
  virtual void swap_coefficients(size_t index);
  /// remove a redundant stored entry prior to combine_coefficients
  /// (default is pop_back)
  virtual void remove_stored_coefficients(size_t index = _NPOS);

  /// combine the current coefficients with a previously stored set
  virtual void combine_coefficients(short combine_type, size_t swap_index);

  /// print the coefficient array computed in compute_coefficients()
  virtual void print_coefficients(std::ostream& s, bool normalized);

  /// return the coefficient array computed by compute_coefficients()
  virtual RealVector approximation_coefficients(bool normalized) const;
  /// set the coefficient array from external sources, rather than
  /// computing with compute_coefficients()
  virtual void approximation_coefficients(const RealVector& approx_coeffs,
					  bool normalized);

  /// retrieve a vector of coefficient label strings, one per expansion term
  virtual void coefficient_labels(std::vector<std::string>& coeff_labels) const;

  //
  //- Heading: Member functions
  //

  /// assign letter or replace existing letter with a new one
  void assign_rep(BasisApproximation* approx_rep, bool ref_count_incr);

  /// returns approxRep for access to derived class member functions
  /// that are not mapped to the top Approximation level
  BasisApproximation* approx_rep() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  BasisApproximation(BaseConstructor, const SharedBasisApproxData& shared_data);

  //
  //- Heading: Data
  //

  /// contains the approximation data that is shared among the response set
  SharedBasisApproxData* sharedDataRep;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// basisApproxRep to the appropriate derived type.
  BasisApproximation*
    get_basis_approx(const SharedBasisApproxData& shared_data);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  BasisApproximation* basisApproxRep;
  /// number of objects sharing basisApproxRep
  int referenceCount;
};


inline BasisApproximation* BasisApproximation::approx_rep() const
{ return basisApproxRep; }

} // namespace Pecos

#endif
