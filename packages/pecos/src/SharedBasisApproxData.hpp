/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SHARED_BASIS_APPROX_DATA_HPP
#define SHARED_BASIS_APPROX_DATA_HPP

#include "pecos_data_types.hpp"


namespace Pecos {

class ExpansionConfigOptions;
class BasisConfigOptions;
class RegressionConfigOptions;


/// Base class for multivariate basis approximations used for
/// projection of random variables through time or space

/** The base class for basis approximations defined from Fourier
    functions, eigenfunctions, or polynomial functions. */

class SharedBasisApproxData
{
  //
  //- Heading: Friends
  //

  friend class BasisApproximation;
  friend class PolynomialApproximation;
  friend class InterpPolyApproximation;
  friend class NodalInterpPolyApproximation;
  friend class HierarchInterpPolyApproximation;
  friend class OrthogPolyApproximation;
  friend class ProjectOrthogPolyApproximation;
  friend class RegressOrthogPolyApproximation;

public:

  //
  //- Heading: Constructors, destructor, operator=
  //

  /// default constructor
  SharedBasisApproxData();
  /// standard constructor for envelope
  SharedBasisApproxData(short basis_type, const UShortArray& approx_order,
			size_t num_vars);
  /// alternate constructor for envelope
  SharedBasisApproxData(short basis_type, const UShortArray& approx_order,
			size_t num_vars,
			const ExpansionConfigOptions& ec_options,
			const BasisConfigOptions& bc_options,
			const RegressionConfigOptions& rc_options);
  /// copy constructor
  SharedBasisApproxData(const SharedBasisApproxData& shared_data);

  /// destructor
  virtual ~SharedBasisApproxData();

  /// assignment operator
  SharedBasisApproxData operator=(const SharedBasisApproxData& shared_data);

  //
  //- Heading: Member functions
  //

  /// assign letter or replace existing letter with a new one
  void assign_rep(SharedBasisApproxData* data_rep, bool ref_count_incr);

  /// returns dataRep for access to derived class member functions
  /// that are not mapped to the top Approximation level
  SharedBasisApproxData* data_rep() const;

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  SharedBasisApproxData(BaseConstructor, short basis_type, size_t num_vars);

  //
  //- Heading: Data members
  //

  /// type of derived instance as well as sub-type for interpolation
  short basisType;

  /// number of variables used in the approximation
  size_t numVars;

private:

  //
  //- Heading: Member functions
  //

  /// used by the standard envelope constructor to initialize dataRep
  /// to the appropriate derived type
  SharedBasisApproxData*
    get_shared_data(short basis_type, const UShortArray& approx_order,
		    size_t num_vars);
  /// used by the alternate envelope constructor to initialize dataRep
  /// to the appropriate derived type
  SharedBasisApproxData*
    get_shared_data(short basis_type, const UShortArray& approx_order,
		    size_t num_vars, const ExpansionConfigOptions& ec_options,
		    const BasisConfigOptions& bc_options,
		    const RegressionConfigOptions& rc_options);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  SharedBasisApproxData* dataRep;
  /// number of objects sharing dataRep
  int referenceCount;
};


inline SharedBasisApproxData* SharedBasisApproxData::data_rep() const
{ return dataRep; }

} // namespace Pecos

#endif
