/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef DATA_TRANSFORMATION_HPP
#define DATA_TRANSFORMATION_HPP

#include "pecos_data_types.hpp"
#include "ProbabilityTransformation.hpp"
//#include "BasisFunction.hpp"


namespace Pecos {


/// Base class for forward/inverse transformations between time and
/// frequency domain data

/** The base class for data transformations based on forward/inverse
    mappings between the time and frequency domain based on
    spectral/FFT, Karhunen-Loeve, and sampling-based approaches. */

class DataTransformation
{
public:

  /// default constructor
  DataTransformation();
  /// standard constructor for envelope
  DataTransformation(const String& data_trans_type);
  /// copy constructor
  DataTransformation(const DataTransformation& data_trans);

  /// destructor
  virtual ~DataTransformation();

  /// assignment operator
  DataTransformation operator=(const DataTransformation& data_trans);

  //
  //- Heading: Virtual functions
  //

  /// set scalar data
  virtual void initialize(const Real& total_t, const Real& w_bar, size_t seed);

  /// set PSD to standard embedded function
  virtual void power_spectral_density(const String& psd_name,
				      const Real& param = 0.);
  // define PSD from a user-defined function
  //virtual void power_spectral_density(fn_ptr);
  /// pass a discretized PSD directly: vector or pairs...
  virtual void power_spectral_density(const RealVector& psd);

  //virtual void correlation_function(const String& fn_name,
  //                                  const Real& param = 0.);
  //virtual void correlation_function(fn_ptr);
  //virtual void correlation_function(const RealRealPairArray& corr_fn);

  /// compute and return InverseTransformation::inverseSample
  virtual const RealVector& compute_sample();
  /// compute and return InverseTransformation::inverseSamples
  virtual const RealMatrix& compute_samples(size_t num_samples);

  //
  //- Heading: Member functions
  //

protected:

  //
  //- Heading: Constructors
  //

  /// constructor initializes the base class part of letter classes
  /// (BaseConstructor overloading avoids infinite recursion in the
  /// derived class constructors - Coplien, p. 139)
  DataTransformation(BaseConstructor);

  //
  //- Heading: Member functions
  //


  //
  //- Heading: Data members
  //

  /// nonlinear variable transformation
  ProbabilityTransformation probTransform;
  // set of Fourier, eigen, or polynomial basis functions
  //BasisFunctionArray basisFns;

private:

  //
  //- Heading: Member functions
  //

  /// Used only by the standard envelope constructor to initialize
  /// transRep to the appropriate derived type.
  DataTransformation* get_data_trans(const String& data_trans_type);

  //
  //- Heading: Data members
  //

  /// pointer to the letter (initialized only for the envelope)
  DataTransformation* dataTransRep;
  /// number of objects sharing dataTransRep
  int referenceCount;
};

} // namespace Pecos

#endif
