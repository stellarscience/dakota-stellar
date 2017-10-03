/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef INVERSE_TRANSFORMATION_HPP
#define INVERSE_TRANSFORMATION_HPP

#include "DataTransformation.hpp"
#include "LHSDriver.hpp"


namespace Pecos {


/// Class for inverse data transformation.

/** The InverseTransformation employs an inverse transform to map from
    the frequency domain to the time domain. */

class InverseTransformation: public DataTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  InverseTransformation();  ///< constructor
  ~InverseTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void initialize(const Real& total_t, const Real& w_bar, size_t seed);

  void power_spectral_density(const String& psd_name, const Real& param = 0.);
  //void power_spectral_density(fn_ptr);
  void power_spectral_density(const RealRealPairArray& psd);

  //void correlation_function(const String& fn_name, Real param = 0.);
  //void correlation_function(fn_ptr);
  //void correlation_function(const RealRealPairArray& corr_fn);

  // return inverseSample
  //const RealVector& sample() const;
  // return inverseSamples
  //const RealMatrix& samples() const;

  //
  //- Heading: Data
  //

  /// total time window
  Real totalTime;
  /// time increment
  Real deltaTime;
  /// discretized time sequence
  RealVector timeSequence;

  /// cut-off frequency (rad/s)
  Real omegaBar;
  /// frequency increment (rad/s)
  Real deltaOmega;
  /// discretized frequency sequence
  RealVector omegaSequence;

  /// PSD sequence (frequency domain)
  RealVector psdSequence;

  /// LHS wrapper for generating normal or uniform sample sets
  LHSDriver lhsSampler;

  /// a single computed inverse sample (time domain)
  RealVector inverseSample;
  /// a computed set of inverse samples (time domain)
  RealMatrix inverseSamples;

private:

  //
  //- Heading: Utility routines
  //

};


inline InverseTransformation::InverseTransformation():
  DataTransformation(BaseConstructor()), lhsSampler("lhs", IGNORE_RANKS, false)
{ }


inline InverseTransformation::~InverseTransformation()
{ }


//inline const RealVector& InverseTransformation::sample() const
//{ return inverseSample; }


//inline const RealMatrix& InverseTransformation::samples() const
//{ return inverseSamples; }

} // namespace Pecos

#endif
