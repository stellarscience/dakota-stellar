/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef FOURIER_INVERSE_TRANSFORMATION_HPP
#define FOURIER_INVERSE_TRANSFORMATION_HPP

#include "InverseTransformation.hpp"

#ifdef HAVE_FFTW
#include "fftw3.h"
#endif

namespace Pecos {

// values for fourierMethod indicating use of Shinozuka-Deodatis or Grigoriu
enum { IFFT_SD, IFFT_G };


/// Class for iFFT data transformation.

/** The FourierInverseTransformation employs an inverse fast Fourier
    transform (iFFT) to map from the frequency domain to the time domain. */

class FourierInverseTransformation: public InverseTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  FourierInverseTransformation(const String& data_trans_type); ///< constructor
  ~FourierInverseTransformation();                             ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void initialize(const Real& total_t, const Real& w_bar, size_t seed);

  void power_spectral_density(const String& psd_name, const Real& param = 0.);

  const RealVector& compute_sample();

  const RealMatrix& compute_samples(size_t num_ifft_samples);

private:

  //
  //- Heading: Utility routines
  //

  /// perform iFFT using Shinozuka-Deodatis algorithm
  void compute_sample_shinozuka_deodatis();
  /// perform iFFT using Grigoriu algorithm
  void compute_sample_grigoriu();
  /// use DFFTPACK or FFTW to map B vector into the i-th inverseSamples vector
  void compute_ifft_sample_set(ComplexVector& ifft_vector);

  /// deallocate data allocated in initialize()
  void finalize();

  //
  //- Heading: Data
  //

  /// iFFT approach: IFFT_SD (Shinozuka-Deodatis) or IFFT_G (Grigoriu)
  short fourierMethod;

  /// counter for number of IFFT time domain realizations
  size_t ifftSampleCntr;

  /// sequence of standard deviations computed from psdSequence
  RealVector sigmaSequence;

  /// the complex vector passed through FFTW/DFFTPACK for IFFT conversion
  /// from frequency domain to time domain
  ComplexVector ifftVector;

  /// first LHS parameter vector (means for generate_normal_samples(),
  /// lower bounds for generate_uniform_samples())
  RealVector lhsParam1;
  /// second LHS parameter vector (std devs for generate_normal_samples(),
  /// upper bounds for generate_uniform_samples())
  RealVector lhsParam2;
  /// random variable samples used in Grigoriu (num_terms x 2 variables)
  /// and Shinozuka-Deodatis (num_terms x 1 variable) algorithms
  RealMatrix lhsSamples;

#ifdef HAVE_FFTW
  /// Plan cache for FFTW
  fftw_plan fftwPlan;
#endif  // HAVE_FFTW
};


inline FourierInverseTransformation::
FourierInverseTransformation(const String& data_trans_type): ifftSampleCntr(0)
{
  if (data_trans_type == "inverse_fourier_shinozuka_deodatis")
    fourierMethod = IFFT_SD;
  else if (data_trans_type == "inverse_fourier_grigoriu")
    fourierMethod = IFFT_G;
  else {
    PCerr << "Error: bad data transformation type in "
	  << "FourierInverseTransformation." << std::endl;
    abort_handler(-1);
  }
}


inline FourierInverseTransformation::~FourierInverseTransformation()
{ finalize(); }

} // namespace Pecos

#endif
