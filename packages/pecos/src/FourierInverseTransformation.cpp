/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "FourierInverseTransformation.hpp"

#ifdef HAVE_DFFTPACK
#define ZFFTI_F77 F77_FUNC(zffti,ZFFTI)
extern "C" void ZFFTI_F77(int& n, Pecos::Real* wsave);

#define ZFFTB_F77 F77_FUNC(zfftb,ZFFTB)
extern "C" void ZFFTB_F77(int& n, const void* complex_array,Pecos::Real* wsave);
#endif

#include <algorithm>

static const char rcsId[]="@(#) $Id: FourierInverseTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

//#define DEBUG

namespace Pecos {


void FourierInverseTransformation::
initialize(const Real& total_t, const Real& w_bar, size_t seed)
{
  InverseTransformation::initialize(total_t, w_bar, seed);

  size_t i, num_terms = omegaSequence.length();
  ifftVector.sizeUninitialized(num_terms);

  switch (fourierMethod) {
  case IFFT_SD: // Generate num_terms LHS samples for Psi ~ iid U(0, 2.*PI).
    lhsSamples.shapeUninitialized(num_terms, 1);
    lhsParam1.sizeUninitialized(1);
    lhsParam2.sizeUninitialized(1);
    lhsParam1[0] = 0.;    // lower bound
    lhsParam2[0] = 2.*PI; // upper bound
    break;
  case IFFT_G:  // Generate num_terms LHS samples for V, W ~ iid N(0,1).
    lhsSamples.shapeUninitialized(num_terms, 2);
    lhsParam1.sizeUninitialized(2);
    lhsParam2.sizeUninitialized(2);
    lhsParam1[0] = lhsParam1[1] = 0.; // zero means
    lhsParam2[0] = lhsParam2[1] = 1.; // unit std deviations
    break;
  }

#ifdef HAVE_FFTW
  // For in-place transformation, second and third args deref the same pointer
  fftwPlan = fftw_plan_dft_1d(num_terms, (fftw_complex*)ifftVector.values(),
    (fftw_complex*)ifftVector.values(), FFTW_BACKWARD, FFTW_MEASURE);
#endif  // HAVE_FFTW
}


void FourierInverseTransformation::finalize()
{
#ifdef HAVE_FFTW
  fftw_destroy_plan(fftwPlan);
#endif  // HAVE_FFTW
}


/** Augments InverseTransformation::power_spectral_density()
    definition to include local data initialization (sigmaSequence). */
void FourierInverseTransformation::
power_spectral_density(const String& psd_name, const Real& param)
{
  InverseTransformation::power_spectral_density(psd_name, param);

  // sigma_i^2 = psd(omega_i)*deltaOmega_i
  size_t i, num_terms = psdSequence.length();
  sigmaSequence.sizeUninitialized(num_terms);
  for (i=0; i<num_terms; i++)
    sigmaSequence[i] = std::sqrt(psdSequence[i]*deltaOmega);
}


const RealVector& FourierInverseTransformation::compute_sample()
{
  size_t i, num_terms = omegaSequence.length();
  inverseSample.sizeUninitialized(num_terms);

  switch (fourierMethod) {
  case IFFT_SD:
    compute_sample_shinozuka_deodatis(); break;
  case IFFT_G:
    compute_sample_grigoriu();           break;
  }

  // FFTW and DFFTPACK return unnormalized IFFTs
  for (i=0; i<num_terms; i++)
    inverseSample[i] = ifftVector[i].real(); //= num_terms*ifftVector[i].real();
  ifftSampleCntr++;

  return inverseSample;
}


const RealMatrix& FourierInverseTransformation::
compute_samples(size_t num_ifft_samples)
{
  size_t i, num_terms = omegaSequence.length();
  inverseSamples.shapeUninitialized(num_ifft_samples, num_terms);

  for (ifftSampleCntr=0; ifftSampleCntr<num_ifft_samples; ifftSampleCntr++) {
    switch (fourierMethod) {
    case IFFT_SD:
      compute_sample_shinozuka_deodatis(); break;
    case IFFT_G:
      compute_sample_grigoriu();           break;
    }

    // FFTW and DFFTPACK return unnormalized IFFTs
    for (i=0; i<num_terms; i++)
      inverseSamples(ifftSampleCntr,i) = ifftVector[i].real();
                           //= num_terms*ifftVector[i].real();
  }

  return inverseSamples;
}


void FourierInverseTransformation::compute_sample_shinozuka_deodatis()
{
  // Function to generate num_ifft_samples independent samples of a zero-mean,
  // stationary, real-valued Gaussian process using the FFT algorithm.
  // The model is (Shinozuka and Deodatis):
  //
  //          m
  // Xm(t) = SUM sqrt(2) * s_k * cos(v_k*t + Psi_k)
  //         k=1
  //
  // where  Psi_k ~ iid U(0,2 pi) are random variables
  //        v_k   = (k-1)*delv is the frequency discretization
  //        s_k^2 = g(v_k)*delv is a discretization of one-sided PSD g(v)

  /* MATLAB code:
  // sample Psi ~ U(0,2*pi)
  rand('seed',seed);
  Psi=2*pi*rand(m,num_ifft_samples);

  // form num_ifft_samples samples of B vector
  IMAG=std::sqrt(-1);
  for i=1:num_ifft_samples,
    B(:,i)=std::sqrt(2)*s.*exp(IMAG*Psi(:,i));
  end

  // use ifft to get samples of X
  X=m*real(ifft(B));
  */

  size_t i, num_terms = omegaSequence.length();
  if (ifftSampleCntr)
    lhsSampler.advance_seed_sequence();
  lhsSampler.generate_uniform_samples(lhsParam1, lhsParam2, num_terms,
				      lhsSamples);

  for (i=0; i<num_terms; i++) {
    //Real A = sigmaSequence[i]*std::sqrt(2.);
    //ifftVector[i] = std::complex<Real>(A*cos(Psi_i), A*sin(Psi_i)); // Euler
    ifftVector[i] = std::polar(sigmaSequence[i]*std::sqrt(2.), lhsSamples(i,0));
  }

  compute_ifft_sample_set(ifftVector); // ifftVector: freq -> time domain
}


void FourierInverseTransformation::compute_sample_grigoriu()
{
  // Function to generate num_ifft_samples independent samples of a zero-mean,
  // stationary, real-valued Gaussian process using the FFT algorithm.
  // The model is (Grigoriu):
  //
  //          m
  // Xm(t) = SUM s_k * [V_k * cos(v_k*t) + W_k * sin(v_k*t) ]
  //         k=1
  //
  // where  V_k,W_k ~ iid N(0,1) are random variables
  //        v_k     = (k-1)*delv is the frequency discretization
  //        s_k^2   = g(v_k)*delv is a discretization of one-sided PSD g(v)

  /* MATLAB code:
  // form num_ifft_samples samples of B vector
  randn('seed',seed);
  IMAG=std::sqrt(-1);
  for i=1:num_ifft_samples,
    V=randn(m,1);W=randn(m,1);  // V, W ~ iid N(0,1)
    A=s.*std::sqrt(V.^2 + W.^2);// A ~ Rayleigh
    Psi=-atan2(W,V);            // Psi ~ U(-pi,pi)
    B(:,i)=A.*exp(IMAG*Psi);
  end
  // use ifft to get samples of X
  X=m*real(ifft(B));
  */

  size_t i, num_terms = omegaSequence.length();
  RealVector empty_rv;
  RealSymMatrix empty_correl;
  if (ifftSampleCntr)
    lhsSampler.advance_seed_sequence();
  lhsSampler.generate_normal_samples(lhsParam1, lhsParam2, empty_rv, empty_rv,
				     num_terms, empty_correl, lhsSamples);

  for (i=0; i<num_terms; i++) {
    const Real& v_i = lhsSamples(i, 0);
    const Real& w_i = lhsSamples(i, 1);
    //Real A = sigmaSequence[i]*std::sqrt(v_i*v_i + w_i*w_i); // A ~ Rayleigh
    //Real Psi = -std::atan2(w_i, v_i);                       // Psi ~ U(-pi,pi)
    //ifftVector[i] = std::complex<Real>(A*cos(Psi), A*sin(Psi)); // Euler
    ifftVector[i] = std::polar(sigmaSequence[i]*std::sqrt(v_i*v_i + w_i*w_i),
			      -std::atan2(w_i, v_i));
  }

  compute_ifft_sample_set(ifftVector); // ifftVector: freq -> time domain
}


void FourierInverseTransformation::
compute_ifft_sample_set(ComplexVector& ifft_vector)
{
  int num_terms = omegaSequence.length();
#ifdef DEBUG
  for (size_t i=0; i<num_terms; i++)
    PCout << "Freq ifft_vector[" << i << "] = (" << ifft_vector[i].real()
	  << ", " << ifft_vector[i].imag() << ")\n";
#endif // DEBUG

#ifdef HAVE_FFTW
  // default FFT package
  fftw_execute(fftwPlan);
#elif HAVE_DFFTPACK
  // fallback FFT package
  Real* wsave = new Real [4*num_terms+15];
  ZFFTI_F77(num_terms, wsave);
  ZFFTB_F77(num_terms, ifft_vector.values(), wsave); // transforms in place
  delete [] wsave;
#else
  PCerr << "Error: FFTW or DFFTPACK required for inverse FFT." << std::endl;
  abort_handler(-1);
#endif

#ifdef DEBUG
  for (size_t i=0; i<num_terms; i++)
    PCout << "Time ifft_vector[" << i << "] = (" << ifft_vector[i].real()
	  << ", " << ifft_vector[i].imag() << ")\n";
#endif // DEBUG
}

} // namespace Pecos
