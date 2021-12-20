/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HistogramBinRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HISTOGRAM_BIN_RANDOM_VARIABLE_HPP
#define HISTOGRAM_BIN_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for continuous histogram random variables.

/** A skyline PDF is defined using valueProbPairs (bins can be of unequal
    width and probability densities are employed rather than counts). */

class HistogramBinRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HistogramBinRandomVariable();
  /// alternate constructor
  HistogramBinRandomVariable(const RealRealMap& bin_prs);
  /// destructor
  ~HistogramBinRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real cdf(Real x) const;
  Real ccdf(Real x) const;
  Real inverse_cdf(Real p_cdf) const;
  Real inverse_ccdf(Real p_ccdf) const;

  Real pdf(Real x) const;
  Real pdf_gradient(Real x) const;
  Real pdf_hessian(Real x) const;

  Real mean() const;
  //Real median() const;
  Real mode() const;
  //Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  RealRealPair distribution_bounds() const;

  Real coefficient_of_variation() const;

  void pull_parameter(short dist_param, RealRealMap& val) const;
  void push_parameter(short dist_param, const RealRealMap& val);

  void copy_parameters(const RandomVariable& rv);

  //
  //- Heading: Member functions
  //

  void update(const RealRealMap& bin_prs);

  //
  //- Heading: Static member functions (global utilities)
  //

  static Real pdf(Real x, const RealRealMap& bin_prs);
  static Real cdf(Real x, const RealRealMap& bin_prs);
  static Real ccdf(Real x, const RealRealMap& bin_prs);
  static Real inverse_cdf(Real p_cdf, const RealRealMap& bin_prs);
  static Real inverse_ccdf(Real p_ccdf, const RealRealMap& bin_prs);

  static Real pdf(Real x, const RealVector& bin_prs);
  //static Real cdf(Real x, const RealVector& bin_prs);
  //static Real inverse_cdf(Real cdf, const RealVector& bin_prs);

  static Real mode(const RealRealMap& bin_prs);

  static void moments_from_params(const RealRealMap& bin_prs,
				  Real& mean, Real& std_dev);
  static void central_moments_from_params(const RealRealMap& bin_prs,
					  Real& mean, Real& var);

protected:

  //
  //- Heading: Data
  //

  /// pairings of abscissas to bin counts
  /// Note: assigned counts are assumed to be normalized (e.g., by Dakota NIDR)
  RealRealMap valueProbPairs;
};


inline HistogramBinRandomVariable::HistogramBinRandomVariable():
  RandomVariable(BaseConstructor())
{ ranVarType = HISTOGRAM_BIN; }


inline HistogramBinRandomVariable::
HistogramBinRandomVariable(const RealRealMap& bin_prs):
  RandomVariable(BaseConstructor()), valueProbPairs(bin_prs)
{ ranVarType = HISTOGRAM_BIN; }


inline HistogramBinRandomVariable::~HistogramBinRandomVariable()
{ }


inline void HistogramBinRandomVariable::
pull_parameter(short dist_param, RealRealMap& val) const
{
  switch (dist_param) {
  case H_BIN_PAIRS:  val = valueProbPairs;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HistogramBinRandomVariable::pull_parameter(RRM)."<< std::endl;
    abort_handler(-1); break;
  }
}


inline void HistogramBinRandomVariable::
push_parameter(short dist_param, const RealRealMap& val)
{
  switch (dist_param) {
  case H_BIN_PAIRS:  valueProbPairs = val;  break;
  default:
    PCerr << "Error: update failure for distribution parameter " << dist_param
	  << " in HistogramBinRandomVariable::push_parameter(RRM)."<< std::endl;
    abort_handler(-1); break;
  }
}


inline void HistogramBinRandomVariable::
copy_parameters(const RandomVariable& rv)
{ rv.pull_parameter(H_BIN_PAIRS, valueProbPairs); }


inline Real HistogramBinRandomVariable::cdf(Real x, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x <= cit->first)
    return 0.;
  else if (x >= (--bin_prs.end())->first)
    return 1.;
  else {
    Real p_cdf = 0., prob, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; prob = cit->second; ++cit;
      upr = cit->first;
      if (x <= upr) {
	p_cdf += prob * (x - lwr);
	break;
      }
      else
	p_cdf += prob * (upr - lwr);
    }
    return p_cdf;
  }
}


inline Real HistogramBinRandomVariable::cdf(Real x) const
{ return cdf(x, valueProbPairs); }


inline Real HistogramBinRandomVariable::ccdf(Real x, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x <= cit->first)
    return 1.;
  else if (x >= (--bin_prs.end())->first)
    return 0.;
  else {
    Real p_ccdf = 1., prob, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; prob = cit->second; ++cit;
      upr = cit->first;
      if (x < upr) {
	p_ccdf -= prob * (x - lwr);
	break;
      }
      else
	p_ccdf -= prob * (upr - lwr);
    }
    return p_ccdf;
  }
}


inline Real HistogramBinRandomVariable::ccdf(Real x) const
{ return ccdf(x, valueProbPairs); }


inline Real HistogramBinRandomVariable::
inverse_cdf(Real p_cdf, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (p_cdf <= 0.)
    return cit->first;                 // lower bound abscissa
  else if (p_cdf >= 1.)
    return (--bin_prs.end())->first; // upper bound abscissa
  else {
    Real upr_cdf = 0., prob, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; prob = cit->second; ++cit;
      upr = cit->first;
      upr_cdf += prob * (upr - lwr);
      if (p_cdf <= upr_cdf)
	return upr - (upr_cdf - p_cdf) / prob;
    }
    // If still no return due to numerical roundoff, return upper bound
    return (--bin_prs.end())->first; // upper bound abscissa
  }
}


inline Real HistogramBinRandomVariable::inverse_cdf(Real p_cdf) const
{ return inverse_cdf(p_cdf, valueProbPairs); }


inline Real HistogramBinRandomVariable::
inverse_ccdf(Real p_ccdf, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (p_ccdf >= 1.)
    return cit->first;                // lower bound abscissa
  else if (p_ccdf <= 0.)
    return (--bin_prs.end())->first; // upper bound abscissa
  else {
    Real upr_ccdf = 1., prob, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; prob = cit->second; ++cit;
      upr = cit->first;
      upr_ccdf -= prob * (upr - lwr);
      if (p_ccdf > upr_ccdf)
	return upr - (p_ccdf - upr_ccdf) / prob;
    }
    // If still no return due to numerical roundoff, return upper bound
    return (--bin_prs.end())->first; // upper bound abscissa
  }
}


inline Real HistogramBinRandomVariable::inverse_ccdf(Real p_ccdf) const
{ return inverse_ccdf(p_ccdf, valueProbPairs); }


inline Real HistogramBinRandomVariable::pdf(Real x, const RealRealMap& bin_prs)
{
  // The PDF is discontinuous at the bin bounds.  To resolve this, this
  // function employs an (arbitrary) convention: closed/inclusive lower bound
  // and open/exclusive upper bound.
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x < cit->first || x >= (--bin_prs.end())->first) // outside support
    return 0.;
  else {
    Real upr, prob;
    for (i=0; i<num_bins; ++i) {
      prob = cit->second;  ++cit;  upr = cit->first;
      if (x < upr) // return density
	return prob;
    }
    return 0.;
  }
}


inline Real HistogramBinRandomVariable::pdf(Real x) const
{ return pdf(x, valueProbPairs); }


inline Real HistogramBinRandomVariable::pdf_gradient(Real x) const
{ return 0.; }


inline Real HistogramBinRandomVariable::pdf_hessian(Real x) const
{ return 0.; }


inline void HistogramBinRandomVariable::
central_moments_from_params(const RealRealMap& bin_prs, Real& mean, Real& var)
{
  // in bin case, (x,y) and (x,c) are not equivalent since bins have non-zero
  // width -> assume (x,c) count/frequency-based (already converted from (x,y)
  // skyline/density-based) with normalization (counts sum to 1.)
  Real sum1 = 0., sum2 = 0., lwr, prob, count, upr, clu;
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  for (i=0; i<num_bins; ++i) {
    lwr   = cit->first;  prob = cit->second;  ++cit;  upr = cit->first;
    count = prob * (upr - lwr);  clu = count * (lwr + upr);
    sum1 += clu;
    sum2 += upr * clu + count * lwr * lwr;
  }
  mean = sum1 / 2.;
  var  = sum2 / 3. - mean * mean;
}


inline void HistogramBinRandomVariable::
moments_from_params(const RealRealMap& bin_prs, Real& mean, Real& std_dev)
{
  Real var;
  central_moments_from_params(bin_prs, mean, var);
  std_dev = std::sqrt(var);
}


inline Real HistogramBinRandomVariable::mean() const
{
  Real sum = 0., lwr, prob, upr;
  size_t i, num_bins = valueProbPairs.size() - 1;
  RRMCIter cit = valueProbPairs.begin();
  for (i=0; i<num_bins; ++i) {
    lwr  = cit->first;  prob = cit->second;  ++cit;
    upr  = cit->first;
    sum += prob * (upr * upr - lwr * lwr); // = count * (lwr + upr)
  }
  return sum / 2.;
}


inline Real HistogramBinRandomVariable::mode(const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  Real mode = cit->first, mode_p = 0., prob, lwr, count, upr;
  for (i=0; i<num_bins; ++i) {
    lwr = cit->first;  prob = cit->second;  ++cit;
    upr = cit->first;
    if (prob > mode_p)
      { mode_p = prob; mode = (lwr + upr) / 2.; }
  }
  return mode;
}


inline Real HistogramBinRandomVariable::mode() const
{ return mode(valueProbPairs); }


inline Real HistogramBinRandomVariable::variance() const
{
  Real mean, var;
  central_moments_from_params(valueProbPairs, mean, var);
  return var;
}


inline RealRealPair HistogramBinRandomVariable::moments() const
{
  Real mean, std_dev;
  moments_from_params(valueProbPairs, mean, std_dev);
  return RealRealPair(mean, std_dev);
}


inline RealRealPair HistogramBinRandomVariable::distribution_bounds() const
{
  return RealRealPair(valueProbPairs.begin()->first,
		      (--valueProbPairs.end())->first);
}


inline Real HistogramBinRandomVariable::coefficient_of_variation() const
{
  Real sum1 = 0., sum2 = 0., lwr, prob, count, upr, clu;
  size_t i, num_bins = valueProbPairs.size() - 1;
  RRMCIter cit = valueProbPairs.begin();
  for (i=0; i<num_bins; ++i) {
    lwr   = cit->first;  prob = cit->second;  ++cit;  upr = cit->first;
    count = prob * (upr - lwr);  clu = count * (lwr + upr);
    sum1 += clu;
    sum2 += upr * clu + count * lwr * lwr;
  }
  // COV^2 = (sum2 / 3. - sum1 * sum1 / 4) / (sum1 * sum1 / 4)
  return std::sqrt(4. * sum2 / (3. * sum1 * sum1) - 1.);
}


inline void HistogramBinRandomVariable::update(const RealRealMap& bin_prs)
{ valueProbPairs = bin_prs; }


inline Real HistogramBinRandomVariable::pdf(Real x, const RealVector& bin_prs)
{
  // Need this case to be fast for usage in NumericGenOrthogPolynomial...
  // Cleaner, but induces unnecessary copy overhead:
  //RealRealMap bin_prs_rrm;
  //copy_data(bin_prs_rv, bin_prs_rrm);
  //return HistogramBinRandomVariable::pdf(x, bin_prs_rrm);

  size_t num_bins = bin_prs.length() / 2 - 1;
  if (x < bin_prs[0] || x >= bin_prs[2*num_bins])
    return 0.;
  else {
    Real upr;
    for (size_t i=0; i<num_bins; ++i) {
      upr = bin_prs[2*i+2];
      if (x < upr) // return density = count / (upr - lwr);
	return bin_prs[2*i+1];
    }
    return 0.;
  }
}


/*  Versions using value-count pairs:

inline Real HistogramBinRandomVariable::pdf(Real x, const RealRealMap& bin_prs)
{
  // The PDF is discontinuous at the bin bounds.  To resolve this, this
  // function employs an (arbitrary) convention: closed/inclusive lower bound
  // and open/exclusive upper bound.
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x < cit->first || x >= (--bin_prs.end())->first) // outside support
    return 0.;
  else {
    Real lwr, upr, count;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      if (x < upr) // return density
	return count / (upr - lwr);
    }
    return 0.;
  }
}


inline Real HistogramBinRandomVariable::cdf(Real x, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x <= cit->first)
    return 0.;
  else if (x >= (--bin_prs.end())->first)
    return 1.;
  else {
    Real p_cdf = 0., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      if (x <= upr) {
	p_cdf += count * (x - lwr) / (upr - lwr);
	break;
      }
      else
	p_cdf += count;
    }
    return p_cdf;
  }
}


inline Real HistogramBinRandomVariable::ccdf(Real x, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (x <= cit->first)
    return 1.;
  else if (x >= (--bin_prs.end())->first)
    return 0.;
  else {
    Real p_ccdf = 1., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      if (x < upr) {
	p_ccdf -= count * (x - lwr) / (upr - lwr);
	break;
      }
      else
	p_ccdf -= count;
    }
    return p_ccdf;
  }
}


inline Real HistogramBinRandomVariable::
inverse_cdf(Real p_cdf, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (p_cdf <= 0.)
    return cit->first;                  // lower bound abscissa
  else if (p_cdf >= 1.)
    return (--bin_prs.end())->first; // upper bound abscissa
  else {
    Real upr_cdf = 0., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      upr_cdf += count;
      if (p_cdf <= upr_cdf)
	return upr - (upr_cdf - p_cdf) / count * (upr - lwr);
    }
    // If still no return due to numerical roundoff, return upper bound
    return (--bin_prs.end())->first; // upper bound abscissa
  }
}


inline Real HistogramBinRandomVariable::
inverse_ccdf(Real p_ccdf, const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  if (p_ccdf >= 1.)
    return cit->first;                // lower bound abscissa
  else if (p_ccdf <= 0.)
    return (--bin_prs.end())->first; // upper bound abscissa
  else {
    Real upr_ccdf = 1., count, upr, lwr;
    for (i=0; i<num_bins; ++i) {
      lwr = cit->first; count = cit->second; ++cit;
      upr = cit->first;
      upr_ccdf -= count;
      if (p_ccdf > upr_ccdf)
	return upr - (p_ccdf - upr_ccdf) / count * (upr - lwr);
    }
    // If still no return due to numerical roundoff, return upper bound
    return (--bin_prs.end())->first; // upper bound abscissa
  }
}


inline Real HistogramBinRandomVariable::cdf(Real x, const RealVector& bin_prs)
{
  // Cleaner, but induces unnecessary copy overhead:
  //RealRealMap bin_prs_rrm;
  //copy_data(bin_prs_rv, bin_prs_rrm);
  //return HistogramBinRandomVariable::cdf(x, bin_prs_rrm);

  size_t num_bins = bin_prs.length() / 2 - 1;
  if (x <= bin_prs[0])
    return 0.;
  else if (x >= bin_prs[2*num_bins])
    return 1.;
  else {
    Real p_cdf = 0., count, upr, lwr;
    for (int i=0; i<num_bins; ++i) {
      count = bin_prs[2*i+1]; upr = bin_prs[2*i+2];
      if (x < upr) {
        lwr = bin_prs[2*i];
	p_cdf += count * (x - lwr) / (upr - lwr);
	break;
      }
      else
	p_cdf += count;
    }
    return p_cdf;
  }
}


inline Real HistogramBinRandomVariable::
inverse_cdf(Real p_cdf, const RealVector& bin_prs)
{
  // Cleaner, but induces unnecessary copy overhead:
  //RealRealMap bin_prs_rrm;
  //copy_data(bin_prs_rv, bin_prs_rrm);
  //return HistogramBinRandomVariable::inverse_cdf(p_cdf, bin_prs_rrm);

  size_t num_bins = bin_prs.length() / 2 - 1;
  if (p_cdf <= 0.)
    return bin_prs[0];          // lower bound abscissa
  else if (p_cdf >= 1.)
    return bin_prs[2*num_bins]; // upper bound abscissa
  else {
    Real upr_cdf = 0., count, upr, lwr;
    for (int i=0; i<num_bins; ++i) {
      count = bin_prs[2*i+1];
      upr_cdf += count;
      if (p_cdf < upr_cdf) {
        upr = bin_prs[2*i+2]; lwr = bin_prs[2*i];
	return upr - (upr_cdf - p_cdf) / count * (upr - lwr);
      }
    }
    // If still no return due to numerical roundoff, return upper bound
    return bin_prs[2*num_bins]; // upper bound abscissa
  }
}


inline Real HistogramBinRandomVariable::mean() const
{
  Real sum = 0., lwr, count, upr;
  size_t i, num_bins = valueProbPairs.size() - 1;
  RRMCIter cit = valueProbPairs.begin();
  for (i=0; i<num_bins; ++i) {
    lwr  = cit->first; count = cit->second; ++cit;
    upr  = cit->first;
    sum += count * (lwr + upr);
  }
  return sum / 2.;
}


inline Real HistogramBinRandomVariable::mode(const RealRealMap& bin_prs)
{
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  Real mode = cit->first, mode_p = 0., density, lwr, count, upr;
  for (i=0; i<num_bins; ++i) {
    lwr  = cit->first; count   = cit->second; ++cit;
    upr  = cit->first; density = count / (upr - lwr);
    if (density > mode_p)
      { mode_p = density; mode = (lwr+upr)/2.; }
  }
  return mode;
}


inline Real HistogramBinRandomVariable::variance() const
{
  Real sum1 = 0., sum2 = 0., lwr, count, upr, clu;
  size_t i, num_bins = valueProbPairs.size() - 1;
  RRMCIter cit = valueProbPairs.begin();
  for (i=0; i<num_bins; ++i) {
    lwr  = cit->first; count = cit->second; ++cit;
    upr  = cit->first; clu   = count * (lwr + upr);
    sum1 += clu;
    sum2 += upr*clu + count*lwr*lwr;
  }
  return sum2 / 3. - sum1 * sum1 / 4.;// raw - mean * mean, for mean = sum1 / 2
}


inline void HistogramBinRandomVariable::
moments_from_params(const RealRealMap& bin_prs, Real& mean, Real& std_dev)
{
  // in bin case, (x,y) and (x,c) are not equivalent since bins have non-zero
  // width -> assume (x,c) count/frequency-based (already converted from (x,y)
  // skyline/density-based) with normalization (counts sum to 1.)
  Real sum1 = 0., sum2 = 0., lwr, count, upr, clu;
  size_t i, num_bins = bin_prs.size() - 1;
  RRMCIter cit = bin_prs.begin();
  for (i=0; i<num_bins; ++i) {
    lwr   = cit->first; count = cit->second; ++cit;
    upr   = cit->first; clu   = count * (lwr + upr);
    sum1 += clu;
    sum2 += upr*clu + count*lwr*lwr;
  }
  mean    = sum1 / 2.;
  std_dev = std::sqrt(sum2 / 3. - mean * mean);
}


inline Real HistogramBinRandomVariable::pdf(Real x, const RealVector& bin_prs)
{
  // Need this case to be fast for usage in NumericGenOrthogPolynomial...
  // Cleaner, but induces unnecessary copy overhead:
  //RealRealMap bin_prs_rrm;
  //copy_data(bin_prs_rv, bin_prs_rrm);
  //return HistogramBinRandomVariable::pdf(x, bin_prs_rrm);

  size_t num_bins = bin_prs.length() / 2 - 1;
  if (x < bin_prs[0] || x >= bin_prs[2*num_bins])
    return 0.;
  else {
    Real upr;
    for (int i=0; i<num_bins; ++i) {
      upr = bin_prs[2*i+2];
      if (x < upr) // return density = count / (upr - lwr);
	return bin_prs[2*i+1] / (upr - bin_prs[2*i]);
    }
    return 0.;
  }
}
*/

} // namespace Pecos

#endif
