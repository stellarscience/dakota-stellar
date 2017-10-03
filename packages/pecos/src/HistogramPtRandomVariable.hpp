/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 HistogramPtRandomVariable
//- Description: Encapsulates random variable data and utilities
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef HISTOGRAM_PT_RANDOM_VARIABLE_HPP
#define HISTOGRAM_PT_RANDOM_VARIABLE_HPP

#include "RandomVariable.hpp"

namespace Pecos {


/// Derived random variable class for gumbel random variables.

/** Manages int, string, or real PtPairs mappings.  At most, one
    mapping is active at a time. */

class HistogramPtRandomVariable: public RandomVariable
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  HistogramPtRandomVariable();
  /// alternate constructor
  HistogramPtRandomVariable(const IntRealMap& i_prs);
  /// alternate constructor
  HistogramPtRandomVariable(const StringRealMap& s_prs);
  /// alternate constructor
  HistogramPtRandomVariable(const RealRealMap& r_prs);
  /// destructor
  ~HistogramPtRandomVariable();

  //
  //- Heading: Virtual function redefinitions
  //

  Real mean() const;
  //Real median() const;
  Real mode() const;
  Real standard_deviation() const;
  Real variance() const;
  
  RealRealPair moments() const;
  RealRealPair bounds() const;

  Real coefficient_of_variation() const;

  //
  //- Heading: Member functions
  //

  void update(const IntRealMap& i_prs);
  void update(const StringRealMap& s_prs);
  void update(const RealRealMap& r_prs);

  //
  //- Heading: Static member functions (global utilities)
  //

  static void moments_from_params(const IntRealMap& i_prs,
				  Real& mean, Real& std_dev);
  static void moments_from_params(const StringRealMap& s_prs,
				  Real& mean, Real& std_dev);
  static void moments_from_params(const RealRealMap& r_prs,
				  Real& mean, Real& std_dev);

protected:

  //
  //- Heading: Data
  //

  /// value-count pairs for int values within a set
  IntRealMap intPtPairs;
  /// value-count pairs for string values within a set
  StringRealMap stringPtPairs;
  /// value-count pairs for real values within a set
  RealRealMap realPtPairs;
};


inline HistogramPtRandomVariable::HistogramPtRandomVariable():
  RandomVariable(BaseConstructor())
{ }


inline HistogramPtRandomVariable::
HistogramPtRandomVariable(const IntRealMap& i_prs):
  RandomVariable(BaseConstructor()), intPtPairs(i_prs)
{ ranVarType = HISTOGRAM_PT_INT; }


inline HistogramPtRandomVariable::
HistogramPtRandomVariable(const StringRealMap& s_prs):
  RandomVariable(BaseConstructor()), stringPtPairs(s_prs)
{ ranVarType = HISTOGRAM_PT_STRING; }


inline HistogramPtRandomVariable::
HistogramPtRandomVariable(const RealRealMap& r_prs):
  RandomVariable(BaseConstructor()), realPtPairs(r_prs)
{ ranVarType = HISTOGRAM_PT_REAL; }


inline HistogramPtRandomVariable::~HistogramPtRandomVariable()
{ }


inline Real HistogramPtRandomVariable::mean() const
{ return moments().first; }


//inline Real HistogramPtRandomVariable::median() const
//{ return inverse_cdf(.5); } // default


inline Real HistogramPtRandomVariable::mode() const
{
  Real mode, mode_cnt;
  switch (ranVarType) {
  case HISTOGRAM_PT_INT: {
    IRMCIter cit = intPtPairs.begin();
    mode = (Real)cit->first; mode_cnt = cit->second; ++cit;
    for (; cit!=intPtPairs.end(); ++cit)
      if (cit->second > mode_cnt)
        { mode = (Real)cit->first; mode_cnt = cit->second; }
    break;
  }
  case HISTOGRAM_PT_STRING: {
    SRMCIter cit = stringPtPairs.begin();
    mode = 0.; mode_cnt = cit->second; ++cit;
    for (size_t index=1; cit!=stringPtPairs.end(); ++cit, ++index)
      if (cit->second > mode_cnt)
        { mode = (Real)index; mode_cnt = cit->second; }
    break;
  }
  case HISTOGRAM_PT_REAL: {
    RRMCIter cit = realPtPairs.begin();
    mode = cit->first; mode_cnt = cit->second; ++cit;
    for (; cit!=realPtPairs.end(); ++cit)
      if (cit->second > mode_cnt)
        { mode = cit->first; mode_cnt = cit->second; }
    break;
  }
  }
  return mode;
}


inline Real HistogramPtRandomVariable::standard_deviation() const
{ return moments().second; }


inline Real HistogramPtRandomVariable::variance() const
{ Real std_dev = moments().second; return std_dev * std_dev; }


inline RealRealPair HistogramPtRandomVariable::moments() const
{
  Real mean, std_dev;
  switch (ranVarType) {
  case HISTOGRAM_PT_INT:
    moments_from_params(intPtPairs,    mean, std_dev); break;
  case HISTOGRAM_PT_STRING:
    moments_from_params(stringPtPairs, mean, std_dev); break;
  case HISTOGRAM_PT_REAL:
    moments_from_params(realPtPairs,   mean, std_dev); break;
  }    
  return RealRealPair(mean, std_dev);
}


inline RealRealPair HistogramPtRandomVariable::bounds() const
{
  Real l_bnd, u_bnd;
  switch (ranVarType) {
  case HISTOGRAM_PT_INT:    // moments and bounds based on values
    l_bnd = (Real)intPtPairs.begin()->first;
    u_bnd = (Real)(--intPtPairs.end())->first;    break;
  case HISTOGRAM_PT_STRING: // moments and bounds based on indices
    l_bnd = 0.; u_bnd = (Real)(stringPtPairs.size() - 1); break;
  case HISTOGRAM_PT_REAL:   // moments and bounds based on values
    l_bnd = realPtPairs.begin()->first;
    u_bnd = (--realPtPairs.end())->first;         break;
  }    
  return RealRealPair(l_bnd, u_bnd);
}


inline Real HistogramPtRandomVariable::coefficient_of_variation() const
{ RealRealPair mom = moments(); return mom.second / mom.first; }


inline void HistogramPtRandomVariable::update(const IntRealMap& i_prs)
{
  intPtPairs = i_prs; stringPtPairs.clear(); realPtPairs.clear();
  ranVarType = HISTOGRAM_PT_INT;
}


inline void HistogramPtRandomVariable::update(const StringRealMap& s_prs)
{
  stringPtPairs = s_prs; intPtPairs.clear(); realPtPairs.clear();
  ranVarType = HISTOGRAM_PT_STRING;
}


inline void HistogramPtRandomVariable::update(const RealRealMap& r_prs)
{
  realPtPairs = r_prs; intPtPairs.clear(); stringPtPairs.clear(); 
  ranVarType = HISTOGRAM_PT_REAL;
}


/// for integer-valued histogram, return a real-valued mean and std dev
inline void HistogramPtRandomVariable::
moments_from_params(const IntRealMap& i_prs, Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = std_dev = 0.;
  Real val, count, prod;
  IRMCIter cit = i_prs.begin(), cit_end = i_prs.end();
  for ( ; cit != cit_end; ++cit) {
    val = cit->first; count = cit->second; prod = count * val;
    mean    += prod;
    std_dev += prod * val;
  }
  std_dev = std::sqrt(std_dev - mean * mean);
}


/// for string variables, define the mean as the count-weighted mean
/// of a zero-based index
inline void HistogramPtRandomVariable::
moments_from_params(const StringRealMap& s_prs, Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = std_dev = 0.;
  Real val, count, prod;
  size_t index = 0;
  SRMCIter cit = s_prs.begin(), cit_end = s_prs.end();
  for ( ; cit != cit_end; ++cit, ++index) {
    val = index; count = cit->second; prod = count * val;
    mean    += prod;
    std_dev += prod * val;
  }
  std_dev = std::sqrt(std_dev - mean * mean);
}


/// return the mean and standard deviation of a real-valued point histogram
inline void HistogramPtRandomVariable::
moments_from_params(const RealRealMap& r_prs, Real& mean, Real& std_dev)
{
  // in point case, (x,y) and (x,c) are equivalent since bins have zero-width.
  // assume normalization (counts sum to 1.).
  mean = std_dev = 0.;
  Real val, count, prod;
  RRMCIter cit = r_prs.begin(), cit_end = r_prs.end();
  for ( ; cit != cit_end; ++cit) {
    val = cit->first; count = cit->second; prod = count * val;
    mean    += prod;
    std_dev += prod * val;
  }
  std_dev = std::sqrt(std_dev - mean * mean);
}

} // namespace Pecos

#endif
