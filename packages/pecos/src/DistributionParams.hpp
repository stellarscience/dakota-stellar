/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef DISTRIBUTION_PARAMS_HPP
#define DISTRIBUTION_PARAMS_HPP

#include "pecos_data_types.hpp"
#include "RandomVariable.hpp"


namespace Pecos {

/// The representation of a set of aleatory distribution parameters.
/// This representation, or body, may be shared by multiple
/// AleatoryDistParams handle instances.

/** The AleatoryDistParams/AleatoryDistParamsRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class AleatoryDistParamsRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class AleatoryDistParams;

private:

  //
  //- Heading: Private member functions
  //

  /// default constructor
  AleatoryDistParamsRep();
  /// constructor
  AleatoryDistParamsRep(const RealVector& nuv_means,
    const RealVector& nuv_std_devs,     const RealVector& nuv_l_bnds,
    const RealVector& nuv_u_bnds,       const RealVector& lnuv_means,
    const RealVector& lnuv_std_devs,    const RealVector& lnuv_lambdas,
    const RealVector& lnuv_zetas,       const RealVector& lnuv_err_facts,
    const RealVector& lnuv_l_bnds,      const RealVector& lnuv_u_bnds,
    const RealVector& uuv_l_bnds,       const RealVector& uuv_u_bnds,
    const RealVector& luuv_l_bnds,      const RealVector& luuv_u_bnds,
    const RealVector& tuv_modes,        const RealVector& tuv_l_bnds,
    const RealVector& tuv_u_bnds,       const RealVector& euv_betas,
    const RealVector& beuv_alphas,      const RealVector& beuv_betas,
    const RealVector& beuv_l_bnds,      const RealVector& beuv_u_bnds,
    const RealVector& gauv_alphas,      const RealVector& gauv_betas,
    const RealVector& guuv_alphas,      const RealVector& guuv_betas,
    const RealVector& fuv_alphas,       const RealVector& fuv_betas,
    const RealVector& wuv_alphas,       const RealVector& wuv_betas,
    const RealRealMapArray& hbuv_prs,   const RealVector& puv_lambdas,
    const RealVector& biuv_p_per_tr,    const IntVector& biuv_num_trials, 
    const RealVector& nbuv_p_per_tr,    const IntVector& nbuv_num_trials, 
    const RealVector& geuv_p_per_tr,    const IntVector& hguv_tot_pop,
    const IntVector& hguv_sel_pop,      const IntVector& hguv_num_drawn,
    const IntRealMapArray& hpuiv_prs,   const StringRealMapArray& hpusv_prs,
    const RealRealMapArray& hpurv_prs,  const RealSymMatrix& uv_corr);
  /// destructor
  ~AleatoryDistParamsRep();

  //
  //- Heading: Private data members
  //

  /// normal uncertain variable means
  RealVector normalMeans;
  /// normal uncertain variable standard deviations
  RealVector normalStdDevs;
  /// normal uncertain variable lower bounds
  RealVector normalLowerBnds;
  /// normal uncertain variable upper bounds
  RealVector normalUpperBnds;
  /// lognormal uncertain variable means
  RealVector lognormalMeans;
  /// lognormal uncertain variable standard deviations
  RealVector lognormalStdDevs;
  /// lognormal uncertain variable lambdas
  RealVector lognormalLambdas;
  /// lognormal uncertain variable zetas
  RealVector lognormalZetas;
  /// lognormal uncertain variable error factors
  RealVector lognormalErrFacts;
  /// lognormal uncertain variable lower bounds
  RealVector lognormalLowerBnds;
  /// lognormal uncertain variable upper bounds
  RealVector lognormalUpperBnds;
  /// uniform uncertain variable lower bounds
  RealVector uniformLowerBnds;
  /// uniform uncertain variable upper bounds
  RealVector uniformUpperBnds;
  /// loguniform uncertain variable lower bounds
  RealVector loguniformLowerBnds;
  /// loguniform uncertain variable upper bounds
  RealVector loguniformUpperBnds;
  /// triangular uncertain variable modes
  RealVector triangularModes;
  /// triangular uncertain variable lower bounds
  RealVector triangularLowerBnds;
  /// triangular uncertain variable upper bounds
  RealVector triangularUpperBnds;
  /// exponential uncertain variable betas
  RealVector exponentialBetas;
  /// beta uncertain variable alphas
  RealVector betaAlphas;
  /// beta uncertain variable betas
  RealVector betaBetas;
  /// beta uncertain variable lower bounds
  RealVector betaLowerBnds;
  /// beta uncertain variable upper bounds
  RealVector betaUpperBnds;
  /// gamma uncertain variable alphas
  RealVector gammaAlphas;
  /// gamma uncertain variable betas
  RealVector gammaBetas;
  /// gumbel uncertain variable alphas
  RealVector gumbelAlphas;
  /// gumbel uncertain variable betas
  RealVector gumbelBetas;
  /// frechet uncertain variable alphas
  RealVector frechetAlphas;
  /// frechet uncertain variable betas
  RealVector frechetBetas;
  /// weibull uncertain variable alphas
  RealVector weibullAlphas;
  /// weibull uncertain variable betas
  RealVector weibullBetas;
  /// histogram uncertain (x,y) bin pairs (continuous linear histogram)
  RealRealMapArray histogramBinPairs;

  /// poisson uncertain variable lambdas
  RealVector poissonLambdas;
  /// binomial uncertain variable probability per trial
  RealVector binomialProbPerTrial;
  /// binomial uncertain variable numbers of trials
  IntVector binomialNumTrials;
  /// negative binomial uncertain variable probability per trial
  RealVector negBinomialProbPerTrial;
  /// negative binomial uncertain variable numbers of trials
  IntVector negBinomialNumTrials;
  /// geometric uncertain variable probability per trial
  RealVector geometricProbPerTrial;
  /// hypergeometric uncertain variable numbers in total population
  IntVector hyperGeomTotalPopulation;
  /// hypergeometric uncertain variable numbers in selected population
  IntVector hyperGeomSelectedPopulation;
  /// hypergeometric uncertain variable numbers failed in population
  IntVector hyperGeomNumDrawn;
  /// histogram uncertain (i,y) point pairs (discrete histogram)
  IntRealMapArray histogramPointIntPairs;
  /// histogram uncertain (s,y) point pairs (discrete histogram)
  StringRealMapArray histogramPointStringPairs;
  /// histogram uncertain (x,y) point pairs (discrete histogram)
  RealRealMapArray histogramPointRealPairs;

  /// uncertain variable correlation matrix (rank correlations for sampling
  /// and correlation coefficients for reliability)
  RealSymMatrix uncertainCorrelations;

  /// number of handle objects sharing adpRep
  int referenceCount;
};


inline AleatoryDistParamsRep::AleatoryDistParamsRep() : referenceCount(1)
{ }


inline AleatoryDistParamsRep::
AleatoryDistParamsRep(const RealVector& nuv_means,
  const RealVector& nuv_std_devs,     const RealVector& nuv_l_bnds,
  const RealVector& nuv_u_bnds,       const RealVector& lnuv_means,
  const RealVector& lnuv_std_devs,    const RealVector& lnuv_lambdas,
  const RealVector& lnuv_zetas,       const RealVector& lnuv_err_facts,
  const RealVector& lnuv_l_bnds,      const RealVector& lnuv_u_bnds,
  const RealVector& uuv_l_bnds,       const RealVector& uuv_u_bnds,
  const RealVector& luuv_l_bnds,      const RealVector& luuv_u_bnds,
  const RealVector& tuv_modes,        const RealVector& tuv_l_bnds,
  const RealVector& tuv_u_bnds,       const RealVector& euv_betas,
  const RealVector& beuv_alphas,      const RealVector& beuv_betas,
  const RealVector& beuv_l_bnds,      const RealVector& beuv_u_bnds,
  const RealVector& gauv_alphas,      const RealVector& gauv_betas,
  const RealVector& guuv_alphas,      const RealVector& guuv_betas,
  const RealVector& fuv_alphas,       const RealVector& fuv_betas,
  const RealVector& wuv_alphas,       const RealVector& wuv_betas,
  const RealRealMapArray& hbuv_prs,   const RealVector& puv_lambdas,
  const RealVector& biuv_p_per_tr,    const IntVector& biuv_num_trials, 
  const RealVector& nbuv_p_per_tr,    const IntVector& nbuv_num_trials, 
  const RealVector& geuv_p_per_tr,    const IntVector& hguv_tot_pop,
  const IntVector& hguv_sel_pop,      const IntVector& hguv_num_drawn,
  const IntRealMapArray& hpuiv_prs,   const StringRealMapArray& hpusv_prs,
  const RealRealMapArray& hpurv_prs,  const RealSymMatrix& uv_corr):
  normalMeans(nuv_means), normalStdDevs(nuv_std_devs),
  normalLowerBnds(nuv_l_bnds), normalUpperBnds(nuv_u_bnds),
  lognormalMeans(lnuv_means), lognormalStdDevs(lnuv_std_devs),
  lognormalLambdas(lnuv_lambdas), lognormalZetas(lnuv_zetas),
  lognormalErrFacts(lnuv_err_facts), lognormalLowerBnds(lnuv_l_bnds),
  lognormalUpperBnds(lnuv_u_bnds), uniformLowerBnds(uuv_l_bnds),
  uniformUpperBnds(uuv_u_bnds), loguniformLowerBnds(luuv_l_bnds),
  loguniformUpperBnds(luuv_u_bnds), triangularModes(tuv_modes),
  triangularLowerBnds(tuv_l_bnds), triangularUpperBnds(tuv_u_bnds),
  exponentialBetas(euv_betas), betaAlphas(beuv_alphas), betaBetas(beuv_betas),
  betaLowerBnds(beuv_l_bnds), betaUpperBnds(beuv_u_bnds),
  gammaAlphas(gauv_alphas), gammaBetas(gauv_betas), gumbelAlphas(guuv_alphas),
  gumbelBetas(guuv_betas), frechetAlphas(fuv_alphas), frechetBetas(fuv_betas),
  weibullAlphas(wuv_alphas), weibullBetas(wuv_betas),
  histogramBinPairs(hbuv_prs), poissonLambdas(puv_lambdas),
  binomialProbPerTrial(biuv_p_per_tr), binomialNumTrials(biuv_num_trials), 
  negBinomialProbPerTrial(nbuv_p_per_tr), negBinomialNumTrials(nbuv_num_trials),
  geometricProbPerTrial(geuv_p_per_tr), hyperGeomTotalPopulation(hguv_tot_pop),
  hyperGeomSelectedPopulation(hguv_sel_pop),  hyperGeomNumDrawn(hguv_num_drawn),
  histogramPointIntPairs(hpuiv_prs), histogramPointStringPairs(hpusv_prs),
  histogramPointRealPairs(hpurv_prs), uncertainCorrelations(uv_corr),
  referenceCount(1)
{ }


inline AleatoryDistParamsRep::~AleatoryDistParamsRep()
{ }


/// The representation of a set of epistemic distribution parameters.
/// This representation, or body, may be shared by multiple
/// EpistemicDistParams handle instances.

/** The EpistemicDistParams/EpistemicDistParamsRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class EpistemicDistParamsRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class EpistemicDistParams;

private:

  //
  //- Heading: Private member functions
  //

  /// default constructor
  EpistemicDistParamsRep();
  /// constructor
  EpistemicDistParamsRep(const RealRealPairRealMapArray& ceuv_bpas,
    const IntIntPairRealMapArray& deuiv_bpas,
    const IntRealMapArray&    deusiv_vals_probs,
    const StringRealMapArray& deussv_vals_probs,
    const RealRealMapArray&   deusrv_vals_probs);
  /// destructor
  ~EpistemicDistParamsRep();

  //
  //- Heading: Private data members
  //

  /*
  /// basic probability values for continuous interval uncertain variables
  RealVectorArray contIntervalProbs;
  /// lower bounds for continuous interval uncertain variables
  RealVectorArray contIntervalLowerBnds;
  /// upper bounds for continuous interval uncertain variables
  RealVectorArray contIntervalUpperBnds;
  */
  /// basic probability assignments (continuous interval bounds + probability)
  /// for continuous interval uncertain variables
  RealRealPairRealMapArray contIntervalBPA; // consider multimap from exp elicit

  /*
  /// basic probability values for discrete interval uncertain variables
  RealVectorArray discIntervalProbs;
  /// lower bounds for discrete interval uncertain variables
  IntVectorArray discIntervalLowerBnds;
  /// upper bounds for discrete interval uncertain variables
  IntVectorArray discIntervalUpperBnds;
  */
  /// basic probability assignments (discrete range bounds + probability)
  /// for discrete interval uncertain variables
  IntIntPairRealMapArray discIntervalBPA; // consider multimap from exp elicit

  /// admissible values and basic probability assignments for discrete
  /// uncertain set integer variables
  IntRealMapArray discSetIntValsProbs;
  /// admissible values and basic probability assignments for discrete
  /// uncertain set string variables
  StringRealMapArray discSetStringValsProbs;
  /// admissible values and basic probability assignments for discrete
  /// uncertain set real variables
  RealRealMapArray discSetRealValsProbs;

  /// number of handle objects sharing edpRep
  int referenceCount;
};


inline EpistemicDistParamsRep::EpistemicDistParamsRep() : referenceCount(1)
{ }


inline EpistemicDistParamsRep::
EpistemicDistParamsRep(const RealRealPairRealMapArray& ceuv_bpas,
		       const IntIntPairRealMapArray& deuiv_bpas,
		       const IntRealMapArray& deusiv_vals_probs,
		       const StringRealMapArray& deussv_vals_probs,
		       const RealRealMapArray& deusrv_vals_probs):
  contIntervalBPA(ceuv_bpas), discIntervalBPA(deuiv_bpas),
  discSetIntValsProbs(deusiv_vals_probs),
  discSetStringValsProbs(deussv_vals_probs),
  discSetRealValsProbs(deusrv_vals_probs), referenceCount(1)
{ }


inline EpistemicDistParamsRep::~EpistemicDistParamsRep()
{ }


/// Container class encapsulating distribution parameters for aleatory
/// random variables.

/** This class consolidates aleatory distribution data and simplifies the
    APIs that require distribution parameters. */

class AleatoryDistParams
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  AleatoryDistParams();
  /// standard constructor
  AleatoryDistParams(const RealVector& nuv_means,
    const RealVector& nuv_std_devs,     const RealVector& nuv_l_bnds,
    const RealVector& nuv_u_bnds,       const RealVector& lnuv_means,
    const RealVector& lnuv_std_devs,    const RealVector& lnuv_lambdas,
    const RealVector& lnuv_zetas,       const RealVector& lnuv_err_facts,
    const RealVector& lnuv_l_bnds,      const RealVector& lnuv_u_bnds,
    const RealVector& uuv_l_bnds,       const RealVector& uuv_u_bnds,
    const RealVector& luuv_l_bnds,      const RealVector& luuv_u_bnds,
    const RealVector& tuv_modes,        const RealVector& tuv_l_bnds,
    const RealVector& tuv_u_bnds,       const RealVector& euv_betas,
    const RealVector& beuv_alphas,      const RealVector& beuv_betas,
    const RealVector& beuv_l_bnds,      const RealVector& beuv_u_bnds,
    const RealVector& gauv_alphas,      const RealVector& gauv_betas,
    const RealVector& guuv_alphas,      const RealVector& guuv_betas,
    const RealVector& fuv_alphas,       const RealVector& fuv_betas,
    const RealVector& wuv_alphas,       const RealVector& wuv_betas,
    const RealRealMapArray& hbuv_prs,   const RealVector& puv_lambdas,
    const RealVector& biuv_p_per_tr,    const IntVector& biuv_num_trials, 
    const RealVector& nbuv_p_per_tr,    const IntVector& nbuv_num_trials, 
    const RealVector& geuv_p_per_tr,    const IntVector& hguv_tot_pop,
    const IntVector& hguv_sel_pop,      const IntVector& hguv_num_drawn,
    const IntRealMapArray& hpuiv_prs,   const StringRealMapArray& hpusv_prs,
    const RealRealMapArray& hpurv_prs,  const RealSymMatrix& uv_corr);
  /// copy constructor
  AleatoryDistParams(const AleatoryDistParams& adp);
  /// destructor
  ~AleatoryDistParams();

  /// assignment operator
  AleatoryDistParams& operator=(const AleatoryDistParams& adp);

  //
  //- Heading: member functions
  //

  /// return total number of continuous aleatory uncertain variables
  size_t cauv()  const;

  /// return total number of discrete aleatory uncertain integer variables
  size_t dauiv() const;
  /// return total number of discrete aleatory uncertain string variables
  size_t dausv() const;
  /// return total number of discrete aleatory uncertain real variables
  size_t daurv() const;
  /// return total number of discrete aleatory uncertain variables
  size_t dauv()   const;

  /// size the number of normal uncertain variables
  void nuv(size_t num_nuv);
  /// size the number of lognormal uncertain variables (mean, std deviation)
  void lnuv_ms(size_t num_lnuv);
  /// size the number of lognormal uncertain variables (lambda, zeta)
  void lnuv_lz(size_t num_lnuv);
  /// size the number of lognormal uncertain variables (mean, error factor)
  void lnuv_me(size_t num_lnuv);
  /// size the number of uniform uncertain variables
  void uuv(size_t num_uuv);
  /// size the number of loguniform uncertain variables
  void luuv(size_t num_luuv);
  /// size the number of triangular uncertain variables
  void tuv(size_t num_tuv);
  /// size the number of exponential uncertain variables
  void euv(size_t num_euv);
  /// size the number of beta uncertain variables
  void buv(size_t num_buv);
  /// size the number of gamma uncertain variables
  void gauv(size_t num_gauv);
  /// size the number of gumbel uncertain variables
  void guuv(size_t num_guuv);
  /// size the number of frechet uncertain variables
  void fuv(size_t num_nuv);
  /// size the number of weibull uncertain variables
  void wuv(size_t num_nuv);
  /// size the number of histogram bin uncertain variables
  void hbuv(size_t num_nuv);
  /// size the number of poisson uncertain variables
  void puv(size_t num_nuv);
  /// size the number of binominal uncertain variables
  void biuv(size_t num_nuv);
  /// size the number of negative binominal uncertain variables
  void nbuv(size_t num_nuv);
  /// size the number of geometric uncertain variables
  void geuv(size_t num_nuv);
  /// size the number of hypergeometric uncertain variables
  void hguv(size_t num_nuv);
  /// size the number of histogram point uncertain variables
  void hpuv(size_t num_hpiuv, size_t num_hpsuv, size_t num_hpruv);

  /// deep copy (as opposed to operator= shallow copy)
  void copy(const AleatoryDistParams& adp);
  /// data update (no changes to representation (unless null))
  void update(const AleatoryDistParams& adp);
  /// partial data update of the distribution data not affected by an
  /// x->u variable transformation
  void update_partial(const AleatoryDistParams& adp_x,
		      const std::vector<RandomVariable>& x_ran_vars,
		      const ShortArray& u_types);

  /// return the normal uncertain variable means
  const RealVector& normal_means() const;
  /// return the ith normal uncertain variable mean
  const Real& normal_mean(size_t i) const;
  /// set the normal uncertain variable means
  void normal_means(const RealVector& n_means);
  /// set the ith normal uncertain variable mean
  void normal_mean(const Real& n_mean, size_t i);
  /// return the normal uncertain variable standard deviations
  const RealVector& normal_std_deviations() const;
  /// return the ith normal uncertain variable standard deviation
  const Real& normal_std_deviation(size_t i) const;
  /// set the normal uncertain variable standard deviations
  void normal_std_deviations(const RealVector& n_std_devs);
  /// set the ith normal uncertain variable standard deviation
  void normal_std_deviation(const Real& n_std_dev, size_t i);
  /// return the normal uncertain variable lower bounds
  const RealVector& normal_lower_bounds() const;
  /// return the ith normal uncertain variable lower bound
  const Real& normal_lower_bound(size_t i) const;
  /// set the normal uncertain variable lower bounds
  void normal_lower_bounds(const RealVector& n_lower_bnds);
  /// set the ith normal uncertain variable lower bound
  void normal_lower_bound(const Real& n_lower_bnd, size_t i);
  /// return the normal uncertain variable upper bounds
  const RealVector& normal_upper_bounds() const;
  /// return the ith normal uncertain variable upper bound
  const Real& normal_upper_bound(size_t i) const;
  /// set the normal uncertain variable upper bounds
  void normal_upper_bounds(const RealVector& n_upper_bnds);
  /// set the ith normal uncertain variable upper bound
  void normal_upper_bound(const Real& n_upper_bnd, size_t i);

  /// return the lognormal uncertain variable means
  const RealVector& lognormal_means() const;
  /// return the ith lognormal uncertain variable mean
  const Real& lognormal_mean(size_t i) const;
  /// set the lognormal uncertain variable means
  void lognormal_means(const RealVector& ln_means);
  /// set the ith lognormal uncertain variable mean
  void lognormal_mean(const Real& ln_mean, size_t i);
  /// return the lognormal uncertain variable standard deviations
  const RealVector& lognormal_std_deviations() const;
  /// return the ith lognormal uncertain variable standard deviation
  const Real& lognormal_std_deviation(size_t i) const;
  /// set the lognormal uncertain variable standard deviations
  void lognormal_std_deviations(const RealVector& ln_std_devs);
  /// set the ith lognormal uncertain variable standard deviation
  void lognormal_std_deviation(const Real& ln_std_dev, size_t i);
  /// return the lognormal uncertain variable lambdas
  const RealVector& lognormal_lambdas() const;
  /// return the ith lognormal uncertain variable lambda
  const Real& lognormal_lambda(size_t i) const;
  /// set the lognormal uncertain variable lambdas
  void lognormal_lambdas(const RealVector& ln_lambdas);
  /// set the ith lognormal uncertain variable lambda
  void lognormal_lambda(const Real& ln_lambda, size_t i);
  /// return the lognormal uncertain variable zetas
  const RealVector& lognormal_zetas() const;
  /// return the ith lognormal uncertain variable zeta
  const Real& lognormal_zeta(size_t i) const;
  /// set the lognormal uncertain variable zetas
  void lognormal_zetas(const RealVector& ln_std_devs);
  /// set the ith lognormal uncertain variable zeta
  void lognormal_zeta(const Real& ln_std_dev, size_t i);
  /// return the lognormal uncertain variable error factors
  const RealVector& lognormal_error_factors() const;
  /// return the ith lognormal uncertain variable error factor
  const Real& lognormal_error_factor(size_t i) const;
  /// set the lognormal uncertain variable error factors
  void lognormal_error_factors(const RealVector& ln_err_facts);
  /// set the ith lognormal uncertain variable error factor
  void lognormal_error_factor(const Real& ln_err_fact, size_t i);
  /// return the lognormal uncertain variable lower bounds
  const RealVector& lognormal_lower_bounds() const;
  /// return the ith lognormal uncertain variable lower bound
  const Real& lognormal_lower_bound(size_t i) const;
  /// set the lognormal uncertain variable lower bounds
  void lognormal_lower_bounds(const RealVector& ln_lower_bnds);
  /// set the ith lognormal uncertain variable lower bound
  void lognormal_lower_bound(const Real& ln_lower_bnd, size_t i);
  /// return the lognormal uncertain variable upper bounds
  const RealVector& lognormal_upper_bounds() const;
  /// return the ith lognormal uncertain variable upper bound
  const Real& lognormal_upper_bound(size_t i) const;
  /// set the lognormal uncertain variable upper bounds
  void lognormal_upper_bounds(const RealVector& ln_upper_bnds);
  /// set the ith lognormal uncertain variable upper bound
  void lognormal_upper_bound(const Real& ln_upper_bnd, size_t i);

  /// return the uniform uncertain variable lower bounds
  const RealVector& uniform_lower_bounds() const;
  /// return the ith uniform uncertain variable lower bound
  const Real& uniform_lower_bound(size_t i) const;
  /// set the uniform uncertain variable lower bounds
  void uniform_lower_bounds(const RealVector& u_lower_bnds);
  /// set the ith uniform uncertain variable lower bound
  void uniform_lower_bound(const Real& u_lower_bnd, size_t i);
  /// return the uniform uncertain variable upper bounds
  const RealVector& uniform_upper_bounds() const;
  /// return the ith uniform uncertain variable upper bound
  const Real& uniform_upper_bound(size_t i) const;
  /// set the uniform uncertain variable upper bounds
  void uniform_upper_bounds(const RealVector& u_upper_bnds);
  /// set the ith uniform uncertain variable upper bound
  void uniform_upper_bound(const Real& u_upper_bnd, size_t i);

  /// return the loguniform uncertain variable lower bounds
  const RealVector& loguniform_lower_bounds() const;
  /// return the ith loguniform uncertain variable lower bound
  const Real& loguniform_lower_bound(size_t i) const;
  /// set the loguniform uncertain variable lower bounds
  void loguniform_lower_bounds(const RealVector& lu_lower_bnds);
  /// set the ith loguniform uncertain variable lower bound
  void loguniform_lower_bound(const Real& lu_lower_bnd, size_t i);
  /// return the loguniform uncertain variable upper bounds
  const RealVector& loguniform_upper_bounds() const;
  /// return the ith loguniform uncertain variable upper bound
  const Real& loguniform_upper_bound(size_t i) const;
  /// set the loguniform uncertain variable upper bounds
  void loguniform_upper_bounds(const RealVector& lu_upper_bnds);
  /// set the ith loguniform uncertain variable upper bound
  void loguniform_upper_bound(const Real& lu_upper_bnd, size_t i);

  /// return the triangular uncertain variable modes
  const RealVector& triangular_modes() const;
  /// return the ith triangular uncertain variable mode
  const Real& triangular_mode(size_t i) const;
  /// set the triangular uncertain variable modes
  void triangular_modes(const RealVector& t_modes);
  /// set the ith triangular uncertain variable mode
  void triangular_mode(const Real& t_mode, size_t i);
  /// return the triangular uncertain variable lower bounds
  const RealVector& triangular_lower_bounds() const;
  /// return the ith triangular uncertain variable lower bound
  const Real& triangular_lower_bound(size_t i) const;
  /// set the triangular uncertain variable lower bounds
  void triangular_lower_bounds(const RealVector& t_lower_bnds);
  /// set the ith triangular uncertain variable lower bound
  void triangular_lower_bound(const Real& t_lower_bnd, size_t i);
  /// return the triangular uncertain variable upper bounds
  const RealVector& triangular_upper_bounds() const;
  /// return the ith triangular uncertain variable upper bound
  const Real& triangular_upper_bound(size_t i) const;
  /// set the triangular uncertain variable upper bounds
  void triangular_upper_bounds(const RealVector& t_upper_bnds);
  /// set the ith triangular uncertain variable upper bound
  void triangular_upper_bound(const Real& t_upper_bnd, size_t i);

  /// return the exponential uncertain variable beta parameters
  const RealVector& exponential_betas() const;
  /// return the ith exponential uncertain variable beta parameter
  const Real& exponential_beta(size_t i) const;
  /// set the exponential uncertain variable beta parameters
  void exponential_betas(const RealVector& e_betas);
  /// set the ith exponential uncertain variable beta parameter
  void exponential_beta(const Real& e_beta, size_t i);

  /// return the beta uncertain variable alphas
  const RealVector& beta_alphas() const;
  /// return the ith beta uncertain variable alpha
  const Real& beta_alpha(size_t i) const;
  /// set the beta uncertain variable alphas
  void beta_alphas(const RealVector& b_alphas);
  /// set the ith beta uncertain variable alpha
  void beta_alpha(const Real& b_alpha, size_t i);
  /// return the beta uncertain variable betas
  const RealVector& beta_betas() const;
  /// return the ith beta uncertain variable beta
  const Real& beta_beta(size_t i) const;
  /// set the beta uncertain variable betas
  void beta_betas(const RealVector& b_betas);
  /// set the ith beta uncertain variable beta
  void beta_beta(const Real& b_beta, size_t i);
  /// return the beta uncertain variable lower bounds
  const RealVector& beta_lower_bounds() const;
  /// return the ith beta uncertain variable lower bound
  const Real& beta_lower_bound(size_t i) const;
  /// set the beta uncertain variable lower bounds
  void beta_lower_bounds(const RealVector& b_lower_bnds);
  /// set the ith beta uncertain variable lower bound
  void beta_lower_bound(const Real& b_lower_bnd, size_t i);
  /// return the beta uncertain variable upper bounds
  const RealVector& beta_upper_bounds() const;
  /// return the ith beta uncertain variable upper bound
  const Real& beta_upper_bound(size_t i) const;
  /// set the beta uncertain variable upper bounds
  void beta_upper_bounds(const RealVector& b_upper_bnds);
  /// set the ith beta uncertain variable upper bound
  void beta_upper_bound(const Real& b_upper_bnd, size_t i);

  /// return the gamma uncertain variable alpha parameters
  const RealVector& gamma_alphas() const;
  /// return the ith gamma uncertain variable alpha parameter
  const Real& gamma_alpha(size_t i) const;
  /// set the gamma uncertain variable alpha parameters
  void gamma_alphas(const RealVector& ga_alphas);
  /// set the ith gamma uncertain variable alpha parameter
  void gamma_alpha(const Real& ga_alpha, size_t i);
  /// return the gamma uncertain variable beta parameters
  const RealVector& gamma_betas() const;
  /// return the ith gamma uncertain variable beta parameter
  const Real& gamma_beta(size_t i) const;
  /// set the gamma uncertain variable beta parameters
  void gamma_betas(const RealVector& ga_betas);
  /// set the ith gamma uncertain variable beta parameter
  void gamma_beta(const Real& ga_beta, size_t i);

  /// return the gumbel uncertain variable alphas
  const RealVector& gumbel_alphas() const;
  /// return the ith gumbel uncertain variable alpha
  const Real& gumbel_alpha(size_t i) const;
  /// set the gumbel uncertain variable alphas
  void gumbel_alphas(const RealVector& gu_alphas);
  /// set the ith gumbel uncertain variable alpha
  void gumbel_alpha(const Real& gu_alpha, size_t i);
  /// return the gumbel uncertain variable betas
  const RealVector& gumbel_betas() const;
  /// return the ith gumbel uncertain variable beta
  const Real& gumbel_beta(size_t i) const;
  /// set the gumbel uncertain variable betas
  void gumbel_betas(const RealVector& gu_betas);
  /// set the ith gumbel uncertain variable beta
  void gumbel_beta(const Real& gu_beta, size_t i);

  /// return the frechet uncertain variable alpha parameters
  const RealVector& frechet_alphas() const;
  /// return the ith frechet uncertain variable alpha parameter
  const Real& frechet_alpha(size_t i) const;
  /// set the frechet uncertain variable alpha parameters
  void frechet_alphas(const RealVector& f_alphas);
  /// set the ith frechet uncertain variable alpha parameter
  void frechet_alpha(const Real& f_alpha, size_t i);
  /// return the frechet uncertain variable beta parameters
  const RealVector& frechet_betas() const;
  /// return the ith frechet uncertain variable beta parameter
  const Real& frechet_beta(size_t i) const;
  /// set the frechet uncertain variable beta parameters
  void frechet_betas(const RealVector& f_betas);
  /// set the ith frechet uncertain variable beta parameter
  void frechet_beta(const Real& f_beta, size_t i);

  /// return the weibull uncertain variable alpha parameters
  const RealVector& weibull_alphas() const;
  /// return the ith weibull uncertain variable alpha parameter
  const Real& weibull_alpha(size_t i) const;
  /// set the weibull uncertain variable alpha parameters
  void weibull_alphas(const RealVector& w_alphas);
  /// set the ith weibull uncertain variable alpha parameter
  void weibull_alpha(const Real& w_alpha, size_t i);
  /// return the weibull uncertain variable beta parameters
  const RealVector& weibull_betas() const;
  /// return the ith weibull uncertain variable beta parameter
  const Real& weibull_beta(size_t i) const;
  /// set the weibull uncertain variable beta parameters
  void weibull_betas(const RealVector& w_betas);
  /// set the ith weibull uncertain variable beta parameter
  void weibull_beta(const Real& w_beta, size_t i);

  /// return the histogram uncertain bin pairs
  const RealRealMapArray& histogram_bin_pairs() const;
  /// return the ith histogram uncertain bin pair
  const RealRealMap& histogram_bin_pairs(size_t i) const;
  /// set the histogram uncertain bin pairs
  void histogram_bin_pairs(const RealRealMapArray& h_bin_pairs);
  /// set the ith histogram uncertain bin pair
  void histogram_bin_pairs(const RealRealMap& h_bin_pairs_i, size_t i);

  /// return the poisson uncertain variable lambda parameters
  const RealVector& poisson_lambdas() const;
  /// return the ith poisson uncertain variable lambda parameter
  const Real& poisson_lambda(size_t i) const;
  /// set the poisson uncertain variable lambda parameters
  void poisson_lambdas(const RealVector& p_lambdas);
  /// set the ith poisson uncertain variable lambda parameter
  void poisson_lambda(const Real& p_lambda, size_t i);

  /// return the binomial probability per each trial (p) 
  const RealVector& binomial_probability_per_trial() const;
  /// return the ith binomial probability per each trial (p) 
  const Real& binomial_probability_per_trial(size_t i) const;
  /// set the binomial probability per each trial (p) 
  void binomial_probability_per_trial(const RealVector& probs_per_trial);
  /// set the ith binomial probability per each trial (p) 
  void binomial_probability_per_trial(const Real& prob_per_trial, size_t i);
  /// return the binomial number of trials (N)
  const IntVector& binomial_num_trials() const;
  /// return the ith binomial number of trials (N)
  int binomial_num_trials(size_t i) const;
  /// set the binomial number of trials (N)
  void binomial_num_trials(const IntVector& num_trials);
  /// set the ith binomial number of trials (N)
  void binomial_num_trials(int num_trials, size_t i);

  /// return the negative binomial probability per each trial (p) 
  const RealVector& negative_binomial_probability_per_trial() const;
  /// return the ith negative binomial probability per each trial (p) 
  const Real& negative_binomial_probability_per_trial(size_t i) const;
  /// set the negative binomial probability per each trial (p) 
  void negative_binomial_probability_per_trial(
    const RealVector& probs_per_trial);
  /// set the ith negative binomial probability per each trial (p) 
  void negative_binomial_probability_per_trial(
    const Real& prob_per_trial, size_t i);
  /// return the negative binomial number of trials (N)
  const IntVector& negative_binomial_num_trials() const;
  /// return the ith negative binomial number of trials (N)
  int negative_binomial_num_trials(size_t i) const;
  /// set the negative binomial number of trials (N)
  void negative_binomial_num_trials(const IntVector& num_trials);
  /// set the ith negative binomial number of trials (N)
  void negative_binomial_num_trials(int num_trials, size_t i);

  /// return the geometric probability per each trial (p) 
  const RealVector& geometric_probability_per_trial() const;
  /// return the ith geometric probability per each trial (p) 
  const Real& geometric_probability_per_trial(size_t i) const;
  /// set the geometric probability per each trial (p) 
  void geometric_probability_per_trial(const RealVector& probs_per_trial);
  /// set the ith geometric probability per each trial (p) 
  void geometric_probability_per_trial(const Real& prob_per_trial, size_t i);

  /// return the hypergeometric number in total population 
  const IntVector& hypergeometric_total_population() const;
  /// return the ith hypergeometric number in total population 
  int hypergeometric_total_population(size_t i) const;
  /// set the hypergeometric number in total population
  void hypergeometric_total_population(const IntVector& total_pop);
  /// set the ith hypergeometric number in total population
  void hypergeometric_total_population(int total_pop, size_t i);
  /// return the hypergeometric number in selected population
  const IntVector& hypergeometric_selected_population() const;
  /// return the ith hypergeometric number in selected population
  int hypergeometric_selected_population(size_t i) const;
  /// set the hypergeometric number in selected population
  void hypergeometric_selected_population(const IntVector& sel_pop);
  /// set the ith hypergeometric number in selected population
  void hypergeometric_selected_population(int sel_pop, size_t i);
  /// return the hypergeometric number failed
  const IntVector& hypergeometric_num_drawn() const;
  /// return the ith hypergeometric number failed
  int hypergeometric_num_drawn(size_t i) const;
  /// set the hypergeometric number in total population
  void hypergeometric_num_drawn(const IntVector& num_drawn);
  /// set the ith hypergeometric number in total population
  void hypergeometric_num_drawn(int num_drawn, size_t i);

  /// return the histogram uncertain point pairs
  const IntRealMapArray& histogram_point_int_pairs() const;
  /// return the ith histogram uncertain point pair
  const IntRealMap& histogram_point_int_pairs(size_t i) const;
  /// set the histogram uncertain point pairs
  void histogram_point_int_pairs(const IntRealMapArray& h_pt_pairs);
  /// set the ith histogram uncertain point pair
  void histogram_point_int_pairs(const IntRealMap& h_pt_pairs_i, size_t i);

  /// return the histogram uncertain point pairs
  const StringRealMapArray& histogram_point_string_pairs() const;
  /// return the ith histogram uncertain point pair
  const StringRealMap& histogram_point_string_pairs(size_t i) const;
  /// set the histogram uncertain point pairs
  void histogram_point_string_pairs(const StringRealMapArray& h_pt_pairs);
  /// set the ith histogram uncertain point pair
  void histogram_point_string_pairs(const StringRealMap& h_pt_pairs_i,
				    size_t i);

  /// return the histogram uncertain point pairs
  const RealRealMapArray& histogram_point_real_pairs() const;
  /// return the ith histogram uncertain point pair
  const RealRealMap& histogram_point_real_pairs(size_t i) const;
  /// set the histogram uncertain point pairs
  void histogram_point_real_pairs(const RealRealMapArray& h_pt_pairs);
  /// set the ith histogram uncertain point pair
  void histogram_point_real_pairs(const RealRealMap& h_pt_pairs_i, size_t i);

  /// return the uncertain variable correlations
  const RealSymMatrix& uncertain_correlations() const;
  /// set the uncertain variable correlations
  void uncertain_correlations(const RealSymMatrix& uncertain_corr);

  /// function to check adpRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //

  /// pointer to the body (handle-body idiom)
  AleatoryDistParamsRep* adpRep;
};


inline AleatoryDistParams::AleatoryDistParams():
  adpRep(new AleatoryDistParamsRep())
{ }


inline AleatoryDistParams::
AleatoryDistParams(const RealVector& nuv_means,
  const RealVector& nuv_std_devs,     const RealVector& nuv_l_bnds,
  const RealVector& nuv_u_bnds,       const RealVector& lnuv_means,
  const RealVector& lnuv_std_devs,    const RealVector& lnuv_lambdas,
  const RealVector& lnuv_zetas,       const RealVector& lnuv_err_facts,
  const RealVector& lnuv_l_bnds,      const RealVector& lnuv_u_bnds,
  const RealVector& uuv_l_bnds,       const RealVector& uuv_u_bnds,
  const RealVector& luuv_l_bnds,      const RealVector& luuv_u_bnds,
  const RealVector& tuv_modes,        const RealVector& tuv_l_bnds,
  const RealVector& tuv_u_bnds,       const RealVector& euv_betas,
  const RealVector& beuv_alphas,      const RealVector& beuv_betas,
  const RealVector& beuv_l_bnds,      const RealVector& beuv_u_bnds,
  const RealVector& gauv_alphas,      const RealVector& gauv_betas,
  const RealVector& guuv_alphas,      const RealVector& guuv_betas,
  const RealVector& fuv_alphas,       const RealVector& fuv_betas,
  const RealVector& wuv_alphas,       const RealVector& wuv_betas,
  const RealRealMapArray& hbuv_prs,   const RealVector& puv_lambdas,
  const RealVector& biuv_p_per_tr,    const IntVector& biuv_num_trials, 
  const RealVector& nbuv_p_per_tr,    const IntVector& nbuv_num_trials, 
  const RealVector& geuv_p_per_tr,    const IntVector& hguv_tot_pop,
  const IntVector& hguv_sel_pop,      const IntVector& hguv_num_drawn,
  const IntRealMapArray& hpuiv_prs,   const StringRealMapArray& hpusv_prs,
  const RealRealMapArray& hpurv_prs,  const RealSymMatrix& uv_corr):
  adpRep(new AleatoryDistParamsRep(nuv_means, nuv_std_devs, nuv_l_bnds,
	 nuv_u_bnds, lnuv_means, lnuv_std_devs, lnuv_lambdas, lnuv_zetas,
	 lnuv_err_facts, lnuv_l_bnds, lnuv_u_bnds, uuv_l_bnds, uuv_u_bnds,
	 luuv_l_bnds, luuv_u_bnds, tuv_modes, tuv_l_bnds, tuv_u_bnds, euv_betas,
	 beuv_alphas, beuv_betas, beuv_l_bnds, beuv_u_bnds, gauv_alphas,
	 gauv_betas, guuv_alphas, guuv_betas, fuv_alphas, fuv_betas, wuv_alphas,
	 wuv_betas, hbuv_prs, puv_lambdas, biuv_p_per_tr, biuv_num_trials,
	 nbuv_p_per_tr, nbuv_num_trials, geuv_p_per_tr, hguv_tot_pop,
	 hguv_sel_pop, hguv_num_drawn, hpuiv_prs, hpusv_prs, hpurv_prs,
	 uv_corr))
{ }


inline AleatoryDistParams::AleatoryDistParams(const AleatoryDistParams& adp)
{
  // Increment new (no old to decrement)
  adpRep = adp.adpRep;
  if (adpRep) // Check for an assignment of NULL
    adpRep->referenceCount++;
}


inline AleatoryDistParams::~AleatoryDistParams()
{
  if (adpRep) { // Check for NULL
    --adpRep->referenceCount; // decrement
    if (adpRep->referenceCount == 0)
      delete adpRep;
  }
}


inline AleatoryDistParams& AleatoryDistParams::
operator=(const AleatoryDistParams& adp)
{
  // Decrement old
  if (adpRep) // Check for NULL
    if ( --adpRep->referenceCount == 0 ) 
      delete adpRep;
  // Increment new
  adpRep = adp.adpRep;
  if (adpRep) // Check for an assignment of NULL
    adpRep->referenceCount++;
  return *this;
}


inline size_t AleatoryDistParams::cauv() const
{
  return adpRep->normalMeans.length() +
    // only one of means / lambdas should be sized
    adpRep->lognormalMeans.length() + adpRep->lognormalLambdas.length() +
    adpRep->uniformLowerBnds.length() + adpRep->loguniformLowerBnds.length() +
    adpRep->triangularModes.length() + adpRep->exponentialBetas.length() + 
    adpRep->betaAlphas.length() + adpRep->gammaAlphas.length() +
    adpRep->gumbelAlphas.length() + adpRep->frechetAlphas.length() +
    adpRep->weibullAlphas.length() + adpRep->histogramBinPairs.size();
}


inline size_t AleatoryDistParams::dauiv() const
{
  return adpRep->poissonLambdas.length() +
    adpRep->binomialProbPerTrial.length() +
    adpRep->negBinomialProbPerTrial.length() +
    adpRep->geometricProbPerTrial.length() +
    adpRep->hyperGeomNumDrawn.length() + adpRep->histogramPointIntPairs.size();
}


inline size_t AleatoryDistParams::dausv() const
{ return adpRep->histogramPointStringPairs.size(); }


inline size_t AleatoryDistParams::daurv() const
{ return adpRep->histogramPointRealPairs.size(); }


inline size_t AleatoryDistParams::dauv()  const
{ return dauiv() + dausv() + daurv(); }


inline void AleatoryDistParams::nuv(size_t num_nuv)
{
  adpRep->normalMeans.sizeUninitialized(num_nuv);
  adpRep->normalStdDevs.sizeUninitialized(num_nuv);
  adpRep->normalLowerBnds.sizeUninitialized(num_nuv);
  adpRep->normalUpperBnds.sizeUninitialized(num_nuv);
}


inline void AleatoryDistParams::lnuv_ms(size_t num_lnuv)
{
  adpRep->lognormalMeans.sizeUninitialized(num_lnuv);
  adpRep->lognormalStdDevs.sizeUninitialized(num_lnuv);
}


inline void AleatoryDistParams::lnuv_lz(size_t num_lnuv)
{
  adpRep->lognormalLambdas.sizeUninitialized(num_lnuv);
  adpRep->lognormalZetas.sizeUninitialized(num_lnuv);
}


inline void AleatoryDistParams::lnuv_me(size_t num_lnuv)
{
  adpRep->lognormalMeans.sizeUninitialized(num_lnuv);
  adpRep->lognormalErrFacts.sizeUninitialized(num_lnuv);
}


inline void AleatoryDistParams::uuv(size_t num_uuv)
{
  adpRep->uniformLowerBnds.sizeUninitialized(num_uuv);
  adpRep->uniformUpperBnds.sizeUninitialized(num_uuv);
}


inline void AleatoryDistParams::luuv(size_t num_luuv)
{
  adpRep->loguniformLowerBnds.sizeUninitialized(num_luuv);
  adpRep->loguniformUpperBnds.sizeUninitialized(num_luuv);
}


inline void AleatoryDistParams::tuv(size_t num_tuv)
{
  adpRep->triangularModes.sizeUninitialized(num_tuv);
  adpRep->triangularLowerBnds.sizeUninitialized(num_tuv);
  adpRep->triangularUpperBnds.sizeUninitialized(num_tuv);
}


inline void AleatoryDistParams::euv(size_t num_euv)
{ adpRep->exponentialBetas.sizeUninitialized(num_euv); }


inline void AleatoryDistParams::buv(size_t num_buv)
{
  adpRep->betaAlphas.sizeUninitialized(num_buv);
  adpRep->betaBetas.sizeUninitialized(num_buv);
  adpRep->betaLowerBnds.sizeUninitialized(num_buv);
  adpRep->betaUpperBnds.sizeUninitialized(num_buv);
}


inline void AleatoryDistParams::gauv(size_t num_gauv)
{
  adpRep->gammaAlphas.sizeUninitialized(num_gauv);
  adpRep->gammaBetas.sizeUninitialized(num_gauv);
}


inline void AleatoryDistParams::guuv(size_t num_guuv)
{
  adpRep->gumbelAlphas.sizeUninitialized(num_guuv);
  adpRep->gumbelBetas.sizeUninitialized(num_guuv);
}


inline void AleatoryDistParams::fuv(size_t num_fuv)
{
  adpRep->frechetAlphas.sizeUninitialized(num_fuv);
  adpRep->frechetBetas.sizeUninitialized(num_fuv);
}


inline void AleatoryDistParams::wuv(size_t num_wuv)
{
  adpRep->weibullAlphas.sizeUninitialized(num_wuv);
  adpRep->weibullBetas.sizeUninitialized(num_wuv);
}


inline void AleatoryDistParams::hbuv(size_t num_hbuv)
{ adpRep->histogramBinPairs.resize(num_hbuv); }


inline void AleatoryDistParams::puv(size_t num_puv)
{ adpRep->poissonLambdas.sizeUninitialized(num_puv); }


inline void AleatoryDistParams::biuv(size_t num_biuv)
{
  adpRep->binomialProbPerTrial.sizeUninitialized(num_biuv);
  adpRep->binomialNumTrials.sizeUninitialized(num_biuv);
}


inline void AleatoryDistParams::nbuv(size_t num_nbuv)
{
  adpRep->negBinomialProbPerTrial.sizeUninitialized(num_nbuv);
  adpRep->negBinomialNumTrials.sizeUninitialized(num_nbuv);
}


inline void AleatoryDistParams::geuv(size_t num_geuv)
{ adpRep->geometricProbPerTrial.sizeUninitialized(num_geuv); }


inline void AleatoryDistParams::hguv(size_t num_hguv)
{
  adpRep->hyperGeomTotalPopulation.sizeUninitialized(num_hguv);
  adpRep->hyperGeomSelectedPopulation.sizeUninitialized(num_hguv);
  adpRep->hyperGeomNumDrawn.sizeUninitialized(num_hguv);
}


inline void AleatoryDistParams::
hpuv(size_t num_hpiuv, size_t num_hpsuv, size_t num_hpruv)
{
  adpRep->histogramPointIntPairs.resize(num_hpiuv);
  adpRep->histogramPointStringPairs.resize(num_hpsuv);
  adpRep->histogramPointRealPairs.resize(num_hpruv);
}


inline void AleatoryDistParams::copy(const AleatoryDistParams& adp)
{ 
  // Decrement old
  if (adpRep) // Check for NULL
    if ( --adpRep->referenceCount == 0 ) 
      delete adpRep;
  // Create new
  adpRep = new AleatoryDistParamsRep(adp.normal_means(),
    adp.normal_std_deviations(), adp.normal_lower_bounds(),
    adp.normal_upper_bounds(), adp.lognormal_means(),
    adp.lognormal_std_deviations(), adp.lognormal_lambdas(),
    adp.lognormal_zetas(), adp.lognormal_error_factors(),
    adp.lognormal_lower_bounds(), adp.lognormal_upper_bounds(),
    adp.uniform_lower_bounds(), adp.uniform_upper_bounds(),
    adp.loguniform_lower_bounds(), adp.loguniform_upper_bounds(),
    adp.triangular_modes(), adp.triangular_lower_bounds(),
    adp.triangular_upper_bounds(), adp.exponential_betas(), adp.beta_alphas(),
    adp.beta_betas(), adp.beta_lower_bounds(), adp.beta_upper_bounds(),
    adp.gamma_alphas(), adp.gamma_betas(), adp.gumbel_alphas(),
    adp.gumbel_betas(), adp.frechet_alphas(), adp.frechet_betas(),
    adp.weibull_alphas(), adp.weibull_betas(), adp.histogram_bin_pairs(),
    adp.poisson_lambdas(), adp.binomial_probability_per_trial(),
    adp.binomial_num_trials(), adp.negative_binomial_probability_per_trial(),
    adp.negative_binomial_num_trials(), adp.geometric_probability_per_trial(),
    adp.hypergeometric_total_population(),
    adp.hypergeometric_selected_population(), adp.hypergeometric_num_drawn(),
    adp.histogram_point_int_pairs(), adp.histogram_point_string_pairs(),
    adp.histogram_point_real_pairs(), adp.uncertain_correlations());
}


inline const RealVector& AleatoryDistParams::normal_means() const
{ return adpRep->normalMeans; }


inline const Real& AleatoryDistParams::normal_mean(size_t i) const
{ return adpRep->normalMeans[i]; }


inline void AleatoryDistParams::normal_means(const RealVector& n_means)
{ adpRep->normalMeans = n_means; }


inline void AleatoryDistParams::normal_mean(const Real& n_mean, size_t i)
{ adpRep->normalMeans[i] = n_mean; }


inline const RealVector& AleatoryDistParams::normal_std_deviations() const
{ return adpRep->normalStdDevs; }


inline const Real& AleatoryDistParams::normal_std_deviation(size_t i) const
{ return adpRep->normalStdDevs[i]; }


inline void AleatoryDistParams::
normal_std_deviations(const RealVector& n_std_devs)
{ adpRep->normalStdDevs = n_std_devs; }


inline void AleatoryDistParams::
normal_std_deviation(const Real& n_std_dev, size_t i)
{ adpRep->normalStdDevs[i] = n_std_dev; }


inline const RealVector& AleatoryDistParams::normal_lower_bounds() const
{ return adpRep->normalLowerBnds; }


inline const Real& AleatoryDistParams::normal_lower_bound(size_t i) const
{ return adpRep->normalLowerBnds[i]; }


inline void AleatoryDistParams::
normal_lower_bounds(const RealVector& n_lower_bnds)
{ adpRep->normalLowerBnds = n_lower_bnds; }


inline void AleatoryDistParams::
normal_lower_bound(const Real& n_lower_bnd, size_t i)
{ adpRep->normalLowerBnds[i] = n_lower_bnd; }


inline const RealVector& AleatoryDistParams::normal_upper_bounds() const
{ return adpRep->normalUpperBnds; }


inline const Real& AleatoryDistParams::normal_upper_bound(size_t i) const
{ return adpRep->normalUpperBnds[i]; }


inline void AleatoryDistParams::
normal_upper_bounds(const RealVector& n_upper_bnds)
{ adpRep->normalUpperBnds = n_upper_bnds; }


inline void AleatoryDistParams::
normal_upper_bound(const Real& n_upper_bnd, size_t i)
{ adpRep->normalUpperBnds[i] = n_upper_bnd; }


inline const RealVector& AleatoryDistParams::lognormal_means() const
{ return adpRep->lognormalMeans; }


inline const Real& AleatoryDistParams::lognormal_mean(size_t i) const
{ return adpRep->lognormalMeans[i]; }


inline void AleatoryDistParams::lognormal_means(const RealVector& ln_means)
{ adpRep->lognormalMeans = ln_means; }


inline void AleatoryDistParams::lognormal_mean(const Real& ln_mean, size_t i)
{ adpRep->lognormalMeans[i] = ln_mean; }


inline const RealVector& AleatoryDistParams::lognormal_std_deviations() const
{ return adpRep->lognormalStdDevs; }


inline const Real& AleatoryDistParams::lognormal_std_deviation(size_t i) const
{ return adpRep->lognormalStdDevs[i]; }


inline void AleatoryDistParams::
lognormal_std_deviations(const RealVector& ln_std_devs)
{ adpRep->lognormalStdDevs = ln_std_devs; }


inline void AleatoryDistParams::
lognormal_std_deviation(const Real& ln_std_dev, size_t i)
{ adpRep->lognormalStdDevs[i] = ln_std_dev; }


inline const RealVector& AleatoryDistParams::lognormal_lambdas() const
{ return adpRep->lognormalLambdas; }


inline const Real& AleatoryDistParams::lognormal_lambda(size_t i) const
{ return adpRep->lognormalLambdas[i]; }


inline void AleatoryDistParams::lognormal_lambdas(const RealVector& ln_lambdas)
{ adpRep->lognormalLambdas = ln_lambdas; }


inline void AleatoryDistParams::
lognormal_lambda(const Real& ln_lambda, size_t i)
{ adpRep->lognormalLambdas[i] = ln_lambda; }


inline const RealVector& AleatoryDistParams::lognormal_zetas() const
{ return adpRep->lognormalZetas; }


inline const Real& AleatoryDistParams::lognormal_zeta(size_t i) const
{ return adpRep->lognormalZetas[i]; }


inline void AleatoryDistParams::lognormal_zetas(const RealVector& ln_zetas)
{ adpRep->lognormalZetas = ln_zetas; }


inline void AleatoryDistParams::lognormal_zeta(const Real& ln_zeta, size_t i)
{ adpRep->lognormalZetas[i] = ln_zeta; }


inline const RealVector& AleatoryDistParams::lognormal_error_factors() const
{ return adpRep->lognormalErrFacts; }


inline const Real& AleatoryDistParams::lognormal_error_factor(size_t i) const
{ return adpRep->lognormalErrFacts[i]; }


inline void AleatoryDistParams::
lognormal_error_factors(const RealVector& ln_err_facts)
{ adpRep->lognormalErrFacts = ln_err_facts; }


inline void AleatoryDistParams::
lognormal_error_factor(const Real& ln_err_fact, size_t i)
{ adpRep->lognormalErrFacts[i] = ln_err_fact; }


inline const RealVector& AleatoryDistParams::lognormal_lower_bounds() const
{ return adpRep->lognormalLowerBnds; }


inline const Real& AleatoryDistParams::lognormal_lower_bound(size_t i) const
{ return adpRep->lognormalLowerBnds[i]; }


inline void AleatoryDistParams::
lognormal_lower_bounds(const RealVector& ln_lower_bnds)
{ adpRep->lognormalLowerBnds = ln_lower_bnds; }


inline void AleatoryDistParams::
lognormal_lower_bound(const Real& ln_lower_bnd, size_t i)
{ adpRep->lognormalLowerBnds[i] = ln_lower_bnd; }


inline const RealVector& AleatoryDistParams::lognormal_upper_bounds() const
{ return adpRep->lognormalUpperBnds; }


inline const Real& AleatoryDistParams::lognormal_upper_bound(size_t i) const
{ return adpRep->lognormalUpperBnds[i]; }


inline void AleatoryDistParams::
lognormal_upper_bounds(const RealVector& ln_upper_bnds)
{ adpRep->lognormalUpperBnds = ln_upper_bnds; }


inline void AleatoryDistParams::
lognormal_upper_bound(const Real& ln_upper_bnd, size_t i)
{ adpRep->lognormalUpperBnds[i] = ln_upper_bnd; }


inline const RealVector& AleatoryDistParams::uniform_lower_bounds() const
{ return adpRep->uniformLowerBnds; }


inline const Real& AleatoryDistParams::uniform_lower_bound(size_t i) const
{ return adpRep->uniformLowerBnds[i]; }


inline void AleatoryDistParams::
uniform_lower_bounds(const RealVector& u_lower_bnds)
{ adpRep->uniformLowerBnds = u_lower_bnds; }


inline void AleatoryDistParams::
uniform_lower_bound(const Real& u_lower_bnd, size_t i)
{ adpRep->uniformLowerBnds[i] = u_lower_bnd; }


inline const RealVector& AleatoryDistParams::uniform_upper_bounds() const
{ return adpRep->uniformUpperBnds; }


inline const Real& AleatoryDistParams::uniform_upper_bound(size_t i) const
{ return adpRep->uniformUpperBnds[i]; }


inline void AleatoryDistParams::
uniform_upper_bounds(const RealVector& u_upper_bnds)
{ adpRep->uniformUpperBnds = u_upper_bnds; }


inline void AleatoryDistParams::
uniform_upper_bound(const Real& u_upper_bnd, size_t i)
{ adpRep->uniformUpperBnds[i] = u_upper_bnd; }


inline const RealVector& AleatoryDistParams::loguniform_lower_bounds() const
{ return adpRep->loguniformLowerBnds; }


inline const Real& AleatoryDistParams::loguniform_lower_bound(size_t i) const
{ return adpRep->loguniformLowerBnds[i]; }


inline void AleatoryDistParams::
loguniform_lower_bounds(const RealVector& lu_lower_bnds)
{ adpRep->loguniformLowerBnds = lu_lower_bnds; }


inline void AleatoryDistParams::
loguniform_lower_bound(const Real& lu_lower_bnd, size_t i)
{ adpRep->loguniformLowerBnds[i] = lu_lower_bnd; }


inline const RealVector& AleatoryDistParams::loguniform_upper_bounds() const
{ return adpRep->loguniformUpperBnds; }


inline const Real& AleatoryDistParams::loguniform_upper_bound(size_t i) const
{ return adpRep->loguniformUpperBnds[i]; }


inline void AleatoryDistParams::
loguniform_upper_bounds(const RealVector& lu_upper_bnds)
{ adpRep->loguniformUpperBnds = lu_upper_bnds; }


inline void AleatoryDistParams::
loguniform_upper_bound(const Real& lu_upper_bnd, size_t i)
{ adpRep->loguniformUpperBnds[i] = lu_upper_bnd; }


inline const RealVector& AleatoryDistParams::triangular_modes() const
{ return adpRep->triangularModes; }


inline const Real& AleatoryDistParams::triangular_mode(size_t i) const
{ return adpRep->triangularModes[i]; }


inline void AleatoryDistParams::triangular_modes(const RealVector& t_modes)
{ adpRep->triangularModes = t_modes; }


inline void AleatoryDistParams::triangular_mode(const Real& t_mode, size_t i)
{ adpRep->triangularModes[i] = t_mode; }


inline const RealVector& AleatoryDistParams::triangular_lower_bounds() const
{ return adpRep->triangularLowerBnds; }


inline const Real& AleatoryDistParams::triangular_lower_bound(size_t i) const
{ return adpRep->triangularLowerBnds[i]; }


inline void AleatoryDistParams::
triangular_lower_bounds(const RealVector& t_lower_bnds)
{ adpRep->triangularLowerBnds = t_lower_bnds; }


inline void AleatoryDistParams::
triangular_lower_bound(const Real& t_lower_bnd, size_t i)
{ adpRep->triangularLowerBnds[i] = t_lower_bnd; }


inline const RealVector& AleatoryDistParams::triangular_upper_bounds() const
{ return adpRep->triangularUpperBnds; }


inline const Real& AleatoryDistParams::triangular_upper_bound(size_t i) const
{ return adpRep->triangularUpperBnds[i]; }


inline void AleatoryDistParams::
triangular_upper_bounds(const RealVector& t_upper_bnds)
{ adpRep->triangularUpperBnds = t_upper_bnds; }


inline void AleatoryDistParams::
triangular_upper_bound(const Real& t_upper_bnd, size_t i)
{ adpRep->triangularUpperBnds[i] = t_upper_bnd; }


inline const RealVector& AleatoryDistParams::exponential_betas() const
{ return adpRep->exponentialBetas; }


inline const Real& AleatoryDistParams::exponential_beta(size_t i) const
{ return adpRep->exponentialBetas[i]; }


inline void AleatoryDistParams::exponential_betas(const RealVector& e_betas)
{ adpRep->exponentialBetas = e_betas; }


inline void AleatoryDistParams::exponential_beta(const Real& e_beta, size_t i)
{ adpRep->exponentialBetas[i] = e_beta; }


inline const RealVector& AleatoryDistParams::beta_alphas() const
{ return adpRep->betaAlphas; }


inline const Real& AleatoryDistParams::beta_alpha(size_t i) const
{ return adpRep->betaAlphas[i]; }


inline void AleatoryDistParams::beta_alphas(const RealVector& b_alphas)
{ adpRep->betaAlphas = b_alphas; }


inline void AleatoryDistParams::beta_alpha(const Real& b_alpha, size_t i)
{ adpRep->betaAlphas[i] = b_alpha; }


inline const RealVector& AleatoryDistParams::beta_betas() const
{ return adpRep->betaBetas; }


inline const Real& AleatoryDistParams::beta_beta(size_t i) const
{ return adpRep->betaBetas[i]; }


inline void AleatoryDistParams::beta_betas(const RealVector& b_betas)
{ adpRep->betaBetas = b_betas; }


inline void AleatoryDistParams::beta_beta(const Real& b_beta, size_t i)
{ adpRep->betaBetas[i] = b_beta; }


inline const RealVector& AleatoryDistParams::beta_lower_bounds() const
{ return adpRep->betaLowerBnds; }


inline const Real& AleatoryDistParams::beta_lower_bound(size_t i) const
{ return adpRep->betaLowerBnds[i]; }


inline void AleatoryDistParams::
beta_lower_bounds(const RealVector& b_lower_bnds)
{ adpRep->betaLowerBnds = b_lower_bnds; }


inline void AleatoryDistParams::
beta_lower_bound(const Real& b_lower_bnd, size_t i)
{ adpRep->betaLowerBnds[i] = b_lower_bnd; }


inline const RealVector& AleatoryDistParams::beta_upper_bounds() const
{ return adpRep->betaUpperBnds; }


inline const Real& AleatoryDistParams::beta_upper_bound(size_t i) const
{ return adpRep->betaUpperBnds[i]; }


inline void AleatoryDistParams::
beta_upper_bounds(const RealVector& b_upper_bnds)
{ adpRep->betaUpperBnds = b_upper_bnds; }


inline void AleatoryDistParams::
beta_upper_bound(const Real& b_upper_bnd, size_t i)
{ adpRep->betaUpperBnds[i] = b_upper_bnd; }


inline const RealVector& AleatoryDistParams::gamma_alphas() const
{ return adpRep->gammaAlphas; }


inline const Real& AleatoryDistParams::gamma_alpha(size_t i) const
{ return adpRep->gammaAlphas[i]; }


inline void AleatoryDistParams::gamma_alphas(const RealVector& ga_alphas)
{ adpRep->gammaAlphas = ga_alphas; }


inline void AleatoryDistParams::gamma_alpha(const Real& ga_alpha, size_t i)
{ adpRep->gammaAlphas[i] = ga_alpha; }


inline const RealVector& AleatoryDistParams::gamma_betas() const
{ return adpRep->gammaBetas; }


inline const Real& AleatoryDistParams::gamma_beta(size_t i) const
{ return adpRep->gammaBetas[i]; }


inline void AleatoryDistParams::gamma_betas(const RealVector& ga_betas)
{ adpRep->gammaBetas = ga_betas; }


inline void AleatoryDistParams::gamma_beta(const Real& ga_beta, size_t i)
{ adpRep->gammaBetas[i] = ga_beta; }


inline const RealVector& AleatoryDistParams::gumbel_alphas() const
{ return adpRep->gumbelAlphas; }


inline const Real& AleatoryDistParams::gumbel_alpha(size_t i) const
{ return adpRep->gumbelAlphas[i]; }


inline void AleatoryDistParams::gumbel_alphas(const RealVector& gu_alphas)
{ adpRep->gumbelAlphas = gu_alphas; }


inline void AleatoryDistParams::gumbel_alpha(const Real& gu_alpha, size_t i)
{ adpRep->gumbelAlphas[i] = gu_alpha; }


inline const RealVector& AleatoryDistParams::gumbel_betas() const
{ return adpRep->gumbelBetas; }


inline const Real& AleatoryDistParams::gumbel_beta(size_t i) const
{ return adpRep->gumbelBetas[i]; }


inline void AleatoryDistParams::gumbel_betas(const RealVector& gu_betas)
{ adpRep->gumbelBetas = gu_betas; }


inline void AleatoryDistParams::gumbel_beta(const Real& gu_beta, size_t i)
{ adpRep->gumbelBetas[i] = gu_beta; }


inline const RealVector& AleatoryDistParams::frechet_alphas() const
{ return adpRep->frechetAlphas; }


inline const Real& AleatoryDistParams::frechet_alpha(size_t i) const
{ return adpRep->frechetAlphas[i]; }


inline void AleatoryDistParams::frechet_alphas(const RealVector& f_alphas)
{ adpRep->frechetAlphas = f_alphas; }


inline void AleatoryDistParams::frechet_alpha(const Real& f_alpha, size_t i)
{ adpRep->frechetAlphas[i] = f_alpha; }


inline const RealVector& AleatoryDistParams::frechet_betas() const
{ return adpRep->frechetBetas; }


inline const Real& AleatoryDistParams::frechet_beta(size_t i) const
{ return adpRep->frechetBetas[i]; }


inline void AleatoryDistParams::frechet_betas(const RealVector& f_betas)
{ adpRep->frechetBetas = f_betas; }


inline void AleatoryDistParams::frechet_beta(const Real& f_beta, size_t i)
{ adpRep->frechetBetas[i] = f_beta; }


inline const RealVector& AleatoryDistParams::weibull_alphas() const
{ return adpRep->weibullAlphas; }


inline const Real& AleatoryDistParams::weibull_alpha(size_t i) const
{ return adpRep->weibullAlphas[i]; }


inline void AleatoryDistParams::weibull_alphas(const RealVector& w_alphas)
{ adpRep->weibullAlphas = w_alphas; }


inline void AleatoryDistParams::weibull_alpha(const Real& alpha, size_t i)
{ adpRep->weibullAlphas[i] = alpha; }


inline const RealVector& AleatoryDistParams::weibull_betas() const
{ return adpRep->weibullBetas; }


inline const Real& AleatoryDistParams::weibull_beta(size_t i) const
{ return adpRep->weibullBetas[i]; }


inline void AleatoryDistParams::weibull_betas(const RealVector& w_betas)
{ adpRep->weibullBetas = w_betas; }


inline void AleatoryDistParams::weibull_beta(const Real& beta, size_t i)
{ adpRep->weibullBetas[i] = beta; }


inline const RealRealMapArray& AleatoryDistParams::histogram_bin_pairs() const
{ return adpRep->histogramBinPairs; }


inline const RealRealMap& AleatoryDistParams::
histogram_bin_pairs(size_t i) const
{ return adpRep->histogramBinPairs[i]; }


inline void AleatoryDistParams::
histogram_bin_pairs(const RealRealMapArray& h_bin_pairs)
{ adpRep->histogramBinPairs = h_bin_pairs; }


inline void AleatoryDistParams::
histogram_bin_pairs(const RealRealMap& h_bin_pr, size_t i)
{ adpRep->histogramBinPairs[i] = h_bin_pr; }


inline const RealVector& AleatoryDistParams::poisson_lambdas() const
{ return adpRep->poissonLambdas; }


inline const Real& AleatoryDistParams::poisson_lambda(size_t i) const
{ return adpRep->poissonLambdas[i]; }


inline void AleatoryDistParams::poisson_lambdas(const RealVector& p_lambdas)
{ adpRep->poissonLambdas = p_lambdas; }


inline void AleatoryDistParams::poisson_lambda(const Real& p_lambda, size_t i)
{ adpRep->poissonLambdas[i] = p_lambda; }


inline const RealVector& AleatoryDistParams::
binomial_probability_per_trial() const
{ return adpRep->binomialProbPerTrial; }


inline const Real& AleatoryDistParams::
binomial_probability_per_trial(size_t i) const
{ return adpRep->binomialProbPerTrial[i]; }


inline void AleatoryDistParams::
binomial_probability_per_trial(const RealVector& probs_per_tr)
{ adpRep->binomialProbPerTrial = probs_per_tr; }


inline void AleatoryDistParams::
binomial_probability_per_trial(const Real& prob_per_tr, size_t i)
{ adpRep->binomialProbPerTrial[i] = prob_per_tr; }


inline const IntVector& AleatoryDistParams::binomial_num_trials() const
{ return adpRep->binomialNumTrials; }


inline int AleatoryDistParams::binomial_num_trials(size_t i) const
{ return adpRep->binomialNumTrials[i]; }


inline void AleatoryDistParams::binomial_num_trials(const IntVector& num_tr)
{ adpRep->binomialNumTrials = num_tr; }


inline void AleatoryDistParams::binomial_num_trials(int num_tr, size_t i)
{ adpRep->binomialNumTrials[i] = num_tr; }


inline const RealVector& AleatoryDistParams::
negative_binomial_probability_per_trial() const
{ return adpRep->negBinomialProbPerTrial; }


inline const Real& AleatoryDistParams::
negative_binomial_probability_per_trial(size_t i) const
{ return adpRep->negBinomialProbPerTrial[i]; }


inline void AleatoryDistParams::
negative_binomial_probability_per_trial(const RealVector& probs_per_tr)
{
  if (adpRep->negBinomialProbPerTrial.empty())//Teuchos operator= doesn't resize
    adpRep->negBinomialProbPerTrial.sizeUninitialized(probs_per_tr.length());
  adpRep->negBinomialProbPerTrial = probs_per_tr; }


inline void AleatoryDistParams::
negative_binomial_probability_per_trial(const Real& prob_per_tr, size_t i)
{ adpRep->negBinomialProbPerTrial[i] = prob_per_tr; }


inline const IntVector& AleatoryDistParams::negative_binomial_num_trials() const
{ return adpRep->negBinomialNumTrials; }


inline int AleatoryDistParams::negative_binomial_num_trials(size_t i) const
{ return adpRep->negBinomialNumTrials[i]; }


inline void AleatoryDistParams::
negative_binomial_num_trials(const IntVector& num_tr)
{ adpRep->negBinomialNumTrials = num_tr; }


inline void AleatoryDistParams::
negative_binomial_num_trials(int num_tr, size_t i)
{ adpRep->negBinomialNumTrials[i] = num_tr; }


inline const RealVector& AleatoryDistParams::
geometric_probability_per_trial() const
{ return adpRep->geometricProbPerTrial; }


inline const Real& AleatoryDistParams::
geometric_probability_per_trial(size_t i) const
{ return adpRep->geometricProbPerTrial[i]; }


inline void AleatoryDistParams::
geometric_probability_per_trial(const RealVector& probs_per_tr)
{ adpRep->geometricProbPerTrial = probs_per_tr; }


inline void AleatoryDistParams::
geometric_probability_per_trial(const Real& prob_per_tr, size_t i)
{ adpRep->geometricProbPerTrial[i] = prob_per_tr; }


inline const IntVector& AleatoryDistParams::
hypergeometric_total_population() const
{ return adpRep->hyperGeomTotalPopulation; }


inline int AleatoryDistParams::hypergeometric_total_population(size_t i) const
{ return adpRep->hyperGeomTotalPopulation[i]; }


inline void AleatoryDistParams::
hypergeometric_total_population(const IntVector& total_pop)
{ adpRep->hyperGeomTotalPopulation = total_pop; }


inline void AleatoryDistParams::
hypergeometric_total_population(int total_pop, size_t i)
{ adpRep->hyperGeomTotalPopulation[i] = total_pop; }


inline const IntVector& AleatoryDistParams::
hypergeometric_selected_population() const
{ return adpRep->hyperGeomSelectedPopulation; }


inline int AleatoryDistParams::
hypergeometric_selected_population(size_t i) const
{ return adpRep->hyperGeomSelectedPopulation[i]; }


inline void AleatoryDistParams::
hypergeometric_selected_population(const IntVector& sel_pop)
{ adpRep->hyperGeomSelectedPopulation = sel_pop; }


inline void AleatoryDistParams::
hypergeometric_selected_population(int sel_pop, size_t i)
{ adpRep->hyperGeomSelectedPopulation[i] = sel_pop; }


inline const IntVector& AleatoryDistParams::hypergeometric_num_drawn() const
{ return adpRep->hyperGeomNumDrawn; }


inline int AleatoryDistParams::hypergeometric_num_drawn(size_t i) const
{ return adpRep->hyperGeomNumDrawn[i]; }


inline void AleatoryDistParams::
hypergeometric_num_drawn(const IntVector& num_drawn)
{ adpRep->hyperGeomNumDrawn = num_drawn; }


inline void AleatoryDistParams::
hypergeometric_num_drawn(int num_drawn, size_t i)
{ adpRep->hyperGeomNumDrawn[i] = num_drawn; }


inline const IntRealMapArray& AleatoryDistParams::
histogram_point_int_pairs() const
{ return adpRep->histogramPointIntPairs; }


inline const IntRealMap& AleatoryDistParams::
histogram_point_int_pairs(size_t i) const
{ return adpRep->histogramPointIntPairs[i]; }


inline void AleatoryDistParams::
histogram_point_int_pairs(const IntRealMapArray& h_pt_int_pairs)
{ adpRep->histogramPointIntPairs = h_pt_int_pairs; }


inline void AleatoryDistParams::
histogram_point_int_pairs(const IntRealMap& h_pt_int_pairs_i, size_t i)
{ adpRep->histogramPointIntPairs[i] = h_pt_int_pairs_i; }


inline const StringRealMapArray& AleatoryDistParams::
histogram_point_string_pairs() const
{ return adpRep->histogramPointStringPairs; }


inline const StringRealMap& AleatoryDistParams::
histogram_point_string_pairs(size_t i) const
{ return adpRep->histogramPointStringPairs[i]; }


inline void AleatoryDistParams::
histogram_point_string_pairs(const StringRealMapArray& h_pt_string_pairs)
{ adpRep->histogramPointStringPairs = h_pt_string_pairs; }


inline void AleatoryDistParams::
histogram_point_string_pairs(const StringRealMap& h_pt_string_pairs_i, size_t i)
{ adpRep->histogramPointStringPairs[i] = h_pt_string_pairs_i; }


inline const RealRealMapArray& AleatoryDistParams::
histogram_point_real_pairs() const
{ return adpRep->histogramPointRealPairs; }


inline const RealRealMap& AleatoryDistParams::
histogram_point_real_pairs(size_t i) const
{ return adpRep->histogramPointRealPairs[i]; }


inline void AleatoryDistParams::
histogram_point_real_pairs(const RealRealMapArray& h_pt_real_pairs)
{ adpRep->histogramPointRealPairs = h_pt_real_pairs; }


inline void AleatoryDistParams::
histogram_point_real_pairs(const RealRealMap& h_pt_real_pairs_i, size_t i)
{ adpRep->histogramPointRealPairs[i] = h_pt_real_pairs_i; }


inline const RealSymMatrix& AleatoryDistParams::uncertain_correlations() const
{ return adpRep->uncertainCorrelations; }


inline void AleatoryDistParams::
uncertain_correlations(const RealSymMatrix& uncertain_corr)
{ adpRep->uncertainCorrelations = uncertain_corr; }


inline void AleatoryDistParams::update(const AleatoryDistParams& adp)
{
  if (!adpRep) // if no rep, create a new instance
    copy(adp);
  else {      // update data of existing instance
    normal_means(adp.normal_means());
    normal_std_deviations(adp.normal_std_deviations());
    normal_lower_bounds(adp.normal_lower_bounds());
    normal_upper_bounds(adp.normal_upper_bounds());
    lognormal_means(adp.lognormal_means());
    lognormal_std_deviations(adp.lognormal_std_deviations());
    lognormal_lambdas(adp.lognormal_lambdas());
    lognormal_zetas(adp.lognormal_zetas());
    lognormal_error_factors(adp.lognormal_error_factors());
    lognormal_lower_bounds(adp.lognormal_lower_bounds());
    lognormal_upper_bounds(adp.lognormal_upper_bounds());
    uniform_lower_bounds(adp.uniform_lower_bounds());
    uniform_upper_bounds(adp.uniform_upper_bounds());
    loguniform_lower_bounds(adp.loguniform_lower_bounds());
    loguniform_upper_bounds(adp.loguniform_upper_bounds());
    triangular_modes(adp.triangular_modes());
    triangular_lower_bounds(adp.triangular_lower_bounds());
    triangular_upper_bounds(adp.triangular_upper_bounds());
    exponential_betas(adp.exponential_betas());
    beta_alphas(adp.beta_alphas());
    beta_betas(adp.beta_betas());
    beta_lower_bounds(adp.beta_lower_bounds());
    beta_upper_bounds(adp.beta_upper_bounds());
    gamma_alphas(adp.gamma_alphas());
    gamma_betas(adp.gamma_betas());
    gumbel_alphas(adp.gumbel_alphas());
    gumbel_betas(adp.gumbel_betas());
    frechet_alphas(adp.frechet_alphas());
    frechet_betas(adp.frechet_betas());
    weibull_alphas(adp.weibull_alphas());
    weibull_betas(adp.weibull_betas());
    histogram_bin_pairs(adp.histogram_bin_pairs());
    poisson_lambdas(adp.poisson_lambdas());
    binomial_probability_per_trial(adp.binomial_probability_per_trial());
    binomial_num_trials(adp.binomial_num_trials());
    negative_binomial_probability_per_trial(
      adp.negative_binomial_probability_per_trial());
    negative_binomial_num_trials(adp.negative_binomial_num_trials());
    geometric_probability_per_trial(adp.geometric_probability_per_trial());
    hypergeometric_total_population(adp.hypergeometric_total_population());
    hypergeometric_selected_population(
      adp.hypergeometric_selected_population());
    hypergeometric_num_drawn(adp.hypergeometric_num_drawn());
    histogram_point_int_pairs(adp.histogram_point_int_pairs());
    histogram_point_string_pairs(adp.histogram_point_string_pairs());
    histogram_point_real_pairs(adp.histogram_point_real_pairs());
    uncertain_correlations(adp.uncertain_correlations());
  }
}


inline void AleatoryDistParams::
update_partial(const AleatoryDistParams& adp_x,
	       const std::vector<RandomVariable>& x_ran_vars,
	       const ShortArray& u_types)
{
  if (!adpRep) { // if no rep, error 
    PCerr << "Error: AleatoryDistParams::update_partial() requires a valid "
	  << "representation." << std::endl;
    abort_handler(-1);
  }
  else { // update data of existing instance
    size_t i, num_vars = x_ran_vars.size(), nuv = 0, lnuv = 0, luuv = 0,
      tuv = 0, buv = 0, gauv = 0, guuv = 0, fuv = 0, wuv = 0, hbuv = 0;
    if (u_types.size() != num_vars) {
      PCerr << "Error: AleatoryDistParams::update_partial() requires "
	    << "transformation variable types." << std::endl;
      abort_handler(-1);
    }
    short x_type, u_type;
    for (i=0; i<num_vars; ++i) {
      x_type = x_ran_vars[i].type(); u_type = u_types[i];
      if (x_type == u_type)
	switch (u_type) {
	case BOUNDED_NORMAL:
	  normal_mean(adp_x.normal_mean(nuv), nuv);
	  normal_std_deviation(adp_x.normal_std_deviation(nuv), nuv);
	  normal_lower_bound(adp_x.normal_lower_bound(nuv), nuv);
	  normal_upper_bound(adp_x.normal_upper_bound(nuv), nuv);
	  ++nuv; break;
	case LOGNORMAL:	case BOUNDED_LOGNORMAL:
	  if (!adp_x.lognormal_means().empty()) {
	    lognormal_mean(adp_x.lognormal_mean(lnuv), lnuv);
	    if (!adp_x.lognormal_std_deviations().empty())
	      lognormal_std_deviation(adp_x.lognormal_std_deviation(lnuv),lnuv);
	    else
	      lognormal_error_factor(adp_x.lognormal_error_factor(lnuv), lnuv);
	  }
	  else if (!adp_x.lognormal_lambdas().empty()) {
	    lognormal_lambda(adp_x.lognormal_lambda(lnuv), lnuv);
	    lognormal_zeta(adp_x.lognormal_zeta(lnuv), lnuv);
	  }
	  if (u_type == BOUNDED_LOGNORMAL) {
	    lognormal_lower_bound(adp_x.lognormal_lower_bound(lnuv), lnuv);
	    lognormal_upper_bound(adp_x.lognormal_upper_bound(lnuv), lnuv);
	  }
	  ++lnuv; break;
	case LOGUNIFORM:
	  loguniform_lower_bound(adp_x.loguniform_lower_bound(luuv), luuv);
	  loguniform_upper_bound(adp_x.loguniform_upper_bound(luuv), luuv);
	  ++luuv; break;
	case TRIANGULAR:
	  triangular_mode(adp_x.triangular_mode(tuv), tuv);
	  triangular_lower_bound(adp_x.triangular_lower_bound(tuv), tuv);
	  triangular_upper_bound(adp_x.triangular_upper_bound(tuv), tuv);
	  ++tuv; break;
	case GUMBEL:
	  gumbel_alpha(adp_x.gumbel_alpha(guuv), guuv);
	  gumbel_beta(adp_x.gumbel_beta(guuv),   guuv);
	  ++guuv; break;
	case FRECHET:
	  frechet_alpha(adp_x.frechet_alpha(fuv), fuv);
	  frechet_beta(adp_x.frechet_beta(fuv),   fuv);
	  ++fuv; break;
	case WEIBULL:
	  weibull_alpha(adp_x.weibull_alpha(wuv), wuv);
	  weibull_beta(adp_x.weibull_beta(wuv),   wuv);
	  ++wuv; break;
	case HISTOGRAM_BIN:
	  histogram_bin_pairs(adp_x.histogram_bin_pairs(hbuv), hbuv);
	  ++hbuv; break;
	//default: no-op
	}
      else if (u_type == STD_BETA && x_type == BETA) {
	// lower,upper bounds are handled in conversion to STD_BETA
	beta_alpha(adp_x.beta_alpha(buv), buv);
	beta_beta(adp_x.beta_beta(buv), buv);
	++buv;
      }
      else if (u_type == STD_GAMMA && x_type == GAMMA) {
	// beta is handled in conversion to STD_GAMMA
	gamma_alpha(adp_x.gamma_alpha(gauv), gauv);
	++gauv;
      }
    }
  }
}


inline bool AleatoryDistParams::is_null() const
{ return (adpRep) ? false : true; }


/// Container class encapsulating distribution parameters for epistemic
/// random variables.

/** This class consolidates epistemic distribution data and simplifies the
    APIs that require distribution parameters. */

class EpistemicDistParams
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  EpistemicDistParams();
  /// standard constructor
  EpistemicDistParams(const RealRealPairRealMapArray& ceuv_bpas,
		      const IntIntPairRealMapArray& deuiv_bpas,
		      const IntRealMapArray& deusiv_vals_probs,
		      const StringRealMapArray& deussv_vals_probs,
		      const RealRealMapArray& deusrv_vals_probs);
  /// copy constructor
  EpistemicDistParams(const EpistemicDistParams& edp);
  /// destructor
  ~EpistemicDistParams();

  /// assignment operator
  EpistemicDistParams& operator=(const EpistemicDistParams& edp);

  //
  //- Heading: member functions
  //

  /// return total number of continuous epistemic uncertain variables
  size_t ceuv()  const;

  /// return total number of discrete epistemic uncertain integer variables
  size_t deuiv() const;
  /// return total number of discrete epistemic uncertain string variables
  size_t deusv() const;
  /// return total number of discrete epistemic uncertain real variables
  size_t deurv() const;
  /// return total number of discrete epistemic uncertain variables
  size_t deuv()   const;

  /// deep copy (as opposed to operator= shallow copy)
  void copy(const EpistemicDistParams& edp);
  /// data update (no changes to representation (unless null))
  void update(const EpistemicDistParams& edp);

  /// return the interval basic probability values
  const RealRealPairRealMapArray&
    continuous_interval_basic_probabilities() const;
  /// return the ith interval basic probability value
  const RealRealPairRealMap&
    continuous_interval_basic_probabilities(size_t i) const;
  /// set the interval basic probability values
  void continuous_interval_basic_probabilities(
    const RealRealPairRealMapArray& ci_probs);
  /// set the ith interval basic probability value
  void continuous_interval_basic_probabilities(
    const RealRealPairRealMap& ci_probs_i, size_t i);

  /*
  /// return the interval bounds
  const RealVectorArray& continuous_interval_lower_bounds() const;
  /// return the ith interval bound
  const RealVector& continuous_interval_lower_bounds(size_t i) const;
  /// set the interval bounds
  void continuous_interval_lower_bounds(const RealVectorArray& ci_l_bnds);
  /// set the ith interval bound
  void continuous_interval_lower_bounds(const RealVector& ci_l_bnds_i,
					size_t i);
  /// return the interval bounds
  const RealVectorArray& continuous_interval_upper_bounds() const;
  /// return the ith interval bound
  const RealVector& continuous_interval_upper_bounds(size_t i) const;
  /// set the interval bounds
  void continuous_interval_upper_bounds(const RealVectorArray& ci_u_bnds);
  /// set the ith interval bound
  void continuous_interval_upper_bounds(const RealVector& ci_u_bnds_i,
					size_t i);
  */

  /// return the interval basic probability values
  const IntIntPairRealMapArray& discrete_interval_basic_probabilities() const;
  /// return the ith interval basic probability value
  const IntIntPairRealMap&
    discrete_interval_basic_probabilities(size_t i) const;
  /// set the interval basic probability values
  void discrete_interval_basic_probabilities(
    const IntIntPairRealMapArray& di_bpa);
  /// set the ith interval basic probability value
  void discrete_interval_basic_probabilities(
    const IntIntPairRealMap& di_bpa_i, size_t i);

  /*
  /// return the interval bounds
  const IntVectorArray& discrete_interval_lower_bounds() const;
  /// return the ith interval bound
  const IntVector& discrete_interval_lower_bounds(size_t i) const;
  /// set the interval bounds
  void discrete_interval_lower_bounds(const IntVectorArray& di_l_bnds);
  /// set the ith interval bound
  void discrete_interval_lower_bounds(const IntVector& di_l_bnds_i, size_t i);
  /// return the interval bounds
  const IntVectorArray& discrete_interval_upper_bounds() const;
  /// return the ith interval bound
  const IntVector& discrete_interval_upper_bounds(size_t i) const;
  /// set the interval bounds
  void discrete_interval_upper_bounds(const IntVectorArray& di_u_bnds);
  /// set the ith interval bound
  void discrete_interval_upper_bounds(const IntVector& di_u_bnds_i, size_t i);
  */

  /// get the discrete integer set values and probabilities
  const IntRealMapArray& discrete_set_int_values_probabilities() const;
  /// set the discrete integer set values and probabilities
  void discrete_set_int_values_probabilities(const IntRealMapArray&
					     dsi_vals_probs);

  /// get the discrete string set values and probabilities
  const StringRealMapArray& discrete_set_string_values_probabilities() const;
  /// set the discrete string set values and probabilities
  void discrete_set_string_values_probabilities(const StringRealMapArray&
						dss_vals_probs);

  /// get the discrete real set values and probabilities
  const RealRealMapArray& discrete_set_real_values_probabilities() const;
  /// set the discrete real set values and probabilities
  void discrete_set_real_values_probabilities(const RealRealMapArray&
					      dsr_vals_probs);

  /// function to check edpRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //

  /// pointer to the body (handle-body idiom)
  EpistemicDistParamsRep* edpRep;
};


inline EpistemicDistParams::EpistemicDistParams():
  edpRep(new EpistemicDistParamsRep())
{ }


inline EpistemicDistParams::
EpistemicDistParams(const RealRealPairRealMapArray& ceuv_bpas,
		    const IntIntPairRealMapArray& deuiv_bpas,
		    const IntRealMapArray& deusiv_vals_probs,
		    const StringRealMapArray& deussv_vals_probs,
		    const RealRealMapArray& deusrv_vals_probs):
  edpRep(new EpistemicDistParamsRep(ceuv_bpas, deuiv_bpas, deusiv_vals_probs,
				    deussv_vals_probs, deusrv_vals_probs))
{ }


inline EpistemicDistParams::EpistemicDistParams(const EpistemicDistParams& edp)
{
  // Increment new (no old to decrement)
  edpRep = edp.edpRep;
  if (edpRep) // Check for an assignment of NULL
    edpRep->referenceCount++;
}


inline EpistemicDistParams::~EpistemicDistParams()
{
  if (edpRep) { // Check for NULL
    --edpRep->referenceCount; // decrement
    if (edpRep->referenceCount == 0)
      delete edpRep;
  }
}


inline EpistemicDistParams& EpistemicDistParams::
operator=(const EpistemicDistParams& edp)
{
  // Decrement old
  if (edpRep) // Check for NULL
    if ( --edpRep->referenceCount == 0 ) 
      delete edpRep;
  // Increment new
  edpRep = edp.edpRep;
  if (edpRep) // Check for an assignment of NULL
    edpRep->referenceCount++;
  return *this;
}


inline size_t EpistemicDistParams::ceuv() const
{ return edpRep->contIntervalBPA.size(); }


inline size_t EpistemicDistParams::deuiv() const
{ return edpRep->discIntervalBPA.size() + edpRep->discSetIntValsProbs.size();}


inline size_t EpistemicDistParams::deusv() const
{ return edpRep->discSetStringValsProbs.size(); }


inline size_t EpistemicDistParams::deurv() const
{ return edpRep->discSetRealValsProbs.size(); }


inline size_t EpistemicDistParams::deuv()  const
{ return deuiv() + deusv() + deurv(); }


inline void EpistemicDistParams::copy(const EpistemicDistParams& edp)
{ 
  // Decrement old
  if (edpRep) // Check for NULL
    if ( --edpRep->referenceCount == 0 ) 
      delete edpRep;
  // Create new
  edpRep = new EpistemicDistParamsRep(
    edp.continuous_interval_basic_probabilities(),
    //edp.continuous_interval_lower_bounds(),
    //edp.continuous_interval_upper_bounds(),
    edp.discrete_interval_basic_probabilities(),
    //edp.discrete_interval_lower_bounds(),
    //edp.discrete_interval_upper_bounds(),
    edp.discrete_set_int_values_probabilities(),
    edp.discrete_set_string_values_probabilities(),
    edp.discrete_set_real_values_probabilities());
}


inline const RealRealPairRealMapArray& EpistemicDistParams::
continuous_interval_basic_probabilities() const
{ return edpRep->contIntervalBPA; }


inline const RealRealPairRealMap& EpistemicDistParams::
continuous_interval_basic_probabilities(size_t i) const
{ return edpRep->contIntervalBPA[i]; }


inline void EpistemicDistParams::
continuous_interval_basic_probabilities(
  const RealRealPairRealMapArray& ci_bpa)
{ edpRep->contIntervalBPA = ci_bpa; }


inline void EpistemicDistParams::
continuous_interval_basic_probabilities(const RealRealPairRealMap& ci_bpa_i,
					size_t i)
{ edpRep->contIntervalBPA[i] = ci_bpa_i; }


/*
inline const RealVectorArray& EpistemicDistParams::
continuous_interval_lower_bounds() const
{ return edpRep->contIntervalLowerBnds; }


inline const RealVector& EpistemicDistParams::
continuous_interval_lower_bounds(size_t i) const
{ return edpRep->contIntervalLowerBnds[i]; }


inline void EpistemicDistParams::
continuous_interval_lower_bounds(const RealVectorArray& ci_l_bnds)
{ edpRep->contIntervalLowerBnds = ci_l_bnds; }


inline void EpistemicDistParams::
continuous_interval_lower_bounds(const RealVector& ci_l_bnds_i, size_t i)
{ edpRep->contIntervalLowerBnds[i] = ci_l_bnds_i; }


inline const RealVectorArray& EpistemicDistParams::
continuous_interval_upper_bounds() const
{ return edpRep->contIntervalUpperBnds; }


inline const RealVector& EpistemicDistParams::
continuous_interval_upper_bounds(size_t i) const
{ return edpRep->contIntervalUpperBnds[i]; }


inline void EpistemicDistParams::
continuous_interval_upper_bounds(const RealVectorArray& ci_u_bnds)
{ edpRep->contIntervalUpperBnds = ci_u_bnds; }


inline void EpistemicDistParams::
continuous_interval_upper_bounds(const RealVector& ci_u_bnds_i, size_t i)
{ edpRep->contIntervalUpperBnds[i] = ci_u_bnds_i; }
*/


inline const IntIntPairRealMapArray& EpistemicDistParams::
discrete_interval_basic_probabilities() const
{ return edpRep->discIntervalBPA; }


inline const IntIntPairRealMap& EpistemicDistParams::
discrete_interval_basic_probabilities(size_t i) const
{ return edpRep->discIntervalBPA[i]; }


inline void EpistemicDistParams::
discrete_interval_basic_probabilities(const IntIntPairRealMapArray& di_bpa)
{ edpRep->discIntervalBPA = di_bpa; }


inline void EpistemicDistParams::
discrete_interval_basic_probabilities(const IntIntPairRealMap& di_bpa_i,
				      size_t i)
{ edpRep->discIntervalBPA[i] = di_bpa_i; }


/*
inline const IntVectorArray& EpistemicDistParams::
discrete_interval_lower_bounds() const
{ return edpRep->discIntervalLowerBnds; }


inline const IntVector& EpistemicDistParams::
discrete_interval_lower_bounds(size_t i) const
{ return edpRep->discIntervalLowerBnds[i]; }


inline void EpistemicDistParams::
discrete_interval_lower_bounds(const IntVectorArray& di_l_bnds)
{ edpRep->discIntervalLowerBnds = di_l_bnds; }


inline void EpistemicDistParams::
discrete_interval_lower_bounds(const IntVector& di_l_bnds_i, size_t i)
{ edpRep->discIntervalLowerBnds[i] = di_l_bnds_i; }


inline const IntVectorArray& EpistemicDistParams::
discrete_interval_upper_bounds() const
{ return edpRep->discIntervalUpperBnds; }


inline const IntVector& EpistemicDistParams::
discrete_interval_upper_bounds(size_t i) const
{ return edpRep->discIntervalUpperBnds[i]; }


inline void EpistemicDistParams::
discrete_interval_upper_bounds(const IntVectorArray& di_u_bnds)
{ edpRep->discIntervalUpperBnds = di_u_bnds; }


inline void EpistemicDistParams::
discrete_interval_upper_bounds(const IntVector& di_u_bnds_i, size_t i)
{ edpRep->discIntervalUpperBnds[i] = di_u_bnds_i; }
*/


inline const IntRealMapArray& EpistemicDistParams::
discrete_set_int_values_probabilities() const
{ return edpRep->discSetIntValsProbs; }


inline void EpistemicDistParams::
discrete_set_int_values_probabilities(const IntRealMapArray& dsi_vals_probs)
{ edpRep->discSetIntValsProbs = dsi_vals_probs; }


inline const StringRealMapArray& EpistemicDistParams::
discrete_set_string_values_probabilities() const
{ return edpRep->discSetStringValsProbs; }


inline void EpistemicDistParams::
discrete_set_string_values_probabilities(
  const StringRealMapArray& dsi_vals_probs)
{ edpRep->discSetStringValsProbs = dsi_vals_probs; }


inline const RealRealMapArray& EpistemicDistParams::
discrete_set_real_values_probabilities() const
{ return edpRep->discSetRealValsProbs; }


inline void EpistemicDistParams::
discrete_set_real_values_probabilities(const RealRealMapArray& dsr_vals_probs)
{ edpRep->discSetRealValsProbs = dsr_vals_probs; }


inline void EpistemicDistParams::update(const EpistemicDistParams& edp)
{
  if (!edpRep) // if no rep, create a new instance
    copy(edp);
  else {      // update data of existing instance
    continuous_interval_basic_probabilities(
      edp.continuous_interval_basic_probabilities());
    //continuous_interval_lower_bounds(edp.continuous_interval_lower_bounds());
    //continuous_interval_upper_bounds(edp.continuous_interval_upper_bounds());
    discrete_interval_basic_probabilities(
      edp.discrete_interval_basic_probabilities());
    //discrete_interval_lower_bounds(edp.discrete_interval_lower_bounds());
    //discrete_interval_upper_bounds(edp.discrete_interval_upper_bounds());
    discrete_set_int_values_probabilities(
      edp.discrete_set_int_values_probabilities());
    discrete_set_string_values_probabilities(
      edp.discrete_set_string_values_probabilities());
    discrete_set_real_values_probabilities(
      edp.discrete_set_real_values_probabilities());
  }
}


inline bool EpistemicDistParams::is_null() const
{ return (edpRep) ? false : true; }

} // namespace Pecos

#endif
