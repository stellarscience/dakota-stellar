/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "boost_test_dist.hpp"
#include "pecos_data_types.hpp"
#include "pecos_global_defs.hpp"
#include "pecos_stat_util.hpp"
//#include <algorithm>
#include <iomanip>


void boost_test_dist::print_comparison()
{
  double cdf_value,quantile_value,this_gamma;
  int i,j;
  PCout.precision(16); 
  PCout.setf(std::ios::scientific);
  PCout.setf(std::ios::showpoint);

  // WJB: note GSL is no longer supported as a TPL

  PCout << "Quantiles for normal" << '\n';
  for (i=0; i<11; i++){
    cdf_value = 0.1*i;
//#ifdef HAVE_BOOST
    Pecos::normal_dist norm(0,1);
    // Pecos::Real q = bmth::quantile(norm, cdf_value);
    PCout << std::setw(16) << "Boost " << bmth::quantile(norm,cdf_value) << '\n';
//#endif
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_ugaussian_Pinv(cdf_value) << '\n';
//#endif
  }

  PCout << std::setw(16) << "CDF values for normal" << '\n';
  for (i=-20; i<20; i++){
    quantile_value = 1.0*i;
//#ifdef HAVE_BOOST
    Pecos::normal_dist norm(0,1);
    PCout << std::setw(16) << "Boost " << bmth::cdf(norm,quantile_value) << '\n';
//#endif
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_ugaussian_P(quantile_value) << '\n';
//#endif // HAVE_GSL
  }
  
  PCout << std::setw(16) << "PDF values for normal" << '\n';
  for (i=-20; i<20; i++){
    quantile_value = 1.0*i;
//#ifdef HAVE_BOOST
    Pecos::normal_dist norm(0,1);
    PCout << std::setw(16) << "Boost " << bmth::pdf(norm,quantile_value) << '\n';
//#endif
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_ran_ugaussian_pdf(quantile_value) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Gamma function values" << '\n';
  for (i=2; i<6; i++){
    quantile_value = 1.0*i+0.5;
//#ifdef HAVE_BOOST
    this_gamma = bmth::tgamma(quantile_value);
    PCout << std::setw(16) << "Boost " << this_gamma << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //this_gamma = gsl_sf_gamma(quantile_value);
    //PCout << std::setw(16) << "  GSL " << this_gamma << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Beta function values" << '\n';
  for (i=1; i<5; i++){
    Pecos::Real alpha1 = 1.0*i;
//#ifdef HAVE_BOOST
    this_gamma = bmth::beta(alpha1, 1.0);
    PCout << std::setw(16) << "Boost " << this_gamma << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //this_gamma = gsl_sf_beta(alpha1, 1.0);
    //PCout << std::setw(16) << "  GSL " << this_gamma << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Quantiles for gamma distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
//#ifdef HAVE_BOOST
    Pecos::gamma_dist gamma1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::quantile(gamma1,cdf_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_gamma_Pinv(cdf_value, 3.,1.) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "CDF values for gamma distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 1.0*i; 
//#ifdef HAVE_BOOST
    Pecos::gamma_dist gamma1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::cdf(gamma1,quantile_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_gamma_P(quantile_value, 3.,1.) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "PDF values for gamma distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 1.0*i; 
//#ifdef HAVE_BOOST
    Pecos::gamma_dist gamma1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(gamma1,quantile_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_ran_gamma_pdf(quantile_value, 3.,1.) << '\n';
//#endif // HAVE_GSL
  }
  
  PCout << std::setw(16) << "PDF values for weibull distribution" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i; 
//#ifdef HAVE_BOOST
    Pecos::weibull_dist weibull1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(weibull1,quantile_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_ran_weibull_pdf(quantile_value, 1.,3.) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "CDF values for weibull distribution" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i; 
//#ifdef HAVE_BOOST
    Pecos::weibull_dist weibull1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::cdf(weibull1,quantile_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_weibull_P(quantile_value, 1.,3.) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "CDF values for beta distribution" << '\n';
  for (i=0; i<11; i++){
    quantile_value = i*0.1; 
//#ifdef HAVE_BOOST
    Pecos::beta_dist beta1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::cdf(beta1,quantile_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_beta_P(quantile_value, 3.,1.) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Quantile values for beta distribution" << '\n';
  for (i=0; i<11; i++){
    cdf_value = 0.1*i; 
//#ifdef HAVE_BOOST
    Pecos::beta_dist beta1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::quantile(beta1,cdf_value) <<'\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //Pecos::ProbabilityTransformation this_trans;
    //Pecos::Real q = this_trans.cdf_beta_Pinv(cdf_value,3.0,1.0); 
    //PCout << std::setw(16) << "  GSL " << q  << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "PDF values for beta distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 0.1*i; 
//#ifdef HAVE_BOOST
    Pecos::beta_dist beta1(3.,1.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(beta1,quantile_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_ran_beta_pdf(quantile_value, 3., 1.) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "PDF values for F distribution" << '\n';
  for (i=1; i<6; i++){
    quantile_value = 0.1*i; 
//#ifdef HAVE_BOOST
    Pecos::fisher_f_dist fisher(10.,4.);
    PCout << std::setw(16) << "Boost " << bmth::pdf(fisher,quantile_value) << '\n';
//#endif //HAVE_BOOST
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_ran_fdist_pdf(quantile_value, 10., 4.) << '\n';
//#endif // HAVE_GSL
  }

  PCout << std::setw(16) << "Quantiles for t-distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
//#ifdef HAVE_BOOST
    Pecos::students_t_dist tdist(10);
    PCout << std::setw(16) << "Boost " << bmth::quantile(tdist,cdf_value) << '\n';
//#endif
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_tdist_Pinv(cdf_value,10) << '\n';
//#endif
  }

  PCout << std::setw(16) << "Quantiles for Chi-squared distribution" << '\n';
  for (i=1; i<10; i++){
    cdf_value = 0.1*i;
//#ifdef HAVE_BOOST
    Pecos::chi_squared_dist chisq(10);
    PCout << std::setw(16) << "Boost " << bmth::quantile(chisq,cdf_value) << '\n';
//#endif
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_chisq_Pinv(cdf_value,10) << '\n';
//#endif
  }

 PCout << std::setw(16) << "CDF values for exponential" << '\n';
  for (i=1; i<5; i++){
    quantile_value = 0.1*i;
//#ifdef HAVE_BOOST
    Pecos::exponential_dist expon1(1.5);
    PCout << std::setw(16) << "Boost " << bmth::cdf(expon1,quantile_value) << '\n';
//#endif
//#ifdef HAVE_GSL
    //PCout << std::setw(16) << "  GSL " << gsl_cdf_exponential_P(quantile_value,(1/1.5)) << '\n';
//#endif // HAVE_GSL
  }
}


int main(int argc, char* argv[])
{
  boost_test_dist bt;
  bt.print_comparison();

  return 0;
}
