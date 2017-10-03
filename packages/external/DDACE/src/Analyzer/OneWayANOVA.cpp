/// Creates a '1-Way Analysis of Variance' table

///
/// Programmed: EJW, 2005
/// Modifiers:  MMC, 1/30/2006    Added comments 
///                               Call to gslib to print out p-value in ANOVA table  
///                               Renamed to OneWayANOVA

#include "OneWayANOVA.h"
#include <iomanip>
#include <iostream>

#ifdef HAVE_BOOST
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/expm1.hpp>
#endif

using std::setw;

namespace DDaceMainEffects {


OneWayANOVA::OneWayANOVA(std::vector<Factor> factors)
	: factors_(factors)
{

	/// Check to make sure that there are factors (i.e. that the vector
	/// is not empty. And if it is, throw a runtime_error
	if(factors_.empty()) 
	  throw std::runtime_error("Error in MainEffects ctor: factors is empty");

	/// Check to make sure that the size of each Factor vector is the
	/// same. If they are not the same, then throw a runtime_error
	int factorSize = factors_[0].getNumberOfObservations();
	
	for(int i = 1; i < (int) factors_.size(); i++)
	{
	  if(factorSize != factors_[i].getNumberOfObservations())
	    throw std::runtime_error("Error in MainEffects ctor: factors are different sizes");				
	}

}

std::string OneWayANOVA::getANOVATable(int factor)
{  	
        std::ostringstream ss;

	/// The 1-Way ANOVA table contains the analysis of variance components
	/// associated with apportioning the variance in a response variable
	/// to just one factor ("Between Groups") and to the remaining error ("Within Groups"):
	///
	/// Source of                  Sum of  
	/// Variation         DoF     Squares     Fdata    P
	/// -----------------------------------------------------------------
	/// Between Groups    N1        XXX       XXX      XXX
	/// Within Groups     N2        XXX             
	/// Total             N1+N2     XXX
	///
	///
	/// The "Fdata" term corresponds to the ratio of these variances and implies importance.
	/// The p-value gives significance of importance.  
	/// If the p-value is smaller than your acceptable Type I error (or risk alpha), then this implies
	/// a (1-p) Confidence that the "Between Groups" effect on the response variable variance 
	/// due to the Factor is significant.  


        Factor currFac = getFactor(factor);
	double SSBG = currFac.sumOfSquaresBetweenGroups();
	double SSWG = currFac.sumOfSquaresWithinGroups();
	int doFBetween = currFac.doFBetween();
	int doFWithin = currFac.doFWithin();
	double Fdata = currFac.Fdata();
#ifdef HAVE_BOOST
        namespace bmth = boost::math;
        namespace bmp  = bmth::policies;

        typedef bmth::
          fisher_f_distribution< double,
             bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
          fisher_f_dist;
 
         fisher_f_dist  fdist(doFBetween,doFWithin);
         double p_value =  1.0 - bmth::cdf(fdist, Fdata);
#endif

//NOTE:  AS OF JULY 5, 2012, we are removing the dependency of 
//DDACE on GSL.  Also, the code below is incorrect:  the p-value 
//should be calculated based on CDF, not PDF.  So previous versions 
//are returning an incorrect number.
//#ifdef HAVE_GSL
//	double p_value = gsl_ran_fdist_pdf(Fdata,doFBetween,doFWithin);
//#endif

        ss.setf(std::ios::scientific);
        ss<<std::setprecision(5);
	ss <<"Source of        "  << "      " << "        Sum of" << "      Mean Sum" << std::endl;
        ss <<"Variation        " << "   DoF"  << "       Squares" << "    of Squares"  << "         Fdata" << "       p-value"  << std::endl;
#ifdef HAVE_BOOST
        ss <<"Between Groups   " << setw(6) << doFBetween << setw(14) << SSBG << setw(14) << currFac.varianceBetweenGroups() << setw(14) << Fdata       << setw(14)    << p_value     << std::endl;
#else
        ss <<"Between Groups   " << setw(6) << doFBetween << setw(14) << SSBG << setw(14) << currFac.varianceBetweenGroups() << setw(14) << Fdata       << setw(14)    << " ************ "     << std::endl;
#endif
        ss  <<"Within Groups    " << setw(6)<< doFWithin  << setw(14) << SSWG << setw(14) <<currFac.varianceWithinGroups() << std::endl;
	ss  << "Total            " << setw(6) << (doFBetween + doFWithin)  << setw(14) << ( SSBG + SSWG ) << std::endl;  
        return(ss.str());

}


std::string OneWayANOVA::getANOVATables()
{
        std::ostringstream ss;
        for(int i = 0; i < (int) factors_.size(); i++)
        {
                
              ss << "\n"  << "ANOVA Table for Factor (Variable) " << i+1 
                 << "\n"  << getANOVATable(i);
        }

	//         ss << getANOVATable(0);
         return(ss.str());                         
        
}

void OneWayANOVA::printANOVATable(int factor)
{
        std::cout << getANOVATable(factor);                               
}
 
void OneWayANOVA::printANOVATables()
{
        std::cout << getANOVATables();
}

}//namespace
