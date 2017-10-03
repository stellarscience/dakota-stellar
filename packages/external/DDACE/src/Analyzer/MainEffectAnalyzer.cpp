// *******************************************************************
// *******************************************************************
// Operators for the class MainEffectAnalyzer  
// Author : Charles Tong
// Date   : December, 1998
// *******************************************************************
// *******************************************************************

#include "MainEffectAnalyzer.h"

#define abs(x) (((x) > 0.0) ? (x) : -(x))

// *******************************************************************
// generate group means 
// -------------------------------------------------------------------

void MainEffectAnalyzer::generateGroupMeans(int Nsamples, int Nintervals, 
                     const Array<double> x, const Array<double> y,
                     Array<double>& x2, Array<double>& gmean)
{
   int    j, k, nreplications, ncount;
   double fixedinput;

   // ---------------------------------------------------------------
   // error checking
   // ---------------------------------------------------------------

   nreplications = Nsamples / Nintervals;
   if ( nreplications*Nintervals != Nsamples )
   {
      cerr << "Error : Nsamples should be multiples of Nsymbols." << endl;
      exit( 1 );
   } 

   // ---------------------------------------------------------------
   // compute group mean
   // ---------------------------------------------------------------

   for ( j = 0; j < Nintervals; j++ )
   {
      gmean[j] = y[j];
      fixedinput = x[j];
      ncount = 1;
      for ( k = Nintervals; k < Nsamples; k++ )
      {
         if ( abs(x[k] - fixedinput ) < 1.0E-8 )
         {
            gmean[j] += y[k];
            ncount++;
         }
      }
      if ( ncount != nreplications )
      {
         cerr << "Error : samples are not replicated LHS." << endl;
         cerr << "      expected = " << nreplications << endl;
         cerr << "      actual   = " << ncount << endl;
         exit( 1 );
      }
      gmean[j] /= (double) nreplications;
   }
   for ( j = 0; j < Nintervals; j++ ) x2[j] = x[j];
   sortDouble(Nintervals, x2, gmean);
}

// ********************************************************************
// Compute agglomerated mean and variance
// --------------------------------------------------------------------
void MainEffectAnalyzer::computeMeanVariance(int Ninputs, int Noutputs, 
                        int Nsamples, const Array<double> y, 
                        Array<double>& amean, Array<double>& avariance)
{
   int    i, m;
   double mean, variance;

   amean.resize(Noutputs);
   avariance.resize(Noutputs);

   for ( m = 0; m < Noutputs; m++ )
   {
      mean = 0.0;
      for ( i = 0; i < Nsamples; i++ ) mean += y[Noutputs*i+m];
      mean    /= (double) Nsamples;
      variance = 0.0;
      for ( i = 0; i < Nsamples; i++ )
      {
          variance += ( ( y[Noutputs*i+m] - mean ) *
                        ( y[Noutputs*i+m] - mean ) );
      }
      variance /= (double) Nsamples;
      amean[m]     = mean;
      avariance[m] = variance;
   }
}

// *******************************************************************
// sortDouble : sort an double array
// -------------------------------------------------------------------
void MainEffectAnalyzer::sortDouble(int leng, Array<double>& array, 
                                    Array<double>& array2)
{
  int    i, j, jm1;
  double dtmp;

  if (leng == 1) return;
  for (i=0; i<leng; i++) {
    for (j=1; j<leng-i; j++) {
      jm1 = j - 1;
      if (array[jm1] > array[j]) {
        dtmp        = array[jm1];
        array[jm1]  = array[j];
        array[j]    = dtmp;
        dtmp        = array2[jm1];
        array2[jm1] = array2[j];
        array2[j]   = dtmp;
      }
    }
  }
}

