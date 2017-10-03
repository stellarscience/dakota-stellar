// ********************************************************************
// This is the main driver for DDACE ANOVA
// Author : Charles Tong
// Data   : August, 1998
// ********************************************************************

// --------------------------------------------------------------------
// internal and external library declarations
// --------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include "DDace.h"
#include "Mars.h"
#include "MainEffectAnalyzer.h"
#include "DDaceArchiveReader.h"
#include "Array.inl"
#include <iostream>
#include <fstream>
using namespace std;

int  ComputeNsymbols(int Ninputs, int Nsamples, const Array<double>& x);
void DDaceAnalysis(char*,int,int*,int*,int*,double**,double**,char***,double**,
             double**, double**, double**, double**);

#define abs(x) (((x) > 0.0) ? (x) : -(x))

void DDaceAnalysis(char *DDace_fname, int DDace_option, int *DDace_axes, 
             int *DDace_nsamples, int *DDace_nrep, double **DDace_lbound,
             double **DDace_ubound, char ***DDace_labels, double **DDace_xdata,
             double **DDace_ydata, double **DDace_zdata, double **DDace_xmean,
             double **DDace_ymean)
{
   int                     i, nSamples, nInputs, nOutputs, nReps=1;
   int                     xindex, yindex, zindex, pts_per_dim=20, length;
   int                     nSymbols;
   double                  xmin, xmax, ymin, ymax;
   const char              *cstr1, *cstr2, *cstr3;
   DDaceSampler            sampler;
   MainEffectAnalyzer      analyzer;
   Array<double>           sampleOne, resultOne, fixedInputs, ygrid;
   Array<double>           xlower, xupper, xmean, ymean;
   Array<Array<double> >   xgrid, results;
   Array<String>           varNames, outputNames;
   Array<DDaceSamplePoint> samplePts;
   Array<DDaceRunStatus>   runStatus;

   //------------------------------------------------------------------
   // read data from file
   //------------------------------------------------------------------

   try {

			 //		 ofstream ls("ddace.log");
			 //ls << "reading filename = " << DDace_fname << endl;
      cerr << "reading filename = [" << DDace_fname << "]" << endl;
			cerr << "string length = " << strlen(DDace_fname) << endl;
			cerr << "address = " << (int*) DDace_fname << endl;
      DDaceReader reader = new DDaceArchiveReader(DDace_fname);
      //ls << "getting result " << endl;
      reader.getArchivedData(samplePts, results, runStatus);
      //ls << samplePts << endl;
      //ls << results << endl;

      nSamples = samplePts.length();
      //ls << "number of samples: " << nSamples << endl;

      // get the input names and number of inputs

      reader.getVarNames(varNames);
      //ls << "variable names: " << varNames << endl;
      cerr << "variable names: " << varNames << endl;
      nInputs = varNames.length();
      //ls << "number of inputs: " << nInputs << endl;

      // get the output names and number of outputs

      reader.getOutputNames(outputNames);
      //ls << "output names " << outputNames << endl;
      nOutputs = outputNames.length();
      //ls << "number of outputs: " << nOutputs << endl;

      // get the sampler

      reader.getSampler(sampler);
			xupper = sampler.upperBounds();
			xlower = sampler.lowerBounds();
      //ls << "sample method: " << sampler.typeName() << endl;

      // if the sampler is Latin hypercube, get the number of replications

      if (sampler.typeName() == "DDaceLHSampler")
      {
         nReps = sampler.getParameter("Replications");
	 //cerr << "number of replications: " << nReps << endl;
      }
   }
   catch(ExceptionBase& e)
   {
      e.print();
      exit(1);
   }
   //------------------------------------------------------------------
   // plot or statistics options
   //------------------------------------------------------------------

   Mars mars(nInputs, nSamples);
   switch ( DDace_option )
   {
      case 1 : xindex = DDace_axes[0];
               yindex = DDace_axes[1];
               zindex = DDace_axes[2];
               resultOne.resize(nSamples);
               for (i=0;i<nSamples;i++) 
                  resultOne[i] = results[i][zindex];
               //resultOne = sliceArray(results,zindex);
               fixedInputs.resize(nInputs);
               for (i=0;i<nInputs;i++) 
                  fixedInputs[i] = 0.5*(xlower[i]+xupper[i]);
               length = pts_per_dim * pts_per_dim;
               xgrid.resize(length);
               for (i=0;i<length;i++) xgrid[i].resize(nInputs);
               ygrid.resize(length);
               mars.setNPtsPerDim(pts_per_dim);
               mars.setBounds(xlower, xupper);
               mars.generate2DGridData(samplePts, resultOne, xindex, yindex, 
                                       fixedInputs, xgrid, ygrid);

               (*DDace_nsamples) = length;
               (*DDace_nrep)  = 1;
               (*DDace_lbound) = new double[2];
               (*DDace_lbound)[0] = xlower[xindex];
               (*DDace_lbound)[1] = xlower[yindex];
               (*DDace_ubound) = new double[2];
               (*DDace_ubound)[0] = xupper[xindex];
               (*DDace_ubound)[1] = xupper[yindex];
               (*DDace_labels) = new char*[3];
               for ( i = 0; i < 3; i++ ) (*DDace_labels)[i] = new char[60];
               cstr1 = varNames[xindex].cString();
               strcpy((*DDace_labels)[0], cstr1); 
               cstr2 = varNames[yindex].cString();
               strcpy((*DDace_labels)[1], cstr2); 
               cstr3 = outputNames[zindex].cString();
               strcpy((*DDace_labels)[2], cstr3); 
               (*DDace_xdata) = new double[length];
               (*DDace_ydata) = new double[length];
               (*DDace_zdata) = new double[length];
               for ( i = 0; i < length; i++ )
               {
                  (*DDace_xdata)[i] = xgrid[i][xindex];
                  (*DDace_ydata)[i] = xgrid[i][yindex];
                  (*DDace_zdata)[i] = ygrid[i];
               }
               break;

      case 2 : xindex = DDace_axes[0];
               yindex = DDace_axes[1];
               sampleOne.resize(nSamples);
               resultOne.resize(nSamples);
               xmax  = -1.0e50;
               ymax  = -1.0e50;
               xmin  =  1.0e50;
               ymin  =  1.0e50;
               for (i=0;i<nSamples;i++) 
                  resultOne[i] = results[i][yindex];
               //resultOne = sliceArray( results, yindex );
               for ( i = 0; i < nSamples; i++ ) 
               {
                  sampleOne[i] = samplePts[i][xindex];
                  if ( sampleOne[i] < xmin ) xmin = sampleOne[i];
                  if ( sampleOne[i] > xmax ) xmax = sampleOne[i];
                  if ( resultOne[i] < ymin ) ymin = resultOne[i];
                  if ( resultOne[i] > ymax ) ymax = resultOne[i];
               }

               nSymbols = nSamples / nReps;
               (*DDace_xmean) = new double[nSymbols];
               (*DDace_ymean) = new double[nSymbols];
               xmean.resize(nSymbols);
               ymean.resize(nSymbols);
               analyzer.generateGroupMeans(nSamples, nSymbols, sampleOne,
                                           resultOne, xmean, ymean);
               (*DDace_xdata) = new double[nSamples];
               (*DDace_ydata) = new double[nSamples];
               for ( i = 0; i < nSamples; i++ )
               {
                  (*DDace_xdata)[i] = sampleOne[i];
                  (*DDace_ydata)[i] = resultOne[i];
               }
               for ( i = 0; i < nSymbols; i++ ) {
                  (*DDace_xmean)[i] = xmean[i];
                  (*DDace_ymean)[i] = ymean[i];
               }
               (*DDace_labels) = new char*[2];
               for ( i = 0; i < 2; i++ ) (*DDace_labels)[i] = new char[30];

               cstr1 = varNames[xindex].cString();
               strcpy((*DDace_labels)[0], cstr1); 
               cstr2 = varNames[yindex].cString();
               strcpy((*DDace_labels)[1], cstr2); 

               (*DDace_lbound) = new double[2];
               (*DDace_lbound)[0] = xmin;
               (*DDace_lbound)[1] = ymin;
               (*DDace_ubound) = new double[2];
               (*DDace_ubound)[0] = xmax;
               (*DDace_ubound)[1] = ymax;
               (*DDace_nsamples)  = nSamples;
               (*DDace_nrep)      = nSymbols;
               break;

      default : 
	ExceptionBase::raise("invalid option in DDaceAnalysis");
   }
   return;
}


// ********************************************************************
// Compute Nsymbols
// --------------------------------------------------------------------

int ComputeNsymbols(int Ninputs, int Nsamples, const Array<double>& x)
{
    int    i, nreplications;
    double fixedpt;

    nreplications = 0;
    fixedpt = x[0];
    for ( i = 0; i < Nsamples; i++ )
    {
       if ( (x[Ninputs*i] - fixedpt) < 1.0E-10 &&      
            (fixedpt - x[Ninputs*i]) < 1.0E-10 )      
          nreplications++;
    }
    return Nsamples/nreplications;
}

