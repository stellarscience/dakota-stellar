// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the terms
// and conditions of the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//*************************************************************************
// * Functions for the class MOATAnalyzer  
// ************************************************************************

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "MOATAnalyzer.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// BMA: adding the following for DAKOTA compatibility
#include <cfloat>
#define PSUADE_UNDEFINED -DBL_MAX

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------

MOATAnalyzer::MOATAnalyzer()
{
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------

MOATAnalyzer::~MOATAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------

double MOATAnalyzer::analyze(aData &adata)
{
   using std::printf;
   int    nInputs, nOutputs, nSamples, outputID, ii, diffIndex;
   int    iD, *indexTrack, index, diffCnt, *counts;
   double *X, *Y, *YY, *YG, *XG, *Xbase;
   double *means, *modifiedMeans, *stds, xtemp1, xtemp2;
   double ytemp1, ytemp2, scale, *xLower, *xUpper;

   // ---------------------------------------------------------------
   // extract test information and data
   // ---------------------------------------------------------------

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   X        = adata.sampleInputs_;
   Y        = adata.sampleOutputs_;
   outputID = adata.outputID_;

   // ---------------------------------------------------------------
   // error checking
   // ---------------------------------------------------------------

   if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0 || 
       outputID < 0 || outputID >= nOutputs)
   {
      printf("MOATAnalyzer:analyze - invalid arguments.\n");
      std::exit(1);
   } 

   // ---------------------------------------------------------------
   // display header 
   // ---------------------------------------------------------------

   printf("\n*************************************************************\n");
   printf("*********************** MOAT Analysis ***********************\n");
   printf("-------------------------------------------------------------\n");

   // ---------------------------------------------------------------
   // set up for calculation
   // ---------------------------------------------------------------

   YY = new double[nSamples];
   YG = new double[nSamples];
   XG = new double[nSamples];
   for (iD = 0; iD < nSamples; iD++) YY[iD] = Y[nOutputs*iD+outputID];
   counts = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) counts[ii] = 0;
   means = new double[nInputs];
   modifiedMeans = new double[nInputs];
   stds = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++)
      means[ii] = modifiedMeans[ii] = stds[ii] = 0.0;
   indexTrack = new int[nSamples];
   for (iD = 0; iD < nSamples; iD++) indexTrack[iD] = -1;
   Xbase = new double[nSamples];
   for (iD = 0; iD < nSamples; iD++) Xbase[iD] = 0.0;

   // ---------------------------------------------------------------
   // first compute the approximate gradients (with filtering)
   // ---------------------------------------------------------------

   indexTrack[0] = -1;
   for (iD = 1; iD < nSamples; iD++)
   {
      Xbase[iD] = 0.0;
      ytemp1 = YY[iD-1]; 
      ytemp2 = YY[iD]; 
      diffCnt = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         xtemp1 = X[(iD-1)*nInputs+ii]; 
         xtemp2 = X[iD*nInputs+ii]; 
         if (xtemp1 != xtemp2 && ytemp1 !=  PSUADE_UNDEFINED &&
             ytemp2 != PSUADE_UNDEFINED)
         {
            diffCnt++;
            diffIndex = ii;
         }
      }
      // this check should be modified for self checking
      if (diffCnt == 1 && ((iD % (nInputs+1)) != 0))
      {
         indexTrack[iD] = diffIndex;
         xtemp1 = X[(iD-1)*nInputs+diffIndex]; 
         xtemp2 = X[iD*nInputs+diffIndex]; 
         scale  = xUpper[diffIndex] - xLower[diffIndex];
         YG[iD] = (ytemp2-ytemp1)/(xtemp2-xtemp1)*scale;
         if (xtemp2 > xtemp1) XG[iD] = xtemp2;
         else                 XG[iD] = xtemp1;
         counts[diffIndex]++;
         Xbase[iD] = xtemp1;
      }
      else
      {
         YG[iD] = PSUADE_UNDEFINED;
         indexTrack[iD] = -1;
      }
   }

   // ---------------------------------------------------------------
   // fix up the indexTrack array
   // ---------------------------------------------------------------

   if (nSamples / (nInputs+1) * (nInputs+1) == nSamples)
   {
      for (iD = 0; iD < nSamples; iD+=(nInputs+1))
         indexTrack[iD] = -1;
   }

   // ---------------------------------------------------------------
   // next compute the basic statistics
   // ---------------------------------------------------------------

   for (iD = 0; iD < nSamples; iD++)
   {
      if (YG[iD] != PSUADE_UNDEFINED)
      {
         index = indexTrack[iD];
         if (index >= 0)
         {
            means[index] += YG[iD];
            modifiedMeans[index] += PABS(YG[iD]);
         }
      }
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      if (counts[ii] > 0)
      {
         means[ii] /= (double) counts[ii];
         modifiedMeans[ii] /= (double) counts[ii];
      }
      else 
      {
         printf("MOATAnalyzer:analyze - zero data points for input %d\n",
                ii+1);
         means[ii] = 0.0;
         modifiedMeans[ii] = 0.0;
      }
   }
   for (iD = 0; iD < nSamples; iD++)
   {
      if (YG[iD] != PSUADE_UNDEFINED)
      {
         index = indexTrack[iD];
         if (index >= 0)
            stds[index] += (YG[iD] - means[index])*(YG[iD] - means[index]);
      }
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      if (counts[ii] > 1)
         stds[ii] /= (double) (counts[ii] - 1);
      else
      {
         printf("MOATAnalyzer:analyze - %d data points for input %d\n",
                counts[ii], ii+1);
         stds[ii] = 0.0;
      }
      if (stds[ii] < 0.0) stds[ii] = -std::sqrt(-stds[ii]);
      else                stds[ii] = std::sqrt(stds[ii]);
   }

   // ---------------------------------------------------------------
   // print the MOAT analysis data
   // ---------------------------------------------------------------

   for (ii = 0; ii < nInputs; ii++) {
     // BMA: modification to print means of elementary effects for comparison
     // printf("Input %3d (mu, mu*, std) = %12.4e %12.4e %12.4e\n",
     //        ii+1, means[ii], modifiedMeans[ii], stds[ii]);
     printf("Input %3d (mod. mean & std) = %12.4e %12.4e \n",
            ii+1, modifiedMeans[ii], stds[ii]);
   }

   // ---------------------------------------------------------------
   // clean up
   // ---------------------------------------------------------------

   delete [] counts;
   delete [] YY;
   delete [] YG;
   delete [] XG;
   delete [] means;
   delete [] modifiedMeans;
   delete [] stds;
   delete [] indexTrack;
   delete [] Xbase;
   return 0.0;
}

