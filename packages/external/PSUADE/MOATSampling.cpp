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
// ************************************************************************
// Functions for the class MOATAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <cstdlib>
#include <cmath>
#include "MOATSampling.h"

// DAKOTA/BMA: for max/min
#include <algorithm>
// DAKOTA/BMA: not sure why this wasn't included
#include <cstdio>
// DAKOTA/BMA: includes random number generator and glue to DAKOTA
#include "DakotaPsuade.H"

//*************************************************************************
//* Constructor
//*------------------------------------------------------------------------

MOATSampling::MOATSampling(): DakotaPsuade()
{
   P_            = 4;
   nSamples_     = 0;
   nInputs_      = 0;
   sampleMatrix_ = NULL;
}

// BMA: alternate constructor with seed for use with DAKOTA
MOATSampling::MOATSampling(int seed, int partitions): DakotaPsuade(seed)
{
   P_            = partitions;
   nSamples_     = 0;
   nInputs_      = 0;
   sampleMatrix_ = NULL;
}

//********************************************************************
//* destructor 
//*-------------------------------------------------------------------

MOATSampling::~MOATSampling()
{
   if (sampleMatrix_ != NULL)
   {
      for (int ii = 0; ii < nSamples_; ii++)
         if (sampleMatrix_[ii] != NULL)
            delete [] sampleMatrix_[ii];
      delete [] sampleMatrix_;
      sampleMatrix_ = NULL;
   }
}

//*************************************************************************
//* initialize the sampling data
//*------------------------------------------------------------------------

int MOATSampling::initialize(sData &sdata)
{
   int    ii, ii2, rr, ss, nReps, maxReps, maxSamples, index, base1, base2, 
          kk1, kk2;
   double **BS, *ranges, ddata, delta, dDist, maxDist;

   // ----------------------------------------------------------------
   // clean up 
   // ----------------------------------------------------------------

   nSamples_ = sdata.nSamples_;
   nInputs_  = sdata.nInputs_;
   if (nSamples_/(nInputs_+1) * (nInputs_+1) != nSamples_) 
   {
      std::printf("MOATSampling: nSamples should be multiples of nInputs+1.\n");
      std::printf("              nSamples reset to be 10*(nInputs+1).\n");
      nSamples_ = 10 * (nInputs_ + 1);
   }

   // ----------------------------------------------------------------
   // preparation for generating the samples
   // ----------------------------------------------------------------

   sampleMatrix_ = new double*[nSamples_];
   for (ii = 0; ii < nSamples_; ii++)
      sampleMatrix_[ii] = new double[nInputs_];

   ranges = new double[nInputs_];
   for (ii = 0;  ii < nInputs_;  ii++) 
      ranges[ii] = sdata.iUpperB_[ii] - sdata.iLowerB_[ii];

   // ----------------------------------------------------------------
   // allocate pattern matrix
   // generate maxReps MOAT samples, then downselect based on distance
   // ----------------------------------------------------------------

   // maxReps = 1000;
   nReps = nSamples_ / (nInputs_ + 1);
   // DAKOTA/BMA: default to 250*nReps, no more than 1000, with a min of nReps
   maxReps = std::max(nReps, std::min(250*nReps, 1000));
   maxSamples = (nInputs_ + 1) * maxReps;
   BS = new double*[maxSamples];
   for (ii = 0; ii < maxSamples; ii++) BS[ii] = new double[nInputs_];

   // ----------------------------------------------------------------
   // generate a MOAT sample
   // ----------------------------------------------------------------

   for (rr = 0; rr < maxReps; rr++) generate(&BS[rr*(nInputs_+1)]);

   // ----------------------------------------------------------------
   // path selection
   // ----------------------------------------------------------------

   for (rr = 1; rr < nReps; rr++) 
   {
      maxDist = 0;
      index = rr;
      base1 = (rr - 1) * (nInputs_ + 1);
      for (ss = rr; ss < maxReps; ss++)
      {
         dDist = 0;
         base2 = ss * (nInputs_ + 1);
         for (kk1 = 0; kk1 <= nInputs_; kk1++) 
         {
            for (kk2 = 0; kk2 <= nInputs_; kk2++) 
            {
               for (ii = 0; ii < nInputs_; ii++) 
               {
                  ddata = BS[base1+kk1][ii] - BS[base2+kk2][ii];
                  dDist += ddata * ddata;
               }
            }
         }
         if (dDist > maxDist)
         {
            maxDist = dDist;
            index = ss;
         }
      }
      if (index != rr)
      { 
         base1 = rr * (nInputs_ + 1);
         base2 = index * (nInputs_ + 1);
         for (kk1 = 0; kk1 <= nInputs_; kk1++) 
         {
            for (ii = 0; ii < nInputs_; ii++)
            {
               ddata = BS[base1+kk1][ii];
               BS[base1+kk1][ii] = BS[base2+kk1][ii];
               BS[base2+kk1][ii] = ddata;
            }
         }
      }
   }

   // ----------------------------------------------------------------
   // convert the sample to application ranges
   // ----------------------------------------------------------------

   for (ss = 0; ss < nSamples_; ss+=(nInputs_+1))
   {
      for (ii = 0; ii <= nInputs_; ii++)
      {
         for (ii2 = 0; ii2 < nInputs_; ii2++)
         {
            ddata = BS[ss+ii][ii2];
            ddata = ddata * ranges[ii2] + sdata.iLowerB_[ii2];
            sampleMatrix_[ss+ii][ii2] = ddata;
         }
      }
   }

   // ----------------------------------------------------------------
   // clean up
   // ----------------------------------------------------------------

   delete [] ranges;
   for (ii = 0;  ii < maxSamples; ii++) delete [] BS[ii];
   delete [] BS;
   return 0;
}

//*************************************************************************
//* generate the BS matrix
//*------------------------------------------------------------------------

int MOATSampling::generate(double **BS)
{
   int    *permute, ss, ii, ii2, idata, imax;
   double **B, *D, *X, **B2, delta, ddata;

   // ----------------------------------------------------------------
   // initialize
   // ----------------------------------------------------------------

   PSUADE_randInit(); // initialize random number generator
   delta = P_ / ((double) (2*P_) - 2.0);

   // ----------------------------------------------------------------
   // build the B matrix, and allocate for the other
   // ----------------------------------------------------------------

   B = new double*[nInputs_+1];
   for (ii = 0; ii <= nInputs_; ii++)
   {
      B[ii] = new double[nInputs_];
      for (ii2 = 0; ii2 < ii; ii2++) B[ii][ii2] = 1.0;
      for (ii2 = ii; ii2 < nInputs_; ii2++) B[ii][ii2] = 0.0;
   }
   D = new double[nInputs_];
   X = new double[nInputs_];
   permute = new int[nInputs_];
   B2   = new double*[nInputs_+1];
   for (ii = 0; ii <= nInputs_; ii++) B2[ii] = new double[nInputs_];

   // ----------------------------------------------------------------
   // initialize the D matrix
   // ----------------------------------------------------------------

   for (ii = 0; ii < nInputs_; ii++)
   {
      D[ii] = PSUADE_drand(); // generate a random double in [0,1]
      if (D[ii] > 0.5) D[ii] = 1.0;
      else             D[ii] = -1.0;
   }

   // -----------------------------------------------------------------
   // initialize the X vector
   // ----------------------------------------------------------------

   imax = (P_ - 1) / 2;
   for (ii = 0; ii < nInputs_; ii++)
   {
       X[ii] = PSUADE_drand(); // generate a random double in [0,1]
       idata = (int) (X[ii] * (imax + 1));
       if (idata > imax) idata--;
       X[ii] = (double) idata / (double) (P_ - 1);
   }
   
   // ----------------------------------------------------------------
   // initialize the permutation matrix
   // ----------------------------------------------------------------

   // put [0:nInputs-1] in permute in random order
   generateRandomIvector(nInputs_, permute);
   
   // ----------------------------------------------------------------
   // form the BS matrix
   // ----------------------------------------------------------------

   for (ii = 0; ii <= nInputs_; ii++)
      for (ii2 = 0; ii2 < nInputs_; ii2++)
         B2[ii][ii2] = X[ii2]+delta/2*((B[ii][ii2]*2-1.0)*D[ii2]+1.0);
   for (ii = 0; ii <= nInputs_; ii++)
      for (ii2 = 0; ii2 < nInputs_; ii2++)
         BS[ii][ii2] = B2[ii][permute[ii2]];

   // ----------------------------------------------------------------
   // clean up
   // ----------------------------------------------------------------

   for (ii = 0;  ii <= nInputs_; ii++)
   {
      delete [] B[ii];
      delete [] B2[ii];
   }
   delete [] B;
   delete [] B2;
   delete [] D;
   delete [] X;
   delete [] permute;
   return 0;
}

