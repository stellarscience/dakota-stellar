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
// * Definition for the Morris one-at-a-time sampling class
// ************************************************************************

#ifndef __MOATSAMPLINGH_
#define __MOATSAMPLINGH_

#include "sData.h"
#include "DakotaPsuade.H"

class MOATSampling: public DakotaPsuade
{
   int    P_;
   int    nSamples_;
   int    nInputs_;

public:

   MOATSampling();
   MOATSampling(int seed, int partitions);
   ~MOATSampling();

   int initialize(sData &);

private:

   int generate(double **);
};

#endif // __MOATSAMPLINGH_

