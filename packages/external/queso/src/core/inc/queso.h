//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// grvy.h: Basic API Definitions
//
//--------------------------------------------------------------------------

#ifndef QUESO_H_
#define QUESO_H_

// Library version/build info

#define QUESO_MAJOR_VERSION  0
#define QUESO_MINOR_VERSION  54
#define QUESO_MICRO_VERSION  0

#define QUESO_BUILD_USER     "briadam"
#define QUESO_BUILD_ARCH     "x86_64-unknown-linux-gnu"
#define QUESO_BUILD_HOST     "rem.sandia.gov"
#define QUESO_BUILD_DATE     "2015-10-08 14:14"
#define QUESO_BUILD_VERSION  "547efb1"

#define QUESO_LIB_VERSION    "0.54.0"
#define QUESO_LIB_RELEASE    "Development Build"

#define QUESO_CXX            "mpic++"
#define QUESO_CXXFLAGS       "-g -O2 -Wall"

// External libraries

#define QUESO_TRILINOS_DIR  ""
#define QUESO_GSL_DIR       "-L/apps/gsl/1.15/lib -lgsl -lgslcblas -lm"
#define QUESO_GRVY_DIR      ""
#define QUESO_GLPK_DIR      ""
#define QUESO_HDF5_DIR      ""

#endif
