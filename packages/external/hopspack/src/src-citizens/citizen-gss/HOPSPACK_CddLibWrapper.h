/*
  $Id: HOPSPACK_CddLibWrapper.h 149 2009-11-12 02:40:41Z tplante $
  $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_CddLibWrapper.h $
*/

/*@HEADER
 * ************************************************************************
 * 
 *         HOPSPACK: Hybrid Optimization Parallel Search Package
 *                 Copyright 2009 Sandia Corporation
 * 
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This file is part of HOPSPACK.
 *
 * HOPSPACK is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *  
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.  If not, see http://www.gnu.org/licenses/.
 * 
 * Questions? Contact Tammy Kolda (tgkolda@sandia.gov) 
 *                 or Todd Plantenga (tplante@sandia.gov)
 * 
 * ************************************************************************
 *@HEADER
 */

/*!
  \file HOPSPACK_CddLibWrapper.h
  \brief No classes - declare C functions for invoking the CDD library.
  \author Robert Michael Lewis, College of William & Mary, 2005
*/

#ifndef HOPSPACK_CDDLIBWRAPPER_H
#define HOPSPACK_CDDLIBWRAPPER_H


#ifdef __cplusplus
extern "C" {
#endif


/*!
  \brief Function to call CDDLIB for computing cone generators.

  \param[out] num_pointy  The number of generators for the "pointy
                          part" of the tangent cone.

  \param[out] P  On exit, *P points to a (num_pointy x n) matrix
   stored in double** format whose rows form a positive spanning set
   for the "pointy part" of the tangent cone.  That is, *P[i] points
   to a double array of length n.  If we define \f$p_i\f$ to be the
   mathematical vector stored in *P[i], then \f$\{p_1,\ldots,p_{\rm
   num\_pointy}\}\f$ generates the pointy part of the space.

  \param[out] num_lineality  The number of vectors in a basis for
                             the lineality space of the tangent cone.

  \param[out] L  On exit, *L points to a (num_lineality x n) matrix
  stored in double** format whose rows form a basis for th lineality
  space of the tangent cone.  That is, *L[i] points to a double array
  of length n.  If we define \f$\ell_i\f$ to be the mathematical
  vector store in *L[i], then \f$\{\ell_1,\ldots,\ell_{\rm
  num\_lineality}\}\f$ forms a basis for the lineality space.

  \param[in] n   The array length, a.k.a. the length of array Eq[i]
                 and Iq[i], a.k.a, the number of columns in Eq and Iq.

  \param[in] num_equalities  Number of rows in Eq.  

  \param[in] Eq  A (num_equalities x n) matrix stored in double** format whose
                 rows span the lineality space of the normal cone.

  \param[in] num_inequalities  Number of rows in Iq.

  \param[in] Iq  A (num_inequalities x n) matrix stored in double** format whose
                 rows generate the "pointy part" of the normal cone. 

  \param[in] append  Used as a logical operator to signal whether or not
                     P and L are already partially formed and being appended to.
  \return 0 if successful.

  Caution:  Memory for arguments P and L is allocated internally 
  within compute_cone_generators()
  "c-style", and must therefore be deallocated "c-style" externally.

  \author Robert Michael Lewis, College of William & Mary, 2005
*/
int compute_cone_generators(int *num_pointy,    double ***P,
                            int *num_lineality, double ***L,
                            int n,
                            int num_equalities, double **Eq,
                            int num_inequalities, double **Iq,
                            int append);

#ifdef __cplusplus
}
#endif

#endif     /*-- HOPSPACK_CDDLIBWRAPPER_H */
