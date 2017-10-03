/*
  $Id: setoper.h 149 2009-11-12 02:40:41Z tplante $
  $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/cddlib/setoper.h $
*/

/*@HEADER
 * ************************************************************************
 * 
 *         HOPSPACK: Hybrid Optimization Parallel Search Package
 *               Copyright (2003-2009) Sandia Corporation
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
  \file setoper.h 
  \brief Part of <A href="http://www.ifor.math.ethz.ch/~fukuda/cdd_home/cdd.html">cdd</a> library.

  This file is from cddlib, version 0.94, and was downloaded from 
  http://www.ifor.math.ethz.ch/~fukuda/cdd_home/cdd.html
  on Sept. 14, 2005. 
  \author <A href="http://www.ifor.math.ethz.ch/~fukuda">Komei Fukuda</a>, August 4, 2005. 
  \sa \ref sectionCddlib.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/* Header file for setoper.c  */

/* setoper.c: 
 * A set operation library 
 * created by Komei Fukuda, Nov.14, 1993
 * last modified on June 1, 2000
 */

#ifndef  __SETOPER_H
#define  __SETOPER_H
#endif  /* __SETOPER_H */

#include <stdio.h>
#include <stdlib.h>

typedef unsigned long *set_type;   /* set type definition */

typedef unsigned char set_card_lut_t;

unsigned long set_blocks(long len);
void set_initialize(set_type *setp,long len);
void set_free(set_type set);
void set_emptyset(set_type set);
void set_copy(set_type setcopy,set_type set);
void set_addelem(set_type set, long elem);
void set_delelem(set_type set, long elem);
void set_int(set_type set,set_type set1,set_type set2);
void set_uni(set_type set,set_type set1,set_type set2);
void set_diff(set_type set,set_type set1,set_type set2);
void set_compl(set_type set,set_type set1);
int set_subset(set_type set1,set_type set2);
int set_member(long elem, set_type set);
long set_card(set_type set);
long set_groundsize(set_type set); /* output the size of the ground set */
void set_write(set_type set);
void set_fwrite(FILE *f,set_type set);
void set_fwrite_compl(FILE *f,set_type set); /* write the complement */
void set_binwrite(set_type set);
void set_fbinwrite(FILE *f,set_type set);



/* End of File: setoper.h */

#endif
