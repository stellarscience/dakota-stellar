/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include <stdio.h>
#include "dot.hpp"

double dot(int len, double* vec1, double* vec2)
{
    int i;
    double d;

    d = 0;
    for(i=0;i<len;i++)
        d += vec1[i]*vec2[i];

    return d;
}
