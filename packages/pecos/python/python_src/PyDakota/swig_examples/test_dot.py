#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

from PyDakota.swig_examples import dot
import numpy
vec1=[1,2,3]
vec2=[4,5,6]
assert dot.dot(vec1,vec2)==numpy.dot(vec1,vec2)
