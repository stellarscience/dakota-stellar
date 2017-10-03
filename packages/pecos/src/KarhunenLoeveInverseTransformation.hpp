/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef KARHUNEN_LOEVE_INVERSE_TRANSFORMATION_HPP
#define KARHUNEN_LOEVE_INVERSE_TRANSFORMATION_HPP

#include "InverseTransformation.hpp"


namespace Pecos {


/// Class for KL data transformation.

/** The KarhunenLoeveInverseTransformation employs an KL decomposition
    to map from the frequency domain to the time domain. */

class KarhunenLoeveInverseTransformation: public InverseTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  KarhunenLoeveInverseTransformation();  ///< constructor
  ~KarhunenLoeveInverseTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //


private:

  //
  //- Heading: Utility routines
  //

};


inline KarhunenLoeveInverseTransformation::KarhunenLoeveInverseTransformation()
{ }


inline KarhunenLoeveInverseTransformation::~KarhunenLoeveInverseTransformation()
{ }

} // namespace Pecos

#endif
