/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SAMPLING_INVERSE_TRANSFORMATION_HPP
#define SAMPLING_INVERSE_TRANSFORMATION_HPP

#include "InverseTransformation.hpp"


namespace Pecos {


/// Class for sampling-based data transformation.

/** The SamplingInverseTransformation employs a sampling-based inverse
    transform to map from the frequency domain to the time domain. */

class SamplingInverseTransformation: public InverseTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  SamplingInverseTransformation();  ///< constructor
  ~SamplingInverseTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //


private:

  //
  //- Heading: Utility routines
  //

};


inline SamplingInverseTransformation::SamplingInverseTransformation()
{ }


inline SamplingInverseTransformation::~SamplingInverseTransformation()
{ }

} // namespace Pecos

#endif
