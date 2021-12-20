/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_DATA_SCALER_HPP
#define PECOS_SURROGATES_DATA_SCALER_HPP

#include "teuchos_data_types.hpp"

namespace Pecos {
namespace surrogates {

/**
 This class contains functions that scale the input data.
*/

class DataScaler {

  public: 
    DataScaler();

    ~DataScaler();

    RealVector scaleFeatures(const RealVector &x);

    RealVector getScalerFeaturesOffsets() {return scalerFeaturesOffsets;}
    RealVector getScalerFeaturesScaleFactors() {return scalerFeaturesScaleFactors;}
    RealMatrix getScaledFeatures() {return scaledFeatures;}

    bool has_scaling = false;

  protected: 

    RealVector scalerFeaturesOffsets;
    RealVector scalerFeaturesScaleFactors;
    RealMatrix scaledFeatures;
};

class NormalizationScaler: public DataScaler {
  public:
    NormalizationScaler();

    ~NormalizationScaler();

    NormalizationScaler(const RealMatrix &features, const bool mean_normalization, 
                        const Real norm_factor = 1.0);
};

class StandardizationScaler: public DataScaler {
  public:
    StandardizationScaler();

    ~StandardizationScaler();

    StandardizationScaler(const RealMatrix &features, const Real norm_factor = 1.0);
};


}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
