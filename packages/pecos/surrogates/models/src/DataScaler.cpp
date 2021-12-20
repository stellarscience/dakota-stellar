/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

// Two scalers:
// - Normalizing scaler - performs mean or max-min normalization
// - Standarization - makes each feature zero mean and unit variance
#include <DataScaler.hpp>
#include "math_tools.hpp"

namespace Pecos {
namespace surrogates {

DataScaler::DataScaler(){}

DataScaler::~DataScaler(){}

RealVector DataScaler::scaleFeatures(const RealVector &unscaled_x) {
  int M = unscaled_x.length();
  RealVector scaledInput(M);
  for (int i = 0; i < M; i++) {
    scaledInput[i] = (unscaled_x(i) - scalerFeaturesOffsets(i))/
                     scalerFeaturesScaleFactors(i);
  }
  return scaledInput;
}

NormalizationScaler::NormalizationScaler(){}

NormalizationScaler::~NormalizationScaler(){}

NormalizationScaler::NormalizationScaler(const RealMatrix &features, 
                                         const bool mean_normalization, const Real norm_factor) {

  int M = features.numRows();
  int K = features.numCols();

  scalerFeaturesOffsets.resize(M);
  scalerFeaturesScaleFactors.resize(M);
  scaledFeatures.reshape(M, K);

  Real min_elem, max_elem, mean_elem;
  
  RealMatrix features_transpose(features, Teuchos::TRANS);

  for (int i = 0; i < M; i++) {
    Real* f = features_transpose[i];
    min_elem = f[util::argmin(K, f)];
    max_elem = f[util::argmax(K, f)];
    mean_elem = util::mean(K, f);
    scalerFeaturesOffsets(i) = (mean_normalization) ? mean_elem : min_elem;
    scalerFeaturesScaleFactors(i) = (max_elem - min_elem)/norm_factor;
    for (int j = 0; j < K; j++) {
      scaledFeatures(i,j) = (features(i,j) - scalerFeaturesOffsets(i))/
                            scalerFeaturesScaleFactors(i);
   }
  }

  has_scaling = true;
}

StandardizationScaler::StandardizationScaler(){}

StandardizationScaler::~StandardizationScaler(){}

StandardizationScaler::StandardizationScaler(const RealMatrix &features, 
                                             const Real norm_factor) {

  int M = features.numRows();
  int K = features.numCols();

  scalerFeaturesOffsets.resize(M);
  scalerFeaturesScaleFactors.resize(M);
  scaledFeatures.reshape(M,K);

  Real mean_val, var_val;
  
  RealMatrix features_transpose(features, Teuchos::TRANS);

  for (int i = 0; i < M; i++) {
    Real* f = features_transpose[i];
    mean_val = util::mean(K, f);
    var_val = util::variance(K, f, 0);
    scalerFeaturesOffsets(i) = mean_val;
    scalerFeaturesScaleFactors(i) = std::sqrt(var_val)/norm_factor;
    for (int j = 0; j < K; j++) {
      scaledFeatures(i,j) = (features(i,j) - scalerFeaturesOffsets(i))/
                            scalerFeaturesScaleFactors(i);
   }
  }

  has_scaling = true;
}


}  // namespace surrogates
}  // namespace Pecos
