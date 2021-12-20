/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_CROSS_VALIDATION_ITERATOR_HPP
#define PECOS_UTIL_CROSS_VALIDATION_ITERATOR_HPP

#include "math_tools.hpp"
#include "OptionsList.hpp"
#include "teuchos_data_types.hpp"

namespace Pecos {
namespace util {

  /**
   * \class CrossValidationIterator
   * \brief Base class for k-folds cross validation.
   */
class CrossValidationIterator{
protected:
  /// The number of folds to use in cross-validation
  int numFolds_;

  /// The total number of points
  int numPts_;

  /// The randomly permuted indices which are used to partition the active
  /// points (points without total failures)
  IntVector indices_;

  /// The index of the first point indices in indices_ for each fold
  IntVector foldStartingIndices_;

  /// Seed for the random number generator
  int seed_;

  /// The amount of data for each point
  int numEquationsPerPoint_;

  /// Compressed storage of what data failed for each point. Base implementation
  /// only considers value failures. 
  IntVector failedRespData_;

public:

  /// Default Constructor
  CrossValidationIterator();

  /// Destructor
  ~CrossValidationIterator();

  /**
   * \brief Set the number of cross-validation folds
   * \param num_folds - The number of folds
   */
  void set_num_folds(int num_folds);
  
  /**
   * \brief Set the total number of points
   * \param num_points - The number of points
   */
  void set_num_points(int num_points);

  /**
   * \brief Set the seed of the random number generator
   * \param seed The random number seed
   */
  void set_seed(int seed);

  /**
   *\brief Return the number of cross-validation folds
   */
  int num_folds() const;

  /**
   *\brief Return the total number of points
   */
  int num_points() const;

  /**
   *\brief Clear the cached data
   */
  void clear();
  
  /**
   * \brief Return the randomly permuted indices which are used to partition
   * the points
   *
   * \return result vector of length (numPts_) -
   *     The randomly permuted indices
   */
  void get_point_indices(IntVector &result)  const;

  /**
   * \brief Return the indices of each point in a fold and the indices of all 
      points not in that fold
   *
   * \param iter -
   *     The fold number 
   *
   * \return result_0 vector of length (numPts-l) - 
   *     The training indices i.e indices not in fold iter
   * \return result_1 vector of length (l) - 
   *     The validation indices i.e. indices in fold iter
   */
  void get_fold_indices(int iter, IntVector &result_0,
			IntVector &result_1)  const;

    /**
   * \brief Return the indices of each point in a fold and the indices of all 
      points not in that fold
   *
   * \param iter -
   *     The fold number 
   *
   * \return result_0 vector of length (numPts-l) - 
   *     The training indices, i.e indices not in fold iter
   */
  void get_fold_training_indices(int iter, IntVector &result_0)  const;


    /**
   * \brief Return the indices of each point in a fold and the indices of all 
      points not in that fold
   *
   * \param iter -
   *     The fold number 
   *
   * \return result_0 vector of length (l) - 
   *     The validation indices i.e. indices in fold iter
   */
  void get_fold_validation_indices(int iter, IntVector &result_0)  const;

  
  /**
   *\brief Copy a cross validation iterator
   */
  void copy(const CrossValidationIterator &source);

  /**
   *\brief Get the number of points in a fold and the total number of points in 
   * all other folds
   *
   * \param iter The fold number
   *
   * \return training_size The number of points NOT in the ith fold
   * \return validation_size The number of points in the ith fold
   */
  void get_fold_size(int iter, int &training_size, int& validation_size) const;

  /**
   *\brief Return the set of values associated with a set of chosen indices
   *
   * \param values matrix of size (numPts_ x num_rhs) - 
   *     The values associated  with each point
   * \param indices vector of length (l<numPts) - 
   *     The selected indices
   *
   * \return result_0 vector of length (l x num_rhs) -
   *     The values at the selected points
   */
  void extract_values(const RealMatrix &B,
		      const IntVector &indices,
		      RealMatrix &result_0) const;

  /**
   *\brief Return the set of points associated with a set of chosen indices
   *
   * \param poitns matrix of length (num_vars x numPts_) - 
   *     The coordinates of each point
   * \param indices vector of length (l<numPts) - 
   *     The selected indices
   *
   * \return result_0 matrix of length (num_varrs x l) -
   *     The coordinates of the selected points
   */
  void extract_points(const RealMatrix &points,
		      const IntVector &training_indices, 
		      RealMatrix &result_0) const;

  /**
   * \brief Compute the cross validation score from a set of fold differences.
   * \param fold_diffs matrix of size (m x n) - 
   *     The differences between the truth value at a validation point and 
   *     the approximation using k-1 folds. n is the number of scenarios
   *     tested by cross validation. Typically n=1.
   *
   * \return result vector of size (n) - 
   *     The cross validation scores of the fold for each of the n 
   *     scenarios tested
   */
  void compute_fold_score(const RealMatrix &fold_diffs, RealVector &result) const;

  void accumulate_fold_scores(const RealMatrix &fold_scores, RealVector &scores)const;

  /**
  * \brief Pass in knowledge about faulty data
  * \param failed_response_data vector of size (numPts_) - 
  *     Non-negative integer values used with a bit wise comparison
  *     to indicate the presence of failed values (1) and failed gradients (2)
  */
  void set_fault_data(const IntVector &failed_response_data);

  /// Create partition into k folds taking into account the presence of
  /// and faulty data
  void create_partitions();

};

}  // namespace util
}  // namespace Pecos

#endif  // include guard
