#ifndef SPLITVECTOR_H_
#define SPLITVECTOR_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {

    /** @class SplitVector
        @brief Provides a mechanism for splitting a vector into multiple pieces.
        @details Imagine you have a vector \f$x\f$ of length $M+N$ that you want to
        into two vectors \f$y_1\f$ and \f$y_2\f$ with lengths \f$M\f$ and \f$N\f$,
        respectively.  This class enables that type of operation.

        Typical usage:
        @code{.cpp}

        // Construct a vector x = [1,2,3,4,5]
        Eigen::VectorXd x(5);
        x << 1, 2, 3, 4, 5;

        // Create a SplitVector instance that will extract y1=[1,2,3] and y2=[4,5]
        Eigen::VectorXi startInds(2);
        startInds << 0, 3; // The indices of x where y1 and y2 start

        // Set the lengths of each segment.  Note that if this was 2,2 then y1=[1,2] and y2=[4,5]
        Eigen::VectorXi sizes(2);
        sizes << 3,2;

        // Create the SplitVector instance
        auto split = std::make_shared<SplitVector>(startInds, sizes, 5);

        // Evaluate the SplitVector instance
        std::vector<Eigen::VectorXd> ys = split->Evaluate(x);

        std::cout << "y1 = [" << ys.at(0).transpose() << "]" << std::endl; // should print [1,2,3]
        std::cout << "y2 = [" << ys.at(1).transpose() << "]" << std::endl; // should print [4,5]

        @endcode

    */
    class SplitVector : public ModPiece {
    public:

      /**
        @param[in] ind The first index of the segment for each output
        @param[in] size The size of each segment
        @param[in] insize The size of the input vector
      */
      SplitVector(Eigen::VectorXi const& ind, Eigen::VectorXi const& size, unsigned int const insize);

      virtual ~SplitVector() = default;

      /** Returns a vector containing the indices of the input vecto where each
          segment returned in the output begins.
      */
      Eigen::VectorXi StartIndices() const{return ind;};

    private:

      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void JacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void GradientImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& sens) override;

      virtual void ApplyJacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& targ) override;

      const Eigen::VectorXi ind;

      const Eigen::VectorXi size;

    };
  } // namespace Modeling
} // namespace muq

#endif
