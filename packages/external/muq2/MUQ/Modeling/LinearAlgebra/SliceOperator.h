#ifndef SLICEOPERATOR_H
#define SLICEOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"


namespace muq
{
namespace Modeling
{

    /** @class SliceOperator
        @brief Defines a "slice" or "range" of an input vector.
        @details Using numpy notation, for an input vector x, this ModPiece returns
                 a vector containg x[startInd:endInd:skip].  Like numpy, this class
                 supports negative indices to access the last components of the
                 input vector.  For example x[0:-1] will return all components of
                 x except the last component.
    */
    class SliceOperator : public LinearOperator
    {
    public:
        SliceOperator(unsigned int vecSize,
                      int startIndIn,
                      int endIndIn,
                      int skipIn=1);

        virtual ~SliceOperator() = default;

        /** Apply the linear operator to a vector */
        virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

        /** Apply the transpose of the linear operator to a vector. */
        virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

        virtual Eigen::MatrixXd GetMatrix() override;

    private:
        static unsigned int ComputeRows(unsigned int vecSize, int startIndIn, int endIndIn, int skipIn);
        
        unsigned int startInd, endInd;
        int skip;

    }; // class SliceOperator


} // namespace Modeling
} // namespace muq

#endif // #ifndef SUMOPERATOR_H
