#include "MUQ/Modeling/LinearAlgebra/SliceOperator.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

TEST(Utilties_LinearOperator, SliceOperator)
{

    const unsigned int inputDim = 100;
    const int startInd = 4;
    const int endInd = 85;
    const int skip = 3;

    Eigen::VectorXd x = Eigen::VectorXd::Random(inputDim);

    auto op = std::make_shared<SliceOperator>(inputDim, startInd, endInd, skip);

    Eigen::VectorXd y = op->Apply(x);
    ASSERT_EQ(std::ceil(double(endInd-startInd)/double(skip)), y.rows());

    for(unsigned int i=0; i<y.rows(); ++i)
      EXPECT_EQ(x(skip*i+startInd), y(i));


    Eigen::VectorXd x2 = op->ApplyTranspose(y);
    ASSERT_EQ(x.rows(), x2.rows());

    for(unsigned int i=startInd; i<endInd; i+=skip)
      EXPECT_EQ(x(i), x2(i));
}

TEST(Utilties_LinearOperator, SliceOperatorReverse)
{

    const unsigned int inputDim = 100;
    const int startInd = 80;
    const int endInd = 21;
    const int skip = -3;

    Eigen::VectorXd x = Eigen::VectorXd::Random(inputDim);

    auto op = std::make_shared<SliceOperator>(inputDim, startInd, endInd, skip);

    Eigen::VectorXd y = op->Apply(x);
    ASSERT_EQ(std::ceil(double(endInd-startInd)/double(skip)), y.rows());

    for(unsigned int i=0; i<y.rows(); ++i)
      EXPECT_EQ(x(skip*i+startInd), y(i));


    Eigen::VectorXd x2 = op->ApplyTranspose(y);
    ASSERT_EQ(x.rows(), x2.rows());

    for(unsigned int i=startInd; i>endInd; i+=skip)
      EXPECT_EQ(x(i), x2(i));
}

TEST(Utilties_LinearOperator, SliceOperatorNegativeIndex)
{

    const unsigned int inputDim = 100;
    const int startInd = 80;
    const int endInd = -10;
    const int skip = 3;

    Eigen::VectorXd x = Eigen::VectorXd::Random(inputDim);

    auto op = std::make_shared<SliceOperator>(inputDim, startInd, endInd, skip);

    Eigen::VectorXd y = op->Apply(x);
    ASSERT_EQ(std::ceil(double(inputDim+endInd-startInd)/double(skip)), y.rows());

    for(unsigned int i=0; i<y.rows(); ++i)
      EXPECT_EQ(x(skip*i+startInd), y(i));


    Eigen::VectorXd x2 = op->ApplyTranspose(y);
    ASSERT_EQ(x.rows(), x2.rows());

    for(unsigned int i=startInd; i>inputDim+endInd; i+=skip)
      EXPECT_EQ(x(i), x2(i));
}
