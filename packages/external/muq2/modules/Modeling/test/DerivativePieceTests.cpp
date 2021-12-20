#include "gtest/gtest.h"

#include "MUQ/Modeling/GradientPiece.h"
#include "MUQ/Modeling/JacobianPiece.h"
#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"

using namespace muq::Modeling;

TEST(DerivativePiecesTest, GradientPiece)
{
  unsigned int dim = 10;

  auto basePiece = std::make_shared<SinOperator>(dim);
  auto gradPiece = std::make_shared<GradientPiece>(basePiece,0,0);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);

  Eigen::VectorXd sens = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd trueGrad = basePiece->Gradient(0,0,input,sens);
  Eigen::VectorXd testGrad = gradPiece->Evaluate(input,sens).at(0);


  EXPECT_EQ(0,basePiece->GetNumCalls("HessianAction"));
  Eigen::VectorXd v = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd hessApply = gradPiece->ApplyJacobian(0,0,input,sens,v);

  EXPECT_EQ(1,basePiece->GetNumCalls("HessianAction"));
  EXPECT_EQ(0,basePiece->GetNumCalls("HessianActionFD"));

  for(unsigned int i=0; i<dim; ++i){
    EXPECT_DOUBLE_EQ(trueGrad(i), testGrad(i));
  }
}

TEST(DerivativePiecesTest, JacobianPiece)
{
  unsigned int dim = 10;

  auto basePiece = std::make_shared<SinOperator>(dim);
  auto jacPiece = std::make_shared<JacobianPiece>(basePiece,0,0);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);

  Eigen::VectorXd v = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd trueJac = basePiece->ApplyJacobian(0,0,input,v);

  Eigen::VectorXd testJac = jacPiece->Evaluate(input,v).at(0);


  for(unsigned int i=0; i<dim; ++i){
    EXPECT_DOUBLE_EQ(trueJac(i), testJac(i));
  }
}
