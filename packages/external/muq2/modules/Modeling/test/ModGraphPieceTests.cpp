#include "gtest/gtest.h"

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Modeling/IdentityPiece.h"

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"
#include "MUQ/Modeling/SplitVector.h"
#include "MUQ/Modeling/SumPiece.h"

#include "WorkPieceTestClasses.h"

using namespace muq::Modeling;
using namespace std;

// define a simple single input forward model
class SquareMod : public ModPiece {
public:

  /** Constructor taking vector dimension and resizing the State.*/
  SquareMod(int dim) : ModPiece(dim*Eigen::VectorXi::Ones(1), dim*Eigen::VectorXi::Ones(1)) {}

  virtual ~SquareMod() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override
  {
    outputs.resize(1);
    outputs.at(0) = input.at(0).get().array().square();
  }

  virtual void GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override
  {
    gradient = 2.0 * sensitivity.array() * input.at(0).get().array();
  }

  virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input) override
  {
    jacobian = 2.0*input.at(0).get().asDiagonal();
  }

  virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec) override
  {
    jacobianAction = 2.0*input.at(0).get().array()*vec.array();
  };

  virtual void ApplyHessianImpl(unsigned int                const  outWrt,
                                unsigned int                const  inWrt1,
                                unsigned int                const  inWrt2,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd             const& sens,
                                Eigen::VectorXd             const& vec) override
  {
    // If the hessian is with respect to the sensitivity
    if(inWrt2==inputSizes.size()){
      hessAction = 2.0*input.at(0).get().array() * vec.array();

    // Otherwise, the hessian is with respect to the input
    }else{
      hessAction = 2.0 * sens.array() * vec.array();
    }

  }

};


// define a simple single input forward model
class SinSumMod : public ModPiece {
public:

  /** Constructor taking vector dimension and resizing the State.*/
  SinSumMod(int dim) : ModPiece(dim*Eigen::VectorXi::Ones(2), dim*Eigen::VectorXi::Ones(1)) {}

  virtual ~SinSumMod() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override
  {
    Eigen::VectorXd out = Eigen::VectorXd::Zero(outputSizes(0));

    for (unsigned int j = 0; j < input.size(); ++j) {
      for (int i = 0; i < outputSizes(0); ++i) {
        out[i] += sin(input.at(j)(i));
      }
    }

    outputs.resize(1);
    outputs.at(0) = out;
  }

  virtual void GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity) override
  {
    gradient = Eigen::VectorXd::Zero(outputSizes(0));

    for (int i = 0; i < outputSizes(0); ++i)
      gradient(i) = sensitivity(i) * cos(input.at(inputDimWrt)(i));
  }

  virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input) override
  {
    jacobian = Eigen::MatrixXd::Zero(outputSizes(0), inputSizes(inputDimWrt));
    for (int i = 0; i < outputSizes(0); ++i)
      jacobian(i, i) = cos(input.at(inputDimWrt)(i));
  }

  virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec) override
  {
    jacobianAction = Eigen::VectorXd::Zero(outputSizes(0));
    for(int i=0; i<outputSizes(0); ++i)
      jacobianAction(i) = cos(input.at(inputDimWrt)(i))*vec(i);
  };

  virtual void ApplyHessianImpl(unsigned int                const  outWrt,
                                unsigned int                const  inWrt1,
                                unsigned int                const  inWrt2,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd             const& sens,
                                Eigen::VectorXd             const& vec) override
  {
    hessAction = Eigen::VectorXd::Zero(vec.size());

    // If the second derivative is with respect to the sensitivity
    if(inWrt2==inputSizes.size()){

      for (int i = 0; i < outputSizes(0); ++i)
        hessAction(i) = cos(input.at(inWrt1)(i))*vec(i);

    // Otherwise, if the derivatives are with respect to the same input
    }if(inWrt1==inWrt2){
      for (int i = 0; i < outputSizes(0); ++i)
        hessAction(i) = -1.0*sens(i) * sin(input.at(inWrt1)(i))*vec(i);
    }
  }

};


TEST(Modeling_ModGraphPiece, MatchInputs)
{
  auto myGraph = make_shared<WorkGraph>();

  // add nodes
  myGraph->AddNode(make_shared<SinSumMod>(2), "f2");
  myGraph->AddNode(make_shared<SinSumMod>(2), "f1");
  myGraph->AddNode(make_shared<SquareMod>(2), "x");
  myGraph->AddNode(make_shared<SquareMod>(2), "y1");
  myGraph->AddNode(make_shared<SquareMod>(2), "y2");

  // add connectivity
  myGraph->AddEdge("x", 0, "f2", 0);
  myGraph->AddEdge("y1", 0, "f1", 0);
  myGraph->AddEdge("y2", 0, "f1", 1);
  myGraph->AddEdge("f1", 0, "f2", 1);

  auto piece1 = myGraph->CreateModPiece("f1");
  auto piece2 = myGraph->CreateModPiece("f2");

  EXPECT_EQ(2, piece1->inputSizes.size());
  EXPECT_EQ(3, piece2->inputSizes.size());

  std::vector<int> sharedIns = piece1->MatchInputs(piece2);
  EXPECT_EQ(3, sharedIns.size());
  EXPECT_EQ(0, sharedIns.at(0));
  EXPECT_EQ(1, sharedIns.at(1));
  EXPECT_EQ(-1, sharedIns.at(2));

  sharedIns = piece2->MatchInputs(piece1);
  EXPECT_EQ(2, sharedIns.size());
  EXPECT_EQ(0, sharedIns.at(0));
  EXPECT_EQ(1, sharedIns.at(1));
}


TEST(Modeling_ModGraphPiece, BasicTest)
{
  auto myGraph = make_shared<WorkGraph>();

  // add nodes
  auto f1 = make_shared<SinSumMod>(2);
  auto f2 = make_shared<SinSumMod>(2);
  myGraph->AddNode(f2, "f2");
  myGraph->AddNode(f1, "f1");
  myGraph->AddNode(make_shared<SquareMod>(2), "x");
  myGraph->AddNode(make_shared<SquareMod>(2), "y1");
  myGraph->AddNode(make_shared<SquareMod>(2), "y2");

  // add connectivity
  myGraph->AddEdge("x", 0, "f2", 0);
  myGraph->AddEdge("y1", 0, "f1", 0);
  myGraph->AddEdge("y2", 0, "f1", 1);
  myGraph->AddEdge("f1", 0, "f2", 1);

  auto graphMod = myGraph->CreateModPiece("f2");
  auto gradPiece = graphMod->GradientGraph(0,0);
  auto jacPiece = graphMod->JacobianGraph(0,0);

  EXPECT_EQ(3, graphMod->inputSizes.size());
  EXPECT_EQ(2, graphMod->inputSizes(0));
  EXPECT_EQ(2, graphMod->outputSizes(0));

  boost::any anyVal1 = Eigen::VectorXd::Constant(2,0.1).eval();
  boost::any anyVal2 = Eigen::VectorXd::Constant(2,0.2).eval();

  myGraph->BindNode("y1", {anyVal1});
  myGraph->BindNode("y2", {anyVal2});

  graphMod = myGraph->CreateModPiece("f2");
  //myGraph->Visualize("BasicTest2.pdf");
  gradPiece = graphMod->GradientGraph(0,0);
  jacPiece = graphMod->JacobianGraph(0,0);
  //gradPiece->GetGraph()->Visualize("BasicTestGrad2.pdf");

  // make sure this modpiece is the size we expect
  EXPECT_EQ(1, graphMod->inputSizes.size());
  EXPECT_EQ(2, graphMod->inputSizes(0));
  EXPECT_EQ(2, graphMod->outputSizes(0));

  // evaluation testing
  Eigen::VectorXd input  = Eigen::VectorXd::Ones(2);
  Eigen::VectorXd output = graphMod->Evaluate(input).at(0);

  EXPECT_DOUBLE_EQ(sin(input[0] * input[0]) + sin(sin(0.1) + sin(0.2)), output[0]);
  EXPECT_DOUBLE_EQ(sin(input[1] * input[1]) + sin(sin(0.1) + sin(0.2)), output[1]);

  // gradient testing (same as J^T*x)
  Eigen::VectorXd vec = Eigen::VectorXd::Ones(2);
  Eigen::VectorXd jacAct = graphMod->ApplyJacobian(0,0,input,vec);
  Eigen::VectorXd grad = graphMod->Gradient(0,0, input, Eigen::VectorXd::Ones(2).eval());

  EXPECT_DOUBLE_EQ(2.0 * input[0] * cos(input[0] * input[0]), grad[0]);
  EXPECT_DOUBLE_EQ(2.0 * input[1] * cos(input[1] * input[1]), grad[1]);

  Eigen::VectorXd jacAct2 = jacPiece->Evaluate(input,vec).at(0);
  Eigen::VectorXd grad2 = gradPiece->Evaluate(input, Eigen::VectorXd::Ones(2).eval()).at(0);
  EXPECT_DOUBLE_EQ(2.0 * input[0] * cos(input[0] * input[0]), grad2[0]);
  EXPECT_DOUBLE_EQ(2.0 * input[1] * cos(input[1] * input[1]), grad2[1]);

  EXPECT_NEAR(jacAct(0), jacAct2(0), 5e-8);
  EXPECT_NEAR(jacAct(1), jacAct2(1), 5e-8);

}

TEST(Modeling_ModGraphPiece, Caching)
{

  auto graph = make_shared<WorkGraph>();
  auto x = make_shared<IdentityOperator>(2);
  auto f = make_shared<ExpOperator>(2);
  auto g = make_shared<SinOperator>(2);

  graph->AddNode(x, "x");
  graph->AddNode(f, "f");
  graph->AddNode(g, "g");
  graph->AddEdge("x",0,"f",0);
  graph->AddEdge("f",0,"g",0);

  auto graphMod = graph->CreateModPiece("g");
  auto gradPiece = graphMod->GradientGraph(0,0);

  Eigen::VectorXd input(2);
  input << 0.5, 0.75;
  Eigen::VectorXd sens(2);
  sens << 1.0, 1.0;

  EXPECT_EQ(0,f->GetNumCalls("Evaluate"));
  graphMod->Evaluate(input);
  EXPECT_EQ(1,f->GetNumCalls("Evaluate"));
  graphMod->Gradient(0,0,input,sens);
  EXPECT_EQ(2,f->GetNumCalls("Evaluate"));
  graphMod->Gradient(0,0,input,sens);
  EXPECT_EQ(3,f->GetNumCalls("Evaluate"));

  f->EnableCache();
  graphMod->Evaluate(input);
  EXPECT_EQ(4,f->GetNumCalls("Evaluate"));
  graphMod->Gradient(0,0,input,sens);
  EXPECT_EQ(4,f->GetNumCalls("Evaluate"));
}


TEST(Modeling_ModGraphPiece, MultiSplit)
{
  auto myGraph = make_shared<WorkGraph>();

  // add nodes
  myGraph->AddNode(std::make_shared<IdentityOperator>(4), "x");
  myGraph->AddNode(std::make_shared<SplitVector>(Eigen::Vector2i{0,2}, Eigen::Vector2i{2,2}, 4), "x1,x2");
  myGraph->AddNode(std::make_shared<ExpOperator>(2), "exp(x1)");
  myGraph->AddNode(std::make_shared<ExpOperator>(2), "exp(x2)");
  myGraph->AddNode(std::make_shared<CosOperator>(2), "cos(x1)");
  myGraph->AddNode(std::make_shared<CosOperator>(2), "cos(x2)");
  myGraph->AddNode(std::make_shared<SumPiece>(2,4), "sum");

  // add connectivity
  myGraph->AddEdge("x", 0, "x1,x2", 0);
  myGraph->AddEdge("x1,x2", 0, "exp(x1)", 0);
  myGraph->AddEdge("x1,x2", 0, "cos(x1)", 0);
  myGraph->AddEdge("x1,x2", 1, "exp(x2)", 0);
  myGraph->AddEdge("x1,x2", 1, "cos(x2)", 0);
  myGraph->AddEdge("exp(x1)",0,"sum",0);
  myGraph->AddEdge("cos(x1)",0,"sum",1);
  myGraph->AddEdge("exp(x2)",0,"sum",2);
  myGraph->AddEdge("cos(x2)",0,"sum",3);

  auto graphMod = myGraph->CreateModPiece("sum");
  auto gradPiece = graphMod->GradientGraph(0,0);
  auto jacPiece = graphMod->JacobianGraph(0,0);

  EXPECT_EQ(1, graphMod->inputSizes.size());
  EXPECT_EQ(4, graphMod->inputSizes(0));
  EXPECT_EQ(2, graphMod->outputSizes(0));


  gradPiece->GetGraph()->Visualize("MulitSplitGraph.pdf");
}


// TEST(Modeling_ModGraphPiece, NodeOrdering)
// {
//   auto myGraph = make_shared<ModGraph>();
//
//   // add nodes
//   myGraph->AddNode(make_shared<sinSumMod>(2), "f1");
//   myGraph->AddNode(make_shared<squareMod>(2), "z");
//   myGraph->AddNode(make_shared<squareMod>(2), "x1");
//
//   // add connectivity
//   myGraph->AddEdge("z", "f1", 0);
//   myGraph->AddEdge("x1", "f1", 1);
//
//   myGraph->writeGraphViz("results/tests/GraphViz/NodeOrderTest.pdf");
//
//   auto GraphMod = ModGraphPiece::Create(myGraph, "f1");
//   GraphMod->writeGraphViz("results/tests/GraphViz/NodeOrderPieceTest.pdf");
//
//   std::vector<std::string> newOrder(2);
//   newOrder[0] = "z";
//   newOrder[1] = "x1";
//   auto GraphMod2 = ModGraphPiece::Create(myGraph, "f1", newOrder);
//   GraphMod->writeGraphViz("results/tests/GraphViz/NodeOrderPieceTest2.pdf");
// }
//
//
TEST(Modeling_ModGraphPiece, DiamondTest)
{
  auto myGraph = make_shared<WorkGraph>();

  // add nodes
  myGraph->AddNode(make_shared<SinSumMod>(2), "f1");
  myGraph->AddNode(make_shared<SquareMod>(2), "x");
  myGraph->AddNode(make_shared<SquareMod>(2), "y1");
  myGraph->AddNode(make_shared<SquareMod>(2), "y2");

  // add connectivity
  myGraph->AddEdge("x", 0, "y1", 0);
  myGraph->AddEdge("x", 0, "y2", 0);
  myGraph->AddEdge("y1", 0, "f1", 0);
  myGraph->AddEdge("y2", 0, "f1", 1);

  //myGraph->Visualize("DiamondForward.pdf");
  auto graphMod = myGraph->CreateModPiece("f1");
  //graphMod->GetGraph()->Visualize("DiamondPieceTest.pdf");

  auto gradPiece = graphMod->GradientGraph(0,0);
  auto jacPiece = graphMod->JacobianGraph(0,0);
  //gradPiece->GetGraph()->Visualize("DiamondGrad.pdf");
  //gradPiece->JacobianGraph(0,0)->GetGraph()->Visualize("DiamondHess.pdf");

  // make sure this modpiece is the size we expect
  EXPECT_EQ(1, graphMod->inputSizes.size());
  EXPECT_EQ(2, graphMod->inputSizes(0));
  EXPECT_EQ(2, graphMod->outputSizes(0));

  // evaluation testing
  Eigen::VectorXd input  = 0.5 * Eigen::VectorXd::Ones(2);
  Eigen::VectorXd output = graphMod->Evaluate(input).at(0);

  EXPECT_DOUBLE_EQ(2.0 * sin(pow(input[0], 4.0)), output[0]);
  EXPECT_DOUBLE_EQ(2.0 * sin(pow(input[1], 4.0)), output[1]);

  Eigen::VectorXd ones = Eigen::VectorXd::Ones(2);
  Eigen::VectorXd grad1 = graphMod->Gradient(0,0,input,ones);
  Eigen::VectorXd grad2 = gradPiece->Evaluate(input,ones).at(0);

  EXPECT_NEAR(grad1(0), grad2(0), 5e-8);
  EXPECT_NEAR(grad1(1), grad2(1), 5e-8);

  Eigen::VectorXd jac1 = graphMod->ApplyJacobian(0,0,input,ones);
  Eigen::VectorXd jac2 = jacPiece->Evaluate(input,ones).at(0);

  EXPECT_NEAR(jac1(0), jac1(0), 5e-8);
  EXPECT_NEAR(jac1(1), jac1(1), 5e-8);

  // Hessian checks
  Eigen::VectorXd vec = Eigen::VectorXd::Random(input.size());
  Eigen::VectorXd hessAct = graphMod->ApplyHessian(0, 0, 0, std::vector<Eigen::VectorXd>{input}, ones, vec);
  Eigen::VectorXd hessActFD = graphMod->ApplyHessianByFD(0, 0, 0, std::vector<Eigen::VectorXd>{input}, ones, vec);

  EXPECT_NEAR(hessActFD(0), hessAct(0), 1e-7);
  EXPECT_NEAR(hessActFD(1), hessAct(1), 1e-7);
}


TEST(Modeling_ModGraphPiece, LinearGaussPosterior)
{
  const unsigned int numNodes = 200;
  const unsigned int dataDim = 10;
  const double noiseStd = 1e-3;

  auto forwardMod = LinearOperator::Create(Eigen::MatrixXd::Random(dataDim,numNodes).eval());
  auto prior = std::make_shared<Gaussian>(Eigen::VectorXd::Zero(numNodes), Eigen::VectorXd::Ones(numNodes));

  Eigen::VectorXd data = forwardMod->Apply(prior->Sample());
  auto noiseModel = std::make_shared<Gaussian>(data, noiseStd*Eigen::VectorXd::Ones(dataDim));

  WorkGraph graph;
  graph.AddNode(std::make_shared<IdentityOperator>(numNodes), "Parameters");
  graph.AddNode(prior->AsDensity(), "Prior");
  graph.AddNode(forwardMod, "ForwardModel");
  graph.AddNode(noiseModel->AsDensity(), "Likelihood");
  graph.AddNode(std::make_shared<DensityProduct>(2),"Posterior");
  graph.AddEdge("Parameters",0,"Prior",0);
  graph.AddEdge("Parameters",0,"ForwardModel",0);
  graph.AddEdge("ForwardModel",0,"Likelihood",0);
  graph.AddEdge("Prior",0,"Posterior",0);
  graph.AddEdge("Likelihood",0,"Posterior",1);

  auto logLikely = graph.CreateModPiece("Likelihood");
  logLikely->GradientGraph(0,0)->GetGraph()->Visualize("Posterior_GradientGraph.png");
  logLikely->GradientGraph(0,0)->JacobianGraph(0,0)->GetGraph()->Visualize("Posterior_HessianGraph.png");

  Eigen::VectorXd input = prior->Sample();

  // Hessian checks
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(1);
  Eigen::VectorXd vec = Eigen::VectorXd::Random(numNodes);
  Eigen::VectorXd hessAct = logLikely->ApplyHessian(0, 0, 0, std::vector<Eigen::VectorXd>{input}, ones, vec);

  Eigen::VectorXd hessActFD = logLikely->ApplyHessianByFD(0, 0, 0, std::vector<Eigen::VectorXd>{input}, ones, vec);

  for(unsigned int i=0; i<numNodes; ++i)
    EXPECT_NEAR(hessActFD(i), hessAct(i), 1e-2*std::abs(hessActFD(i)));
}


TEST(Modeling_ModGraphPiece, DiamondTestOperator)
{
  auto myGraph = make_shared<WorkGraph>();

  // add nodes
  myGraph->AddNode(make_shared<SinSumMod>(2), "f1");
  myGraph->AddNode(make_shared<SquareOperator>(2), "x");
  myGraph->AddNode(make_shared<SquareOperator>(2), "y1");
  myGraph->AddNode(make_shared<SquareOperator>(2), "y2");

  // add connectivity
  myGraph->AddEdge("x", 0, "y1", 0);
  myGraph->AddEdge("x", 0, "y2", 0);
  myGraph->AddEdge("y1", 0, "f1", 0);
  myGraph->AddEdge("y2", 0, "f1", 1);

  //myGraph->Visualize("DiamondForward.pdf");
  auto graphMod = myGraph->CreateModPiece("f1");
  //graphMod->GetGraph()->Visualize("DiamondPieceTest.pdf");

  auto gradPiece = graphMod->GradientGraph(0,0);
  auto jacPiece = graphMod->JacobianGraph(0,0);
  //gradPiece->GetGraph()->Visualize("DiamondGrad.pdf");
  //gradPiece->JacobianGraph(0,0)->GetGraph()->Visualize("DiamondHess.pdf");

  // make sure this modpiece is the size we expect
  EXPECT_EQ(1, graphMod->inputSizes.size());
  EXPECT_EQ(2, graphMod->inputSizes(0));
  EXPECT_EQ(2, graphMod->outputSizes(0));

  // evaluation testing
  Eigen::VectorXd input  = 0.5 * Eigen::VectorXd::Ones(2);
  Eigen::VectorXd output = graphMod->Evaluate(input).at(0);

  EXPECT_DOUBLE_EQ(2.0 * sin(pow(input[0], 4.0)), output[0]);
  EXPECT_DOUBLE_EQ(2.0 * sin(pow(input[1], 4.0)), output[1]);

  Eigen::VectorXd ones = Eigen::VectorXd::Ones(2);
  Eigen::VectorXd grad1 = graphMod->Gradient(0,0,input,ones);
  Eigen::VectorXd grad2 = gradPiece->Evaluate(input,ones).at(0);

  EXPECT_NEAR(grad1(0), grad2(0), 5e-8);
  EXPECT_NEAR(grad1(1), grad2(1), 5e-8);

  Eigen::VectorXd jac1 = graphMod->ApplyJacobian(0,0,input,ones);
  Eigen::VectorXd jac2 = jacPiece->Evaluate(input,ones).at(0);

  EXPECT_NEAR(jac1(0), jac1(0), 5e-8);
  EXPECT_NEAR(jac1(1), jac1(1), 5e-8);

  // Hessian checks
  Eigen::VectorXd vec = Eigen::VectorXd::Random(input.size());
  Eigen::VectorXd hessAct = graphMod->ApplyHessian(0, 0, 0, std::vector<Eigen::VectorXd>{input}, ones, vec);
  Eigen::VectorXd hessActFD = graphMod->ApplyHessianByFD(0, 0, 0, std::vector<Eigen::VectorXd>{input}, ones, vec);

  EXPECT_NEAR(hessActFD(0), hessAct(0), 1e-7);
  EXPECT_NEAR(hessActFD(1), hessAct(1), 1e-7);
}

//
// TEST(Modeling_ModGraphPiece, UnionTest)
// {
//   auto a = make_shared<ModGraph>();
//   auto b = make_shared<ModGraph>();
//   auto sourceNode = make_shared<squareMod>(2);
//
//   // add nodes
//   a->AddNode(make_shared<sinSumMod>(2), "af2");
//   a->AddNode(make_shared<sinSumMod>(2), "af1");
//   a->AddNode(sourceNode, "x");
//   a->AddNode(make_shared<squareMod>(2), "ay1");
//   a->AddNode(make_shared<squareMod>(2), "ay2");
//
//   // add connectivity
//   a->AddEdge("x", "af2", 0);
//   a->AddEdge("ay1", "af1", 0);
//   a->AddEdge("ay2", "af1", 1);
//   a->AddEdge("af1", "af2", 1);
//
//     // add nodes
//   b->AddNode(make_shared<sinSumMod>(2), "bf2");
//   b->AddNode(make_shared<sinSumMod>(2), "bf1");
//   b->AddNode(sourceNode, "x");
//   b->AddNode(make_shared<squareMod>(2), "by1");
//   b->AddNode(make_shared<squareMod>(2), "by2");
//
//   // add connectivity
//   b->AddEdge("x", "bf2", 0);
//   b->AddEdge("by1", "bf1", 0);
//   b->AddEdge("by2", "bf1", 1);
//   b->AddEdge("bf1", "bf2", 1);
//
//   auto unionGraph = ModGraph::FormUnion(a,b);
//   unionGraph->writeGraphViz("results/tests/UnionTest.pdf");
//
//   EXPECT_EQ(5, unionGraph->NumInputs()); //this tests that the "x" node is not duplicated
//   EXPECT_EQ(2, unionGraph->NumOutputs()); //and we have both model outputs
// }
//
// TEST(Modeling_ModGraphPiece, UnionNameClashDeath)
// {
//   auto a = make_shared<ModGraph>();
//   auto b = make_shared<ModGraph>();
//
//   // add nodes - x is distinct but has the same name
//   a->AddNode(make_shared<sinSumMod>(2), "af2");
//   a->AddNode(make_shared<sinSumMod>(2), "af1");
//   a->AddNode(make_shared<squareMod>(2), "x");
//
//   b->AddNode(make_shared<squareMod>(2), "x");
//   b->AddNode(make_shared<squareMod>(2), "by1");
//   b->AddNode(make_shared<squareMod>(2), "by2");
//
//   ASSERT_DEATH(ModGraph::FormUnion(a,b), "(result->GetNodeModel(currentName) == b->ModelGraph[v]->piece)*");
//
// }
//
// TEST(Modeling_ModGraphPiece, UnionEdgeClashDeath)
// {
// 	//Both graphs share a node correctly, but both try to provide an input, so we don't know how to
// 	//uniquely resolve it and hence assert out
//   auto a = make_shared<ModGraph>();
//   auto b = make_shared<ModGraph>();
//   auto sourceNode = make_shared<squareMod>(2);
//
//   // add nodes
//   a->AddNode(make_shared<sinSumMod>(2), "af2");
//   a->AddNode(make_shared<sinSumMod>(2), "af1");
//   a->AddNode(sourceNode, "x");
//   a->AddNode(make_shared<squareMod>(2), "ay1");
//   a->AddNode(make_shared<squareMod>(2), "ay2");
//
//   // add connectivity
//   a->AddEdge("x", "af2", 0);
//   a->AddEdge("ay1", "af1", 0);
//   a->AddEdge("ay1", "x", 0); //clashes for the first input of x
//   a->AddEdge("ay2", "af1", 1);
//   a->AddEdge("af1", "af2", 1);
//
//     // add nodes
//   b->AddNode(make_shared<sinSumMod>(2), "bf2");
//   b->AddNode(make_shared<sinSumMod>(2), "bf1");
//   b->AddNode(sourceNode, "x");
//   b->AddNode(make_shared<squareMod>(2), "by1");
//   b->AddNode(make_shared<squareMod>(2), "by2");
//
//   // add connectivity
//   b->AddEdge("x", "bf2", 0);
//   b->AddEdge("by1", "bf1", 0);
//   b->AddEdge("by2", "x", 0); //clashes for the first input of x
//   b->AddEdge("by2", "bf1", 1);
//   b->AddEdge("bf1", "bf2", 1);
//
//   ASSERT_DEATH(ModGraph::FormUnion(a,b), "(result->ModelGraph[*e_result]->GetDim() != b->ModelGraph[e]->GetDim())*");
// }
