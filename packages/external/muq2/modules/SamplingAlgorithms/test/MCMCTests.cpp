#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/MALAProposal.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/SamplingAlgorithms/DRKernel.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

TEST(MCMC, MHKernel_ThinScheduler) {
  const unsigned int N = 1e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.BurnIn", 0);
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.ThinIncrement", 10); // How often to save the sample
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 3.0); // the variance of the isotropic MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu + RandomGenerator::GetNormal(2);

  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"),problem);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  EXPECT_EQ(int(std::floor(pt.get<double>("MyMCMC.NumSamples")/pt.get<double>("MyMCMC.ThinIncrement"))), samps->size());

  //boost::any anyMean = samps.Mean();
  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 8e-2);
  EXPECT_NEAR(mu(1), mean(1), 8e-2);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1.2e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1.2e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1.2e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1.2e-1);

}

TEST(MCMC, MHKernel_MHProposal) {
  const unsigned int N = 1e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu;

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"),problem);

  // Make sure the kernel and proposal are correct
  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  ASSERT_TRUE(kernelBase);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  ASSERT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<MHProposal> proposalMH = std::dynamic_pointer_cast<MHProposal>(proposalBase);
  ASSERT_TRUE(proposalMH);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  EXPECT_EQ(pt.get<int>("MyMCMC.NumSamples"), samps->size());

  //boost::any anyMean = samps.Mean();
  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}


TEST(MCMC, MHKernel_RepeatedRuns) {
  const unsigned int batchSize = 1e3;
  const unsigned int numBatches = 10;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", batchSize); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu;

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"),problem);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);
  for(unsigned int batchInd=1; batchInd<numBatches; ++batchInd){
    mcmc->AddNumSamps(batchSize);
    mcmc->Run();
  }

  EXPECT_EQ(batchSize*numBatches, samps->size());

  //boost::any anyMean = samps.Mean();
  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}


TEST(MCMC, MHKernel_MixtureProposal) {
  const unsigned int N = 1e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N); // number of Monte Carlo samples
  pt.put("PrintLevel",0);
  pt.put("KernelList", "Kernel1"); // the transition kernel
  pt.put("Kernel1.Method","MHKernel");
  pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("Kernel1.MyProposal.Method", "MixtureProposal");
  pt.put("Kernel1.MyProposal.Components", "SubProp1,SubProp2");
  pt.put("Kernel1.MyProposal.Weights","0.4,0.6");
  pt.put("Kernel1.MyProposal.SubProp1.Method","MHProposal");
  pt.put("Kernel1.MyProposal.SubProp1.ProposalVariance", 0.5); // the variance of the isotropic MH proposal
  pt.put("Kernel1.MyProposal.SubProp2.Method","MHProposal");
  pt.put("Kernel1.MyProposal.SubProp2.ProposalVariance", 1.0);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu;

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt,problem);
  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  EXPECT_EQ(pt.get<int>("NumSamples"), samps->size());

  //boost::any anyMean = samps.Mean();
  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}


TEST(MCMC, MHKernel_MALAProposal) {
  const unsigned int N = 1e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N); // number of Monte Carlo samples
  pt.put("PrintLevel",0);
  pt.put("KernelList", "Kernel1"); // the transition kernel
  pt.put("Kernel1.Method","MHKernel");
  pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("Kernel1.MyProposal.Method", "MALAProposal");
  pt.put("Kernel1.MyProposal.StepSize", 1.0); // the variance of the isotropic MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu;

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt,problem);

  // Make sure the kernel and proposal are correct
  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  ASSERT_TRUE(kernelBase);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  ASSERT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<MALAProposal> proposalMALA = std::dynamic_pointer_cast<MALAProposal>(proposalBase);
  ASSERT_TRUE(proposalMALA);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  EXPECT_GE(N, dist->GetNumCalls("Evaluate"));
  EXPECT_GE(N, dist->GetNumCalls("Gradient"));
  EXPECT_EQ(pt.get<int>("NumSamples"), samps->size());
  
  //boost::any anyMean = samps.Mean();
  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}

TEST(MCMC, MetropolisInGibbs_IsoGauss) {

  const unsigned int N = 1e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1,Kernel2"); // the transition kernel

  // MH Kernel for first block
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 0.5);

  pt.put("MyMCMC.Kernel2.Method","MHKernel");
  pt.put("MyMCMC.Kernel2.Proposal", "MyProposal");
  pt.put("MyMCMC.Kernel2.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel2.MyProposal.ProposalVariance", 0.5);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu1 = 1.0*Eigen::VectorXd::Ones(2);
  const Eigen::VectorXd mu2 = 1.5*Eigen::VectorXd::Ones(2);
  auto dist1 = std::make_shared<Gaussian>(mu1); // standard normal Gaussian
  auto dist2 = std::make_shared<Gaussian>(mu2); // standard normal Gaussian

  auto graph = std::make_shared<WorkGraph>();
  graph->AddNode(dist1->AsDensity(), "Gaussian1");
  graph->AddNode(dist2->AsDensity(), "Gaussian2");

  graph->AddNode(std::make_shared<DensityProduct>(2), "ProductDensity");
  graph->AddEdge("Gaussian1",0,"ProductDensity",0);
  graph->AddEdge("Gaussian2",0,"ProductDensity",1);

  auto dens = graph->CreateModPiece("ProductDensity");

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dens);


  // starting point
  std::vector<Eigen::VectorXd> start(2);
  start.at(0) = mu1;
  start.at(1) = mu2;

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"),problem);

  // Make sure the kernel and proposal are correct
  EXPECT_EQ(2,mcmc->Kernels().size());

  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  EXPECT_TRUE(kernelMH);

  kernelBase = mcmc->Kernels().at(1);
  kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  EXPECT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<MHProposal> proposalMH = std::dynamic_pointer_cast<MHProposal>(proposalBase);
  EXPECT_TRUE(proposalMH);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  Eigen::VectorXd mean = samps->Mean();


  EXPECT_NEAR(mu1(0), mean(0), 1e-1);
  EXPECT_NEAR(mu1(1), mean(1), 1e-1);
  EXPECT_NEAR(mu2(0), mean(2), 1e-1);
  EXPECT_NEAR(mu2(1), mean(3), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  for(int j=0; j<cov.cols(); ++j){
    for(int i=0; i<cov.rows(); ++i){
      if(i!=j)
        EXPECT_NEAR(0.0, cov(i,j), 2e-1);
      else
        EXPECT_NEAR(1.0, cov(i,j), 2e-1);
    }
  }
}


TEST(MCMC, MHKernel_AMProposal) {
  const unsigned int N = 5e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "AMProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.InitialVariance", 1.0); // the variance of the isotropic MH proposal
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptSteps",200);
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptStart",2000);
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptScale",1.0);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem);

  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  EXPECT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<AMProposal> proposalAM = std::dynamic_pointer_cast<AMProposal>(proposalBase);
  EXPECT_TRUE(proposalAM);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  Eigen::VectorXd mean = samps->Mean();

  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}

TEST(MCMC, MHKernel_DRAM) {
  const unsigned int N = 5e4;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N); // number of Monte Carlo samples
  pt.put("PrintLevel",0);
  pt.put("KernelList", "Kernel1"); // the transition kernel
  pt.put("Kernel1.Method","DRKernel");
  pt.put("Kernel1.NumStages",3);
  pt.put("Kernel1.Scale",2.0);
  pt.put("Kernel1.ScaleType","Power");
  pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("Kernel1.MyProposal.Method", "AMProposal");
  pt.put("Kernel1.MyProposal.InitialVariance", 1.0); // the variance of the isotropic MH proposal
  pt.put("Kernel1.MyProposal.AdaptSteps",20);
  pt.put("Kernel1.MyProposal.AdaptStart",2000);
  pt.put("Kernel1.MyProposal.AdaptScale",1.0);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  const double rho = 0.8;
  Eigen::MatrixXd tgtCov(2,2);
  tgtCov << 1.0,                rho*std::sqrt(2.0),
            rho*std::sqrt(2.0), 2.0;

  auto dist = std::make_shared<Gaussian>(mu,tgtCov)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt, problem);

  std::shared_ptr<DRKernel> kernel = std::dynamic_pointer_cast<DRKernel>(mcmc->Kernels().at(0));
  std::shared_ptr<AMProposal> proposal = std::dynamic_pointer_cast<AMProposal>(kernel->Proposals().at(0));

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  Eigen::VectorXd mean = samps->Mean();

  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd propCov = proposal->ProposalCovariance();
  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(tgtCov(0,0), cov(0,0), 1e-1);
  EXPECT_NEAR(tgtCov(0,1), cov(0,1), 1e-1);
  EXPECT_NEAR(tgtCov(1,0), cov(1,0), 1e-1);
  EXPECT_NEAR(tgtCov(1,1), cov(1,1), 1e-1);
}
