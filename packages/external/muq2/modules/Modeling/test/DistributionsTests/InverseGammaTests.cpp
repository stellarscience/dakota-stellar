#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/InverseGamma.h"

using namespace muq::Modeling;


TEST(InverseGammaDistributionTests, EvaluateLogDensity) {

  const double alpha = 1.5;
  const double beta = 0.5;

  auto dist = std::make_shared<InverseGamma>(1.5, 0.5);

  double logScale = alpha * std::log(beta) - std::lgamma(alpha);

  // evalute the log-denstion
  Eigen::VectorXd x(1);
  x << 0.5;

  double logdens = dist->LogDensity(x);
  EXPECT_DOUBLE_EQ(logScale + (-alpha-1.0)*std::log(x(0)) - beta/x(0), logdens);

  x << 1.5;
  logdens = dist->LogDensity(x);
  EXPECT_DOUBLE_EQ(logScale + (-alpha-1.0)*std::log(x(0)) - beta/x(0), logdens);

  x << -1.0;
  logdens = dist->LogDensity(x);
  EXPECT_DOUBLE_EQ(-1.0*std::numeric_limits<double>::infinity(), logdens);
}

TEST(InverseGammaDistributionTests, EvaluateLogDensity_Multivariate) {

  Eigen::VectorXd alpha(2);
  alpha << 1.5, 2.0;
  Eigen::VectorXd beta(2);
  beta << 0.5, 0.5;

  auto dist = std::make_shared<InverseGamma>(alpha, beta);

  double logScale = alpha(0) * std::log(beta(0)) - std::lgamma(alpha(0));
  logScale += alpha(1)*std::log(beta(1)) - std::lgamma(alpha(1));

  // evalute the log-denstion
  Eigen::VectorXd x(2);
  x << 0.5, 0.6;

  double logdens = dist->LogDensity(x);
  EXPECT_DOUBLE_EQ(logScale + (-alpha(0)-1.0)*std::log(x(0)) - beta(0)/x(0)+ (-alpha(1)-1.0)*std::log(x(1)) - beta(1)/x(1), logdens);
}

TEST(InverseGammaDistributionTests, Sample) {

  const double alpha = 2.5;
  const double beta = 1.0;

  auto dist = std::make_shared<InverseGamma>(alpha, beta);

  unsigned int numSamps  = 1e6;
  double mu = 0.0;
  for(int i=0; i<numSamps; ++i){
    double samp  = dist->Sample()(0);
    EXPECT_TRUE(samp>0);
    mu += (1.0/double(numSamps))*samp;
  }

  EXPECT_NEAR(beta/(alpha-1.0), mu, 5e-3);
}


TEST(InverseGammaDistributionTests, Sample_Multivariate) {

  Eigen::VectorXd alpha(2);
  alpha << 1.5, 2.0;
  Eigen::VectorXd beta(2);
  beta << 0.5, 0.5;

  auto dist = std::make_shared<InverseGamma>(alpha, beta);

  unsigned int numSamps  = 1e6;
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(2);
  for(int i=0; i<numSamps; ++i){
    Eigen::VectorXd samp  = dist->Sample();
    EXPECT_TRUE(samp(0)>0);
    EXPECT_TRUE(samp(1)>0);
    mu += (1.0/double(numSamps))*samp;
  }

  EXPECT_NEAR(beta(0)/(alpha(0)-1.0), mu(0), 5e-2);
  EXPECT_NEAR(beta(1)/(alpha(1)-1.0), mu(1), 5e-2);
}
