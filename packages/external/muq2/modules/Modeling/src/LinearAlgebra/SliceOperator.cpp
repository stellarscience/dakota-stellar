#include "MUQ/Modeling/LinearAlgebra/SliceOperator.h"

using namespace muq::Modeling;


SliceOperator::SliceOperator(unsigned int vecSize,
                             int startIndIn,
                             int endIndIn,
                             int skipIn) : LinearOperator(ComputeRows(vecSize,startIndIn,endIndIn,skipIn),
                                                          vecSize),
                                           startInd(startIndIn),
                                                    endInd(endIndIn),
                                                    skip(skipIn)
{
  if(startIndIn<0){
    startInd = vecSize+startIndIn;
  }else{
    startInd = startIndIn;
  }

  if(endIndIn<0){
    endInd = vecSize+endIndIn;
  }else{
    endInd = endIndIn;
  }

  assert(skipIn!=0);
  if(skip<0){
    assert(startInd>endInd);
  }else{
    assert(startInd<endInd);
  }
};

unsigned int SliceOperator::ComputeRows(unsigned int vecSize, int startIndIn, int endIndIn, int skipIn)
{
  int startInd, endInd;

  if(startIndIn<0){
    startInd = vecSize+startIndIn;
  }else{
    startInd = startIndIn;
  }

  if(endIndIn<0){
    endInd = vecSize+endIndIn;
  }else{
    endInd = endIndIn;
  }

  return std::ceil(double(endInd-startInd)/double(skipIn));
}


Eigen::MatrixXd SliceOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  Eigen::MatrixXd output(rows(),x.cols());
  unsigned int outInd = 0;
  const int multiplier = (skip<0) ? -1.0 : 1.0;
  for(unsigned int inInd=startInd; multiplier*inInd<multiplier*endInd; inInd+=skip){
    output.row(outInd) = x.row(inInd);
    outInd++;
  }

  return output;
}

/** Apply the transpose of the linear operator to a vector. */
Eigen::MatrixXd SliceOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(cols(), x.cols());
  unsigned int outInd = 0;
  const int multiplier = (skip<0) ? -1.0 : 1.0;
  for(unsigned int inInd=startInd; multiplier*inInd<multiplier*endInd; inInd+=skip){
    output.row(inInd) = x.row(outInd);
    outInd++;
  }

  return output;
}

Eigen::MatrixXd SliceOperator::GetMatrix()
{
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows(), cols());
  unsigned int outInd = 0;
  for(unsigned int inInd=startInd; inInd!=endInd; inInd+=skip){
    output(outInd,inInd) = 1.0;
    outInd++;
  }

  return output;
}
