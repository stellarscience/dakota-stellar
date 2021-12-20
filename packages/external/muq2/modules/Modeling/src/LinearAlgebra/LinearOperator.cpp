#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"


using namespace muq::Modeling;

LinearOperator::LinearOperator(int rowsIn, int colsIn, int numInputCols) : muq::Modeling::ModPiece(colsIn*numInputCols*Eigen::VectorXi::Ones(1),
                                                                                                   rowsIn*numInputCols*Eigen::VectorXi::Ones(1)),
                                                         ncols(colsIn),
                                                         nrows(rowsIn)
{
};

void LinearOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y)
{
    assert(y.cols()==x.cols());
    y = Apply(x);
};


void LinearOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y)
{
    assert(y.cols()==x.cols());
    y = ApplyTranspose(x);
};

Eigen::MatrixXd LinearOperator::GetMatrix()
{

    Eigen::MatrixXd output(nrows, ncols);
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Identity(ncols, ncols);

    for(int i=0; i<ncols; ++i)
        output.col(i) = Apply(rhs.col(i));

    return output;
}

void LinearOperator::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input)
{
  outputs.resize(1);
  outputs.at(0) = Apply(input.at(0).get()).col(0);
}

void LinearOperator::GradientImpl(unsigned int                const  outputDimWrt,
                          unsigned int                const  inputDimWrt,
                          muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                          Eigen::VectorXd             const& sensitivity)
{
  gradient = ApplyTranspose(sensitivity);
}

void LinearOperator::JacobianImpl(unsigned int                const  outputDimWrt,
                          unsigned int                const  inputDimWrt,
                          muq::Modeling::ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = GetMatrix();
}

void LinearOperator::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                               unsigned int                const  inputDimWrt,
                               muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                               Eigen::VectorXd             const& vec)
{
  jacobianAction = Apply(vec);
}

void LinearOperator::ApplyHessianImpl(unsigned int const outWrt,
                                      unsigned int inWrt1,
                                      unsigned int inWrt2,
                                      muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                                      Eigen::VectorXd             const& sens,
                                      Eigen::VectorXd             const& vec)
{

  // If inWrt1==inWrt2, then the Hessian is zero
  if(inWrt1==inWrt2){
    hessAction = Eigen::VectorXd::Zero(cols());
  }else{

    /* Since there is only one input, inWrt1!=inWrt2 will only be possible if inWrt2
       is wrt to the sensitivity.  Since the gradient is A^Ts, the Jacobian of this
       wrt to the sensitivity is A^T.
    */
    hessAction = ApplyTranspose(vec);
  }
}
