//------------------------------------------------------------------------
// P.J. Williams
// Sandia National Laboratories
// pwillia@sandia.gov
// Last modified 11/16/1999
//------------------------------------------------------------------------

#include "NLP.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

// Constructors 
NLP::NLP(): ptr_(0){;}
NLP::NLP(NLPBase* base): ptr_(base){;}

void NLP::setX(const int i, const double& x) 
{
   ptr_->setX(i,x);
}

void NLP::setX(const SerialDenseVector<int,double>& x) 
{
   ptr_->setX(x);
}

void NLP::setF(const double& fx) 
{
   ptr_->setF(fx);
}

void NLP::setIsExpensive(const int e) 
{
   ptr_->setIsExpensive(e);
}

void NLP::setFcnAccrcy(const int i, const double& accrcy) 
{
   ptr_->setFcnAccrcy(i,accrcy);
}

void NLP::setFcnAccrcy(const SerialDenseVector<int,double>& accrcy) 
{
   ptr_->setFcnAccrcy(accrcy);
}

int NLP::getDim() const
{
   int result = ptr_ -> getDim();
   return result;
}

int NLP::getFevals() const
{
   int result = ptr_ -> getFevals();
   return result;
}

int NLP::getIsExpensive() const
{
   int result = ptr_ -> getIsExpensive();
   return result;
}

double NLP::getF() const
{
   double result = ptr_ -> getF();
   return result;
}

SerialDenseVector<int,double> NLP::getFcnAccrcy() const
{
  SerialDenseVector<int,double> result(ptr_->getFcnAccrcy().length());
   result = ptr_ -> getFcnAccrcy();
   return result;
}

SerialDenseVector<int,double> NLP::getXc() const
{
  SerialDenseVector<int,double> result(ptr_->getXc().length()); 
   result = ptr_ -> getXc();
   return result;
}

double NLP::getFcnTime() const
{
   double result = ptr_ -> getFcnTime();
   return result;
}

int NLP::getNumOfCons() const
{
   int result = ptr_ -> getNumOfCons();
   return result;
}

int NLP::getNumOfNLCons() const
{
   int result = ptr_ -> getNumOfNLCons();
   return result;
}


bool NLP::hasConstraints() 
{
   bool result = ptr_ -> hasConstraints();
   return result;
}

void NLP::printConstraints() 
{
   ptr_ -> printConstraints();
}

void NLP::setDebug() 
{
   ptr_ -> setDebug();
}

bool NLP::getDebug() const 
{
   bool result = ptr_ -> getDebug();
   return result;
}

void NLP::reset()
{
   ptr_->reset();
}

void NLP::initFcn()
{
   ptr_->initFcn();
}

void NLP::eval()
{
   ptr_->eval();
}

double NLP::evalF()
{
   double result = ptr_->evalF();
   return result;
}

double NLP::evalF(const SerialDenseVector<int,double>& x)
{
   double result = ptr_->evalF(x);
   return result;
}

SerialDenseVector<int,double> NLP::evalG()
{
  SerialDenseVector<int,double> result(ptr_->evalG());
  //  SerialDenseVector<int,double> result(ptr_->getDim());
  //   result = ptr_->evalG();
   return result;
}

SerialDenseVector<int,double> NLP::evalG(const SerialDenseVector<int,double>& x)
{
  SerialDenseVector<int,double> result(ptr_->evalG(x));
  //  SerialDenseVector<int,double> result(ptr_->getDim());
  //   result = ptr_->evalG(x);
   return result;
}

SerialSymDenseMatrix<int,double> NLP::evalH()
{
  SerialSymDenseMatrix<int,double> result(ptr_->evalH());
  //  SerialSymDenseMatrix<int,double> result(ptr_->getDim());
  //     result = ptr_->evalH();
   return result;
}

SerialSymDenseMatrix<int,double> NLP::evalH(SerialDenseVector<int,double>& x)
{
  SerialSymDenseMatrix<int,double> result(ptr_->evalH(x));
  //  SerialSymDenseMatrix<int,double> result(ptr_->getDim());
  //     result = ptr_->evalH(x);
   return result;
}

SerialDenseVector<int,double> NLP::evalCF(const SerialDenseVector<int,double>& x)
{
  SerialDenseVector<int,double> result(ptr_->evalCF(x));
  //  SerialDenseVector<int,double> result(ptr_->getNumOfNLCons());
  //     result = ptr_->evalCF(x);
   return result;
}


SerialDenseMatrix<int,double> NLP::evalCG(const SerialDenseVector<int,double>& x)
{
  SerialDenseMatrix<int,double> result(ptr_->evalCG(x));
  //  SerialDenseMatrix<int,double> result(ptr_->getDim(),ptr_->getNumOfNLCons());
  //     result = ptr_->evalCG(x);
   return result;
}

SerialSymDenseMatrix<int,double> NLP::evalCH(SerialDenseVector<int,double>& x)
{
  SerialSymDenseMatrix<int,double> result(ptr_->evalCH(x));
  //  SerialSymDenseMatrix<int,double> result(ptr_->getNumOfNLCons());
  //     result = ptr_->evalCH(x);
   return result;
}

OptppArray<SerialSymDenseMatrix<int,double> > NLP::evalCH(SerialDenseVector<int,double>& x, int darg)
{
  OptppArray<SerialSymDenseMatrix<int,double> > result = ptr_->evalCH(x,darg);
   return result;
}

void NLP::evalC(const SerialDenseVector<int,double>& x)
{
  ptr_->evalC(x);
}

void NLP::printState(const char* s)
{
   ptr_->printState(s);
}

void NLP::fPrintState(std::ostream *nlpout, const char* s)
{
   ptr_->fPrintState(nlpout,s);
}

} // namespace OPTPP
