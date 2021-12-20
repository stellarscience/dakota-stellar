#include "MUQ/Modeling/ODEData.h"

#include "boost/none.hpp"

using namespace muq::Modeling;

// construct basic ode data
ODEData::ODEData(std::shared_ptr<ModPiece> const& rhs,
                 ref_vector<Eigen::VectorXd> const& refinputs,
                 bool const autonomous,
                 int const wrtIn,
                 Eigen::VectorXd const& actionVecIn) : rhs(rhs),
                                                     autonomous(autonomous),
                                                     wrtIn(wrtIn),
                                                     actionVec(actionVecIn),
                                                     isAction(actionVecIn.size()>0)
{
  inputs.reserve(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs.push_back(refinputs[i]); }
}

// construct with root function
ODEData::ODEData(std::shared_ptr<ModPiece> const& rhs,
                 std::shared_ptr<ModPiece> const& root,
                 ref_vector<Eigen::VectorXd> const& refinputs,
                 bool const autonomous,
                 int const wrtIn,
                 Eigen::VectorXd const& actionVecIn) : rhs(rhs),
                                                       root(root),
                                                       autonomous(autonomous),
                                                       wrtIn(wrtIn),
                                                       actionVec(actionVecIn),
                                                       isAction(actionVecIn.size()>0)
{
  inputs.reserve(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs.push_back(refinputs[i]); }
}

void ODEData::UpdateInputs(Eigen::Ref<const Eigen::VectorXd> const& newState, double const newTime)
{

  state = newState;

  if( autonomous ) {
    inputs.at(0) = std::cref(state);
  }else{
    time = Eigen::VectorXd::Constant(1, newTime);
    inputs.at(0) = std::cref(time);
    inputs.at(1) = std::cref(state);
  }

}
