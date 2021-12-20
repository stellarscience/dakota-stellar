#ifndef MULTIINDEXSERIALIZER_H
#define MULTIINDEXSERIALIZER_H


#include "cereal/cereal.hpp"

#include "parcer/Eigen.h"


namespace cereal {

  template<class Archive>
  void save(Archive & ar, MultiIndex const& obj)
  {
    ar(obj.GetVector());
  }

  template<class Archive>
  void load(Archive & ar, MultiIndex & obj)
  {

    Eigen::RowVectorXi vector;
    ar(vector);
    obj = muq::Utilities::MultiIndex(vector);
  }

}

#endif
