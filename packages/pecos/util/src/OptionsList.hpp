/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_UTIL_OPTIONS_LIST_HPP
#define PECOS_UTIL_OPTIONS_LIST_HPP

#include <string>
#include <iostream>
#include <map>
#include <boost/any.hpp>
#include <vector>
#include "Teuchos_SerialDenseVector.hpp"

namespace Pecos {
namespace util {

typedef std::map<std::string, boost::any> opts_map;

// declare functions inline to avoid duplicate symbols warning when linking
inline bool equal_any(const boost::any& lhs, const boost::any& rhs);
inline std::ostream& operator<<(std::ostream& os, const boost::any& item);

template<typename T>
std::ostream& print_vector(const std::vector<T> &vec, std::ostream& os){
  os << "[";
  for (size_t i=0; i<vec.size(); ++i){
    if (i>0) os << ", ";
    os << vec[i];
  }
  os << "]";
  return os;
}


class OptionsList {
  
private:
  opts_map map_;
  
public:
  
  OptionsList(){};

  OptionsList(const OptionsList& opts){copy(opts);};
  
  virtual ~OptionsList(){clear();};

  template<typename T>
  T get(const std::string& name) const{
    opts_map::const_iterator iter = map_.find(name);
    if (iter==map_.end()){
      std::string msg = "Item "+name+" not found in options";
      throw(std::runtime_error(msg));
    }
    try{
      return boost::any_cast<T>(iter->second);
    }catch (boost::bad_any_cast &e){
      std::string msg = e.what() + std::string(" in OptionsList::get(") +
	name+")";
      throw(std::runtime_error(msg));
    }
  };

  void get_keys(std::vector<std::string>& keys) const{
    keys.resize(size());
    opts_map::const_iterator iter;
    int i=0;
    for (iter=map_.begin(), i=0; iter!=map_.end(); ++iter, ++i)
      keys[i]=iter->first;
  }

  bool exists(const std::string &name) const{
    return (map_.find(name)!=map_.end());
  }

  template<typename T>
  bool isType(const std::string &name) const{
    opts_map::const_iterator iter=map_.find(name);
    if (iter==map_.end())
      return false;
    if (iter->second.type()==typeid(T))
      return true;
    else 
      return false;
  }

  size_t size() const{return map_.size();}

  void erase(const std::string &name){
    map_.erase(name);
  }

  const boost::any * get_entry(const std::string &name) const{
    opts_map::const_iterator iter = map_.find(name);
    if (iter!=map_.end())
      return &iter->second;
    else{
      std::string msg = "Item not found in options";
      throw(std::runtime_error(msg));
    }
  }

  template<typename T>
  T get(const std::string& name, T default_item) const{
    opts_map::const_iterator iter = map_.find(name);
    if (iter==map_.end()){
      return default_item;
    }
    try{
      return boost::any_cast<T>(iter->second);
    }catch (boost::bad_any_cast &e){
      throw(std::runtime_error(e.what()));
    }
  };

  template<typename T>
  void set(const std::string& name, const T &item){
    map_[name]=item;
  };

  void set(const std::string& name, char const * item){
    map_[name]=std::string(item);
  };

  void print_keys() const{
    this->print_keys(std::cout);
    std::cout << std::flush;
  }

  std::ostream& print_keys(std::ostream& os) const{
    opts_map::const_iterator iter;
    os << "[";
    for (iter=map_.begin(); iter!=map_.end(); ++iter){
      if (iter!=map_.begin()) os << ", ";
      os << iter->first;
    }
    os << "]";
    return os;
  }

  void print() const{
    this->print(std::cout);
    std::cout << std::flush;
  }

  std::ostream& print(std::ostream& os) const{
    opts_map::const_iterator iter;
    os << "{";
    for (iter=map_.begin(); iter!=map_.end(); ++iter){
      if (iter!=map_.begin()) os << ", ";
      os << iter->first << " : " << iter->second;
    }
    os << "}";
    return os;
  };

  void clear(){
    map_.clear();
  }

  void copy(const OptionsList& opts){
    clear();
    opts_map::const_iterator iter;
    for (iter=opts.map_.begin(); iter!=opts.map_.end(); ++iter){
      set(iter->first,iter->second);
    }
  }

  friend inline std::ostream& operator<<(std::ostream& os, const OptionsList& opts){
    return opts.print(os);
  }
  
  friend bool operator==(const OptionsList& list1, const OptionsList& list2 ){
    // Check that the two lists are the same length:
    if (list1.size()!=list2.size()) {
      return false;
    }
    
    opts_map::const_iterator iter1, iter2;
    for (iter1=list1.map_.begin(), iter2=list2.map_.begin();
         iter1!=list1.map_.end() && iter2!=list2.map_.end();
         ++iter1, ++iter2){
      if (iter1->first!=iter2->first)
        return false;
      //if (iter1->second!=iter2->second)//does not work with boost-1.49
      //if (!(iter1->second==iter2->second))
      if (!equal_any(iter1->second,iter2->second))
        return false;
    }
    return true;
  };

  friend bool operator!=(const OptionsList& list1, const OptionsList& list2 ){
    return !(list1==list2);
  }
};

bool equal_any(const boost::any& lhs, const boost::any& rhs){
  if (lhs.type()!=rhs.type())
    return false;
  else if (lhs.type()==typeid(int))
    return boost::any_cast<int>(lhs)==boost::any_cast<int>(rhs);
  else if (lhs.type()==typeid(double))
    return boost::any_cast<double>(lhs)==boost::any_cast<double>(rhs);
  else if (lhs.type()==typeid(std::string))
    return boost::any_cast<std::string>(lhs)==boost::any_cast<std::string>(rhs);
  else if (lhs.type()==typeid(OptionsList))
    return boost::any_cast<OptionsList>(lhs)==boost::any_cast<OptionsList>(rhs);
  else if (lhs.type()==typeid(Teuchos::SerialDenseVector<int,int>))
    return boost::any_cast< Teuchos::SerialDenseVector<int,int> >(lhs)==boost::any_cast< Teuchos::SerialDenseVector<int,int> >(rhs);
  else if (lhs.type()==typeid(Teuchos::SerialDenseVector<int,double>))
    return boost::any_cast< Teuchos::SerialDenseVector<int,double> >(lhs)==boost::any_cast< Teuchos::SerialDenseVector<int,double> >(rhs);
  else if (lhs.type()==typeid(std::vector<OptionsList>)){
    return boost::any_cast< std::vector<OptionsList> >(lhs)==boost::any_cast< std::vector<OptionsList> >(rhs);
  }else{
    std::string msg = "Comparion of any not implemented for specified type";
    throw(std::runtime_error(msg));
  }
  return false;
}

std::ostream& operator<<(std::ostream& os, const boost::any& item){
  if (item.type()==typeid(int))
    os << boost::any_cast<int>(item);
  else if (item.type()==typeid(double))
    os << boost::any_cast<double>(item);
  else if (item.type()==typeid(std::string))
    os << boost::any_cast<std::string>(item);
  else if (item.type()==typeid(OptionsList))
    os << boost::any_cast<OptionsList>(item);
  else if (item.type()==typeid(Teuchos::SerialDenseVector<int,int>))
    os << boost::any_cast< Teuchos::SerialDenseVector<int,int> >(item);
  else if (item.type()==typeid(Teuchos::SerialDenseVector<int,double>))
    os << boost::any_cast< Teuchos::SerialDenseVector<int,double> >(item);
  else if (item.type()==typeid(Teuchos::SerialDenseMatrix<int,int>))
    os << boost::any_cast< Teuchos::SerialDenseMatrix<int,int> >(item);
  else if (item.type()==typeid(Teuchos::SerialDenseMatrix<int,double>))
    os << boost::any_cast< Teuchos::SerialDenseMatrix<int,double> >(item);
  else if (item.type()==typeid(std::vector<OptionsList>)){
    std::vector<OptionsList> vec=
      boost::any_cast< std::vector<OptionsList> >(item);
    print_vector(vec,os);
  }
  else{
    std::string msg = "Print of any not implemented for specified type";
    throw(std::runtime_error(msg));
  }
  return os;
}


}  // namespace util
}  // namespace Pecos

#endif  // include guard
