/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "OptionsList.hpp"
#define BOOST_TEST_MODULE options
#include <boost/test/included/unit_test.hpp>

using Pecos::util::OptionsList;

void test_set(){
  // test basic set and get
  
  OptionsList opts;
  opts.set("a","ab");
  BOOST_CHECK( opts.get<std::string>("a")=="ab");

  //cannot use opts.set("b",1); because of use & in second argument of of set(name, T &)
  int one=1;
  opts.set("b",one);
  BOOST_CHECK(opts.get<int>("b")==1);

  double two=2.;
  opts.set("b",two);
  BOOST_CHECK( opts.get<double>("b")==2.0);

  OptionsList sub_opts;
  sub_opts.set("a","abc");
  opts.set("subopts",sub_opts);
  BOOST_CHECK( opts.get<OptionsList>("subopts")==sub_opts);


  int zero=0;
  int three=3;
  std::vector<OptionsList> list_of_opts(3);
  list_of_opts[0].set("a",zero);
  list_of_opts[1].set("a",one);
  list_of_opts[2].set("a",three);
  opts.set("list_of_opts",list_of_opts);
  std::vector<OptionsList> result=opts.get< std::vector<OptionsList> >("list_of_opts");
  BOOST_CHECK( result.size()==list_of_opts.size());
  for (int i=0; i<result.size(); ++i){
    BOOST_CHECK(result[i].get<int>("a")==list_of_opts[i].get<int>("a"));
    std::cout << result[i].get<int>("a") << std::endl;
  }
  std::cout << opts << std::endl;
}

void test_get_missing_option(){
  // Check that error is thrown if a value is requested that does not exist

  OptionsList opts;
  bool error_caught=false;
  try{
    opts.get<double>("b");
  }catch(std::runtime_error &e){
    error_caught = true;
  }
  BOOST_CHECK(error_caught);
}

void test_get_incorrect_type(){
  // Check that error is thrown if a value exists but wrong type is requested

  int two=2;
  OptionsList opts;
  opts.set("b",two);
  bool error_caught=false;
  try{
    opts.get<double>("b");
  }catch(std::runtime_error &e){
    error_caught = true;
  }
  BOOST_CHECK(error_caught);
}


void test_set_to_default(){
  // Check default is applied correctly

  OptionsList opts;
  
  BOOST_CHECK(opts.get<double>("c",2.0)==2.0);

  BOOST_CHECK(opts.get<int>("c",2)==2);

  BOOST_CHECK(opts.get<std::string>("c","ab")=="ab");
}

void test_equality(){
  int one=1;
  double two=2.0;
  OptionsList opts1;
  opts1.set("a","ab");
  opts1.set("b",one);
  opts1.set("c",two);
  OptionsList sub_opts1;
  sub_opts1.set("a","abc");
  opts1.set("subopts",sub_opts1);

  OptionsList opts2;
  opts2.set("a","ab");
  opts2.set("b",one);
  opts2.set("c",two);
  OptionsList sub_opts2;
  sub_opts2.set("a","abc");
  opts2.set("subopts",sub_opts2);

  BOOST_CHECK(opts1==opts2);  

  opts2.set("a","cd");
  BOOST_CHECK(opts1!=opts2);  
}

void test_copy(){
  int one=1.;
  double two=2.0;
  OptionsList opts;
  opts.set("a","ab");
  opts.set("b",one);
  opts.set("c",two);
  OptionsList sub_opts;
  sub_opts.set("a","abc");
  opts.set("subopts",sub_opts);

  OptionsList opts_copy(opts);
  BOOST_CHECK(opts_copy==opts);
  std::cout << opts_copy << std::endl;
}


BOOST_AUTO_TEST_CASE( test_main ){
  test_set();
  test_get_missing_option();
  test_get_incorrect_type();
  test_set_to_default();
  test_equality();
  test_copy();
  
  int run_result = 0;
  BOOST_CHECK( run_result == 0 || run_result == boost::exit_success );
}

 //cannot make the following work
 // opts_variants v(1);
  // int x = boost::apply_visitor(GetOption(), v);
  // I get
  // error: no viable conversion from 'typename
  //      GetOption::result_type' (aka 'boost::variant
  // but this works
  // opts_variants s = boost::apply_visitor(GetOption(), v);
