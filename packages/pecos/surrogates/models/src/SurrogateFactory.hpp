/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_SURROGATES_SURROGATE_FACTORY_HPP
#define PECOS_SURROGATES_SURROGATE_FACTORY_HPP

namespace Pecos {
namespace surrogates {
  
enum ApproxType {PCE,GP};

/**
   \class SurrogateFactory
   \brief The object used to build surrogates.
*/
class SurrogateFactory {
  public:

  /** \brief Build a surrogate using the specifications in opts
   *
   * \param[in] opts Options used to build the surrogate
   *
   */  
  build(const OptionsList opts){
    ApproxType approx_type = opts.get<ApproxType>("approximation type");
    switch (approx_type){
      case PCE : {
        PCEFactory(opts, approx);
        break;
      }
      case GP : { 
        GPFactory(opts, approx);
        break;
      }
      default : {
        throw(std::runtime_error("Incorrect approximation type"));
      }
    }
  }

  /// Constructor
  SurrogateFactory(){};
  
  /// Destructor
  virtual ~SurrogateFactory(){};
}

}  // namespace surrogates
}  // namespace Pecos

#endif  // include guard
