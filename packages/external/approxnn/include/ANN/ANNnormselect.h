// Dakota extension to override macros and allow norm selection at
// run-time (with a performance penalty).  Thanks to Dan Turner for
// initial implementation.

/* 
Example of using this in client code to override to l_inf norm:
Conditionally compile ANN.h with -D ANN_NORM_SELECT (see ANN.h)

approxnn::normSelector::instance().method(approxnn::LINF_NORM);
ANN_kdtree(...);

// optional, but recommended if another intervening context might
// call into the library
normSelector::instance().reset();
*/

#ifndef ANN_NORMSELECT_H
#define ANN_NORMSELECT_H

namespace approxnn {

enum {L2_NORM, LINF_NORM};

/// \class normSelector
/// \brief singleton class to keep track of which norm method
/// should be used for NORM macro
class normSelector{
public:
  /// return an instance of the singleton
  static normSelector &instance(){
    static normSelector instance_;
    return instance_;
  }

  void reset() {
    method_ = L2_NORM;
  }
  /// set the method to use
  /// \param method an integer value that represents the norm flag
  void method(const int method){
    method_ = method;
  }

  /// returns the current method
  int method(){
    return method_;
  }

  /// runtime-switched implementation of ANN_POW
  template<typename T>
  T pow(const T& v) {
    switch(normSelector::instance().method()) {
    case L2_NORM:
      return v*v;
      break;
    case LINF_NORM:
      return std::fabs(v);
      break;
    }
  }

  /// runtime-switched implementation of ANN_POW
  template<typename T>
  T root(const T& x) {
    switch(normSelector::instance().method()) {
    case L2_NORM:
      return std::sqrt(x);
      break;
    case LINF_NORM:
      return x;
      break;
    }
  }

  /// runtime-switched implementation of ANN_SUM
  template<typename T>
  T sum(const T& x, const T& y) {
    switch(normSelector::instance().method()) {
    case L2_NORM:
      return x + y;
      break;
    case LINF_NORM:
      return (x > y) ? x : y;
      break;
    }
  }

  /// runtime-switched implementation of ANN_DIFF
  template<typename T>
  T diff(const T& x, const T& y) {
    switch(normSelector::instance().method()) {
    case L2_NORM:
      return (y-x);
      break;
    case LINF_NORM:
      return y;
      break;
    }
  }

private:
  /// constructor
  normSelector():method_(L2_NORM){};
  /// copy constructor
  normSelector(normSelector const&);
  /// asignment operator
  void operator=(normSelector const &);
  /// method to use
  int method_;
};


}  // namespace approxnn

#endif
