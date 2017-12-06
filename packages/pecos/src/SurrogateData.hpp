/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SURROGATE_DATA_HPP
#define SURROGATE_DATA_HPP

#include "pecos_data_types.hpp"
#include <boost/math/special_functions/fpclassify.hpp> //for boostmath::isfinite


namespace Pecos {

/// The representation of a SurrogateDataVars instance.  This representation,
/// or body, may be shared by multiple SurrogateDataVars handle instances.

/** The SurrogateDataVars/SurrogateDataVarsRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class SurrogateDataVarsRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class SurrogateDataVars;

private:

  //
  //- Heading: Private member functions
  //

  /// constructor
  SurrogateDataVarsRep(const RealVector& c_vars, const IntVector& di_vars,
		       const RealVector& dr_vars, short mode);
  /// alternate constructor (data sizing only)
  SurrogateDataVarsRep(size_t num_c_vars, size_t num_di_vars,
		       size_t num_dr_vars);
  /// destructor
  ~SurrogateDataVarsRep();

  //
  //- Heading: Private data members
  //

  RealVector continuousVars;   ///< continuous variables
  IntVector  discreteIntVars;  ///< discrete integer variables
  RealVector discreteRealVars; ///< discrete real variables

  int referenceCount;        ///< number of handle objects sharing sdvRep
};


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(const RealVector& c_vars, const IntVector& di_vars,
		     const RealVector& dr_vars, short mode): referenceCount(1)
{
  // Note: provided a way to query DataAccess mode for c_vars, could make
  // greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (mode == DEEP_COPY) {         // enforce deep vector copy
    if (c_vars.length())  copy_data( c_vars, continuousVars);
    if (di_vars.length()) copy_data(di_vars, discreteIntVars);
    if (dr_vars.length()) copy_data(dr_vars, discreteRealVars);
  }
  else if (mode == SHALLOW_COPY) { // enforce shallow vector copy
    if (c_vars.length())
      continuousVars
	= RealVector(Teuchos::View,  c_vars.values(),  c_vars.length());
    if (di_vars.length())
      discreteIntVars
	= IntVector(Teuchos::View,  di_vars.values(), di_vars.length());
    if (dr_vars.length())
      discreteRealVars
	= RealVector(Teuchos::View, dr_vars.values(), dr_vars.length());
  }
  else {                           // default: assume existing Copy/View state
    if (c_vars.length())    continuousVars = c_vars;
    if (di_vars.length())  discreteIntVars = di_vars;
    if (dr_vars.length()) discreteRealVars = dr_vars;
  }
}


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars):
  referenceCount(1)
{
  continuousVars.sizeUninitialized(num_c_vars);
  discreteIntVars.sizeUninitialized(num_di_vars);
  discreteRealVars.sizeUninitialized(num_dr_vars);
}


inline SurrogateDataVarsRep::~SurrogateDataVarsRep()
{ }


/// Container class encapsulating basic parameter data for defining a
/// "truth" data point.

/** A set of these input data points is contained in SurrogateData and
    provides the data to build the approximation.  A handle-body idiom
    is used to avoid excessive data copying overhead. */

class SurrogateDataVars
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  SurrogateDataVars();
  /// standard constructor
  SurrogateDataVars(const RealVector& c_vars, const IntVector& di_vars,
		    const RealVector& dr_vars, short mode = DEFAULT_COPY);
  /// alternate constructor (data sizing only)
  SurrogateDataVars(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars);
  /// copy constructor
  SurrogateDataVars(const SurrogateDataVars& sdv);
  /// destructor
  ~SurrogateDataVars();

  /// assignment operator
  SurrogateDataVars& operator=(const SurrogateDataVars& sdv);
  // equality operator
  //bool operator==(const SurrogateDataVars& sdv) const;

  //
  //- Heading: member functions
  //

  /// return deep copy of SurrogateDataVars instance
  SurrogateDataVars copy() const;

  /// set i^{th} entry within continuousVars
  void continuous_variable(Real c_var, size_t i);
  /// set continuousVars
  void continuous_variables(const RealVector& c_vars,
			    short mode = DEFAULT_COPY);
  /// get continuousVars
  const RealVector& continuous_variables() const;
  /// get view of continuousVars for updating in place
  RealVector continuous_variables_view();

  /// set i^{th} entry within discreteIntVars
  void discrete_int_variable(int di_var, size_t i);
  /// set discreteIntVars
  void discrete_int_variables(const IntVector& di_vars,
			      short mode = DEFAULT_COPY);
  /// get discreteIntVars
  const IntVector& discrete_int_variables() const;
  /// get view of discreteIntVars for updating in place
  IntVector discrete_int_variables_view();

  /// set i^{th} entry within discreteRealVars
  void discrete_real_variable(Real dr_var, size_t i);
  /// set discreteRealVars
  void discrete_real_variables(const RealVector& dr_vars,
			       short mode = DEFAULT_COPY);
  /// get discreteRealVars
  const RealVector& discrete_real_variables() const;
  /// get view of discreteRealVars for updating in place
  RealVector discrete_real_variables_view();

  /// function to check sdvRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  SurrogateDataVarsRep* sdvRep;
};


inline SurrogateDataVars::SurrogateDataVars(): sdvRep(NULL)
{ }


inline SurrogateDataVars::
SurrogateDataVars(const RealVector& c_vars, const IntVector& di_vars,
		  const RealVector& dr_vars, short mode):
  sdvRep(new SurrogateDataVarsRep(c_vars, di_vars, dr_vars, mode))
{ }


inline SurrogateDataVars::
SurrogateDataVars(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars):
  sdvRep(new SurrogateDataVarsRep(num_c_vars, num_di_vars, num_dr_vars))
{ }


inline SurrogateDataVars::SurrogateDataVars(const SurrogateDataVars& sdv)
{
  // Increment new (no old to decrement)
  sdvRep = sdv.sdvRep;
  if (sdvRep) // Check for an assignment of NULL
    ++sdvRep->referenceCount;
}


inline SurrogateDataVars::~SurrogateDataVars()
{
  if (sdvRep) { // Check for NULL
    --sdvRep->referenceCount; // decrement
    if (sdvRep->referenceCount == 0)
      delete sdvRep;
  }
}


inline SurrogateDataVars& SurrogateDataVars::
operator=(const SurrogateDataVars& sdv)
{
  if (sdvRep != sdv.sdvRep) { // prevent re-assignment of same rep
    // Decrement old
    if (sdvRep) // Check for NULL
      if ( --sdvRep->referenceCount == 0 ) 
	delete sdvRep;
    // Increment new
    sdvRep = sdv.sdvRep;
    if (sdvRep) // Check for an assignment of NULL
      ++sdvRep->referenceCount;
  }
  // else if assigning same rep, then leave referenceCount as is

  return *this;
}


//inline bool SurrogateDataVars::operator==(const SurrogateDataVars& sdv) const
//{
//  return (sdvRep->continuousVars   == sdv.sdvRep->continuousVars  &&
//          sdvRep->discreteIntVars  == sdv.sdvRep->discreteIntVars &&
//          sdvRep->discreteRealVars == sdv.sdvRep->discreteRealVars) ?
//    true : false;
//}


/// deep copy of SurrogateDataVars instance
inline SurrogateDataVars SurrogateDataVars::copy() const
{
  SurrogateDataVars sdv(sdvRep->continuousVars,   sdvRep->discreteIntVars,
			sdvRep->discreteRealVars, DEEP_COPY);
  return sdv;
}


inline void SurrogateDataVars::continuous_variable(Real c_var, size_t i)
{ sdvRep->continuousVars[i] = c_var; }


inline void SurrogateDataVars::
continuous_variables(const RealVector& c_vars, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(c_vars, sdvRep->continuousVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    sdvRep->continuousVars
      = RealVector(Teuchos::View, c_vars.values(), c_vars.length());
  else                           // default: assume existing Copy/View state
    sdvRep->continuousVars = c_vars;
}


inline const RealVector& SurrogateDataVars::continuous_variables() const
{ return sdvRep->continuousVars; }


inline RealVector SurrogateDataVars::continuous_variables_view()
{
  return RealVector(Teuchos::View, sdvRep->continuousVars.values(),
		    sdvRep->continuousVars.length());
}


inline void SurrogateDataVars::discrete_int_variable(int di_var, size_t i)
{ sdvRep->discreteIntVars[i] = di_var; }


inline void SurrogateDataVars::
discrete_int_variables(const IntVector& di_vars, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(di_vars, sdvRep->discreteIntVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    sdvRep->discreteIntVars
      = IntVector(Teuchos::View, di_vars.values(), di_vars.length());
  else                           // default: assume existing Copy/View state
    sdvRep->discreteIntVars = di_vars;
}


inline const IntVector& SurrogateDataVars::discrete_int_variables() const
{ return sdvRep->discreteIntVars; }


inline IntVector SurrogateDataVars::discrete_int_variables_view()
{
  return IntVector(Teuchos::View, sdvRep->discreteIntVars.values(),
		   sdvRep->discreteIntVars.length());
}


inline void SurrogateDataVars::discrete_real_variable(Real dr_var, size_t i)
{ sdvRep->discreteRealVars[i] = dr_var; }


inline void SurrogateDataVars::
discrete_real_variables(const RealVector& dr_vars, short mode)
{
  if (mode == DEEP_COPY)         // enforce deep vector copy
    copy_data(dr_vars, sdvRep->discreteRealVars);
  else if (mode == SHALLOW_COPY) // enforce shallow vector copy
    sdvRep->discreteRealVars
      = RealVector(Teuchos::View, dr_vars.values(), dr_vars.length());
  else                           // default: assume existing Copy/View state
    sdvRep->discreteRealVars = dr_vars;
}


inline const RealVector& SurrogateDataVars::discrete_real_variables() const
{ return sdvRep->discreteRealVars; }


inline RealVector SurrogateDataVars::discrete_real_variables_view()
{
  return RealVector(Teuchos::View, sdvRep->discreteRealVars.values(),
		    sdvRep->discreteRealVars.length());
}


inline bool SurrogateDataVars::is_null() const
{ return (sdvRep) ? false : true; }


/// The representation of a surrogate data response.  This representation,
/// or body, may be shared by multiple SurrogateDataResp handle instances.

/** The SurrogateDataResp/SurrogateDataRespRep pairs utilize a
    handle-body idiom (Coplien, Advanced C++). */

class SurrogateDataRespRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class SurrogateDataResp;

private:

  //
  //- Heading: Private member functions
  //

  /// constructor
  SurrogateDataRespRep(Real fn_val, const RealVector& fn_grad,
		       const RealSymMatrix& fn_hess, short bits, short mode);
  /// alternate constructor (data sizing only)
  SurrogateDataRespRep(short data_order, size_t num_vars);
  /// destructor
  ~SurrogateDataRespRep();

  //
  //- Heading: Private data members
  //

  short           activeBits; ///< active data bits: 1 (fn), 2 (grad), 4 (hess)
  Real            responseFn; ///< truth response function value
  RealVector    responseGrad; ///< truth response function gradient
  RealSymMatrix responseHess; ///< truth response function Hessian
  int         referenceCount; ///< number of handle objects sharing sdrRep

};


inline SurrogateDataRespRep::
SurrogateDataRespRep(Real fn_val, const RealVector& fn_grad,
		     const RealSymMatrix& fn_hess, short bits, short mode):
  responseFn(fn_val), // always deep copy for scalars
  activeBits(bits), referenceCount(1)
{
  // Note: provided a way to query incoming grad/hess DataAccess modes,
  // could make greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (mode == DEEP_COPY) {          // enforce vector/matrix deep copy
    if (activeBits & 2) copy_data(fn_grad, responseGrad);
    if (activeBits & 4) copy_data(fn_hess, responseHess);
  }
  else if (mode == SHALLOW_COPY) {  // enforce vector/matrix shallow copy
    if (activeBits & 2)
      responseGrad = RealVector(Teuchos::View, fn_grad.values(),
				fn_grad.length());
    if (activeBits & 4)
      responseHess = RealSymMatrix(Teuchos::View, fn_hess, fn_hess.numRows());
  }
  else {                            // default: assume existing Copy/View state
    if (activeBits & 2) responseGrad = fn_grad;
    if (activeBits & 4) responseHess = fn_hess;
  }
}


inline SurrogateDataRespRep::
SurrogateDataRespRep(short data_order, size_t num_vars):
  activeBits(data_order), referenceCount(1)
{
  if (data_order & 2)
    responseGrad.sizeUninitialized(num_vars);
  if (data_order & 4)
    responseHess.shapeUninitialized(num_vars);
}


inline SurrogateDataRespRep::~SurrogateDataRespRep()
{ }


/// Container class encapsulating basic parameter and response data
/// for defining a "truth" data point.

/** A list of these data points is contained in Approximation instances
    (e.g., Dakota::Approximation::currentPoints) and provides the data
    to build the approximation.  A handle-body idiom is used to avoid
    excessive data copying overhead. */

class SurrogateDataResp
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  /// default constructor
  SurrogateDataResp();
  /// standard constructor
  SurrogateDataResp(Real fn_val, const RealVector& fn_grad,
		    const RealSymMatrix& fn_hess, short bits,
		    short mode = DEFAULT_COPY);
  /// alternate constructor (data sizing only)
  SurrogateDataResp(short data_order, size_t num_vars);
  /// copy constructor
  SurrogateDataResp(const SurrogateDataResp& sdr);
  /// destructor
  ~SurrogateDataResp();

  /// assignment operator
  SurrogateDataResp& operator=(const SurrogateDataResp& sdr);
  // equality operator
  //bool operator==(const SurrogateDataResp& sdr) const;

  //
  //- Heading: member functions
  //

  /// return deep copy of SurrogateDataResp instance
  SurrogateDataResp copy() const;

  /// set activeBits
  void active_bits(short val);
  /// get activeBits
  short active_bits() const;

  /// set responseFn
  void response_function(Real fn);
  /// get responseFn
  Real response_function() const;
  /// get "view" of responseFn for updating in place
  Real& response_function_view();

  /// set i^{th} entry within responseGrad
  void response_gradient(Real grad_i, size_t i);
  /// set responseGrad
  void response_gradient(const RealVector& grad, short mode = DEFAULT_COPY);
  /// get responseGrad
  const RealVector& response_gradient() const;
  /// get view of responseGrad for updating in place
  RealVector response_gradient_view();

  /// set i-j^{th} entry within responseHess
  void response_hessian(Real hess_ij, size_t i, size_t j);
  /// set responseHess
  void response_hessian(const RealSymMatrix& hess, short mode = DEFAULT_COPY);
  /// get responseHess
  const RealSymMatrix& response_hessian() const;
  /// get view of responseHess for updating in place
  RealSymMatrix response_hessian_view();

  /// function to check sdrRep (does this handle contain a body)
  bool is_null() const;
  /// output response function, gradient, and Hessian data
  void write(std::ostream& s) const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  SurrogateDataRespRep* sdrRep;
};


inline SurrogateDataResp::SurrogateDataResp(): sdrRep(NULL)
{ }


inline SurrogateDataResp::
SurrogateDataResp(Real fn_val, const RealVector& fn_grad,
		  const RealSymMatrix& fn_hess, short bits, short mode):
  sdrRep(new SurrogateDataRespRep(fn_val, fn_grad, fn_hess, bits, mode))
{ }


inline SurrogateDataResp::
SurrogateDataResp(short data_order, size_t num_vars):
  sdrRep(new SurrogateDataRespRep(data_order, num_vars))
{ }


inline SurrogateDataResp::SurrogateDataResp(const SurrogateDataResp& sdr)
{
  // Increment new (no old to decrement)
  sdrRep = sdr.sdrRep;
  if (sdrRep) // Check for an assignment of NULL
    ++sdrRep->referenceCount;
}


inline SurrogateDataResp::~SurrogateDataResp()
{
  if (sdrRep) { // Check for NULL
    --sdrRep->referenceCount; // decrement
    if (sdrRep->referenceCount == 0)
      delete sdrRep;
  }
}


inline SurrogateDataResp& SurrogateDataResp::
operator=(const SurrogateDataResp& sdr)
{
  if (sdrRep != sdr.sdrRep) { // prevent re-assignment of same rep
    // Decrement old
    if (sdrRep) // Check for NULL
      if ( --sdrRep->referenceCount == 0 ) 
	delete sdrRep;
    // Increment new
    sdrRep = sdr.sdrRep;
    if (sdrRep) // Check for an assignment of NULL
      ++sdrRep->referenceCount;
  }
  // else if assigning same rep, then leave referenceCount as is

  return *this;
}


//inline bool SurrogateDataResp::operator==(const SurrogateDataResp& sdr) const
//{
//  return ( sdrRep->responseFn   == sdr.sdrRep->responseFn   &&
//	     sdrRep->responseGrad == sdr.sdrRep->responseGrad &&
//	     sdrRep->responseHess == sdr.sdrRep->responseHess ) ? true : false;
//}


/// deep copy of SurrogateDataResp instance
inline SurrogateDataResp SurrogateDataResp::copy() const
{
  SurrogateDataResp sdr(sdrRep->responseFn,   sdrRep->responseGrad,
			sdrRep->responseHess, sdrRep->activeBits, DEEP_COPY);
  return sdr;
}


inline void SurrogateDataResp::active_bits(short val)
{ sdrRep->activeBits = val; }


inline short SurrogateDataResp::active_bits() const
{ return sdrRep->activeBits; }


inline void SurrogateDataResp::response_function(Real fn)
{ sdrRep->responseFn = fn; }


inline Real SurrogateDataResp::response_function() const
{ return sdrRep->responseFn; }


inline Real& SurrogateDataResp::response_function_view()
{ return sdrRep->responseFn; }


inline void SurrogateDataResp::response_gradient(Real grad_i, size_t i)
{ sdrRep->responseGrad[i] = grad_i; }


inline void SurrogateDataResp::
response_gradient(const RealVector& grad, short mode)
{
  if (mode == DEEP_COPY)          // enforce vector deep copy
    copy_data(grad, sdrRep->responseGrad);
  else if (mode == SHALLOW_COPY)  // enforce vector shallow copy
    sdrRep->responseGrad
      = RealVector(Teuchos::View, grad.values(), grad.length());
  else                            // default: assume existing Copy/View state
    sdrRep->responseGrad = grad;
}


inline const RealVector& SurrogateDataResp::response_gradient() const
{ return sdrRep->responseGrad; }


inline RealVector SurrogateDataResp::response_gradient_view()
{
  return RealVector(Teuchos::View, sdrRep->responseGrad.values(),
		    sdrRep->responseGrad.length());
}


inline void SurrogateDataResp::
response_hessian(Real hess_ij, size_t i, size_t j)
{ sdrRep->responseHess(i,j) = hess_ij; }


inline void SurrogateDataResp::
response_hessian(const RealSymMatrix& hess, short mode)
{
  if (mode == DEEP_COPY)          // enforce matrix deep copy
    copy_data(hess, sdrRep->responseHess);
  else if (mode == SHALLOW_COPY)  // enforce matrix shallow copy
    sdrRep->responseHess = RealSymMatrix(Teuchos::View, hess, hess.numRows());
  else                            // default: assume existing Copy/View state
    sdrRep->responseHess = hess;
}


inline const RealSymMatrix& SurrogateDataResp::response_hessian() const
{ return sdrRep->responseHess; }


inline RealSymMatrix SurrogateDataResp::response_hessian_view()
{
  return RealSymMatrix(Teuchos::View, sdrRep->responseHess,
		       sdrRep->responseHess.numRows());
}


inline bool SurrogateDataResp::is_null() const
{ return (sdrRep) ? false : true; }


inline void SurrogateDataResp::write(std::ostream& s) const
{
  if (sdrRep->activeBits & 1)
    s << "function value    =  " << std::setw(WRITE_PRECISION+7)
      << sdrRep->responseFn << '\n';
  if (sdrRep->activeBits & 2) {
    s << "function gradient =\n";
    write_data_trans(s, sdrRep->responseGrad, true, true, true);
  }
  if (sdrRep->activeBits & 4) {
    s << "function Hessian  =\n";
    write_data(s, sdrRep->responseHess, true, true, true);
  }
}


/// std::ostream insertion operator for SurrogateDataResp
inline std::ostream& operator<<(std::ostream& s, const SurrogateDataResp& sdr)
{ sdr.write(s); return s; }


/// Representation of management class for surrogate data defined from
/// input variable data and output response data.

/** Response sets are generally unique for each approximated response
    function, but variable sets are often (but not always) replicated.
    Thus the design accommodates the case where inputs and/or outputs
    involve either shallow or deep copies. */

class SurrogateDataRep
{
  //
  //- Heading: Friends
  //

  /// the handle class can access attributes of the body class directly
  friend class SurrogateData;

public:

private:

  //
  //- Heading: Constructors and destructor
  //

  SurrogateDataRep();  ///< constructor
  ~SurrogateDataRep(); ///< destructor

  //
  //- Heading: Private data members
  //

  /// a special variables sample (often at the center of the
  /// approximation region) for which exact matching is enforced
  /// (e.g., using equality-constrained least squares regression).
  SurrogateDataVars anchorVars;
  /// a set of variables samples used to build the approximation.
  /// These sample points are fit approximately (e.g., using least
  /// squares regression); exact matching is not enforced.
  SDVArray varsData;
  /// sets of variables samples that have been popped off varsData but
  /// which are available for future restoration.  The granularity of
  /// this 2D array corresponds to multiple trial sets, each
  /// contributing an SDVArray.
  SDV2DArray poppedVarsTrials;
  /// a set of variables samples that have been stored for future
  /// restoration.  The granularity of this 2D array corresponds to a
  /// wholesale caching of the current varsData state, where each of
  /// multiple cachings contributes an SDVArray.
  SDV2DArray storedVarsData;

  /// a special response sample (often at the center of the
  /// approximation region) for which exact matching is enforced
  /// (e.g., using equality-constrained least squares regression).
  SurrogateDataResp anchorResp;
  /// a set of response samples used to build the approximation.  These
  /// sample points are fit approximately (e.g., using least squares
  /// regression); exact matching is not enforced.
  SDRArray respData;
  /// sets of response samples that have been popped off respData but
  /// which are available for future restoration.  The granularity of
  /// this 2D array corresponds to multiple trial sets, each
  /// contributing a SDRArray.
  SDR2DArray poppedRespTrials;
  /// a set of response samples that have been stored for future
  /// restoration.  The granularity of this 2D array corresponds to a
  /// wholesale caching of the current respData state, where each of
  /// multiple cachings contributes an SDRArray.
  SDR2DArray storedRespData;

  /// failed anchor data bits; defined in sample_checks() and used for
  /// fault tolerance in regression() and expectation()
  short failedAnchorData;
  /// map from failed respData indices to failed data bits; defined
  /// in sample_checks() and used for fault tolerance
  SizetShortMap failedRespData;

  /// a stack managing the number of points previously added by calls
  /// to append() that can be removed by calls to pop()
  SizetArray popCountStack;

  /// number of handle objects sharing sdRep
  int referenceCount;
};


inline SurrogateDataRep::SurrogateDataRep():
  failedAnchorData(0), referenceCount(1)
{ }


inline SurrogateDataRep::~SurrogateDataRep()
{ }


/// Management class for surrogate data defined from input variable
/// data and output response data.

/** Response sets are generally unique for each approximated response
    function, but variable sets are often (but not always) replicated.
    Thus the design accommodates the case where inputs and/or outputs
    involve either shallow or deep copies. */

class SurrogateData
{
public:

  //
  //- Heading: Constructors, destructor, and operators
  //

  SurrogateData();                        ///< default constructor
  SurrogateData(const SurrogateData& sd); ///< copy constructor
  ~SurrogateData();                       ///< destructor

  /// assignment operator
  SurrogateData& operator=(const SurrogateData& sdv);

  //
  //- Heading: Member functions
  //

  /// deep copy of SurrogateData instance with options for shallow copies
  /// of the SurrogateData{Vars,Resp} objects
  SurrogateData copy(short sdv_mode = DEEP_COPY,
		     short sdr_mode = DEEP_COPY) const;

  /// set anchor{Vars,Resp}
  void anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);
  /// set {vars,resp}Data
  void data_points(const SDVArray& sdv_array, const SDRArray& sdr_array);

  /// set anchorVars
  void anchor_variables(const SurrogateDataVars& sdv);
  /// get anchorVars
  const SurrogateDataVars& anchor_variables() const;
  /// set anchorResp
  void anchor_response(const SurrogateDataResp& sdr);
  /// get anchorResp
  const SurrogateDataResp& anchor_response() const;

  /// set varsData
  void variables_data(const SDVArray& sdv_array);
  /// get varsData
  const SDVArray& variables_data() const;
  /// set respData
  void response_data(const SDRArray& sdr_array);
  /// get respData
  const SDRArray& response_data() const;

  /// get anchorVars.continuous_variables()
  const RealVector& anchor_continuous_variables() const;
  /// get anchorVars.discrete_int_variables()
  const IntVector& anchor_discrete_int_variables() const;
  /// get anchorVars.discrete_real_variables()
  const RealVector& anchor_discrete_real_variables() const;
  /// get anchorResp.active_bits()
  short anchor_active_bits() const;
  /// get anchorResp.response_function()
  Real anchor_function() const;
  /// get anchorResp.response_gradient()
  const RealVector& anchor_gradient() const;
  /// get anchorResp.response_hessian()
  const RealSymMatrix& anchor_hessian() const;

  /// get varsData[i].continuous_variables()
  const RealVector& continuous_variables(size_t i) const;
  /// get varsData[i].discrete_int_variables()
  const IntVector& discrete_int_variables(size_t i) const;
  /// get varsData[i].discrete_real_variables()
  const RealVector& discrete_real_variables(size_t i) const;
  /// get respData[i].active_bits()
  short response_active_bits(size_t i) const;

  /// set respData[i].response_function()
  void response_function(Real fn_val, size_t i);
  /// get respData[i].response_function()
  Real response_function(size_t i) const;
  /// set respData[i].response_gradient()
  void response_gradient(const RealVector& fn_grad, size_t i);
  /// get respData[i].response_gradient()
  const RealVector& response_gradient(size_t i) const;
  /// set respData[i].response_hessian()
  void response_hessian(const RealSymMatrix& fn_hess, size_t i);
  /// get respData[i].response_hessian()
  const RealSymMatrix& response_hessian(size_t i) const;

  /// push sdv onto end of varsData
  void push_back(const SurrogateDataVars& sdv);
  /// push sdr onto end of respData
  void push_back(const SurrogateDataResp& sdr);
  /// push {sdv,sdr} onto ends of {vars,resp}Data
  void push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);

  /// remove num_pop_pts entries from ends of {vars,resp}Data
  void pop(bool save_data = true);
  /// return a previously popped data set (identified by index) to the
  /// ends of {vars,resp}Data
  void push(size_t index, bool erase_popped = true);

  /// append count to popCountStack
  void pop_count(size_t count);
  /// return popCountStack.back()
  size_t pop_count() const;

  /// move all entries from {vars,resp}Data to stored{Vars,Resp}Data
  /// (default is push_back)
  void store(size_t index = _NPOS);
  /// return an entry set from stored{Vars,Resp}Data to {vars,resp}Data
  void restore(size_t index = _NPOS);
  /// remove an entry from stored{Vars,Resp}Data (default is pop_back)
  void remove_stored(size_t index = _NPOS);
  /// swap stored and active data sets (variables and response)
  void swap(size_t index);

  /// query presence of anchor{Vars,Resp}
  bool anchor() const;
  /// return size of {vars,resp}Data arrays (neglecting anchor point)
  size_t points() const;

  /// return total number of available data components
  size_t response_size() const;
  /// return number of failed data components
  size_t failed_response_size() const;
  /// return net number of active data components (total minus failed)
  size_t active_response_size() const;

  /// return number of 1D arrays within stored{Vars,Resp}Data 2D arrays
  size_t stored_sets() const;
  /// return number of 1D arrays within popped{Vars,Resp}Trials 2D arrays
  size_t popped_sets() const;
  /// return number of derivative variables as indicated by size of
  /// gradient/Hessian arrays
  size_t num_derivative_variables() const;

  /// convenience function used by data_checks() for anchorResp and respData
  void response_check(const SurrogateDataResp& sdr, short& failed_data);
  /// screen data sets for samples with Inf/Nan that should be excluded;
  /// defines failedAnchorData and failedRespData
  void data_checks();
  /// return failedAnchorData
  short failed_anchor_data() const;
  /// return failedRespData
  const SizetShortMap& failed_response_data() const;

  /// clear anchor{Vars,Resp}
  void clear_anchor();
  /// clear {vars,resp}Data
  void clear_data();
  /// clear popped{Vars,Resp}Trials
  void clear_popped();
  /// clear stored{Vars,Resp}Data
  void clear_stored();

  /// return sdRep
  SurrogateDataRep* data_rep() const;

private:

  //
  //- Heading: Member functions
  //

  /// set failedAnchorData
  void failed_anchor_data(short fail_anchor);
  /// set failedRespData
  void failed_response_data(const SizetShortMap& fail_resp);

  /// get storedVarsData
  const SDV2DArray& stored_variables_data() const;
  /// set storedVarsData
  void stored_variables_data(const SDV2DArray& stored_vars);
  /// get storedRespData
  const SDR2DArray& stored_response_data() const;
  /// set storedRespData
  void stored_response_data(const SDR2DArray& stored_resp);

  /// get poppedVarsTrials
  const SDV2DArray& popped_variables_trials() const;
  /// set poppedVarsTrials
  void popped_variables_trials(const SDV2DArray& popped_vars);
  /// get poppedRespTrials
  const SDR2DArray& popped_response_trials() const;
  /// set poppedRespTrials
  void popped_response_trials(const SDR2DArray& popped_resp);

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  SurrogateDataRep* sdRep;
};


inline SurrogateData::SurrogateData(): sdRep(new SurrogateDataRep())
{ }


inline SurrogateData::SurrogateData(const SurrogateData& sd)
{
  // Increment new (no old to decrement)
  sdRep = sd.sdRep;
  if (sdRep) // Check for an assignment of NULL
    ++sdRep->referenceCount;
}


inline SurrogateData::~SurrogateData()
{
  if (sdRep) { // Check for NULL
    --sdRep->referenceCount; // decrement
    if (sdRep->referenceCount == 0)
      delete sdRep;
  }
}


inline SurrogateData& SurrogateData::operator=(const SurrogateData& sd)
{
  if (sdRep != sd.sdRep) { // prevent re-assignment of same rep
    // Decrement old
    if (sdRep) // Check for NULL
      if ( --sdRep->referenceCount == 0 ) 
	delete sdRep;
    // Increment new
    sdRep = sd.sdRep;
    if (sdRep) // Check for an assignment of NULL
      ++sdRep->referenceCount;
  }
  // else if assigning same rep, then leave referenceCount as is

  return *this;
}


inline void SurrogateData::
anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{ sdRep->anchorVars = sdv; sdRep->anchorResp = sdr; }


inline void SurrogateData::
data_points(const SDVArray& sdv_array, const SDRArray& sdr_array)
{ sdRep->varsData = sdv_array; sdRep->respData = sdr_array; }


inline void SurrogateData::anchor_variables(const SurrogateDataVars& sdv)
{ sdRep->anchorVars = sdv; }


inline const SurrogateDataVars& SurrogateData::anchor_variables() const
{ return sdRep->anchorVars; }


inline void SurrogateData::anchor_response(const SurrogateDataResp& sdr)
{ sdRep->anchorResp = sdr; }


inline const SurrogateDataResp& SurrogateData::anchor_response() const
{ return sdRep->anchorResp; }


inline void SurrogateData::variables_data(const SDVArray& sdv_array)
{ sdRep->varsData = sdv_array; }


inline const SDVArray& SurrogateData::variables_data() const
{ return sdRep->varsData; }


inline void SurrogateData::response_data(const SDRArray& sdr_array)
{ sdRep->respData = sdr_array; }


inline const SDRArray& SurrogateData::response_data() const
{ return sdRep->respData; }


inline const RealVector& SurrogateData::anchor_continuous_variables() const
{ return sdRep->anchorVars.continuous_variables(); }


inline const IntVector& SurrogateData::anchor_discrete_int_variables() const
{ return sdRep->anchorVars.discrete_int_variables(); }


inline const RealVector& SurrogateData::anchor_discrete_real_variables() const
{ return sdRep->anchorVars.discrete_real_variables(); }


inline const RealVector& SurrogateData::continuous_variables(size_t i) const
{ return sdRep->varsData[i].continuous_variables(); }


inline const IntVector& SurrogateData::discrete_int_variables(size_t i) const
{ return sdRep->varsData[i].discrete_int_variables(); }


inline const RealVector& SurrogateData::discrete_real_variables(size_t i) const
{ return sdRep->varsData[i].discrete_real_variables(); }


inline short SurrogateData::anchor_active_bits() const
{ return sdRep->anchorResp.active_bits(); }


inline short SurrogateData::response_active_bits(size_t i) const
{ return sdRep->respData[i].active_bits(); }


inline Real SurrogateData::anchor_function() const
{ return sdRep->anchorResp.response_function(); }


inline const RealVector& SurrogateData::anchor_gradient() const
{ return sdRep->anchorResp.response_gradient(); }


inline const RealSymMatrix& SurrogateData::anchor_hessian() const
{ return sdRep->anchorResp.response_hessian(); }


inline void SurrogateData::response_function(Real fn_val, size_t i)
{ sdRep->respData[i].response_function(fn_val); }


inline Real SurrogateData::response_function(size_t i) const
{ return sdRep->respData[i].response_function(); }


inline void SurrogateData::
response_gradient(const RealVector& fn_grad, size_t i)
{ sdRep->respData[i].response_gradient(fn_grad); }


inline const RealVector& SurrogateData::response_gradient(size_t i) const
{ return sdRep->respData[i].response_gradient(); }


inline void SurrogateData::
response_hessian(const RealSymMatrix& fn_hess, size_t i)
{ sdRep->respData[i].response_hessian(fn_hess); }


inline const RealSymMatrix& SurrogateData::response_hessian(size_t i) const
{ return sdRep->respData[i].response_hessian(); }


inline void SurrogateData::push_back(const SurrogateDataVars& sdv)
{ sdRep->varsData.push_back(sdv); }


inline void SurrogateData::push_back(const SurrogateDataResp& sdr)
{ sdRep->respData.push_back(sdr); }


inline void SurrogateData::
push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{ sdRep->varsData.push_back(sdv); sdRep->respData.push_back(sdr); }


inline void SurrogateData::pop(bool save_data)
{
  if (sdRep->popCountStack.empty()) {
    PCerr << "\nError: empty count stack in SurrogateData::pop()." << std::endl;
    abort_handler(-1);
  }
  size_t num_pop_pts = sdRep->popCountStack.back();
  if (num_pop_pts) {
    size_t data_size = points();
    if (data_size < num_pop_pts) {
      PCerr << "Error: pop count (" << num_pop_pts << ") exceeds data size ("
	    << data_size << ") in SurrogateData::pop(size_t)." << std::endl;
      abort_handler(-1);
    }
    if (save_data) {
      // append empty arrays and then update them in place
      SDVArray sdv_array; SDRArray sdr_array;
      sdRep->poppedVarsTrials.push_back(sdv_array);
      sdRep->poppedRespTrials.push_back(sdr_array);
      SDVArray& last_sdv_array = sdRep->poppedVarsTrials.back();
      SDRArray& last_sdr_array = sdRep->poppedRespTrials.back();
      /*
      // prevent underflow portability issue w/ compiler coercion of -num_pop
      SDVArray::difference_type reverse_adv_vars = num_pop_pts;
      SDRArray::difference_type reverse_adv_resp = num_pop_pts;
      SDVIter vit = varsData.end(); std::advance(vit, -reverse_adv_vars);
      SDRIter rit = respData.end(); std::advance(rit, -reverse_adv_resp);
      */
      last_sdv_array.insert(last_sdv_array.begin(), //vit,
			    sdRep->varsData.end() - num_pop_pts,
			    sdRep->varsData.end());
      last_sdr_array.insert(last_sdr_array.begin(), //rit,
			    sdRep->respData.end() - num_pop_pts,
			    sdRep->respData.end());
    }
    size_t new_size = data_size - num_pop_pts;
    sdRep->varsData.resize(new_size); sdRep->respData.resize(new_size);
  }
  sdRep->popCountStack.pop_back();
}


inline void SurrogateData::push(size_t index, bool erase_popped)
{
  SDV2DArray::iterator vit = sdRep->poppedVarsTrials.begin();
  SDR2DArray::iterator rit = sdRep->poppedRespTrials.begin();
  std::advance(vit, index); std::advance(rit, index);
  size_t num_pts = std::min(vit->size(), rit->size());
  sdRep->varsData.insert(sdRep->varsData.end(), vit->begin(), vit->end());
  sdRep->respData.insert(sdRep->respData.end(), rit->begin(), rit->end());

  if (erase_popped)
    { sdRep->poppedVarsTrials.erase(vit); sdRep->poppedRespTrials.erase(rit); }

  sdRep->popCountStack.push_back(num_pts);
}


inline void SurrogateData::pop_count(size_t count)
{ sdRep->popCountStack.push_back(count); }


inline size_t SurrogateData::pop_count() const
{ return (sdRep->popCountStack.empty()) ? _NPOS : sdRep->popCountStack.back(); }


inline void SurrogateData::store(size_t index)
{
  size_t stored_len = std::min(sdRep->storedVarsData.size(),
			       sdRep->storedRespData.size());
  if (index == _NPOS || index == stored_len) { // append
    sdRep->storedVarsData.push_back(sdRep->varsData); // shallow copies
    sdRep->storedRespData.push_back(sdRep->respData); // shallow copies
  }
  else if (index < stored_len) { // replace
    sdRep->storedVarsData[index] = sdRep->varsData; // shallow copies
    sdRep->storedRespData[index] = sdRep->respData; // shallow copies
  }
  else {
    PCerr << "Error: bad index (" << index
	  << ") passed in SurrogateData::store()" << std::endl;
    abort_handler(-1);
  }
  clear_data();
}


inline void SurrogateData::restore(size_t index)
{
  size_t stored_len = std::min(sdRep->storedVarsData.size(),
			       sdRep->storedRespData.size());
  if (index == _NPOS) {
    sdRep->varsData = sdRep->storedVarsData.back(); // shallow copies
    sdRep->respData = sdRep->storedRespData.back(); // shallow copies
  }
  else if (index < stored_len) { // replace
    sdRep->varsData = sdRep->storedVarsData[index]; // shallow copies
    sdRep->respData = sdRep->storedRespData[index]; // shallow copies
  }
  else {
    PCerr << "Error: bad index (" << index
	  << ") passed in SurrogateData::store()" << std::endl;
    abort_handler(-1);
  }
}


inline void SurrogateData::remove_stored(size_t index)
{
  size_t stored_len = std::min(sdRep->storedVarsData.size(),
			       sdRep->storedRespData.size());
  if (index == _NPOS || index == stored_len - 1) {
    sdRep->storedVarsData.pop_back(); // shallow copies
    sdRep->storedRespData.pop_back(); // shallow copies
  }
  else if (index < stored_len) { // replace
    SDV2DArray::iterator vit = sdRep->storedVarsData.begin();
    std::advance(vit, index); sdRep->storedVarsData.erase(vit);
    SDR2DArray::iterator rit = sdRep->storedRespData.begin();
    std::advance(rit, index); sdRep->storedRespData.erase(rit);
  }
  else {
    PCerr << "Error: bad index (" << index
	  << ") passed in SurrogateData::remove_stored()" << std::endl;
    abort_handler(-1);
  }
  clear_data();
}


inline void SurrogateData::swap(size_t index)
{
  if (index == _NPOS)
    return;
  else if (index >= sdRep->storedVarsData.size()) {
    PCerr << "Error: index out of range in SurrogateData::swap()" << std::endl;
    abort_handler(-1);
  }

  // swap stored and active using shallow copies
  SDVArray tmp_vars_data = sdRep->varsData;
  SDRArray tmp_resp_data = sdRep->respData;

  sdRep->varsData = sdRep->storedVarsData[index];
  sdRep->respData = sdRep->storedRespData[index];

  sdRep->storedVarsData[index] = tmp_vars_data;
  sdRep->storedRespData[index] = tmp_resp_data;
}


inline bool SurrogateData::anchor() const
{ return (!sdRep->anchorVars.is_null() && !sdRep->anchorResp.is_null()); }


inline size_t SurrogateData::points() const
{ return std::min(sdRep->varsData.size(), sdRep->respData.size()); }


inline size_t SurrogateData::response_size() const
{
  size_t i, data_size = 0, num_resp = sdRep->respData.size(), nh;
  short active_bits;
  if (!anchor_response().is_null()) {
    active_bits = anchor_active_bits();
    if (active_bits & 1) ++data_size;
    if (active_bits & 2) data_size += anchor_gradient().length();
    if (active_bits & 4) {
      nh = anchor_hessian().numRows();
      if (nh) data_size += nh * (nh + 1) / 2;
    }
  }
  for (i=0; i<num_resp; ++i) {
    active_bits = response_active_bits(i);
    if (active_bits & 1) ++data_size;
    if (active_bits & 2) data_size += response_gradient(i).length();
    if (active_bits & 4) {
      nh = response_hessian(i).numRows();
      if (nh) data_size += nh * (nh + 1) / 2;
    }
  }
  return data_size;
}


inline size_t SurrogateData::failed_response_size() const
{
  size_t failed_size = 0, nh; short fail_bits;
  if (!anchor_response().is_null()) {
    fail_bits = failed_anchor_data();
    if (fail_bits & 1) ++failed_size;
    if (fail_bits & 2) failed_size += anchor_gradient().length();
    if (fail_bits & 4) {
      nh = anchor_hessian().numRows();
      if (nh) failed_size += nh * (nh + 1) / 2;
    }
  }
  const SizetShortMap& failed_resp = failed_response_data();
  for (StShMCIter cit=failed_resp.begin(); cit!=failed_resp.end(); ++cit) {
    fail_bits = cit->second;
    if (fail_bits & 1) ++failed_size;
    if (fail_bits & 2) failed_size += response_gradient(cit->first).length();
    if (fail_bits & 4) {
      nh = response_hessian(cit->first).numRows();
      if (nh) failed_size += nh * (nh + 1) / 2;
    }
  }
  return failed_size;
}


inline size_t SurrogateData::active_response_size() const
{ return response_size() - failed_response_size(); }


inline size_t SurrogateData::stored_sets() const
{ return std::min(sdRep->storedVarsData.size(), sdRep->storedRespData.size()); }


inline size_t SurrogateData::popped_sets() const
{
  return std::min(sdRep->poppedVarsTrials.size(),
		  sdRep->poppedRespTrials.size());
}


inline size_t SurrogateData::num_derivative_variables() const
{
  return (sdRep->anchorResp.is_null()) ?
    sdRep->respData[0].response_gradient().length() :
    sdRep->anchorResp.response_gradient().length();
}


inline void SurrogateData::
response_check(const SurrogateDataResp& sdr, short& failed_data)
{
  // We take a conservative approach of rejecting all data of derivative
  // order greater than or equal to a detected failure:

  short resp_bits = sdr.active_bits();
  failed_data = 0;
  if (resp_bits & 1) {
    if (!boost::math::isfinite<Real>(sdr.response_function()))
      failed_data = resp_bits;       // all data for this & higher deriv orders
  }
  if ( (resp_bits & 2) && !failed_data ) {
    const RealVector& grad = sdr.response_gradient();
    size_t j, num_deriv_vars = grad.length();
    for (j=0; j<num_deriv_vars; ++j)
      if (!boost::math::isfinite<Real>(grad[j]))
	{ failed_data = (resp_bits & 6); break; } // this & higher deriv orders
  }
  if ( (resp_bits & 4) && !failed_data ) {
    const RealSymMatrix& hess = sdr.response_hessian();
    size_t j, k, num_deriv_vars = hess.numRows();
    for (j=0; j<num_deriv_vars; ++j)
      for (k=0; k<=j; ++k)
	if (!boost::math::isfinite<Real>(hess(j,k)))
	  { failed_data = 4; break; }             // this & higher deriv orders
  }
}


inline void SurrogateData::data_checks()
{
  if (anchor())
    response_check(anchor_response(), sdRep->failedAnchorData);

  sdRep->failedRespData.clear();
  const SDRArray& resp_data = response_data();
  size_t i, num_resp = resp_data.size(); short failed_data;
  for (i=0; i<num_resp; ++i) {
    response_check(resp_data[i], failed_data);
    if (failed_data)
      sdRep->failedRespData[i] = failed_data;
  }

#ifdef DEBUG
  if (sdRep->failedAnchorData) {
    PCout << "failedAnchorData = " << sdRep->failedAnchorData << '\n';
  if (!sdRep->failedRespData.empty()) {
    PCout << "failedRespData:\n";
    for (SizetShortMap::iterator it=sdRep->failedRespData.begin();
	 it!=sdRep->failedRespData.end(); ++it)
      PCout << "index: " << std::setw(6) << it->first
	    << " data: " << it->second << '\n';
  }
#endif // DEBUG
}


inline short SurrogateData::failed_anchor_data() const
{ return sdRep->failedAnchorData; }


inline const SizetShortMap& SurrogateData::failed_response_data() const
{ return sdRep->failedRespData; }


inline void SurrogateData::clear_anchor()
{
  sdRep->anchorVars = SurrogateDataVars();
  sdRep->anchorResp = SurrogateDataResp();
}


inline void SurrogateData::clear_data()
{ sdRep->varsData.clear(); sdRep->respData.clear(); }


inline void SurrogateData::clear_popped()
{
  sdRep->poppedVarsTrials.clear(); sdRep->poppedRespTrials.clear();
  sdRep->popCountStack.clear();
}


inline void SurrogateData::clear_stored()
{ sdRep->storedVarsData.clear(); sdRep->storedRespData.clear(); }


inline SurrogateDataRep* SurrogateData::data_rep() const
{ return sdRep; }


inline void SurrogateData::failed_anchor_data(short fail_anchor)
{ sdRep->failedAnchorData = fail_anchor; }


inline void SurrogateData::failed_response_data(const SizetShortMap& fail_resp)
{ sdRep->failedRespData = fail_resp; }


inline const SDV2DArray& SurrogateData::stored_variables_data() const
{ return sdRep->storedVarsData; }


inline void SurrogateData::stored_variables_data(const SDV2DArray& stored_vars)
{ sdRep->storedVarsData = stored_vars; }


inline const SDR2DArray& SurrogateData::stored_response_data() const
{ return sdRep->storedRespData; }


inline void SurrogateData::stored_response_data(const SDR2DArray& stored_resp)
{ sdRep->storedRespData = stored_resp; }


inline const SDV2DArray& SurrogateData::popped_variables_trials() const
{ return sdRep->poppedVarsTrials; }


inline void SurrogateData::
popped_variables_trials(const SDV2DArray& popped_vars)
{ sdRep->poppedVarsTrials = popped_vars; }


inline const SDR2DArray& SurrogateData::popped_response_trials() const
{ return sdRep->poppedRespTrials; }


inline void SurrogateData::
popped_response_trials(const SDR2DArray& popped_resp)
{ sdRep->poppedRespTrials = popped_resp; }


inline SurrogateData SurrogateData::copy(short sdv_mode, short sdr_mode) const
{
  SurrogateData sd;
  bool anchor_pt = anchor();

  if (sdv_mode == DEEP_COPY) {
    if (anchor_pt) sd.anchor_variables(sdRep->anchorVars.copy());

    size_t i, num_pts = sdRep->varsData.size();
    SDVArray sdv_array(num_pts);
    for (i=0; i<num_pts; ++i)
      sdv_array[i] = sdRep->varsData[i].copy();
    sd.variables_data(sdv_array);

    size_t j, num_sdva = sdRep->storedVarsData.size();
    SDV2DArray stored_sdv(num_sdva);
    for (i=0; i<num_sdva; ++i) {
      const SDVArray& rep_stored_i = sdRep->storedVarsData[i];
      num_pts = rep_stored_i.size();
      stored_sdv[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	stored_sdv[i][j] = rep_stored_i[j].copy();
    }
    sd.stored_variables_data(stored_sdv);

    num_sdva = sdRep->poppedVarsTrials.size();
    SDV2DArray popped_sdv(num_sdva);
    for (i=0; i<num_sdva; ++i) {
      const SDVArray& rep_popped_i = sdRep->poppedVarsTrials[i];
      num_pts = rep_popped_i.size();
      popped_sdv[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	popped_sdv[i][j] = rep_popped_i[j].copy();
    }
    sd.popped_variables_trials(popped_sdv);
  }
  else { // shallow SDV copies based on operator=
    if (anchor_pt) sd.anchor_variables(sdRep->anchorVars);
    sd.variables_data(sdRep->varsData);
    sd.stored_variables_data(sdRep->storedVarsData);
    sd.popped_variables_trials(sdRep->poppedVarsTrials);
  }

  if (sdr_mode == DEEP_COPY) {
    if (anchor_pt) sd.anchor_response(sdRep->anchorResp.copy());

    size_t i, num_pts = sdRep->respData.size();
    SDRArray sdr_array(num_pts);
    for (i=0; i<num_pts; ++i)
      sdr_array[i] = sdRep->respData[i].copy();
    sd.response_data(sdr_array);

    size_t j, num_sdra = sdRep->storedRespData.size();
    SDR2DArray stored_sdr(num_sdra);
    for (i=0; i<num_sdra; ++i) {
      const SDRArray& rep_stored_i = sdRep->storedRespData[i];
      num_pts = rep_stored_i.size();
      stored_sdr[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	stored_sdr[i][j] = rep_stored_i[j].copy();
    }
    sd.stored_response_data(stored_sdr);

    num_sdra = sdRep->poppedRespTrials.size();
    SDR2DArray popped_sdr(num_sdra);
    for (i=0; i<num_sdra; ++i) {
      const SDRArray& rep_popped_i = sdRep->poppedRespTrials[i];
      num_pts = rep_popped_i.size();
      popped_sdr[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	popped_sdr[i][j] = rep_popped_i[j].copy();
    }
    sd.popped_response_trials(popped_sdr);
  }
  else { // shallow SDR copies based on operator=
    if (anchor_pt) sd.anchor_response(sdRep->anchorResp);
    sd.response_data(sdRep->respData);
    sd.stored_response_data(sdRep->storedRespData);
    sd.popped_response_trials(sdRep->poppedRespTrials);
  }

  if (anchor_pt) sd.failed_anchor_data(sd.failed_anchor_data());
  sd.failed_response_data(sd.failed_response_data());

  return sd;
}

} // namespace Pecos

#endif
