/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef SURROGATE_DATA_HPP
#define SURROGATE_DATA_HPP

#include "pecos_data_types.hpp"
#include "ActiveKey.hpp"
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

public:
  /// destructor
  ~SurrogateDataVarsRep();

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

  /// lightweight constructor (common use case: continuous variables)
  SurrogateDataVarsRep(const RealVector& c_vars, short mode);
  /// alternate lightweight constructor (data sizing only)
  SurrogateDataVarsRep(size_t num_c_vars);


  //
  //- Heading: Private data members
  //

  RealVector continuousVars;   ///< continuous variables
  IntVector  discreteIntVars;  ///< discrete integer variables
  RealVector discreteRealVars; ///< discrete real variables
};


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(const RealVector& c_vars, const IntVector& di_vars,
		     const RealVector& dr_vars, short mode)
{
  // Note: provided a way to query DataAccess mode for c_vars, could make
  // greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (mode == DEEP_COPY) {         // enforce deep vector copy
    if (!c_vars.empty())  copy_data( c_vars, continuousVars);
    if (!di_vars.empty()) copy_data(di_vars, discreteIntVars);
    if (!dr_vars.empty()) copy_data(dr_vars, discreteRealVars);
  }
  else if (mode == SHALLOW_COPY) { // enforce shallow vector copy
    if (!c_vars.empty())
      continuousVars
	= RealVector(Teuchos::View,  c_vars.values(),  c_vars.length());
    if (!di_vars.empty())
      discreteIntVars
	= IntVector(Teuchos::View,  di_vars.values(), di_vars.length());
    if (!dr_vars.empty())
      discreteRealVars
	= RealVector(Teuchos::View, dr_vars.values(), dr_vars.length());
  }
  else {                           // default: assume existing Copy/View state
    if (!c_vars.empty())    continuousVars = c_vars;
    if (!di_vars.empty())  discreteIntVars = di_vars;
    if (!dr_vars.empty()) discreteRealVars = dr_vars;
  }
}


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars)
{
  continuousVars.sizeUninitialized(num_c_vars);
  discreteIntVars.sizeUninitialized(num_di_vars);
  discreteRealVars.sizeUninitialized(num_dr_vars);
}


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(const RealVector& c_vars, short mode)
{
  // Note: provided a way to query DataAccess mode for c_vars, could make
  // greater use of operator= for {DEEP,SHALLOW}_COPY modes
  if (!c_vars.empty())
    switch (mode) {
    case DEEP_COPY:
      copy_data(c_vars, continuousVars); break;
    case SHALLOW_COPY:
      continuousVars
	= RealVector(Teuchos::View, c_vars.values(), c_vars.length());
      break;
    default: // assume existing Copy/View state
      continuousVars = c_vars; break;
    }
}


inline SurrogateDataVarsRep::
SurrogateDataVarsRep(size_t num_c_vars)
{ continuousVars.sizeUninitialized(num_c_vars); }


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

  /// lightweight constructor
  SurrogateDataVars(const RealVector& c_vars, short mode = DEFAULT_COPY);
  /// alternate lightweight constructor (data sizing only)
  SurrogateDataVars(size_t num_c_vars);

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

  /// instantiate a sdvRep
  void create_rep(size_t num_vars);

  /// return number of continuous variables
  size_t cv() const;
  /// return number of discrete integer variables
  size_t div() const;
  /// return number of discrete real variables
  size_t drv() const;

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
  /// output response function, gradient, and Hessian data
  void write(std::ostream& s) const;

private:

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  std::shared_ptr<SurrogateDataVarsRep> sdvRep;
};


inline SurrogateDataVars::SurrogateDataVars()
{ }


// BMA NOTE: The following don't use make_shared<SurrogateDataVarsRep>()
// due to private ctors

inline SurrogateDataVars::
SurrogateDataVars(const RealVector& c_vars, const IntVector& di_vars,
		  const RealVector& dr_vars, short mode):
  sdvRep(new SurrogateDataVarsRep(c_vars, di_vars, dr_vars, mode))
{ }


inline SurrogateDataVars::
SurrogateDataVars(size_t num_c_vars, size_t num_di_vars, size_t num_dr_vars):
  sdvRep(new SurrogateDataVarsRep(num_c_vars, num_di_vars, num_dr_vars))
{ }


inline SurrogateDataVars::
SurrogateDataVars(const RealVector& c_vars, short mode):
  sdvRep(new SurrogateDataVarsRep(c_vars, mode))
{ }


inline SurrogateDataVars::SurrogateDataVars(size_t num_c_vars):
  sdvRep(new SurrogateDataVarsRep(num_c_vars))
{ }


inline SurrogateDataVars::SurrogateDataVars(const SurrogateDataVars& sdv):
  sdvRep(sdv.sdvRep)  
{ }


inline SurrogateDataVars::~SurrogateDataVars()
{ }


inline SurrogateDataVars& SurrogateDataVars::
operator=(const SurrogateDataVars& sdv)
{
  sdvRep = sdv.sdvRep;
  return *this;
}


//inline bool SurrogateDataVars::operator==(const SurrogateDataVars& sdv) const
//{
//  return (sdvRep->continuousVars   == sdv.sdvRep->continuousVars  &&
//          sdvRep->discreteIntVars  == sdv.sdvRep->discreteIntVars &&
//          sdvRep->discreteRealVars == sdv.sdvRep->discreteRealVars) ?
//    true : false;
//}


inline void SurrogateDataVars::create_rep(size_t num_vars)
{ sdvRep.reset(new SurrogateDataVarsRep(num_vars)); }


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


inline size_t SurrogateDataVars::cv() const
{ return sdvRep->continuousVars.length(); }


inline size_t SurrogateDataVars::div() const
{ return sdvRep->discreteIntVars.length(); }


inline size_t SurrogateDataVars::drv() const
{ return sdvRep->discreteRealVars.length(); }


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


inline void SurrogateDataVars::write(std::ostream& s) const
{
  s << "SDV continuous variables:\n"    << sdvRep->continuousVars
    << "SDV discrete int variables:\n"  << sdvRep->discreteIntVars
    << "SDV discrete real variables:\n" << sdvRep->discreteRealVars << '\n';
}


/// std::ostream insertion operator for SurrogateDataVars
inline std::ostream& operator<<(std::ostream& s, const SurrogateDataVars& sdv)
{ sdv.write(s); return s; }


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

public:
  /// destructor
  ~SurrogateDataRespRep();

private:

  //
  //- Heading: Private member functions
  //

  /// constructor
  SurrogateDataRespRep(Real fn_val, const RealVector& fn_grad,
		       const RealSymMatrix& fn_hess, short bits, short mode);
  /// lightweight constructor (common use case)
  SurrogateDataRespRep(Real fn_val);
  /// alternate constructor (data sizing only)
  SurrogateDataRespRep(short bits, size_t num_vars);

  //
  //- Heading: Private data members
  //

  short           activeBits; ///< active data bits: 1 (fn), 2 (grad), 4 (hess)
  Real            responseFn; ///< truth response function value
  RealVector    responseGrad; ///< truth response function gradient
  RealSymMatrix responseHess; ///< truth response function Hessian
};


inline SurrogateDataRespRep::
SurrogateDataRespRep(Real fn_val, const RealVector& fn_grad,
		     const RealSymMatrix& fn_hess, short bits, short mode):
  responseFn(fn_val), // always deep copy for scalars
  activeBits(bits)
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


inline SurrogateDataRespRep::SurrogateDataRespRep(Real fn_val):
  responseFn(fn_val), // always deep copy for scalars
  activeBits(1)
{ }


inline SurrogateDataRespRep::
SurrogateDataRespRep(short bits, size_t num_vars):
  activeBits(bits)
{
  if (bits & 2)
    responseGrad.sizeUninitialized(num_vars);
  if (bits & 4)
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
  /// lightweight constructor
  SurrogateDataResp(Real fn_val);
  /// alternate constructor (data sizing only)
  SurrogateDataResp(short bits, size_t num_vars);
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

  /// instantiate a sdvRep
  void create_rep(short bits, size_t num_vars);

  /// return deep copy of SurrogateDataResp instance
  SurrogateDataResp copy() const;

  /// set activeBits
  void active_bits(short val);
  /// get activeBits
  short active_bits() const;

  /// size response{Grad,Hess}
  void derivative_variables(size_t num_v);
  /// return size of response{Grad,Hess} (0 if inactive)
  size_t derivative_variables() const;

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
  std::shared_ptr<SurrogateDataRespRep> sdrRep;
};


inline SurrogateDataResp::SurrogateDataResp()
{ }

// BMA NOTE: The following don't use make_shared<SurrogateDataRespRep>()
// due to private ctors

inline SurrogateDataResp::
SurrogateDataResp(Real fn_val, const RealVector& fn_grad,
		  const RealSymMatrix& fn_hess, short bits, short mode):
  sdrRep(new SurrogateDataRespRep(fn_val, fn_grad, fn_hess, bits, mode))
{ }


inline SurrogateDataResp::SurrogateDataResp(Real fn_val):
  sdrRep(new SurrogateDataRespRep(fn_val))
{ }


inline SurrogateDataResp::
SurrogateDataResp(short bits, size_t num_vars):
  sdrRep(new SurrogateDataRespRep(bits, num_vars))
{ }


inline SurrogateDataResp::SurrogateDataResp(const SurrogateDataResp& sdr):
  sdrRep(sdr.sdrRep)
{ }


inline SurrogateDataResp::~SurrogateDataResp()
{ }


inline SurrogateDataResp& SurrogateDataResp::
operator=(const SurrogateDataResp& sdr)
{
  sdrRep = sdr.sdrRep;
  return *this;
}


//inline bool SurrogateDataResp::operator==(const SurrogateDataResp& sdr) const
//{
//  return ( sdrRep->responseFn   == sdr.sdrRep->responseFn   &&
//	     sdrRep->responseGrad == sdr.sdrRep->responseGrad &&
//	     sdrRep->responseHess == sdr.sdrRep->responseHess ) ? true : false;
//}


inline void SurrogateDataResp::create_rep(short bits, size_t num_vars)
{ sdrRep.reset(new SurrogateDataRespRep(bits, num_vars)); }


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


inline void SurrogateDataResp::derivative_variables(size_t num_v)
{
  if (sdrRep->activeBits & 2)
    sdrRep->responseGrad.sizeUninitialized(num_v);
  if (sdrRep->activeBits & 4)
    sdrRep->responseHess.shapeUninitialized(num_v);
}


inline size_t SurrogateDataResp::derivative_variables() const
{
  return std::max(sdrRep->responseGrad.length(),
		  sdrRep->responseHess.numRows());
}


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
    s << "SDR function value    =  " << std::setw(WRITE_PRECISION+7)
      << sdrRep->responseFn << '\n';
  if (sdrRep->activeBits & 2) {
    s << "SDR function gradient =\n";
    write_data_trans(s, sdrRep->responseGrad, true, true, true);
  }
  if (sdrRep->activeBits & 4)
    s << "SDR function Hessian  =\n" << sdrRep->responseHess;
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

  ~SurrogateDataRep(); ///< destructor

private:

  //
  //- Heading: Constructors and destructor
  //

  SurrogateDataRep();  ///< constructor

  //
  //- Heading: Member functions
  //

  /// update {varsData,respData}Iter from activeKey
  void update_active_iterators();

  //
  //- Heading: Private data members
  //

  /// database of reference variable data sets, with lookup by model/level index
  std::map<ActiveKey, SDVArray> varsData;
  /// iterator to active entry within varsData
  std::map<ActiveKey, SDVArray>::iterator varsDataIter;
  /// filtered database of variable data sets, a subset drawn from varsData
  std::map<ActiveKey, SDVArray> filteredVarsData;
  
  /// database of reference response data sets, with lookup by model/level index
  std::map<ActiveKey, SDRArray> respData;
  /// iterator to active entry within respData
  std::map<ActiveKey, SDRArray>::iterator respDataIter;
  /// filtered database of response data sets, a subset drawn from respData
  std::map<ActiveKey, SDRArray> filteredRespData;
  /// additive and multiplicative factors for scaling active response
  /// functions to [0,1]:  scaled_fn = (unscaled - min) / range  -->
  /// unscaled = range * scaled_fn + min  where RealRealPair = (min,range)
  /*std::map<ActiveKey,*/ RealRealPair respFnScaling;

  /// optional data identifiers, one per entry in {vars,Resp}Data, allowing
  /// individual data point interactions (e.g., replace(id))
  /** use an un-ordered container (with slow look-ups) since push/pop
      could modify original (sequential) order */
  std::map<ActiveKey, IntArray> dataIdentifiers;
  /// iterator to active entry within dataIdentifiers
  std::map<ActiveKey, IntArray>::iterator dataIdsIter;

  /// sets of popped variables data sets, with lookup by model/level index.
  /// Each popped set is an SDVArray extracted from varsData.
  std::map<ActiveKey, SDVArrayDeque> poppedVarsData;
  /// sets of popped response data sets, with lookup by model/level index.
  /// Each popped set is an SDRArray extracted from respData.
  std::map<ActiveKey, SDRArrayDeque> poppedRespData;
  /// sets of popped data identifiers, with lookup by active key
  std::map<ActiveKey, IntArrayDeque> poppedDataIds;
  /// a stack managing the number of points previously appended that
  /// can be removed by calls to pop()
  std::map<ActiveKey, SizetArray> popCountStack;

  /// database key indicating the currently active {SDV,SDR}Arrays.
  /// the key is a multi-index managing multiple modeling dimensions
  /// such as model form, doscretization level, etc.
  ActiveKey activeKey;

  /// index of anchor point within {vars,resp}Data, _NPOS if none; for now,
  /// we restrict anchor to reference data to simplify bookkeeping (assume
  /// anchor does not migrate within pushed/popped data)
  std::map<ActiveKey, size_t> anchorIndex;

  /// map from failed respData indices to failed data bits; defined
  /// in sample_checks() and used for fault tolerance
  std::map<ActiveKey, SizetShortMap> failedRespData;
};


inline SurrogateDataRep::SurrogateDataRep():
  respFnScaling(0.,0.), dataIdsIter(dataIdentifiers.end())
{ }


inline SurrogateDataRep::~SurrogateDataRep()
{ }


inline void SurrogateDataRep::update_active_iterators()
{
  if (dataIdsIter != dataIdentifiers.end() && dataIdsIter->first == activeKey)
    return;

  varsDataIter = varsData.find(activeKey);
  respDataIter = respData.find(activeKey);
  dataIdsIter  = dataIdentifiers.find(activeKey);

  /* So long as we only create new keys and avoid modifying existing ones,
     this deep copy is not needed.
  ActiveKey active_copy; // share 1 deep copy of current active key
  if (varsDataIter == varsData.end() || respDataIter == respData.end() ||
      dataIdsIter == dataIdentifiers.end())
    active_copy = activeKey.copy();
  */

  if (varsDataIter == varsData.end()) {
    std::pair<ActiveKey, SDVArray> vd_pair(activeKey/*active_copy*/,SDVArray());
    varsDataIter = varsData.insert(vd_pair).first;
  }
  if (respDataIter == respData.end()) {
    std::pair<ActiveKey, SDRArray> vr_pair(activeKey/*active_copy*/,SDRArray());
    respDataIter = respData.insert(vr_pair).first;
  }
  if (dataIdsIter == dataIdentifiers.end()) {
    std::pair<ActiveKey, IntArray> ia_pair(activeKey/*active_copy*/,IntArray());
    dataIdsIter = dataIdentifiers.insert(ia_pair).first;
  }
}


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

  SurrogateData();                        ///< default handle ctor (no body)
  SurrogateData(bool handle);             ///< handle + body constructor
  SurrogateData(const ActiveKey& key);  ///< handle + body constructor
  SurrogateData(const SurrogateData& sd); ///< copy constructor
  ~SurrogateData();                       ///< destructor

  /// assignment operator
  SurrogateData& operator=(const SurrogateData& sdv);

  //
  //- Heading: Member functions
  //

  /// copy incoming sd instance with options for shallow or deep copies
  /// of the SurrogateData{Vars,Resp} objects
  void copy(const SurrogateData& sd, short sdv_mode = DEEP_COPY,
	    short sdr_mode = DEEP_COPY) const;
  /// returns a new SurrogateData instance with options for shallow or
  /// deep copies of the SurrogateData{Vars,Resp} objects
  SurrogateData copy(short sdv_mode = DEEP_COPY,
		     short sdr_mode = DEEP_COPY) const;
  /// copy active data for incoming sd instance with options for
  /// shallow or deep copies of the SurrogateData{Vars,Resp} objects
  void copy_active(const SurrogateData& sd, short sdv_mode,
		   short sdr_mode) const;

  void size_active_sdv(const SurrogateData& sd) const;
  void copy_active_sdv(const SurrogateData& sd, short sdv_mode) const;
  void copy_active_pop_sdv(const SurrogateData& sd, short sdv_mode) const;
  void size_active_sdr(const SDRArray& sdr_array) const;
  void size_active_sdr(const SurrogateData& sd) const;
  void copy_active_sdr(const SurrogateData& sd, short sdr_mode) const;
  void copy_active_pop_sdr(const SurrogateData& sd, short sdr_mode) const;

  /// resize {vars,resp}Data
  void resize(size_t num_pts);
  /// resize {vars,resp}Data and allocate new SDV/SDR
  void resize(size_t num_pts, short bits, size_t num_vars);

  /// augment {vars,resp}Data and define anchorIndex
  void anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr,
		    bool append = true);
  /// set {vars,resp}Data
  void data_points(const SDVArray& sdv_array, const SDRArray& sdr_array);

  // augment varsData and define anchorIndex
  //void anchor_variables(const SurrogateDataVars& sdv);
  // augment respData and define anchorIndex
  //void anchor_response(const SurrogateDataResp& sdr);

  /// get varsData instance corresponding to anchorIndex
  const SurrogateDataVars& anchor_variables() const;
  /// get respData instance corresponding to anchorIndex
  const SurrogateDataResp& anchor_response() const;

  /// set varsData[activeKey]
  void variables_data(const SDVArray& sdv_array);
  /// get varsData[activeKey]
  const SDVArray& variables_data() const;
  /// get varsData[key]
  const SDVArray& variables_data(const ActiveKey& key) const;
  /// get varsData[activeKey]
  SDVArray& variables_data();

  /// return the i-th active continuous variables
  const RealVector& continuous_variables(size_t i) const;

  /// set respData[activeKey]
  void response_data(const SDRArray& sdr_array);
  /// get respData[activeKey]
  const SDRArray& response_data() const;
  /// get respData[key]
  const SDRArray& response_data(const ActiveKey& key) const;
  /// get respData[activeKey]
  SDRArray& response_data();

  /// return the i-th active response function
  Real response_function(size_t i) const;
  /// return the i-th active response function, scaled by active respFnScaling
  Real scaled_response_function(size_t i) const;

  /// return active respFnScaling
  const RealRealPair& response_function_scaling() const;
  /// check for valid response scaling (non-zero range defined from at
  /// least two different points)
  bool valid_response_scaling() const;

  /// set dataIdentifiers
  void data_ids(const IntArray& ids);
  /// get dataIdentifiers
  const IntArray& data_ids() const;

  /// get varsData
  const std::map<ActiveKey, SDVArray>& variables_data_map() const;
  /// build/return filteredVarsData, pulling aggregated/nonaggregated
  /// keys from varsData
  const std::map<ActiveKey, SDVArray>&
    filtered_variables_data_map(short mode) const;
  /// set varsData
  void variables_data_map(const std::map<ActiveKey, SDVArray>& vars_map);

  /// get respData
  const std::map<ActiveKey, SDRArray>& response_data_map() const;
  /// build/return filteredRespData, pulling aggregated/reduced/raw data
  /// sets from respData
  const std::map<ActiveKey, SDRArray>&
    filtered_response_data_map(short mode) const;
  /// set respData
  void response_data_map(const std::map<ActiveKey, SDRArray>& resp_map);

  /// return a particular key instance present within the
  /// aggregated/reduced/raw data subset (respData is used)
  const ActiveKey& filtered_key(short mode, size_t index) const;

  /// push sdv onto end of varsData
  void push_back(const SurrogateDataVars& sdv);
  /// push sdr onto end of respData
  void push_back(const SurrogateDataResp& sdr);
  /// push {sdv,sdr} onto ends of {vars,resp}Data
  void push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr);
  /// push eval_id onto end of dataIdentifiers
  void push_back(int eval_id);

  /// remove the first entry from {vars,resp}Data, managing anchorIndex
  /// Note: inefficient for std::vector's, but needed in rare cases.
  void pop_front(size_t num_pop = 1);
  /// remove the last entry from {vars,resp}Data, managing anchorIndex
  void pop_back(size_t num_pop = 1);

  /// remove num_pop_pts entries from ends of active entry in {vars,resp}Data
  void pop(bool save_data = true);
  /// remove num_pop_pts entries from ends of keyed entries in {vars,resp}Data
  void pop(const ActiveKey& key, bool save_data);
  /// return a previously popped data set (identified by index) to the
  /// ends of active entry in {vars,resp}Data
  void push(size_t index, bool erase_popped = true);
  /// return previously popped data sets (identified by index) to the
  /// ends of keyed entries in {vars,resp}Data
  void push(const ActiveKey& key, size_t index, bool erase_popped = true);

  /// replace the SurrogateDataVars instance identified by data_id
  void replace(const SurrogateDataVars& sdv, int data_id);
  /// replace the SurrogateDataResp instance identified by data_id
  void replace(const SurrogateDataResp& sdr, int data_id);
  /// replace the SurrogateDataVars and SurrogateDataResp instances
  /// identified by data_id
  void replace(const SurrogateDataVars& sdv,
	       const SurrogateDataResp& sdr, int data_id);

  /// append count to popCountStack[activeKey]
  void pop_count(size_t count) const;
  /// return popCountStack[activeKey].back()
  size_t pop_count() const;
  /// assign popCountStack[activeKey]
  void pop_count_stack(const SizetArray& pop_count) const;
  /// return popCountStack[key]
  const SizetArray& pop_count_stack(const ActiveKey& key) const;
  /// return popCountStack[activeKey]
  const SizetArray& pop_count_stack() const;
  /// assign popCountStack
  void pop_count_stack_map(
    const std::map<ActiveKey, SizetArray>& pcs_map) const;
  /// return popCountStack
  const std::map<ActiveKey, SizetArray>& pop_count_stack_map() const;

  /// pop records from front of {vars,resp}Data to achieve target length,
  /// for each key in passed set
  void history_target(size_t target, const ActiveKey& key);

  /// query presence of anchor indexed within {vars,resp}Data
  bool anchor() const;
  /// assign anchorIndex[activeKey] to incoming index
  void anchor_index(size_t index) const;
  /// assign anchorIndex[key] to incoming index
  void anchor_index(size_t index, const ActiveKey& key) const;
  /// return anchorIndex[activeKey], if defined
  size_t anchor_index() const;
  /// return anchorIndex[key], if defined
  size_t anchor_index(const ActiveKey& key) const;
  /// assign anchorIndex
  void anchor_index_map(const std::map<ActiveKey, size_t>& ai_map) const;
  /// return anchorIndex
  const std::map<ActiveKey, size_t>& anchor_index_map() const;
  /// erase anchorIndex[activeKey]
  void clear_anchor_index();
  /// erase anchorIndex[key]
  void clear_anchor_index(const ActiveKey& key);

  /// return size of active key within {vars,resp}Data
  size_t points() const;
  /// return size of {vars,resp}Data instance corresponding to key
  size_t points(const ActiveKey& key) const;

  /// synchonize data size for aggregate reduction key based on data
  /// sizes for embedded keys
  void synchronize_reduction_size();

  /// return total number of available data components
  size_t response_size() const;
  /// return number of failed data components
  size_t failed_response_size() const;
  /// return net number of active data components (total minus failed)
  size_t active_response_size() const;

  /// return number of 1D arrays within popped{Vars,Resp}Data 2D arrays
  /// identified by key
  size_t popped_sets(const ActiveKey& key) const;
  /// return number of 1D arrays within popped{Vars,Resp}Data 2D arrays
  /// identified by activeKey
  size_t popped_sets() const;

  /// return number of gradient variables from size of gradient arrays
  size_t num_gradient_variables() const;
  /// return number of Hessian variables from size of Hessian arrays
  size_t num_hessian_variables() const;
  /// return number of derivative variables as indicated by size of
  /// gradient/Hessian arrays
  size_t num_derivative_variables() const;

  void clear_response_function_scaling() const;
  /// compute a simple bounds-based scaling for response functions
  void compute_response_function_scaling() const;
  /// assign a scalar shift for the set of response function data
  void assign_response_function_shift(Real shift) const;

  /// convenience function used by data_checks() for respData
  void response_check(const SurrogateDataResp& sdr, short& failed_data) const;
  /// screen data sets for samples with Inf/Nan that should be excluded;
  /// defines failedRespData
  void data_checks() const;
  /// return failedRespData corresponding to active anchorIndex
  short failed_anchor_data() const;
  /// assign active failedRespData
  void failed_response_data(const SizetShortMap& fail_data) const;
  /// return failedRespData[key]
  const SizetShortMap& failed_response_data(const ActiveKey& key) const;
  /// return failedRespData[activeKey]
  const SizetShortMap& failed_response_data() const;

  /// assign activeKey and update active iterators
  void active_key(const ActiveKey& key) const;
  /// return activeKey
  const ActiveKey& active_key() const;
  /// searches for key and updates {vars,resp}DataIter only if found
  bool contains(const ActiveKey& key);

  /// clear active key within {vars,resp}Data
  void clear_active_data();
  /// clear a set of keys within {vars,resp}Data
  void clear_active_data(const ActiveKey& key);
  /// clear all inactive data within {vars,resp}Data
  void clear_inactive_data();

  /// clear filtered{Vars,Resp}Data (shallow copy subsets of {vars,resp}Data)
  void clear_filtered();

  /// clear active key within popped{Vars,Resp}Data
  void clear_active_popped();
  /// clear a set of keys within popped{Vars,Resp}Data
  void clear_active_popped(const ActiveKey& key);
  /// clear popped{Vars,Resp}Data and restore to initial popped state
  void clear_popped();

  /// clear {vars,resp}Data and restore to initial data state
  void clear_data(bool initialize = true);

  /// clear all keys for all maps and optionally restore to initial state
  void clear_all(bool initialize = true);
  /// invokes both clear_active_data() and clear_active_popped()
  void clear_all_active();
  /// invokes both clear_active_data(keys) and clear_active_popped(keys)
  void clear_all_active(const ActiveKey& key);

  /// return sdRep
  std::shared_ptr<SurrogateDataRep> data_rep() const;

  /// function to check sdRep (does this handle contain a body)
  bool is_null() const;

private:

  //
  //- Heading: Member functions
  //

  /// define or update anchorIndex[activeKey]
  size_t assign_anchor_index();
  /// retrieve anchorIndex[activeKey]
  size_t retrieve_anchor_index(bool hard_fail) const;
  /// retrieve anchorIndex[key]
  size_t retrieve_anchor_index(const ActiveKey& key, bool hard_fail) const;

  /// helper function for an individual key
  void anchor_index(size_t index, std::map<ActiveKey, size_t>& anchor_index_map,
		    const ActiveKey& key) const;

  /// assign sdv within varsData[activeKey] at indicated index
  void assign_variables(const SurrogateDataVars& sdv, size_t index);
  /// assign sdr within respData[activeKey] at indicated index
  void assign_response(const SurrogateDataResp& sdr, size_t index);

  /// set failedRespData
  void failed_response_data_map(
    const std::map<ActiveKey, SizetShortMap>&	fail_resp) const;
  /// get failedRespData
  const std::map<ActiveKey, SizetShortMap>& failed_response_data_map() const;

  /// set active poppedVarsData
  void popped_variables(const SDVArrayDeque& popped_vars);
  /// get active poppedVarsData
  const SDVArrayDeque& popped_variables() const;
  /// set active poppedRespData
  void popped_response(const SDRArrayDeque& popped_resp);
  /// get active poppedRespData
  const SDRArrayDeque& popped_response() const;

  /// set poppedVarsData
  void popped_variables_map(
    const std::map<ActiveKey, SDVArrayDeque>& popped_vars);
  /// get poppedVarsData
  const std::map<ActiveKey, SDVArrayDeque>& popped_variables_map() const;
  /// set poppedRespData
  void popped_response_map(
    const std::map<ActiveKey, SDRArrayDeque>& popped_resp);
  /// get poppedRespData
  const std::map<ActiveKey, SDRArrayDeque>& popped_response_map() const;

  /// lower level helper for data associated with a particular ActiveKey
  void history_target(size_t target, SDVArray& sdv_array, SDRArray& sdr_array,
		      std::map<ActiveKey, size_t>::iterator anchor_it);

  /// helper function for an individual key
  void pop_back(size_t num_pop, SDVArray& sdv_array, SDRArray& sdr_array);
  /// helper function for an individual key
  void pop_front(size_t num_pop, SDVArray& sdv_array, SDRArray& sdr_array);
  /// helper function for an individual key
  void pop(SDVArray& sdv_array, SDRArray& sdr_array, IntArray& data_ids,
	   std::map<ActiveKey, SizetArray>::iterator pop_cnt_it,
	   SDVArrayDeque& popped_sdv_arrays, SDRArrayDeque& popped_sdr_arrays,
	   IntArrayDeque& popped_ids, SizetShortMap& failed_resp,
	   bool save_data);
  /// helper function for an individual key
  void push(SDVArray& sdv_array, SDRArray& sdr_array, IntArray& data_ids,
	    SizetArray& pop_count_stack,
	    std::map<ActiveKey, SDVArrayDeque>::iterator pvd_it,
	    std::map<ActiveKey, SDRArrayDeque>::iterator prd_it,
	    std::map<ActiveKey, IntArrayDeque>::iterator pid_it,
	    SizetShortMap& failed_resp, size_t index, bool erase_popped);

  /// screen resp_data for samples with Inf/Nan that should be excluded;
  /// defines failed_resp
  void data_checks(const SDRArray& resp_data, SizetShortMap& failed_resp) const;

  //
  //- Heading: Private data members
  //
 
  /// pointer to the body (handle-body idiom)
  std::shared_ptr<SurrogateDataRep> sdRep;
};


inline SurrogateData::SurrogateData()
{ }


// BMA NOTE: The following don't use make_shared<SurrogateDataRep>()
// due to private ctors

inline SurrogateData::SurrogateData(bool handle):
  sdRep(new SurrogateDataRep())
{ sdRep->update_active_iterators(); } // default activeKey is empty array


inline SurrogateData::SurrogateData(const ActiveKey& key):
  sdRep(new SurrogateDataRep())
{ active_key(key); }


inline SurrogateData::SurrogateData(const SurrogateData& sd):
  sdRep(sd.sdRep)
{ }


inline SurrogateData::~SurrogateData()
{ }


inline SurrogateData& SurrogateData::operator=(const SurrogateData& sd)
{
  sdRep = sd.sdRep;
  return *this;
}


inline void SurrogateData::active_key(const ActiveKey& key) const
{
  if (sdRep->activeKey != key) {
    sdRep->activeKey = key;
    sdRep->update_active_iterators();

    // Seems a bit of overkill; for now, prefer synchronization calls from
    // contexts that append/pop data
    //if (key.raw_with_reduction_data())
    //  synchronize_reduction_size(); // reduction always derived from embedded
  }
}


inline const ActiveKey& SurrogateData::active_key() const
{ return sdRep->activeKey; }


inline size_t SurrogateData::points() const
{
  return std::min(sdRep->varsDataIter->second.size(),
		  sdRep->respDataIter->second.size());
}


inline size_t SurrogateData::points(const ActiveKey& key) const
{
  std::map<ActiveKey, SDVArray>::const_iterator sdv_it
    = sdRep->varsData.find(key);
  std::map<ActiveKey, SDRArray>::const_iterator sdr_it
    = sdRep->respData.find(key);
  return (sdv_it == sdRep->varsData.end() || sdr_it == sdRep->respData.end()) ?
    0 : std::min(sdv_it->second.size(), sdr_it->second.size());
}


inline void SurrogateData::synchronize_reduction_size()
{
  const ActiveKey& key = sdRep->activeKey;
  if (!key.aggregated() || !key.raw_with_reduction_data()) return;

  // Could streamline using only 1st embedded key, but retain generality for now
  std::vector<ActiveKey> embedded_keys;
  key.extract_keys(embedded_keys);
  size_t k, num_k = embedded_keys.size(), num_pts_k,
    num_pts = 0;//std::numeric_limits<size_t>::max();
  for (k=0; k<num_k; ++k) {
    num_pts_k = points(embedded_keys[k]);
    //if (num_pts_k < num_pts) num_pts = num_pts_k; // min points
    // Due to use of 2 embedded keys for synthetic surrogate data, which is
    // not yet defined, we use a max operation over these keys:
    if (num_pts_k > num_pts) num_pts = num_pts_k; // max points
  }
  if (num_pts != points()) resize(num_pts);
}


inline bool SurrogateData::contains(const ActiveKey& key)
{
  return (sdRep->varsData.find(key) != sdRep->varsData.end() &&
	  sdRep->respData.find(key) != sdRep->respData.end());
}


inline void SurrogateData::
data_points(const SDVArray& sdv_array, const SDRArray& sdr_array)
{
  sdRep->varsDataIter->second = sdv_array;
  sdRep->respDataIter->second = sdr_array;
}


inline void SurrogateData::
anchor_index(size_t index, std::map<ActiveKey, size_t>& anchor_index_map,
	     const ActiveKey& key) const
{
  if (index == _NPOS) { // remove entry, if present
    std::map<ActiveKey, size_t>::iterator it = anchor_index_map.find(key);
    if (it != anchor_index_map.end())
      anchor_index_map.erase(it);
  }
  else // add or overwrite entry
    anchor_index_map[key] = index;
}


inline void SurrogateData::
anchor_index(size_t index, const ActiveKey& key) const
{
  std::map<ActiveKey, size_t>& anchor_index_map = sdRep->anchorIndex;
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    anchor_index(index, anchor_index_map, key);
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k)
      anchor_index(index, anchor_index_map, embedded_keys[k]);
  }
}


inline void SurrogateData::anchor_index(size_t index) const
{ anchor_index(index, sdRep->activeKey); }


inline size_t SurrogateData::anchor_index(const ActiveKey& key) const
{ return retrieve_anchor_index(key, false); }


inline size_t SurrogateData::anchor_index() const
{ return retrieve_anchor_index(sdRep->activeKey, false); }


inline const std::map<ActiveKey, size_t>& SurrogateData::
anchor_index_map() const
{ return sdRep->anchorIndex; }


inline void SurrogateData::
anchor_index_map(const std::map<ActiveKey, size_t>& ai_map) const
{ sdRep->anchorIndex = ai_map; }


inline void SurrogateData::clear_anchor_index(const ActiveKey& key)
{
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    sdRep->anchorIndex.erase(key);
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k)
      sdRep->anchorIndex.erase(embedded_keys[k]);
  }
}


inline void SurrogateData::clear_anchor_index()
{ clear_anchor_index(sdRep->activeKey); }


inline size_t SurrogateData::assign_anchor_index()
{
  // This is often called in sequence of assign_anchor_variables() and
  // assign_anchor_response() --> use points() for consistent indexing
  size_t index = points(); // push_back() to follow
  // This approach reassigns an existing anchor index to freshly appended
  // anchor data --> previous anchor data is preserved but demoted
  anchor_index(index);

  /*
  // This approach preserves a previously assigned anchor index, which is good
  // for a pair of assign_anchor_{variables,response}() calls, but bad if a
  // previous sdv/sdr anchor assignment has not been cleared --> a subsequent
  // assign_{variables,response}() will overwrite the previous data, corrupting
  // the history for multipoint approximations.
  std::map<ActiveKey, size_t>& anchor_index = sdRep->anchorIndex;
  const ActiveKey& key = sdRep->activeKey;
  std::map<ActiveKey, size_t>::iterator anchor_it = anchor_index.find(key);
  size_t index;
  if (anchor_it == anchor_index.end()) // no anchor defined
    anchor_index[key] = index = points();
  else {
    index = anchor_it->second;
    if (index == _NPOS) // reassign only if undefined, else preserve
      anchor_it->second = index = points();
  }
  */

  return index;
}


inline size_t SurrogateData::
retrieve_anchor_index(const ActiveKey& key, bool hard_fail) const
{
  std::map<ActiveKey, size_t>& anchor_index = sdRep->anchorIndex;
  std::map<ActiveKey, size_t>::iterator anchor_it = anchor_index.find(key);
  if (anchor_it == anchor_index.end() || anchor_it->second == _NPOS) {
    if (hard_fail) {
      PCerr << "Error: lookup failure in SurrogateData::retrieve_anchor_index"
	    << "()." << std::endl;
      abort_handler(-1);
    }
    return _NPOS;
  }
  else
    return anchor_it->second;
}


inline size_t SurrogateData::retrieve_anchor_index(bool hard_fail) const
{ return retrieve_anchor_index(sdRep->activeKey, hard_fail); }


inline void SurrogateData::
assign_variables(const SurrogateDataVars& sdv, size_t index)
{
  SDVArray& sdv_array = sdRep->varsDataIter->second;
  if (index == sdv_array.size() || index == _NPOS) sdv_array.push_back(sdv);
  else sdv_array[index] = sdv;
}


inline void SurrogateData::
assign_response(const SurrogateDataResp& sdr, size_t index)
{
  SDRArray& sdr_array = sdRep->respDataIter->second;
  if (index == sdr_array.size() || index == _NPOS) sdr_array.push_back(sdr);
  else sdr_array[index] = sdr;
}


inline void SurrogateData::
anchor_point(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr,
	     bool append)
{
  size_t new_index;
  if (append)
    new_index = assign_anchor_index();
  else {
    size_t curr_index = retrieve_anchor_index(false); // no hard error
    new_index = (curr_index == _NPOS) ? assign_anchor_index() : curr_index;
  }
  assign_variables(sdv, new_index);
  assign_response(sdr,  new_index);
}


/*
  *** FRAGILE: index assignment relies on points() using min vars,resp size
  *** Migrate to anchor_point() with append = true: default for backwards
      compat --> DiscrepancyCorrection could use false or clear_active_data()
inline void SurrogateData::anchor_variables(const SurrogateDataVars& sdv)
{
  size_t index = assign_anchor_index();
  assign_variables(sdv, index);
}


inline void SurrogateData::anchor_response(const SurrogateDataResp& sdr)
{
  size_t index = assign_anchor_index();
  assign_response(sdr, index);
}
*/


inline const SurrogateDataVars& SurrogateData::anchor_variables() const
{
  size_t index = retrieve_anchor_index(true); // abort on index error
  return sdRep->varsDataIter->second[index];
}


inline const SurrogateDataResp& SurrogateData::anchor_response() const
{
  size_t index = retrieve_anchor_index(true); // abort on index error
  return sdRep->respDataIter->second[index];
}


inline void SurrogateData::variables_data(const SDVArray& sdv_array)
{ sdRep->varsDataIter->second = sdv_array; }


inline const SDVArray& SurrogateData::variables_data() const
{ return sdRep->varsDataIter->second; }


inline const SDVArray& SurrogateData::
variables_data(const ActiveKey& key) const
{ return sdRep->varsData[key]; }


inline SDVArray& SurrogateData::variables_data()
{ return sdRep->varsDataIter->second; }


inline const RealVector& SurrogateData::continuous_variables(size_t i) const
{ return sdRep->varsDataIter->second[i].continuous_variables(); }


inline void SurrogateData::response_data(const SDRArray& sdr_array)
{ sdRep->respDataIter->second = sdr_array; }


inline const SDRArray& SurrogateData::response_data() const
{ return sdRep->respDataIter->second; }


inline const SDRArray& SurrogateData::
response_data(const ActiveKey& key) const
{ return sdRep->respData[key]; }


inline SDRArray& SurrogateData::response_data()
{ return sdRep->respDataIter->second; }


inline Real SurrogateData::response_function(size_t i) const
{ return sdRep->respDataIter->second[i].response_function(); }


inline Real SurrogateData::scaled_response_function(size_t i) const
{
  // return (fn - min) / range
  Real range = sdRep->respFnScaling.second, fn = response_function(i);
  return (range > 0.) ? (fn - sdRep->respFnScaling.first) / range : fn;
}


inline const RealRealPair& SurrogateData::response_function_scaling() const
{ return sdRep->respFnScaling; }


inline bool SurrogateData::valid_response_scaling() const
{ return (sdRep->respFnScaling.second > 0.); } // non-zero range


inline const std::map<ActiveKey, SDVArray>& SurrogateData::
variables_data_map() const
{ return sdRep->varsData; }


inline const std::map<ActiveKey, SDVArray>& SurrogateData::
filtered_variables_data_map(short mode) const
{
  const std::map<ActiveKey, SDVArray>& vars_map = sdRep->varsData;
  std::map<ActiveKey, SDVArray>&  filt_vars_map = sdRep->filteredVarsData;
  std::map<ActiveKey, SDVArray>::const_iterator cit;

  filt_vars_map.clear();
  switch (mode) {
  case SINGLETON_FILTER:
    for (cit=vars_map.begin(); cit!=vars_map.end(); ++cit)
      if (cit->first.singleton())
	filt_vars_map.insert(*cit);
    break;
  case AGGREGATED_FILTER: // less restrictive
    for (cit=vars_map.begin(); cit!=vars_map.end(); ++cit)
      if (cit->first.aggregated())
	filt_vars_map.insert(*cit);
    break;
  case RAW_DATA_FILTER:       // more restrictive
    for (cit=vars_map.begin(); cit!=vars_map.end(); ++cit)
      if (cit->first.raw_data())
	filt_vars_map.insert(*cit);
    break;
  case REDUCTION_DATA_FILTER:       // more restrictive
    for (cit=vars_map.begin(); cit!=vars_map.end(); ++cit)
      if (cit->first.reduction_data())
	filt_vars_map.insert(*cit);
    break;
  case RAW_WITH_REDUCTION_DATA_FILTER:       // more restrictive
    for (cit=vars_map.begin(); cit!=vars_map.end(); ++cit)
      if (cit->first.raw_with_reduction_data())
	filt_vars_map.insert(*cit);
    break;
  case NO_FILTER:
    filt_vars_map = vars_map;
    break;
  }
  return filt_vars_map;
}


inline void SurrogateData::
variables_data_map(const std::map<ActiveKey, SDVArray>& vars_map)
{ sdRep->varsData = vars_map; }


inline const std::map<ActiveKey, SDRArray>& SurrogateData::
response_data_map() const
{ return sdRep->respData; }


inline const std::map<ActiveKey, SDRArray>& SurrogateData::
filtered_response_data_map(short mode) const
{
  const std::map<ActiveKey, SDRArray>& resp_map = sdRep->respData;
  std::map<ActiveKey, SDRArray>&  filt_resp_map = sdRep->filteredRespData;
  std::map<ActiveKey, SDRArray>::const_iterator cit;

  filt_resp_map.clear();
  switch (mode) {
  case SINGLETON_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit)
      if (cit->first.singleton())
	filt_resp_map.insert(*cit);
    break;
  case AGGREGATED_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit)
      if (cit->first.aggregated())
	filt_resp_map.insert(*cit);
    break;
  case RAW_DATA_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit)
      if (cit->first.raw_data())
	filt_resp_map.insert(*cit);
    break;
  case REDUCTION_DATA_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit)
      if (cit->first.reduction_data())
	filt_resp_map.insert(*cit);
    break;
  case RAW_WITH_REDUCTION_DATA_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit)
      if (cit->first.raw_with_reduction_data())
	filt_resp_map.insert(*cit);
    break;
  case NO_FILTER:
    filt_resp_map = resp_map;
    break;
  }
  return filt_resp_map;
}


inline void SurrogateData::
response_data_map(const std::map<ActiveKey, SDRArray>& resp_map)
{ sdRep->respData = resp_map; }


inline const ActiveKey& SurrogateData::
filtered_key(short mode, size_t index) const
{
  const std::map<ActiveKey, SDRArray>& resp_map = sdRep->respData;
  std::map<ActiveKey, SDRArray>::const_iterator cit;

  size_t cntr = 0;
  switch (mode) {
  case SINGLETON_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit) {
      const ActiveKey& key = cit->first;
      if (key.singleton()) {
	if (cntr == index) return key;
	else               ++cntr;
      }
    }
    break;
  case AGGREGATED_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit) {
      const ActiveKey& key = cit->first;
      if (key.aggregated()) {
	if (cntr == index) return key;
	else               ++cntr;
      }
    }
    break;
  case RAW_DATA_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit) {
      const ActiveKey& key = cit->first;
      if (key.raw_data()) {
	if (cntr == index) return key;
	else               ++cntr;
      }
    }
    break;
  case REDUCTION_DATA_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit) {
      const ActiveKey& key = cit->first;
      if (key.reduction_data()) {
	if (cntr == index) return key;
	else               ++cntr;
      }
    }
    break;
  case RAW_WITH_REDUCTION_DATA_FILTER:
    for (cit=resp_map.begin(); cit!=resp_map.end(); ++cit) {
      const ActiveKey& key = cit->first;
      if (key.raw_with_reduction_data()) {
	if (cntr == index) return key;
	else               ++cntr;
      }
    }
    break;
  case NO_FILTER:
    cit = resp_map.begin();
    std::advance(cit, index);
    return cit->first;
    break;
  }

  PCerr << "Error: index not reached in SurrogateData::filtered_key(mode,index)"
	<< std::endl;
  abort_handler(-1);
  return resp_map.begin()->first; // dummy return for compiler
}


inline void SurrogateData::push_back(const SurrogateDataVars& sdv)
{ sdRep->varsDataIter->second.push_back(sdv); }


inline void SurrogateData::push_back(const SurrogateDataResp& sdr)
{ sdRep->respDataIter->second.push_back(sdr); }


inline void SurrogateData::
push_back(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr)
{
  sdRep->varsDataIter->second.push_back(sdv);
  sdRep->respDataIter->second.push_back(sdr);
}


inline void SurrogateData::push_back(int eval_id)
{ sdRep->dataIdsIter->second.push_back(eval_id); }


inline void SurrogateData::
pop_back(size_t num_pop, SDVArray& sdv_array, SDRArray& sdr_array)
{
  size_t len = std::min(sdv_array.size(), sdr_array.size());
  if (len < num_pop) {
    PCerr << "Error: insufficient size (" << len << ") for pop_back("
	  << num_pop << ")." << std::endl;
    abort_handler(-1);
  }

  // could also just use resize(), but retain symmetry with pop_front()
  size_t start = len - num_pop;
  SDVArray::iterator v_it = sdv_array.begin() + start;
  SDRArray::iterator r_it = sdr_array.begin() + start;
  sdv_array.erase(v_it, sdv_array.end());
  sdr_array.erase(r_it, sdr_array.end());
}


inline void SurrogateData::pop_back(size_t num_pop)
{
  size_t start_pts = points(); // count prior to pop

  const ActiveKey& key = sdRep->activeKey;
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    pop_back(num_pop, sdRep->varsDataIter->second, sdRep->respDataIter->second);
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      pop_back(num_pop, sdRep->varsData[key_k], sdRep->respData[key_k]);
    }
  }

  size_t a_index = retrieve_anchor_index(key, false);
  if (a_index != _NPOS && a_index >= start_pts - num_pop)// popped pt was anchor
    clear_anchor_index(key);
}


inline void SurrogateData::
pop_front(size_t num_pop, SDVArray& sdv_array, SDRArray& sdr_array)
{
  size_t len = std::min(sdv_array.size(), sdr_array.size());
  if (len < num_pop) {
    PCerr << "Error: insufficient size (" << len << ") for pop_front("
	  << num_pop << ")." << std::endl;
    abort_handler(-1);
  }

  SDVArray::iterator v_it = sdv_array.begin();
  SDRArray::iterator r_it = sdr_array.begin();
  sdv_array.erase(v_it, v_it + num_pop);
  sdr_array.erase(r_it, r_it + num_pop);
}


inline void SurrogateData::pop_front(size_t num_pop)
{
  const ActiveKey& key = sdRep->activeKey;
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    pop_front(num_pop, sdRep->varsDataIter->second,
	      sdRep->respDataIter->second);
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      pop_front(num_pop, sdRep->varsData[key_k], sdRep->respData[key_k]);
    }
  }

  size_t a_index = retrieve_anchor_index(key, false);
  if (a_index < num_pop)     // anchor has been popped
    clear_anchor_index(key);
  else if (a_index != _NPOS) // anchor (still) exists, decrement its index
    anchor_index(a_index - num_pop, key);
}


inline void SurrogateData::
history_target(size_t target, SDVArray& sdv_array, SDRArray& sdr_array,
	       std::map<ActiveKey, size_t>::iterator anchor_it)
{
  size_t len = std::min(sdv_array.size(), sdr_array.size());
  if (len > target) {
    // erase oldest data (pop from front of array)
    size_t num_pop = len - target;
    pop_front(num_pop, sdv_array, sdr_array);

    // Manage anchorIndex if anchor is popped
    std::map<ActiveKey, size_t>& anchor_index = sdRep->anchorIndex;
    if (anchor_it != anchor_index.end() && anchor_it->second != _NPOS) {
      if (anchor_it->second < num_pop) // a popped point was anchor
	anchor_index.erase(anchor_it);//(key_k);
      else // anchor still exists, decrement its index
	anchor_it->second -= num_pop;
    }
  }
}


inline void SurrogateData::
history_target(size_t target, const ActiveKey& key)
{
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    history_target(target, sdRep->varsData[key], sdRep->respData[key],
		   sdRep->anchorIndex.find(key));
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      history_target(target, sdRep->varsData[key_k], sdRep->respData[key_k],
		     sdRep->anchorIndex.find(key_k));
    }
  }
}


inline void SurrogateData::
pop(SDVArray& sdv_array, SDRArray& sdr_array, IntArray& data_ids,
    std::map<ActiveKey, SizetArray>::iterator pop_cnt_it,
    SDVArrayDeque& popped_sdv_arrays, SDRArrayDeque& popped_sdr_arrays,
    IntArrayDeque& popped_ids, SizetShortMap& failed_resp, bool save_data)
{
  size_t num_pts = std::min(sdv_array.size(), sdr_array.size());

  // harden logic for case of an empty SurrogateData (e.g.,
  // distinct discrepancy at level 0)
  if (pop_cnt_it == sdRep->popCountStack.end()) {
    if (num_pts == 0)
      return; // assume inactive SurrogateData -> ignore pop request
    else {
      PCerr << "\nError: active count stack not found in SurrogateData::pop()"
	    << std::endl;
      abort_handler(-1);
    }
  }

  SizetArray& pop_count_stack = pop_cnt_it->second;
  if (pop_count_stack.empty()) {
    PCerr << "\nError: empty count stack in SurrogateData::pop()" << std::endl;
    abort_handler(-1);
  }
  size_t num_pop_pts = pop_count_stack.back();
  if (num_pop_pts) {
    if (num_pts < num_pop_pts) {
      PCerr << "Error: pop count (" << num_pop_pts << ") exceeds data size ("
	    << num_pts << ") in SurrogateData::pop(size_t)" << std::endl;
      abort_handler(-1);
    }
    if (save_data) {
      // append empty arrays and then update them in place
      popped_sdv_arrays.push_back(SDVArray());
      popped_sdr_arrays.push_back(SDRArray());
      SDVArray& last_popped_sdv_array = popped_sdv_arrays.back();
      SDRArray& last_popped_sdr_array = popped_sdr_arrays.back();
      SDVArray::iterator v_end = sdv_array.end();
      SDRArray::iterator r_end = sdr_array.end();
      last_popped_sdv_array.insert(last_popped_sdv_array.begin(),
				   v_end - num_pop_pts, v_end);
      last_popped_sdr_array.insert(last_popped_sdr_array.begin(),
				   r_end - num_pop_pts, r_end);
    }
    size_t new_size = num_pts - num_pop_pts;
    sdv_array.resize(new_size); sdr_array.resize(new_size);

    // TO DO: prune failedRespData[key] or leave in map ?
    data_checks(sdr_array, failed_resp); // from scratch for now...

    if (!data_ids.empty()) {
      if (save_data) {
	// append empty arrays and then update them in place
	popped_ids.push_back(IntArray());
	IntArray& last_popped_int_array = popped_ids.back();
	IntArray::iterator r_end = data_ids.end();
	last_popped_int_array.insert(last_popped_int_array.begin(),
				     r_end - num_pop_pts, r_end);
      }
      data_ids.resize(new_size);
    }
  }

  pop_count_stack.pop_back();
}


inline void SurrogateData::pop(bool save_data)
{
  const ActiveKey& key = sdRep->activeKey;
  bool agg_key = key.aggregated();
  SDVArrayDeque empty_sdva; SDRArrayDeque empty_sdra; IntArrayDeque empty_ia;
  if (!agg_key || key.reduction_data()) { // process original key
    IntArray& data_ids = sdRep->dataIdsIter->second;
    SDVArrayDeque& popped_vars = (save_data) ?
      sdRep->poppedVarsData[key] : empty_sdva;
    SDRArrayDeque& popped_resp = (save_data) ?
      sdRep->poppedRespData[key] : empty_sdra;
    IntArrayDeque& popped_ids = (save_data && !data_ids.empty()) ?
      sdRep->poppedDataIds[key] : empty_ia;
    pop(sdRep->varsDataIter->second, sdRep->respDataIter->second, data_ids,
	sdRep->popCountStack.find(key), popped_vars, popped_resp, popped_ids,
	sdRep->failedRespData[key], save_data);
  }

  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      IntArray& data_ids = sdRep->dataIdentifiers[key_k];
      SDVArrayDeque& popped_vars = (save_data) ?
	sdRep->poppedVarsData[key_k] : empty_sdva;
      SDRArrayDeque& popped_resp = (save_data) ?
	sdRep->poppedRespData[key_k] : empty_sdra;
      IntArrayDeque& popped_ids = (save_data && !data_ids.empty()) ?
	sdRep->poppedDataIds[key_k] : empty_ia;
      pop(sdRep->varsData[key_k], sdRep->respData[key_k], data_ids,
	  sdRep->popCountStack.find(key_k), popped_vars, popped_resp,
	  popped_ids, sdRep->failedRespData[key_k], save_data);
    }
  }
}


inline void SurrogateData::pop(const ActiveKey& key, bool save_data)
{
  bool agg_key = key.aggregated();
  SDVArrayDeque empty_sdva; SDRArrayDeque empty_sdra; IntArrayDeque empty_ia;
  if (!agg_key || key.reduction_data()) { // process original key
    IntArray& data_ids = sdRep->dataIdentifiers[key];
    SDVArrayDeque& popped_vars = (save_data) ?
      sdRep->poppedVarsData[key] : empty_sdva;
    SDRArrayDeque& popped_resp = (save_data) ?
      sdRep->poppedRespData[key] : empty_sdra;
    IntArrayDeque& popped_ids = (save_data && !data_ids.empty()) ?
      sdRep->poppedDataIds[key] : empty_ia;
    pop(sdRep->varsData[key], sdRep->respData[key], data_ids,
	sdRep->popCountStack.find(key), popped_vars, popped_resp, popped_ids,
	sdRep->failedRespData[key], save_data);
  }

  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      IntArray& data_ids = sdRep->dataIdentifiers[key_k];
      SDVArrayDeque& popped_vars = (save_data) ?
	sdRep->poppedVarsData[key_k] : empty_sdva;
      SDRArrayDeque& popped_resp = (save_data) ?
	sdRep->poppedRespData[key_k] : empty_sdra;
      IntArrayDeque& popped_ids = (save_data && !data_ids.empty()) ?
	sdRep->poppedDataIds[key_k] : empty_ia;
      pop(sdRep->varsData[key_k], sdRep->respData[key_k], data_ids,
	  sdRep->popCountStack.find(key_k), popped_vars, popped_resp,
	  popped_ids, sdRep->failedRespData[key_k], save_data);
    }
  }
}


inline void SurrogateData::
push(SDVArray& sdv_array, SDRArray& sdr_array,
     IntArray& data_ids,  SizetArray& pop_count_stack,
     std::map<ActiveKey, SDVArrayDeque>::iterator pvd_it,
     std::map<ActiveKey, SDRArrayDeque>::iterator prd_it,
     std::map<ActiveKey, IntArrayDeque>::iterator pid_it,
     SizetShortMap& failed_resp, size_t index, bool erase_popped)
{
  size_t num_pts, num_popped
    = (pvd_it != sdRep->poppedVarsData.end() &&
       prd_it != sdRep->poppedRespData.end()) ?
    std::min(pvd_it->second.size(), prd_it->second.size()) : 0;

  // harden logic for case of an empty SurrogateData (e.g.,
  // distinct discrepancy at level 0)
  if (num_popped > index) {
    SDVArrayDeque& popped_sdv_arrays = pvd_it->second;
    SDRArrayDeque& popped_sdr_arrays = prd_it->second;
    SDVArrayDeque::iterator vit = popped_sdv_arrays.begin() + index;
    SDRArrayDeque::iterator rit = popped_sdr_arrays.begin() + index;
    num_pts = std::min(vit->size(), rit->size());

    sdv_array.insert(sdv_array.end(), vit->begin(), vit->end());
    sdr_array.insert(sdr_array.end(), rit->begin(), rit->end());

    // TO DO: update failedRespData[activeKey] ?
    data_checks(sdr_array, failed_resp); // from scratch for now...

    if (erase_popped)
      { popped_sdv_arrays.erase(vit); popped_sdr_arrays.erase(rit); }

    if (pid_it != sdRep->poppedDataIds.end()) { // no error if not tracking ids
      IntArrayDeque& popped_ids = pid_it->second;
      if (index < popped_ids.size()) {
	IntArrayDeque::iterator iit = popped_ids.begin() + index;
	data_ids.insert(data_ids.end(), iit->begin(), iit->end());
	if (erase_popped)
	  popped_ids.erase(iit);
      }
      else { // this is an error
	PCerr << "Error: index (" << index << ") out of bounds (size = "
	      << popped_ids.size() << ") for evaluation id in SurrogateData"
	      << "::push()" << std::endl;
	abort_handler(-1);
      }
    }

    pop_count_stack.push_back(num_pts);
  }
  else if (num_popped) { // not empty
    PCerr << "Error: index out of range for active popped arrays in "
	  << "SurrogateData::push()." << std::endl;
    abort_handler(-1);
  }
  // else ignore push request for empty popped (SurrogateData assumed inactive)
}


inline void SurrogateData::push(size_t index, bool erase_popped)
{
  const ActiveKey& key = sdRep->activeKey;
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    push(sdRep->varsDataIter->second,     sdRep->respDataIter->second,
	 sdRep->dataIdsIter->second,      sdRep->popCountStack[key],
	 sdRep->poppedVarsData.find(key), sdRep->poppedRespData.find(key),
	 sdRep->poppedDataIds.find(key),  sdRep->failedRespData[key],
	 index, erase_popped);
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      ActiveKey& key_k = embedded_keys[k];
      push(sdRep->varsData[key_k],            sdRep->respData[key_k],
	   sdRep->dataIdentifiers[key_k],     sdRep->popCountStack[key_k],
	   sdRep->poppedVarsData.find(key_k), sdRep->poppedRespData.find(key_k),
	   sdRep->poppedDataIds.find(key_k),  sdRep->failedRespData[key_k],
	   index, erase_popped);
    }
  }
}


inline void SurrogateData::
push(const ActiveKey& key, size_t index, bool erase_popped)
{
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    push(sdRep->varsData[key],            sdRep->respData[key],
	 sdRep->dataIdentifiers[key],     sdRep->popCountStack[key],
	 sdRep->poppedVarsData.find(key), sdRep->poppedRespData.find(key),
	 sdRep->poppedDataIds.find(key),  sdRep->failedRespData[key],
	 index, erase_popped);
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      ActiveKey& key_k = embedded_keys[k];
      push(sdRep->varsData[key_k],            sdRep->respData[key_k],
	   sdRep->dataIdentifiers[key_k],     sdRep->popCountStack[key_k],
	   sdRep->poppedVarsData.find(key_k), sdRep->poppedRespData.find(key_k),
	   sdRep->poppedDataIds.find(key_k),  sdRep->failedRespData[key_k],
	   index, erase_popped);
    }
  }
}


inline void SurrogateData::
replace(const SurrogateDataVars& sdv, int id)
{
  std::map<ActiveKey, IntArray>::iterator it
    = sdRep->dataIdentifiers.find(sdRep->activeKey);
  size_t index = (it == sdRep->dataIdentifiers.end()) ? _NPOS :
    find_index(it->second, id);
  // Note: the following logic differs from assign_variables(sdv, index):
  if (index == _NPOS) {
    PCerr << "Error: id lookup failure in SurrogateData::replace()."<<std::endl;
    abort_handler(-1);
  }
  else {
    SDVArray& sdv_array = sdRep->varsDataIter->second;
    if (index >= sdv_array.size()) {
      PCerr << "Error: index out of range in SurrogateData::replace()."
	    << std::endl;
      abort_handler(-1);
    }
    else
      sdv_array[index] = sdv;
  }
}


inline void SurrogateData::
replace(const SurrogateDataResp& sdr, int id)
{
  std::map<ActiveKey, IntArray>::iterator it
    = sdRep->dataIdentifiers.find(sdRep->activeKey);
  size_t index = (it == sdRep->dataIdentifiers.end()) ? _NPOS :
    find_index(it->second, id);
  // Note: the following logic differs from assign_response(sdr, index):
  if (index == _NPOS) {
    PCerr << "Error: id lookup failure in SurrogateData::replace()."<<std::endl;
    abort_handler(-1);
  }
  else {
    SDRArray& sdr_array = sdRep->respDataIter->second;
    if (index >= sdr_array.size()) {
      PCerr << "Error: index out of range in SurrogateData::replace()."
	    << std::endl;
      abort_handler(-1);
    }
    else
      sdr_array[index] = sdr;
  }
}


inline void SurrogateData::
replace(const SurrogateDataVars& sdv, const SurrogateDataResp& sdr, int id)
{
  std::map<ActiveKey, IntArray>::iterator it
    = sdRep->dataIdentifiers.find(sdRep->activeKey);
  size_t index = (it == sdRep->dataIdentifiers.end()) ? _NPOS :
    find_index(it->second, id);
  // Note: the following logic differs from assign_{variables,response}:
  if (index == _NPOS) {
    PCerr << "Error: id lookup failure in SurrogateData::replace()."<<std::endl;
    abort_handler(-1);
  }
  else {
    SDVArray& sdv_array = sdRep->varsDataIter->second;
    SDRArray& sdr_array = sdRep->respDataIter->second;
    if (index >= sdv_array.size() || index >= sdr_array.size()) {
      PCerr << "Error: index out of range in SurrogateData::replace()."
	    << std::endl;
      abort_handler(-1);
    }
    else {
      sdv_array[index] = sdv;
      sdr_array[index] = sdr;
    }
  }
}


inline void SurrogateData::pop_count(size_t count) const
{
  const ActiveKey& key = sdRep->activeKey;
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) // process original key
    sdRep->popCountStack[key].push_back(count);
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k)
      sdRep->popCountStack[embedded_keys[k]].push_back(count);
  }
}


inline size_t SurrogateData::pop_count() const
{
  std::map<ActiveKey, SizetArray>::iterator pop_it
    = sdRep->popCountStack.find(sdRep->activeKey);
  return (pop_it == sdRep->popCountStack.end() || pop_it->second.empty()) ?
    _NPOS : pop_it->second.back();
}


inline const SizetArray& SurrogateData::
pop_count_stack(const ActiveKey& key) const
{ return sdRep->popCountStack[key]; }


inline const SizetArray& SurrogateData::pop_count_stack() const
{ return sdRep->popCountStack[sdRep->activeKey]; }


inline void SurrogateData::pop_count_stack(const SizetArray& pop_count) const
{ sdRep->popCountStack[sdRep->activeKey] = pop_count; }


inline void SurrogateData::
pop_count_stack_map(const std::map<ActiveKey, SizetArray>& pcs_map) const
{ sdRep->popCountStack = pcs_map; }


inline const std::map<ActiveKey, SizetArray>& SurrogateData::
pop_count_stack_map() const
{ return sdRep->popCountStack; }


inline bool SurrogateData::anchor() const
{
  std::map<ActiveKey, size_t>::iterator anchor_it
    = sdRep->anchorIndex.find(sdRep->activeKey);
  return (anchor_it != sdRep->anchorIndex.end() && anchor_it->second != _NPOS);
}


inline size_t SurrogateData::response_size() const
{
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  size_t i, data_size = 0, num_resp = sdr_array.size(), nh;
  short active_bits;
  for (i=0; i<num_resp; ++i) {
    const SurrogateDataResp& sdr = sdr_array[i];
    active_bits = sdr.active_bits();
    if (active_bits & 1) ++data_size;
    if (active_bits & 2) data_size += sdr.response_gradient().length();
    if (active_bits & 4) {
      nh = sdr.response_hessian().numRows();
      if (nh) data_size += nh * (nh + 1) / 2;
    }
  }
  return data_size;
}


inline size_t SurrogateData::failed_response_size() const
{
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  size_t failed_size = 0, nh; short fail_bits;
  const SizetShortMap& failed_resp = failed_response_data();
  for (StShMCIter cit=failed_resp.begin(); cit!=failed_resp.end(); ++cit) {
    fail_bits = cit->second;
    const SurrogateDataResp& sdr = sdr_array[cit->first];
    if (fail_bits & 1) ++failed_size;
    if (fail_bits & 2) failed_size += sdr.response_gradient().length();
    if (fail_bits & 4) {
      nh = sdr.response_hessian().numRows();
      if (nh) failed_size += nh * (nh + 1) / 2;
    }
  }
  return failed_size;
}


inline size_t SurrogateData::active_response_size() const
{ return response_size() - failed_response_size(); }


inline size_t SurrogateData::popped_sets(const ActiveKey& key) const
{
  return std::min(sdRep->poppedVarsData[key].size(),
		  sdRep->poppedRespData[key].size());
}


inline size_t SurrogateData::popped_sets() const
{ return popped_sets(sdRep->activeKey); }


inline size_t SurrogateData::num_gradient_variables() const
{
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  return (sdr_array.empty()) ? 0 : sdr_array[0].response_gradient().length();
}


inline size_t SurrogateData::num_hessian_variables() const
{
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  return (sdr_array.empty()) ? 0 : sdr_array[0].response_hessian().numRows();
}


inline size_t SurrogateData::num_derivative_variables() const
{
  size_t num_grad_vars = num_gradient_variables();
  if (num_grad_vars) return num_grad_vars;           // precedence
  else               return num_hessian_variables(); // fall-back
}


inline void SurrogateData::clear_response_function_scaling() const
{ sdRep->respFnScaling/*[activeKey]*/ = RealRealPair(0., 0.); }


inline void SurrogateData::compute_response_function_scaling() const
{
  // compute scale factors for mapping active response functions to [0,1]
  // current uses do not require key management
  const SDRArray& sdr_array = sdRep->respDataIter->second;
  size_t i, pts = sdr_array.size();
  if (pts <= 1) { clear_response_function_scaling(); return; }

  Real fn, min_fn, max_fn, range;
  min_fn = max_fn = sdr_array[0].response_function();
  for (i=1; i<pts; ++i) {
    fn = sdr_array[i].response_function();
    if      (fn > max_fn) max_fn = fn;
    else if (fn < min_fn) min_fn = fn;
  }
  range = max_fn - min_fn;
  sdRep->respFnScaling/*[activeKey]*/ = RealRealPair(min_fn, range);
}


inline void SurrogateData::assign_response_function_shift(Real shift) const
{ sdRep->respFnScaling/*[activeKey]*/ = RealRealPair(-shift, 1.); }


inline void SurrogateData::
response_check(const SurrogateDataResp& sdr, short& failed_data) const
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


inline void SurrogateData::
data_checks(const SDRArray& resp_data, SizetShortMap& failed_resp) const
{
  failed_resp.clear();
  size_t i, num_resp = resp_data.size();  short failed_data;
  for (i=0; i<num_resp; ++i) {
    response_check(resp_data[i], failed_data);
    if (failed_data)
      failed_resp[i] = failed_data; // include in map
  }

#ifdef DEBUG
  PCout << "failedRespData:\n";
  for (SizetShortMap::iterator it=failed_resp.begin();
       it!=failed_resp.end(); ++it)
    PCout << "index: " << std::setw(6) << it->first
	  << " data: " << it->second << '\n';
#endif // DEBUG
}


inline void SurrogateData::data_checks() const
{
  data_checks(sdRep->respDataIter->second,
	      sdRep->failedRespData[sdRep->activeKey]);
}


inline short SurrogateData::failed_anchor_data() const
{
  // returns 3-bit representation (0-7) of failed response data
  std::map<ActiveKey, SizetShortMap>::const_iterator cit1
    = sdRep->failedRespData.find(sdRep->activeKey);
  if (cit1 == sdRep->failedRespData.end()) return 0;
  else {
    const SizetShortMap& failed_resp_data = cit1->second;
    SizetShortMap::const_iterator cit2 = failed_resp_data.find(anchor_index());
    return (cit2 == failed_resp_data.end()) ? 0 : cit2->second;
  }
}


inline void SurrogateData::
failed_response_data(const SizetShortMap& fail_data) const
{ sdRep->failedRespData[sdRep->activeKey] = fail_data; }


inline const SizetShortMap& SurrogateData::
failed_response_data(const ActiveKey& key) const
{ return sdRep->failedRespData[key]; }


inline const SizetShortMap& SurrogateData::failed_response_data() const
{ return failed_response_data(sdRep->activeKey); }


inline void SurrogateData::
failed_response_data_map(const std::map<ActiveKey, SizetShortMap>& fail_resp)
  const
{ sdRep->failedRespData = fail_resp; }


inline const std::map<ActiveKey, SizetShortMap>& SurrogateData::
failed_response_data_map() const
{ return sdRep->failedRespData; }


inline const SDVArrayDeque& SurrogateData::popped_variables() const
{
  std::map<ActiveKey, SDVArrayDeque>& pop_vars = sdRep->poppedVarsData;
  const ActiveKey& key = sdRep->activeKey;
  std::map<ActiveKey, SDVArrayDeque>::iterator cit = pop_vars.find(key);
  if (cit == pop_vars.end()) {
    std::pair<ActiveKey, SDVArrayDeque> sdv_pair(key, SDVArrayDeque());
    cit = pop_vars.insert(sdv_pair).first; // create empty array for key
  }
  return cit->second;
}


inline void SurrogateData::popped_variables(const SDVArrayDeque& popped_vars)
{ sdRep->poppedVarsData[sdRep->activeKey] = popped_vars; }


inline const SDRArrayDeque& SurrogateData::popped_response() const
{
  std::map<ActiveKey, SDRArrayDeque>& pop_resp = sdRep->poppedRespData;
  const ActiveKey& key = sdRep->activeKey;
  std::map<ActiveKey, SDRArrayDeque>::iterator cit = pop_resp.find(key);
  if (cit == pop_resp.end()) {
    std::pair<ActiveKey, SDRArrayDeque> sdr_pair(key, SDRArrayDeque());
    cit = pop_resp.insert(sdr_pair).first; // create empty array for key
  }
  return cit->second;
}


inline void SurrogateData::
popped_response(const SDRArrayDeque& popped_resp)
{ sdRep->poppedRespData[sdRep->activeKey] = popped_resp; }


inline const std::map<ActiveKey, SDVArrayDeque>& SurrogateData::
popped_variables_map() const
{ return sdRep->poppedVarsData; }


inline void SurrogateData::
popped_variables_map(const std::map<ActiveKey, SDVArrayDeque>& popped_vars)
{ sdRep->poppedVarsData = popped_vars; }


inline const std::map<ActiveKey, SDRArrayDeque>& SurrogateData::
popped_response_map() const
{ return sdRep->poppedRespData; }


inline void SurrogateData::
popped_response_map(const std::map<ActiveKey, SDRArrayDeque>& popped_resp)
{ sdRep->poppedRespData = popped_resp; }


inline void SurrogateData::
copy(const SurrogateData& sd, short sdv_mode, short sdr_mode) const
{
  if (sdv_mode == DEEP_COPY) {
    size_t i, j, num_pts, num_sdva;
    const std::map<ActiveKey, SDVArray>& vars_map = sd.variables_data_map();
    std::map<ActiveKey, SDVArray>::const_iterator v_cit;
    std::map<ActiveKey, SDVArray>& new_vars_map = sdRep->varsData;
    new_vars_map.clear();
    for (v_cit = vars_map.begin(); v_cit != vars_map.end(); ++v_cit) {
      const SDVArray& sdv_array = v_cit->second;
      num_pts = sdv_array.size();
      SDVArray new_sdv_array(num_pts);
      for (i=0; i<num_pts; ++i)
	new_sdv_array[i] = sdv_array[i].copy();
      new_vars_map[v_cit->first] = sdv_array;
    }

    const std::map<ActiveKey, SDVArrayDeque>& popped_vars_map
      = sd.popped_variables_map();
    std::map<ActiveKey, SDVArrayDeque>::const_iterator v2_cit;
    std::map<ActiveKey, SDVArrayDeque>& new_popped_vars_map
      = sdRep->poppedVarsData;
    new_popped_vars_map.clear();
    for (v2_cit = popped_vars_map.begin(); v2_cit != popped_vars_map.end();
	 ++v2_cit) {
      const SDVArrayDeque& sdv_2d_array = v2_cit->second;
      num_sdva = sdv_2d_array.size();
      SDVArrayDeque new_popped_sdv_2d(num_sdva);
      for (i=0; i<num_sdva; ++i) {
	const SDVArray& sdva_i = sdv_2d_array[i];
	num_pts = sdva_i.size();
	new_popped_sdv_2d[i].resize(num_pts);
	for (j=0; j<num_pts; ++j)
	  new_popped_sdv_2d[i][j] = sdva_i[j].copy();
      }
      new_popped_vars_map[v2_cit->first] = new_popped_sdv_2d;
    }
  }
  else { // shallow SDV copies based on operator=
    sdRep->varsData       = sd.variables_data_map();
    sdRep->poppedVarsData = sd.popped_variables_map();
  }

  if (sdr_mode == DEEP_COPY) {
    size_t i, j, num_pts, num_sdra;
    const std::map<ActiveKey, SDRArray>& resp_map = sd.response_data_map();
    std::map<ActiveKey, SDRArray>::const_iterator r_cit;
    std::map<ActiveKey, SDRArray>& new_resp_map = sdRep->respData;
    new_resp_map.clear();
    for (r_cit = resp_map.begin(); r_cit != resp_map.end(); ++r_cit) {
      const SDRArray& sdr_array = r_cit->second;
      num_pts = sdr_array.size();
      SDRArray new_sdr_array(num_pts);
      for (i=0; i<num_pts; ++i)
	new_sdr_array[i] = sdr_array[i].copy();
      new_resp_map[r_cit->first] = new_sdr_array;
    }

    const std::map<ActiveKey, SDRArrayDeque>& popped_resp_map
      = sd.popped_response_map();
    std::map<ActiveKey, SDRArrayDeque>::const_iterator r2_cit;
    std::map<ActiveKey, SDRArrayDeque>& new_popped_resp_map
      = sdRep->poppedRespData;
    new_popped_resp_map.clear();
    for (r2_cit = popped_resp_map.begin(); r2_cit != popped_resp_map.end();
	 ++r2_cit) {
      const SDRArrayDeque& sdr_2d_array = r2_cit->second;
      num_sdra = sdr_2d_array.size();
      SDRArrayDeque new_popped_sdr_2d(num_sdra);
      for (i=0; i<num_sdra; ++i) {
	const SDRArray& sdra_i = sdr_2d_array[i];
	num_pts = sdra_i.size();
	new_popped_sdr_2d[i].resize(num_pts);
	for (j=0; j<num_pts; ++j)
	  new_popped_sdr_2d[i][j] = sdra_i[j].copy();
      }
      new_popped_resp_map[r2_cit->first] = new_popped_sdr_2d;
    }
  }
  else { // shallow SDR copies based on operator=
    sdRep->respData       = sd.response_data_map();
    sdRep->poppedRespData = sd.popped_response_map();
  }

  active_key(sd.active_key());
  pop_count_stack_map(sd.pop_count_stack_map());
  anchor_index_map(sd.anchor_index_map());
  failed_response_data_map(sd.failed_response_data_map());
}


inline SurrogateData SurrogateData::copy(short sdv_mode, short sdr_mode) const
{
  SurrogateData sd;
  sd.copy(*this, sdv_mode, sdr_mode);
  return sd;
}


inline void SurrogateData::
copy_active(const SurrogateData& sd, short sdv_mode, short sdr_mode) const
{
  const ActiveKey& key = sd.active_key();
  if (sdRep->activeKey != key) active_key(key);
  copy_active_sdv(sd, sdv_mode);
  copy_active_pop_sdv(sd, sdv_mode);
  copy_active_sdr(sd, sdr_mode);
  copy_active_pop_sdr(sd, sdr_mode);

  // deep copies of bookkeeping arrays
  anchor_index(sd.anchor_index());
  pop_count_stack(sd.pop_count_stack());
  failed_response_data(sd.failed_response_data());
}


inline void SurrogateData::size_active_sdv(const SurrogateData& sd) const
{
  const ActiveKey& key = sd.active_key();
  if (sdRep->activeKey != key) active_key(key);

  const SDVArray& sdv_array = sd.variables_data();
  size_t num_pts = sdv_array.size();
  if (num_pts) {
    const SurrogateDataVars& sdv0 = sdv_array[0];
    size_t i, num_cv = sdv0.cv(), num_div = sdv0.div(), num_drv = sdv0.drv();
    SDVArray& new_sdv_array = sdRep->varsDataIter->second;
    new_sdv_array.resize(num_pts);
    for (i=0; i<num_pts; ++i)
      new_sdv_array[i] = SurrogateDataVars(num_cv, num_div, num_drv);
  }
}


inline void SurrogateData::
copy_active_sdv(const SurrogateData& sd, short sdv_mode) const
{
  const ActiveKey& key = sd.active_key();
  if (sdRep->activeKey != key) active_key(key);

  SDVArray& new_sdv_array = sdRep->varsDataIter->second;
  if (sdv_mode == DEEP_COPY) {
    const SDVArray& sdv_array = sd.variables_data();
    size_t i, num_pts = sdv_array.size();
    new_sdv_array.resize(num_pts);
    for (i=0; i<num_pts; ++i)
      new_sdv_array[i] = sdv_array[i].copy();
  }
  else // shallow SDV copies based on operator=
    new_sdv_array = sd.variables_data();
}


inline void SurrogateData::
copy_active_pop_sdv(const SurrogateData& sd, short sdv_mode) const
{
  const ActiveKey& key = sd.active_key();
  if (sdRep->activeKey != key) active_key(key);

  SDVArrayDeque& new_pop_sdv_array = sdRep->poppedVarsData[key];
  if (sdv_mode == DEEP_COPY) {
    const SDVArrayDeque& pop_sdv_array = sd.popped_variables();
    size_t i, j, num_pts, num_sdva = pop_sdv_array.size();
    new_pop_sdv_array.resize(num_sdva);
    for (i=0; i<num_sdva; ++i) {
      const SDVArray& sdva_i = pop_sdv_array[i];
      num_pts = sdva_i.size();
      new_pop_sdv_array[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	new_pop_sdv_array[i][j] = sdva_i[j].copy();
    }
  }
  else // shallow SDV copies based on operator=
    new_pop_sdv_array = sd.popped_variables();
}


inline void SurrogateData::size_active_sdr(const SDRArray& sdr_array) const
{
  size_t num_pts = sdr_array.size();
  SDRArray& new_sdr_array = sdRep->respDataIter->second;
  new_sdr_array.resize(num_pts);
  if (num_pts) {
    const SurrogateDataResp& sdr0 = sdr_array[0];
    short bits = sdr0.active_bits(); // assume homogeneity in deriv data
    size_t i, num_deriv_v = sdr0.derivative_variables();
    for (i=0; i<num_pts; ++i)
      new_sdr_array[i] = SurrogateDataResp(bits, num_deriv_v);
  }
}


inline void SurrogateData::size_active_sdr(const SurrogateData& sd) const
{
  const ActiveKey& key = sd.active_key();
  if (sdRep->activeKey != key) active_key(key);
  size_active_sdr(sd.response_data());
}


inline void SurrogateData::
copy_active_sdr(const SurrogateData& sd, short sdr_mode) const
{
  const ActiveKey& key = sd.active_key();
  if (sdRep->activeKey != key) active_key(key);

  SDRArray& new_sdr_array = sdRep->respDataIter->second;
  if (sdr_mode == DEEP_COPY) {
    const SDRArray& sdr_array = sd.response_data();
    size_t i, num_pts = sdr_array.size();
    new_sdr_array.resize(num_pts);
    for (i=0; i<num_pts; ++i)
      new_sdr_array[i] = sdr_array[i].copy();
  }
  else // shallow SDR copies based on operator=
    new_sdr_array = sd.response_data();
}


inline void SurrogateData::
copy_active_pop_sdr(const SurrogateData& sd, short sdr_mode) const
{
  const ActiveKey& key = sd.active_key();
  if (sdRep->activeKey != key) active_key(key);

  SDRArrayDeque& new_pop_sdr_array = sdRep->poppedRespData[key];
  if (sdr_mode == DEEP_COPY) {
    const SDRArrayDeque& pop_sdr_array = sd.popped_response();
    size_t i, j, num_pts, num_sdra = pop_sdr_array.size();
    new_pop_sdr_array.resize(num_sdra);
    for (i=0; i<num_sdra; ++i) {
      const SDRArray& sdra_i = pop_sdr_array[i];
      num_pts = sdra_i.size();
      new_pop_sdr_array[i].resize(num_pts);
      for (j=0; j<num_pts; ++j)
	new_pop_sdr_array[i][j] = sdra_i[j].copy();
    }
  }
  else // shallow SDR copies based on operator=
    new_pop_sdr_array = sd.popped_response();
}


inline void SurrogateData::resize(size_t new_pts)
{
  // new SDV/SDR are empty (no rep)
  sdRep->varsDataIter->second.resize(new_pts);
  sdRep->respDataIter->second.resize(new_pts);
}


inline void SurrogateData::resize(size_t new_pts, short bits, size_t num_vars)
{
  size_t i, pts = points();
  SDVArray& sdv_array = sdRep->varsDataIter->second;
  SDRArray& sdr_array = sdRep->respDataIter->second;

  sdv_array.resize(new_pts);  sdr_array.resize(new_pts);
  for (i=pts; i<new_pts; ++i) {
    // minimal ctors that define a rep and size arrays
    sdv_array[i].create_rep(num_vars);
    sdr_array[i].create_rep(bits, num_vars);
  }
}


inline void SurrogateData::clear_active_data()
{
  const ActiveKey& key = sdRep->activeKey;
  /*
  // Too aggressive due to DataFitSurrModel::build_approximation() call to
  // approxInterface.clear_current_active_data();
  sdRep->varsData.erase(key);
  sdRep->varsDataIter = sdRep->varsData.end();
  sdRep->respData.erase(key);
  sdRep->respDataIter = sdRep->respData.end();
  sdRep->anchorIndex.erase(key);
  sdRep->failedRespData.erase(key);
  */

  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) { // process original key
    // Retain valid {varsData,RespData,dataIds}Iter when clearing data:
    sdRep->varsDataIter->second.clear();
    sdRep->respDataIter->second.clear();
    sdRep->dataIdsIter->second.clear();
    // anchorIndex and failedRespData can be pruned:
    sdRep->anchorIndex.erase(key);   //sdRep->anchorIndex[key] = _NPOS;
    sdRep->failedRespData.erase(key);//sdRep->failedRespData[key].clear();
  }

  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      // clear instead of erase
      sdRep->varsData[key_k].clear();  sdRep->respData[key_k].clear();
      sdRep->dataIdentifiers[key_k].clear();
      // can erase
      sdRep->anchorIndex.erase(key_k); sdRep->failedRespData.erase(key_k);
    }
  }
}


inline void SurrogateData::clear_active_data(const ActiveKey& key)
{
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) { // process original key
    // clear instead of erase
    sdRep->varsData[key].clear();  sdRep->respData[key].clear();
    sdRep->dataIdentifiers[key].clear();
    // can erase
    sdRep->anchorIndex.erase(key); sdRep->failedRespData.erase(key);
  }
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      // clear instead of erase
      sdRep->varsData[key_k].clear();  sdRep->respData[key_k].clear();
      sdRep->dataIdentifiers[key_k].clear();
      // can erase
      sdRep->anchorIndex.erase(key_k); sdRep->failedRespData.erase(key_k);
    }
  }
}


inline void SurrogateData::clear_inactive_data()
{
  // Checking each of the traversed keys against aggregate + embedded keys is
  // inefficient, so instead rebuild maps with only active + embedded sets
  std::map<ActiveKey, SDVArray> new_vd;
  std::map<ActiveKey, SDRArray> new_rd;
  std::map<ActiveKey, IntArray> new_di;
  std::map<ActiveKey, size_t>   new_ai;
  std::map<ActiveKey, SizetShortMap> new_frd;
  std::map<ActiveKey, SDVArray>::iterator vit;
  std::map<ActiveKey, SDRArray>::iterator rit;
  std::map<ActiveKey, IntArray>::iterator dit;
  std::map<ActiveKey, size_t>::iterator ait;
  std::map<ActiveKey, SizetShortMap>::iterator fit;
  const ActiveKey& key = sdRep->activeKey;
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) { // process original key
    new_vd.insert(*sdRep->varsDataIter);
    new_rd.insert(*sdRep->respDataIter);
    new_di.insert(*sdRep->dataIdsIter);
    ait = sdRep->anchorIndex.find(key);
    if (ait != sdRep->anchorIndex.end()) new_ai.insert(*ait);
    fit = sdRep->failedRespData.find(key);
    if (fit != sdRep->failedRespData.end()) new_frd.insert(*fit);
  }
  if (agg_key && key.raw_data()) {
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      vit = sdRep->varsData.find(key_k);
      if (vit != sdRep->varsData.end()) new_vd.insert(*vit);
      rit = sdRep->respData.find(key_k);
      if (rit != sdRep->respData.end()) new_rd.insert(*rit);
      dit = sdRep->dataIdentifiers.find(key_k);
      if (dit != sdRep->dataIdentifiers.end()) new_di.insert(*dit);
      ait = sdRep->anchorIndex.find(key_k);
      if (ait != sdRep->anchorIndex.end()) new_ai.insert(*ait);
      fit = sdRep->failedRespData.find(key_k);
      if (fit != sdRep->failedRespData.end()) new_frd.insert(*fit);
    }
  }
  sdRep->varsData        = new_vd;  sdRep->respData    = new_rd;
  sdRep->dataIdentifiers = new_di;  sdRep->anchorIndex = new_ai;
  sdRep->failedRespData  = new_frd;

  /*
  // Preserves active key but not embedded keys extracted from active key,
  // so is not currently the complement of clear_active_data().

  std::map<ActiveKey, SDVArray>::iterator vd_it = sdRep->varsData.begin();
  std::map<ActiveKey, SDRArray>::iterator rd_it = sdRep->respData.begin();
  while (vd_it != sdRep->varsData.end())
    if (vd_it == sdRep->varsDataIter) // preserve active
      { ++vd_it; ++rd_it; }
    else {                            //  clear inactive
      const ActiveKey& key = vd_it->first;
      // anchorIndex and failedRespData can be pruned:
      sdRep->anchorIndex.erase(key);    // if it exists
      sdRep->failedRespData.erase(key); // if it exists
      // Be more conservative with clearing {vars,resp}Data:
      vd_it->second.clear(); ++vd_it;
      rd_it->second.clear(); ++rd_it;
      // Too aggressive. Note: postfix increments manage iterator invalidations
      //sdRep->varsData.erase(vd_it++); sdRep->respData.erase(rd_it++);
    }
  */
}


inline void SurrogateData::clear_filtered()
{
  sdRep->filteredVarsData.clear();
  sdRep->filteredRespData.clear();
}


inline void SurrogateData::clear_active_popped(const ActiveKey& key)
{
  bool agg_key = key.aggregated();
  if (!agg_key || key.reduction_data()) {// process original key
    // can erase as will be recreated in pop() if needed
    sdRep->poppedVarsData.erase(key);//sdRep->poppedVarsData[key].clear();
    sdRep->poppedRespData.erase(key);//sdRep->poppedRespData[key].clear();
    sdRep->poppedDataIds.erase(key); //sdRep->poppedDataIds[key].clear();
    sdRep->popCountStack.erase(key); //sdRep->popCountStack[key].clear();
  }
  if (agg_key && key.raw_data()) { // enumerate embedded keys
    std::vector<ActiveKey> embedded_keys;
    key.extract_keys(embedded_keys);
    size_t k, num_k = embedded_keys.size();
    for (k=0; k<num_k; ++k) {
      const ActiveKey& key_k = embedded_keys[k];
      sdRep->poppedVarsData.erase(key_k);
      sdRep->poppedRespData.erase(key_k);
      sdRep->poppedDataIds.erase(key_k);
      sdRep->popCountStack.erase(key_k);
    }
  }
}


inline void SurrogateData::clear_active_popped()
{ clear_active_popped(sdRep->activeKey); }


inline void SurrogateData::clear_popped()
{
  sdRep->poppedVarsData.clear();
  sdRep->poppedRespData.clear();
  sdRep->poppedDataIds.clear();
  sdRep->popCountStack.clear();
}


inline void SurrogateData::clear_data(bool initialize)
{
  sdRep->varsData.clear();
  sdRep->respData.clear();
  sdRep->dataIdentifiers.clear();
  clear_filtered();

  sdRep->anchorIndex.clear();
  sdRep->failedRespData.clear();

  // dataIdsIter must be reset in either case below (used at beginning of
  // update_active_iterators())
  sdRep->dataIdsIter = sdRep->dataIdentifiers.end();

  if (initialize) // preserve activeKey and restore to initialization state
    sdRep->update_active_iterators();
  else {
    sdRep->activeKey.clear();
    sdRep->varsDataIter = sdRep->varsData.end();
    sdRep->respDataIter = sdRep->respData.end();
  }
}


inline void SurrogateData::clear_all(bool initialize)
{ clear_data(initialize); clear_popped(); }


inline void SurrogateData::clear_all_active()
{ clear_active_data(); clear_active_popped(); }


inline void SurrogateData::clear_all_active(const ActiveKey& key)
{ clear_active_data(key); clear_active_popped(key); }


inline std::shared_ptr<SurrogateDataRep> SurrogateData::data_rep() const
{ return sdRep; }


inline bool SurrogateData::is_null() const
{ return (sdRep) ? false : true; }

} // namespace Pecos

#endif
