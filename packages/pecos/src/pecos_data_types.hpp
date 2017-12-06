/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef PECOS_DATA_TYPES_H
#define PECOS_DATA_TYPES_H

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  // HAVE_CONFIG_H is STILL set in Dakota/src (EVEN IN THE CMAKE BUILD!) so
  // use a "disable config header" conditional to help manage the transition
  #include "pecos_config.h"
#endif // HAVE_CONFIG_H

#include "pecos_global_defs.hpp"

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"

#include "boost/multi_array.hpp"
#include "boost/dynamic_bitset.hpp"

#include <algorithm>  // for std::find
#include <complex>
#include <map>
#include <set>
#include <string>
#include <vector>


namespace Pecos {

// avoid problems with circular dependencies by using fwd declarations
//class BasisFunction;
class SurrogateDataVars;
class SurrogateDataResp;

// -----------------------------------
// Aliases for fundamental data types:
// -----------------------------------
typedef double Real;

// --------
// Strings:
// --------
typedef std::string String;


// --------------------------------
// Numerical arrays (serial dense):
// --------------------------------
typedef Teuchos::SerialDenseVector<int, Real>                RealVector;
typedef Teuchos::SerialDenseVector<int, int>                 IntVector;
typedef Teuchos::SerialDenseVector<int, std::complex<Real> > ComplexVector;
typedef Teuchos::SerialDenseMatrix<int, Real>                RealMatrix;
typedef Teuchos::SerialDenseMatrix<int, int>                 IntMatrix;
typedef Teuchos::SerialSymDenseMatrix<int, Real>             RealSymMatrix;


// ---------------------------------
// Numerical solvers (serial dense):
// ---------------------------------
typedef Teuchos::SerialDenseSolver<int, Real>    RealSolver;
typedef Teuchos::SerialSpdDenseSolver<int, Real> RealSpdSolver;


// ---------------------------------------
// Admin/bookkeeping arrays (serial only):
// ---------------------------------------
typedef boost::dynamic_bitset<unsigned long> BitArray;

typedef std::pair<unsigned short, unsigned short> UShortUShortPair;
typedef std::pair<int, int>                       IntIntPair;
typedef std::pair<Real, Real>                     RealRealPair;
typedef std::pair<Real, RealVector>               RealRealVectorPair;

typedef std::list<unsigned short>    UShortList;
typedef std::list<size_t>            SizetList;

typedef std::vector<Real>              RealArray;
typedef std::vector<RealArray>         Real2DArray;
typedef std::vector<Real2DArray>       Real3DArray;
typedef std::vector<int>               IntArray;
typedef std::vector<IntArray>          Int2DArray;
typedef std::vector<short>             ShortArray;
typedef std::vector<unsigned short>    UShortArray;
typedef std::vector<unsigned long>     ULongArray;
typedef std::vector<UShortArray>       UShort2DArray;
typedef std::vector<UShort2DArray>     UShort3DArray;
typedef std::vector<UShort3DArray>     UShort4DArray;
typedef std::vector<UShort4DArray>     UShort5DArray;
typedef std::vector<size_t>            SizetArray;
typedef std::vector<SizetArray>        Sizet2DArray;
typedef std::vector<Sizet2DArray>      Sizet3DArray;
typedef std::vector<std::complex<Real> > ComplexArray;
typedef std::vector<RealRealPair>      RealRealPairArray;
typedef std::vector<String>            StringArray;
typedef std::vector<RealVector>        RealVectorArray;
typedef std::vector<RealVectorArray>   RealVector2DArray;
typedef std::vector<RealVector2DArray> RealVector3DArray;
typedef std::vector<IntVector>         IntVectorArray;
typedef std::vector<RealMatrix>        RealMatrixArray;
typedef std::vector<RealMatrixArray>   RealMatrix2DArray;
typedef std::vector<RealMatrix2DArray> RealMatrix3DArray;
typedef std::vector<RealSymMatrix>     RealSymMatrixArray;

//typedef std::vector<BasisFunction>  BasisFunctionArray;
typedef std::vector<SurrogateDataVars> SDVArray;
typedef std::vector<SurrogateDataResp> SDRArray;
typedef std::vector<SDVArray>          SDV2DArray;
typedef std::vector<SDRArray>          SDR2DArray;

typedef std::set<size_t>                  SizetSet;
typedef std::set<int>                     IntSet;
typedef std::set<String>                  StringSet;
typedef std::set<Real>                    RealSet;
typedef std::set<BitArray>                BitArraySet;
typedef std::set<UShortArray>             UShortArraySet;
typedef std::multiset<unsigned short>     UShortMultiSet;
typedef std::multiset<UShortMultiSet>     UShort2DMultiSet;
typedef std::vector<IntSet>               IntSetArray;
typedef std::vector<StringSet>            StringSetArray;
typedef std::vector<RealSet>              RealSetArray;
typedef std::map<int, short>              IntShortMap;
typedef std::map<size_t, short>           SizetShortMap;
typedef std::map<BitArray, unsigned long> BitArrayULongMap;
typedef std::map<unsigned long, unsigned long> ULongULongMap;
typedef std::map<int, int>                IntIntMap;
typedef std::map<int, Real>               IntRealMap;
typedef std::map<String, Real>            StringRealMap;
typedef std::map<Real, Real>              RealRealMap;
typedef std::map<RealRealPair, Real>      RealRealPairRealMap;
typedef std::map<IntIntPair, Real>        IntIntPairRealMap;
typedef std::vector<SizetSet>             SizetSetArray;
typedef std::vector<UShortArraySet>       UShortArraySetArray;
typedef std::vector<IntRealMap>           IntRealMapArray;
typedef std::vector<StringRealMap>        StringRealMapArray;
typedef std::vector<RealRealMap>          RealRealMapArray;
typedef std::vector<RealRealPairRealMap>  RealRealPairRealMapArray;
typedef std::vector<IntIntPairRealMap>    IntIntPairRealMapArray;
typedef std::map<int, RealVector>         IntRealVectorMap;
typedef std::map<UShortMultiSet,   Real>  UShortMultiSetRealMap;
typedef std::map<UShort2DMultiSet, Real>  UShort2DMultiSetRealMap;

typedef boost::multi_array_types::index_range      idx_range;
typedef boost::multi_array<size_t, 1>              SizetMultiArray;
typedef SizetMultiArray::array_view<1>::type       SizetMultiArrayView;
typedef SizetMultiArray::const_array_view<1>::type SizetMultiArrayConstView;

// ---------
// Iterators
// ---------
typedef IntSet::iterator                    ISIter;
typedef IntSet::const_iterator              ISCIter;
typedef SizetSet::iterator                  StSIter;
typedef SizetSet::const_iterator            StSCIter;
typedef BitArraySet::iterator               BASIter;
typedef BitArraySet::const_iterator         BASCIter;
typedef RealSet::iterator                   RSIter;
typedef RealSet::const_iterator             RSCIter;
typedef IntShortMap::iterator               IShMIter;
typedef IntShortMap::const_iterator         IShMCIter;
typedef SizetShortMap::iterator             StShMIter;
typedef SizetShortMap::const_iterator       StShMCIter;
typedef IntIntMap::iterator                 IIMIter;
typedef IntIntMap::const_iterator           IIMCIter;
typedef BitArrayULongMap::iterator          BAULMIter;
typedef BitArrayULongMap::const_iterator    BAULMCIter;
typedef ULongULongMap::iterator             ULULMIter;
typedef ULongULongMap::const_iterator       ULULMCIter;
typedef IntRealMap::iterator                IRMIter;
typedef IntRealMap::const_iterator          IRMCIter;
typedef StringRealMap::iterator             SRMIter;
typedef StringRealMap::const_iterator       SRMCIter;
typedef RealRealMap::iterator               RRMIter;
typedef RealRealMap::const_iterator         RRMCIter;
typedef IntIntPairRealMap::iterator         IIPRMIter;
typedef IntIntPairRealMap::const_iterator   IIPRMCIter;
typedef RealRealPairRealMap::iterator       RRPRMIter;
typedef RealRealPairRealMap::const_iterator RRPRMCIter;


/// equality operator for SizetArray and SizetMultiArrayConstView
inline bool operator==(const SizetArray& sa, SizetMultiArrayConstView smav)
{
  // Check for equality in array lengths
  size_t i, len = sa.size();
  if ( smav.size() != len )
    return false;

  // Check each size_t
  for (i=0; i<len; i++)
    if ( smav[i] != sa[i] )
      return false;

  return true;
}


/// equality operator for SizetArray and SizetMultiArrayConstView
template <typename OrdinalType, typename ScalarType>
bool equivalent(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& v,
		const std::map<ScalarType, ScalarType>& m)
{
  // Check for equality in array lengths
  size_t i, m_len = m.size();
  OrdinalType v_len = v.length();
  if ( m_len * 2 != v_len )
    return false;

  // Check each key,value pair
  typename std::map<ScalarType, ScalarType>::const_iterator cit;
  OrdinalType cntr = 0;
  for (cit=m.begin(); cit!=m.end(); ++cit) {
    if (v[cntr] != cit->first)  return false; ++cntr;
    if (v[cntr] != cit->second) return false; ++cntr;
  }

  return true;
}


template <typename PecosContainerType>
inline size_t find_index(const PecosContainerType& c,
			 const typename PecosContainerType::value_type& val)
{
  // For a default container, employ one traversal

  typename PecosContainerType::const_iterator cit;
  size_t cntr; // force size_t to ensure that _NPOS is valid
  for (cit=c.begin(), cntr=0; cit!=c.end(); ++cit, ++cntr)
    if (*cit == val)
      return cntr;
  return _NPOS;
}


template <typename ValueType>
inline size_t set_value_to_index(const std::set<ValueType>& s,
				 const ValueType& val)
{
  // For a sorted container, use fast lookup + distance()

  typename std::set<ValueType>::const_iterator cit = s.find(val);
  return (cit == s.end()) ? _NPOS : std::distance(s.begin(), cit);
}


template <typename KeyType, typename ValueType>
inline size_t map_key_to_index(const std::map<KeyType, ValueType>& m,
			       const KeyType& key)
{
  // For a sorted container, use fast lookup + distance()

  typename std::map<KeyType, ValueType>::const_iterator cit = m.find(key);
  return (cit == m.end()) ? _NPOS : std::distance(m.begin(), cit);
}


/// compare two Real values using DBL_EPSILON relative tolerance:
/// return true if same, false if different
inline bool real_compare(Real r1, Real r2)
{
  if (r2 >= DBL_MAX || r2 <= -DBL_MAX) // +/-DBL_MAX or +/-dbl_inf
    return (r1 == r2);
  else if  (std::abs(r2) < DBL_MIN)    // +/-0 (see IEEE 754)
    return (std::abs(r1) < DBL_MIN);
  else
    return (std::abs(1. - r1/r2) <= DBL_EPSILON); // relative
}


/// compute 1-norm |i| (sum of indices) for the given index_set
template <typename OrdinalType> 
inline size_t l1_norm(const std::vector<OrdinalType>& index_set)
{
  // hardwire to size_t instead of OrdinalType since could be a large sum
  // of smaller types
  size_t i, norm = 0, len = index_set.size();
  // assume unsigned types since std::abs(index_set[i]) will not compile
  // for unsigned types on some platforms
  for (i=0; i<len; ++i)
    norm += index_set[i];//std::abs(index_set[i]);
  return norm;
}


/// compare two Real values using DBL_EPSILON relative tolerance
template <typename OrdinalType, typename ScalarType> 
void inflate_scalar(std::vector<ScalarType>& v, OrdinalType num_v)
{
  OrdinalType v_len = v.size();
  if (v_len != num_v) {
    if (v_len == 1) {
      ScalarType v0 = v[0];
      v.assign(num_v, v0);
    }
    else {
      PCerr << "Error: specification length (" << v_len
	    << ") does not match target length (" << num_v
	    << ") in Pecos::inflate_scalar()." << std::endl;
      abort_handler(-1);
    }
  }
}


/// copy Teuchos::SerialDenseVector<OrdinalType, ScalarType> to
/// std::vector<ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_data(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv,
	       std::vector<ScalarType>& v)
{
  OrdinalType size_sdv = sdv.length();
  if (size_sdv != v.size())
    v.resize(size_sdv);
  for (OrdinalType i=0; i<size_sdv; ++i)
    v[i] = sdv[i];
}


/// copy std::vector<ScalarType> to
/// Teuchos::SerialDenseVector<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_data(const std::vector<ScalarType>& v,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv)
{
  size_t size_v = v.size();
  if (sdv.length() != size_v)
    sdv.sizeUninitialized(size_v);
  for (OrdinalType i=0; i<size_v; ++i)
    sdv[i] = v[i];
}


/// copy Teuchos::SerialDenseVector<OrdinalType, ScalarType> to same
/// (used in place of operator= when a deep copy is required regardless
/// of Teuchos DataAccess mode)
template <typename OrdinalType, typename ScalarType> 
void copy_data(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv1,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv2)
{
  OrdinalType size_sdv1 = sdv1.length();
  if (size_sdv1 != sdv2.length())
    sdv2.sizeUninitialized(size_sdv1);
  for (OrdinalType i=0; i<size_sdv1; ++i)
    sdv2[i] = sdv1[i];
}


/// copy Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType> to same
/// (used in place of operator= when a deep copy is required regardless
/// of Teuchos DataAccess mode)
template <typename OrdinalType, typename ScalarType> 
void copy_data(const
	       Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType>& ssdm1,
	       Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType>& ssdm2)
{
  OrdinalType size_ssdm1 = ssdm1.numRows();
  if (size_ssdm1 != ssdm2.numRows())
    ssdm2.shapeUninitialized(size_ssdm1);
  ssdm2.assign(ssdm1); // copies values
}


/// copy Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType> to same
/// (used in place of operator= when a deep copy is required regardless
/// of Teuchos DataAccess mode)
template <typename OrdinalType, typename ScalarType1, typename ScalarType2> 
void copy_data(const
	       Teuchos::SerialDenseMatrix<OrdinalType, ScalarType1>& sdm1,
	       Teuchos::SerialDenseMatrix<OrdinalType, ScalarType2>& sdm2)
{
  OrdinalType r, c, nr_sdm1 = sdm1.numRows(), nc_sdm1 = sdm1.numCols();
  if (nr_sdm1 != sdm2.numRows() || nc_sdm1 != sdm2.numCols())
    sdm2.shapeUninitialized(nr_sdm1, nc_sdm1);
  for (r=0; r<nr_sdm1; ++r)
    for (c=0; c<nc_sdm1; ++c)
      sdm2(r,c) = static_cast<ScalarType2>(sdm1(r,c));
}


/// copy ScalarType* to ScalarType*
template <typename OrdinalType, typename ScalarType> 
void copy_data(const ScalarType* ptr1, const OrdinalType ptr_len,
	       ScalarType* ptr2)
{
  for (OrdinalType i=0; i<ptr_len; ++i)
    ptr2[i] = ptr1[i];
}


/// copy ScalarType* to Teuchos::SerialDenseVector<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_data(const ScalarType* ptr, const OrdinalType ptr_len,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv)
{
  if (sdv.length() != ptr_len)
    sdv.sizeUninitialized(ptr_len);
  for (OrdinalType i=0; i<ptr_len; ++i)
    sdv[i] = ptr[i];
}


/// copy Array<T> to T*
template <typename T>
void copy_data(const std::vector<T>& vec, T* ptr, const size_t ptr_len)
{
  if (ptr_len != vec.size()) { // could use <, but enforce exact match
    PCerr << "Error: bad ptr_len in copy_data(std::vector<T>, T* ptr)."
	  << std::endl;
    abort_handler(-1);
  }
  for (size_t i=0; i<ptr_len; ++i)
    ptr[i] = vec[i];
}


/// copy T* to Array<T>
template <typename T>
void copy_data(const T* ptr, const size_t ptr_len, std::vector<T>& vec)
{
  if (ptr_len != vec.size())
    vec.resize(ptr_len);
  for (size_t i=0; i<ptr_len; ++i)
    vec[i] = ptr[i];
}


/// copy ScalarType1* to Array<ScalarType2> with data cast
template <typename ScalarType1, typename ScalarType2>
void copy_data(const ScalarType1* ptr, const size_t ptr_len,
	       std::vector<ScalarType2>& vec)
{
  if (ptr_len != vec.size())
    vec.resize(ptr_len);
  for (size_t i=0; i<ptr_len; ++i)
    vec[i] = static_cast<ScalarType2>(ptr[i]);
}


/// copy ScalarType* to BitArray
template <typename OrdinalType, typename ScalarType> 
void copy_data(const ScalarType* ptr, const OrdinalType ptr_len,
	       BitArray& ba)
{
  if (ba.size() != ptr_len)
    ba.resize(ptr_len);
  for (OrdinalType i=0; i<ptr_len; ++i)
    ba[i] = (ptr[i]) ? true : false;
}


/// copy BitArray to ScalarType*
template <typename OrdinalType> //, typename ScalarType>
void copy_data(const BitArray& ba, bool* ptr, //ScalarType* ptr,
	       const OrdinalType ptr_len)
{
  for (OrdinalType i=0; i<ptr_len; ++i)
    ptr[i] = ba[i]; // or static_cast<ScalarType>(ba[i])
}


template <typename OrdinalType, typename ScalarType>
void copy_data(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& v,
	       std::map<ScalarType, ScalarType>& m)
{
  OrdinalType v_len = v.length();
  if (v_len % 2) {
    PCerr << "Error: vector length (" << v_len << ") must be multiple of 2 in "
	  << "Pecos::copy_data(std::vector, std::map)." << std::endl;
    abort_handler(-1);
  }
  // convert vector of concatenated (x,y) pairs to map[x] = y
  OrdinalType i, j, m_len = v_len / 2;
  m.clear();
  for (i=0; i<m_len; ++i)
    { j = 2*i; m[v[j]] = v[j+1]; }
}


template <typename OrdinalType, typename ScalarType>
void copy_data(const std::map<ScalarType, ScalarType>& m,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& v)
{
  // convert map[x] = y to vector of concatenated (x,y) pairs
  v.sizeUninitialized(m.size() * 2);
  typename std::map<ScalarType, ScalarType>::const_iterator cit;
  OrdinalType cntr = 0;
  for (cit=m.begin(); cit!=m.end(); ++cit) {
    v[cntr] = cit->first;  ++cntr;
    v[cntr] = cit->second; ++cntr;
  }
}

/// specialization for IntScalarTypeMap 
template <typename OrdinalType, typename ScalarType>
void copy_data(const std::map<int, ScalarType>& m,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& v)
{
  // convert map[x] = y to vector of concatenated (x,y) pairs
  v.sizeUninitialized(m.size() * 2);
  typename std::map<int, ScalarType>::const_iterator cit;
  OrdinalType cntr = 0;
  for (cit=m.begin(); cit!=m.end(); ++cit) {
    // possibly promote int to double
    v[cntr] = cit->first;  ++cntr;
    v[cntr] = cit->second; ++cntr;
  }
}

/// specialization for StringScalarTypeMap: store index of string value 
template <typename OrdinalType, typename ScalarType>
void copy_data(const std::map<String, ScalarType>& m,
	       Teuchos::SerialDenseVector<OrdinalType, ScalarType>& v)
{
  // convert map[x] = y to vector of concatenated (x,y) pairs
  v.sizeUninitialized(m.size() * 2);
  typename std::map<String, ScalarType>::const_iterator cit;
  size_t str_index = 0;
  OrdinalType cntr = 0;
  for (cit=m.begin(); cit!=m.end(); ++cit) {
    // possibly promote string index to double for dist params
    v[cntr] = str_index;   ++cntr; ++str_index;
    v[cntr] = cit->second; ++cntr;
  }
}

/// copy Teuchos::SerialDenseVector<OrdinalType, ScalarType> to
/// ith row of Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_row(const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv,
	      Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>& sdm,
	      OrdinalType row)
{
  OrdinalType i, len = sdv.length();
  //if (sdm.numCols() != len)
  //  PCerr << std::endl;
  for (i=0; i<len; ++i)
    sdm(row, i) = sdv[i];
}


/// copy ith row of Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>
/// to Teuchos::SerialDenseVector<OrdinalType, ScalarType>
template <typename OrdinalType, typename ScalarType> 
void copy_row(const Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>& sdm,
	      OrdinalType row,
	      const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv)
{
  OrdinalType i, len = sdm.numCols();
  if (sdv.length() != len)
    sdv.sizeUninitialized(len);
  for (OrdinalType i=0; i<len; ++i)
    sdv[i] = sdm(row, i);
}


/// std::ostream write for Teuchos::SerialDenseVector
template <typename OrdinalType, typename ScalarType>
void write_data(std::ostream& s,
		const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& v)
{
  s << std::scientific << std::setprecision(WRITE_PRECISION);
  OrdinalType len = v.length();
  for (OrdinalType i=0; i<len; i++)
    s << "                     " << std::setw(WRITE_PRECISION+7)
      << v[i] << '\n';
}


/// formatted ostream insertion operator for SerialDenseMatrix
template <typename OrdinalType, typename ScalarType>
void write_data(std::ostream& s,
                const Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>& m,
                bool brackets, bool row_rtn, bool final_rtn)
{
  OrdinalType i, j, nrows = m.numRows(), ncols = m.numCols();
  s << std::scientific << std::setprecision(WRITE_PRECISION);
  if (brackets)  s << "[[ ";
  for (i=0; i<nrows; ++i) {
    for (j=0; j<ncols; ++j)
      s << std::setw(WRITE_PRECISION+7) << m(i,j) << ' ';
    if (row_rtn && i!=m.numRows()-1)
      s << "\n";
  }
  if (brackets)  s << "]] ";
  if (final_rtn) s << '\n';
}


/// formatted ostream insertion operator for SerialSymDenseMatrix
template <typename OrdinalType, typename ScalarType>
void write_data(std::ostream& s,
                const Teuchos::SerialSymDenseMatrix<OrdinalType, ScalarType>& m,
                bool brackets, bool row_rtn, bool final_rtn)
{
  OrdinalType i, j, nrows = m.numRows();
  s << std::scientific << std::setprecision(WRITE_PRECISION);
  if (brackets)  s << "[[ ";
  for (i=0; i<nrows; ++i) {
    for (j=0; j<nrows; ++j)
      s << std::setw(WRITE_PRECISION+7) << m(i,j) << ' ';
    if (row_rtn && i!=m.numRows()-1)
      s << "\n   ";
  }
  if (brackets)  s << "]] ";
  if (final_rtn) s << '\n';
}


/// ostream insertion operator for a column vector from a SerialDenseMatrix
template <typename OrdinalType, typename ScalarType>
void write_data_trans(std::ostream& s, 
  const Teuchos::SerialDenseVector<OrdinalType, ScalarType>& sdv,
  bool brackets, bool row_rtn, bool final_rtn)
{
  OrdinalType i, num_items = sdv.length();
  s << std::scientific << std::setprecision(WRITE_PRECISION);
  if (brackets) s << " [ ";
  else          s << "   ";
  for (i=0; i<num_items; ++i) {
    s << std::setw(WRITE_PRECISION+7) << sdv[i] << ' ';
    if (row_rtn && (i+1)%4 == 0)
      s << "\n   "; // Output 4 gradient components per line
  }
  if (brackets)  s << "] ";
  if (final_rtn) s << '\n';
}


/// global std::ostream insertion operator for std::vector
template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& data)
{
  s << std::scientific << std::setprecision(WRITE_PRECISION);
  size_t i=0, len = data.size();
  for (i=0; i<len; ++i)
    s << "                     " << std::setw(WRITE_PRECISION+7)
      << data[i] << '\n';
  return s;
}


/// global std::ostream insertion operator for std::set
template <typename T>
std::ostream& operator<<(std::ostream& s, const std::set<T>& data)
{
  for (typename std::set<T>::const_iterator cit = data.begin();
       cit != data.end(); ++cit)
    s << "                     " << std::setw(WRITE_PRECISION+7)
      << *cit << '\n';
  return s;
}

} // namespace Pecos

#endif // PECOS_DATA_TYPES_H
