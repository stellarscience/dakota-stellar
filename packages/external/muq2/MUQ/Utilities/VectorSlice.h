#ifndef VECTORSLICE_H
#define VECTORSLICE_H

#include <vector>
#include <Eigen/Core>

namespace muq {
namespace Utilities {

    /** @class VectorSlice
        @brief Enables a subset of a vector to be easily accessed or reversed without copying memory.
        @details Given a vector v, this class provides a way for defining slices
        of the vector of the form v[startInd:endInd:skip] (in python notation).

        This class is designed for use with the std::vector and Eigen::VectorXd classes,
        but can work with any type that exposes the [] operator and has a size()
        function, which includes other VectorViews.

        IMPORTANT:
        This class does not copy data, so segfaults could occur if
        the underlying vector goes out of scope before the view is destructed.
    */
    template<typename VecType, typename ScalarType>
    class VectorSlice
    {
    public:

      VectorSlice(VectorSlice<VecType, ScalarType> const& vec) : VectorSlice(vec.data, vec.startInd, vec.endInd, vec.skip){};

      VectorSlice<VecType, ScalarType>& operator=(VectorSlice<VecType, ScalarType> const& vec){
          data = vec.data;
          startInd = vec.startInd;
          endInd = vec.endInd;
          skip = vec.skip;
          return *this;
      }

      /** Construct the view with a reference to vector of doubles.  Note that
      endInd-startInd should have the same size as skip.  Thus, reversing the
      order of std::vector would be accomplished with
      @code{.cpp}
      std::vector v{1,2,3,4};
      VectorSlice vRev(v, v.size()-1,-1, -1);
      @endcode
      Note that the second argument, endInd, is negative because the end index is
      not inclusive.  Thus, endInd must be -1 to capture v[0].

      @param[in] dataIn The base vector that we want to slice.
      @param[in] startIndIn non-negative integer specifying the starting index in the data vector
      @param[in] endIndIn non-negative integer specifying the ending index in the data vector.  Numpy indexing rules are followed, so this index is not-inclusive.
      @parm[in] skip.  The number of
      */
      VectorSlice(VecType& dataIn,
                  int startIndIn,
                  int endIndIn,
                  int skipIn=1) : data(dataIn),
                                  startInd(startIndIn),
                                  endInd(endIndIn),
                                  skip(skipIn)
      {
        assert(skip!=0);
        assert(startInd>=0);
        assert(endInd>=-1);
        if(skip>0){
          assert(endInd>startInd);
        }else{
          assert(endInd<startInd);
        }
      }

      ScalarType& operator()(int i) {
        CheckBounds(i);
        return data[startInd + i*skip];
      }

      ScalarType operator()(int i) const {
        CheckBounds(i);
        return data[startInd + i*skip];
      }

      ScalarType& operator[](int i) {
        return data[startInd + i*skip];
      }

      ScalarType operator[](int i) const {
        return data[startInd + i*skip];
      }

      unsigned int size() const{
        return std::ceil( double(std::abs(endInd-startInd)) / double(std::abs(skip)) );
      };


      VecType& data;
      int startInd, endInd, skip;

      void CheckBounds(int i) const {
        if(skip>0){
          assert(startInd + i*skip < endInd);
        }else{
          assert(startInd + i*skip > endInd);
        }
      }

    };


    template<typename ScalarType>
    VectorSlice<std::vector<ScalarType>, ScalarType> GetSlice(std::vector<ScalarType>& dataIn,
                         int                            startIndIn,
                         int                            endIndIn,
                         int                            skipIn=1)
    {
      return VectorSlice<std::vector<ScalarType>,ScalarType>(dataIn, startIndIn, endIndIn, skipIn);
    }

    template<typename ScalarType>
    VectorSlice<Eigen::Matrix<ScalarType,Eigen::Dynamic,1>, ScalarType> GetSlice(Eigen::Matrix<ScalarType,Eigen::Dynamic,1>& dataIn,
                         int                                               startIndIn,
                         int                                               endIndIn,
                         int                                               skipIn=1)
    {
      return VectorSlice<Eigen::Matrix<ScalarType,Eigen::Dynamic,1>,ScalarType>(dataIn, startIndIn, endIndIn, skipIn);
    }

    template<typename VectorType, typename ScalarType>
    VectorSlice<VectorType, ScalarType> GetSlice(VectorSlice<VectorType, ScalarType>& dataIn,
                         int                                        startIndIn,
                         int                                        endIndIn,
                         int                                        skipIn=1)
    {
        return VectorSlice<VectorType,ScalarType>(dataIn.data,
                                                  dataIn.startInd + dataIn.skip*startIndIn,
                                                  std::max(dataIn.startInd + dataIn.skip*endIndIn,-1),
                                                  dataIn.skip * skipIn);

    }

}
}


#endif // #ifndef VECTORSLICE_H
