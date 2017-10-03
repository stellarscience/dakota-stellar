// *******************************************************************
// *******************************************************************
// Definition for the class Mars
// Date         : November, 1998
//
// Modified by Kevin Long, April 21 1999.
// *******************************************************************
// *******************************************************************

#ifndef MARS_H
#define MARS_H

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  #include "ddace_config.h"
  #define MARS_F77 F77_FUNC(mars,MARS)
  #define FMODM_F77 F77_FUNC(fmodm,FMODM)
#else
  #include "ddace_fortran.h"
  #define MARS_F77 DDACE_FORTRAN_GLOBAL(mars,MARS)
  #define FMODM_F77 DDACE_FORTRAN_GLOBAL(fmodm,FMODM)
#endif /* HAVE_CONFIG_H */

#include "FuncApproxBase.h"
#include <vector>
#include <iostream>
#include <fstream>

/**
 * Mars is a handle class to FuncApproxBase.
 * This class is an interface to the MARS package for generating
 * response surfaces.
 * <p>
 * EXAMPLE:  Given some (x,y,z) values, MARS can compute
 * the surface in 3D space.  i.e. MARS will divide the surface
 * into equal-area rectangles.  MARS will return the (x,y,z)
 * value of each and every rectangle.
 *
 * Author       : Charles Tong
 * Modified by Kevin Long, April 21 1999.
 */

#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void MARS_F77(int&, int&, float*, float*, float*, int&,
	      int&, int*, float*, int*, float*, double*, int*);

#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void FMODM_F77(int&, int&, float*, float*, int*, float*, float*);

class Mars : public FuncApproxBase 
{
 public:

  /**
   * Constructor
   */
  Mars();

  /**
   * Constructor
   * @param nInputs The number of input, or independent, variables.
   * For example, for 3D space, the number of input variables is 2
   * (the x and y values).
   * @param nSamples The number input data points.
   * For example, for 3D space, the number of (x,y,z) tuples
   * that are given to MARS.
   *
   */
  Mars(int nInputs, int nSamples);


  Mars(int nInputs, int nSamples, int nBasis, 
			 int maxVarPerBasis, const std::vector<int>& varFlags);

  virtual ~Mars(){;}

  
  
  virtual void generateGridData(const std::vector<DDaceSamplePoint>& samplePts,
				const std::vector<double>& sampleResults,
				std::vector<std::vector<double> >& gridPts,
				std::vector<double>& gridResults) ;

  /**
   * Given some (x,y,z) data values, compute the
   * surface of the 3D space.  MARS will divide the surface
   * into equal-area rectangle.  MARS will return the (x,y,z)
   * value of each and every rectangle.  
   * <P>
   * SIDE EFFECT:  The parameter, gridpts, is populated
   * with the (x,y) coordiantes of the computer surface.
   * SIDE EFFECT:  the parameter, gridResults, is populated
   * with the (z) coordiantes of the computer surface.
   * @param filename The name of the output text file.
   * @param samplePts Contains (x1, x2, ...) tuples of data
   * values.  The input variable var1 will select one of
   * these variables to represent the x axis.  the input
   * variable var2 will select one of these variables to
   * represent the y axis. 
   * @param sampleResults The z coordinates of the data that
   * MARS will process.
   * @param var1 The index value of the x-axis variable.  
   * If the index value is set to 0, then the first axis
   * in the input parameter, samplePts, will become the x axis.
   * @param var2 The index value of the y-axis variable.  
   * If the index value is set to 1, then the second axis
   * in the input parameter, samplePts, will become the y axis.
   * @settings The first element contains the number of dimensions
   * of the response surface (which should be set to 2); the second
   * element contains the initial z value of the response surface.
   * @param gridPts Will contain the (x,y) pairs of every
   * rectangle on the computed surface. <br>
   * Example, if the x-y space is divided into 9 rectangles, then <br>
   * &nbsp;&nbsp;girdPts[0] contains (y[0],x[0]) <br>
   * &nbsp;&nbsp;girdPts[1] contains (y[0],x[1]) <br>
   * &nbsp;&nbsp;girdPts[2] contains (y[0],x[2]) <br>
   * &nbsp;&nbsp;girdPts[3] contains (y[1],x[0]) <br>
   * &nbsp;&nbsp;girdPts[4] contains (y[1],x[1]) <br>
   * &nbsp;&nbsp;girdPts[5] contains (y[1],x[2]) <br>
   * &nbsp;&nbsp;girdPts[6] contains (y[2],x[0]) <br>
   * &nbsp;&nbsp;girdPts[7] contains (y[2],x[1]) <br>
   * &nbsp;&nbsp;girdPts[8] contains (y[2],x[2]) <br>
   *
   * @param gridResults Will contain the (z) coordinate of every
   * rectangle on the computer surface. <br>
   * Example, if the x-y space is divided into 9 rectangles, then <br>
   * &nbsp;&nbsp;girdResults[0] contains z value at (y[0],x[0]) <br>
   * &nbsp;&nbsp;girdResults[1] contains z value at (y[0],x[1]) <br>
   * &nbsp;&nbsp;girdResults[2] contains z value at (y[0],x[2]) <br>
   * &nbsp;&nbsp;girdResults[3] contains z value at (y[1],x[0]) <br>
   * &nbsp;&nbsp;girdResults[4] contains z value at (y[1],x[1]) <br>
   * &nbsp;&nbsp;girdResults[5] contains z value at (y[1],x[2]) <br>
   * &nbsp;&nbsp;girdResults[6] contains z value at (y[2],x[0]) <br>
   * &nbsp;&nbsp;girdResults[7] contains z value at (y[2],x[1]) <br>
   * &nbsp;&nbsp;girdResults[8] contains z value at (y[2],x[2]) <br>
   */
  virtual void generate2DGridData(const std::vector<DDaceSamplePoint>& samplePts,
				  const std::vector<double>& sampleResults,
				  int var1, 
				  int var2,
				  const std::vector<double>& settings,
				  std::vector<std::vector<double> >& gridPts,
				  std::vector<double>& gridResults) ;
  
  /**
   * Given some (x,y,z) data values, compute the
   * surface of the 3D space.  MARS will divide the surface
   * into equal-area rectangle.  MARS will return the (x,y,z)
   * value of each and every rectangle.  The (x,y,z) values
   * are written to a text file; the format of the file is: <br>
   * &nbsp;&nbsp;n, number of coordinates per axis <br>
   * &nbsp;&nbsp;x[0], value of 1st coordinate along x axis <br>
   * &nbsp;&nbsp;x[1], value of 2nd coordinate along x axis <br>
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;x[n], value of last coordinate along x axis <br>
   * &nbsp;&nbsp;y[0], value of 1st coordinate along x axis <br>
   * &nbsp;&nbsp;y[1], value of 2nd coordinate along x axis <br>
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;y[n], value of last coordinate along y axis <br>
   * &nbsp;&nbsp;z[0][0], value of z at y[0],x[0] <br>
   * &nbsp;&nbsp;z[0][1], value of z at y[0],x[1] <br>
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;z[0][n], value of z at y[0],x[n] <br>
   * &nbsp;&nbsp;z[1][0], value of z at y[1],x[0] <br>
   * &nbsp;&nbsp;z[1][1], value of z at y[1],x[1] <br>
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;z[1][n], value of z at y[1],x[n] <br>
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;z[n][0], value of z at y[n],x[0] <br>
   * &nbsp;&nbsp;z[n][1], value of z at y[n],x[1] <br>
   * &nbsp;&nbsp;:
   * &nbsp;&nbsp;z[n][n], value of z at y[n],x[n] <br>
   *
   * @param filename The name of the output text file.
   * @param samplePts Contains (x1, x2, ...) tuples of data
   * values.  The input variable var1 will select one of
   * these variables to represent the x axis.  the input
   * variable var2 will select one of these variables to
   * represent the y axis. 
   * @param sampleResults The z coordinates of the data that
   * MARS will process.
   * @param var1 The index value of the x-axis variable.  
   * If the index value is set to 0, then the first axis
   * in the input parameter, samplePts, will become the x axis.
   * @param var2 The index value of the y-axis variable.  
   * If the index value is set to 1, then the second axis
   * in the input parameter, samplePts, will become the y axis.
   * @settings The first element contains the number of dimensions
   * of the response surface (which should be set to 2); the second
   * element contains the initial z value of the response surface.
   */
  virtual void write2DGridData(const std::string& filename,
			       const std::vector<DDaceSamplePoint>& samplePts,
			       const std::vector<double>& sampleResults,
			       int var1, 
			       int var2,
			       const std::vector<double>& settings);
  
  virtual double evaluatePoint(const std::vector<double> & x) const ;
  virtual void evaluatePoints(const std::vector<std::vector<double> >& pts, 
			      std::vector<double>& y) const ;
  
  // preprocess: call fortran mars() to get the regression coeffs.
  void preProcess(const std::vector<DDaceSamplePoint>& samplePts,
		  const std::vector<double>& sampleResults);
  
  virtual FuncApproxBase* clone() const ;
 private:
  
  
  
  std::vector<float> weights_;
  int nBasis_;
  int maxVarPerBasis_;
  std::vector<int> varFlags_;
  std::vector<float> fm_;
  std::vector<int> im_;
  
  // set default values using static variables, so that they're 
  // available during constructor calls.
  static int defaultNBasis_;
  static int defaultMaxVarPerBasis_;
};

#endif // MARS_H






