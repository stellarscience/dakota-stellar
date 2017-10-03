//- File:        aData.h
//- Description: This header defines data type aData for PSUADE
//-              (eventually to be replaced with PSUADE's header)
//- Owner:       Brian M. Adams, Sandia National Laboratories

#ifndef __ADATAH_
#define __ADATAH_

// definition of data type aData (Analyzer Data) needed for use with PSUADE
typedef struct {
  int nInputs_;
  int nOutputs_;
  int nSamples_;
  double *iLowerB_;
  double *iUpperB_;
  double *sampleInputs_;
  double *sampleOutputs_;
  int outputID_;
} aData;

#endif
