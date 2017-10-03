//- File:        aData.h
//- Description: This header defines data type aData for PSUADE
//-              (eventually to be replaced with PSUADE's header)
//- Owner:       Brian M. Adams, Sandia National Laboratories

#ifndef __SDATAH_
#define __SDATAH_

// definition of data type sDate (Sampling Data) needed for use with PSUADE
typedef struct {
  int nInputs_;
  int nSamples_;
  double *iLowerB_;
  double *iUpperB_;
} sData;


#endif
