/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#include "FourierInverseTransformation.hpp"
#include "KarhunenLoeveInverseTransformation.hpp"
#include "SamplingInverseTransformation.hpp"

static const char rcsId[]="@(#) $Id: DataTransformation.cpp 4768 2007-12-17 17:49:32Z mseldre $";

namespace Pecos {


/** This constructor is the one which must build the base class data
    for all derived classes.  get_data_trans() instantiates a derived
    class letter and the derived constructor selects this base class
    constructor in its initialization list (to avoid recursion in the
    base class constructor calling get_data_trans() again).  Since the
    letter IS the representation, its rep pointer is set to NULL (an
    uninitialized pointer causes problems in ~DataTransformation). */
DataTransformation::DataTransformation(BaseConstructor):
  dataTransRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "DataTransformation::DataTransformation(BaseConstructor) called to "
        << "build base class for letter." << std::endl;
#endif
}


/** The default constructor: dataTransRep is NULL in this case.  This
    makes it necessary to check for NULL in the copy constructor,
    assignment operator, and destructor. */
DataTransformation::DataTransformation(): dataTransRep(NULL), referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "DataTransformation::DataTransformation() called to build empty "
        << "envelope." << std::endl;
#endif
}


/** Envelope constructor only needs to extract enough data to properly
    execute get_data_trans, since DataTransformation(BaseConstructor)
    builds the actual base class data for the derived transformations. */
DataTransformation::DataTransformation(const String& data_trans_type):
  referenceCount(1)
{
#ifdef REFCOUNT_DEBUG
  PCout << "DataTransformation::DataTransformation(string&) called to "
        << "instantiate envelope." << std::endl;
#endif

  // Set the rep pointer to the appropriate derived type
  dataTransRep = get_data_trans(data_trans_type);
  if ( !dataTransRep ) // bad type or insufficient memory
    abort_handler(-1);
}


/** Used only by the envelope constructor to initialize dataTransRep to the 
    appropriate derived type. */
DataTransformation* DataTransformation::
get_data_trans(const String& data_trans_type)
{
#ifdef REFCOUNT_DEBUG
  PCout << "Envelope instantiating letter in get_data_trans(string&)."
        << std::endl;
#endif

  if (data_trans_type == "inverse_fourier_shinozuka_deodatis" ||
      data_trans_type == "inverse_fourier_grigoriu")
    return new FourierInverseTransformation(data_trans_type);
  else if (data_trans_type == "inverse_kl")
    return new KarhunenLoeveInverseTransformation();
  else if (data_trans_type == "inverse_sampling")
    return new SamplingInverseTransformation();
  else {
    PCerr << "Error: DataTransformation type " << data_trans_type
	  << " not available." << std::endl;
    return NULL;
  }
}


/** Copy constructor manages sharing of dataTransRep and incrementing
    of referenceCount. */
DataTransformation::DataTransformation(const DataTransformation& data_trans)
{
  // Increment new (no old to decrement)
  dataTransRep = data_trans.dataTransRep;
  if (dataTransRep) // Check for an assignment of NULL
    dataTransRep->referenceCount++;

#ifdef REFCOUNT_DEBUG
  PCout << "DataTransformation::DataTransformation(DataTransformation&)"
        << std::endl;
  if (dataTransRep)
    PCout << "dataTransRep referenceCount = " << dataTransRep->referenceCount
	  << std::endl;
#endif
}


/** Assignment operator decrements referenceCount for old dataTransRep,
    assigns new dataTransRep, and increments referenceCount for new
    dataTransRep. */
DataTransformation DataTransformation::
operator=(const DataTransformation& data_trans)
{
  if (dataTransRep != data_trans.dataTransRep) { // normal case: old != new
    // Decrement old
    if (dataTransRep) // Check for null pointer
      if (--dataTransRep->referenceCount == 0) 
	delete dataTransRep;
    // Assign and increment new
    dataTransRep = data_trans.dataTransRep;
    if (dataTransRep) // Check for an assignment of NULL
      dataTransRep->referenceCount++;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  PCout << "DataTransformation::operator=(DataTransformation&)" << std::endl;
  if (dataTransRep)
    PCout << "dataTransRep referenceCount = " << dataTransRep->referenceCount
	  << std::endl;
#endif

  return *this; // calls copy constructor since returned by value
}


/** Destructor decrements referenceCount and only deletes dataTransRep
    when referenceCount reaches zero. */
DataTransformation::~DataTransformation()
{ 
  // Check for NULL pointer 
  if (dataTransRep) {
    --dataTransRep->referenceCount;
#ifdef REFCOUNT_DEBUG
    PCout << "dataTransRep referenceCount decremented to " 
	  << dataTransRep->referenceCount << std::endl;
#endif
    if (dataTransRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      PCout << "deleting dataTransRep" << std::endl;
#endif
      delete dataTransRep;
    }
  }
}


void DataTransformation::
initialize(const Real& total_t, const Real& w_bar, size_t seed)
{
  if (dataTransRep) // envelope fwd to letter
    dataTransRep->initialize(total_t, w_bar, seed);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine initialize() virtual fn.\n"
          << "       No default defined at DataTransformation base class.\n"
	  << std::endl;
    abort_handler(-1);
  }
}


void DataTransformation::
power_spectral_density(const String& psd_name, const Real& param)
{
  if (dataTransRep) // envelope fwd to letter
    dataTransRep->power_spectral_density(psd_name, param);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine power_spectral_density() "
	  << "virtual fn.\n       No default defined at DataTransformation "
	  << "base class.\n" << std::endl;
    abort_handler(-1);
  }
}


/*
void DataTransformation::power_spectral_density(fn_ptr)
{
  if (dataTransRep) // envelope fwd to letter
    dataTransRep->power_spectral_density(fn_ptr);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine power_spectral_density() "
	  << "virtual fn.\n       No default defined at DataTransformation "
	  << "base class.\n" << std::endl;
    abort_handler(-1);
  }
}
*/


void DataTransformation::power_spectral_density(const RealVector& psd)
{
  if (dataTransRep) // envelope fwd to letter
    dataTransRep->power_spectral_density(psd);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine power_spectral_density() "
	  << "virtual fn.\n       No default defined at DataTransformation "
	  << "base class.\n" << std::endl;
    abort_handler(-1);
  }
}


const RealVector& DataTransformation::compute_sample()
{
  if (dataTransRep) // envelope fwd to letter
    return dataTransRep->compute_sample();
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine compute_sample() virtual "
          << "fn.\nNo default defined at DataTransformation base class.\n"
	  << std::endl;
    abort_handler(-1);
  }
}


const RealMatrix& DataTransformation::compute_samples(size_t num_samples)
{
  if (dataTransRep) // envelope fwd to letter
    return dataTransRep->compute_samples(num_samples);
  else { // letter lacking redefinition of virtual fn
    PCerr << "Error: derived class does not redefine compute_samples() virtual "
          << "fn.\nNo default defined at DataTransformation base class.\n"
	  << std::endl;
    abort_handler(-1);
  }
}

} // namespace Pecos
