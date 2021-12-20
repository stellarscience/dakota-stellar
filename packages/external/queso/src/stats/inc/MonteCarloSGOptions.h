//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef UQ_MOC_SG_OPTIONS_H
#define UQ_MOC_SG_OPTIONS_H

#include <queso/Environment.h>
#include <queso/SequenceStatisticalOptions.h>

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/BoostInputOptionsParser.h>
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

#define UQ_MOC_SG_FILENAME_FOR_NO_FILE "."

// _ODV = option default value
#define UQ_MOC_SG_HELP                             ""
#define UQ_MOC_SG_DATA_OUTPUT_FILE_NAME_ODV        UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_DATA_OUTPUT_ALLOWED_SET_ODV      ""

#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_PERIOD_ODV      0
#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_NAME_ODV   UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_FILE_TYPE_ODV   UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MOC_SG_PSEQ_DATA_OUTPUT_ALLOWED_SET_ODV ""
#define UQ_MOC_SG_PSEQ_COMPUTE_STATS_ODV           0

#define UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_NAME_ODV    UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_QSEQ_DATA_INPUT_FILE_TYPE_ODV    UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MOC_SG_QSEQ_SIZE_ODV                    100
#define UQ_MOC_SG_QSEQ_DISPLAY_PERIOD_ODV          500
#define UQ_MOC_SG_QSEQ_MEASURE_RUN_TIMES_ODV       0
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_PERIOD_ODV      0
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_NAME_ODV   UQ_MOC_SG_FILENAME_FOR_NO_FILE
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_FILE_TYPE_ODV   UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT
#define UQ_MOC_SG_QSEQ_DATA_OUTPUT_ALLOWED_SET_ODV ""
#define UQ_MOC_SG_QSEQ_COMPUTE_STATS_ODV           0

#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
namespace boost {
  namespace program_options {
    class options_description;
  }
}
#endif

namespace QUESO {

/*! \file MonteCarloSGOptions.h
    \brief Classes to allow options to be passed to a Monte Carlo sequence generator.
*/

/*! \class McOptionsValues
 *  \brief This class provides options for the  Monte Carlo sequence generator if no input file is available.
 *
 *  Monte Carlo sequence generator expects options for its methods. This class provides default
 * values for such options if no input file is available. */

class McOptionsValues
{
public:
  //! @name Constructor/Destructor methods
  //@{

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  McOptionsValues            (const SsOptionsValues* alternativePSsOptionsValues,
                                     const SsOptionsValues* alternativeQSsOptionsValues);
  McOptionsValues            (const SsOptionsValues* alternativePSsOptionsValues,
                                     const SsOptionsValues* alternativeQSsOptionsValues,
                                     const BaseEnvironment * env, const char * prefix);
#else
  //! Default constructor.
  /*! Assigns the default suite of options to the Monte Carlo sequence generator.*/
  McOptionsValues            ();

  //! Prefix constructor for reading input options from a file
  McOptionsValues(const BaseEnvironment * env, const char * prefix);
#endif
  //! Copy constructor.
  /*! It assigns the same options values from  \c src to \c this.*/
  McOptionsValues            (const McOptionsValues& src);

  //! Set default values for parameter options
  void set_defaults();

  //! Set parameter option names to begin with prefix
  void set_prefix(const std::string& prefix);

  //! Given prefix, read the input file for parameters named "prefix"+*
  void parse(const BaseEnvironment& env, const std::string& prefix);

  //! Destructor
  virtual ~McOptionsValues            ();
  //@}

  //! @name Set methods
  //@{
  //! Assignment operator; it copies \c rhs to \c this.
  McOptionsValues& operator= (const McOptionsValues& rhs);
  //@}

  std::string                        m_prefix;

  //! If non-empty string, print options and values to output file
  std::string                        m_help;

  std::string                        m_dataOutputFileName;
  std::set<unsigned int>             m_dataOutputAllowedSet;

  unsigned int                       m_pseqDataOutputPeriod;
  std::string                        m_pseqDataOutputFileName;
  std::string                        m_pseqDataOutputFileType;
  std::set<unsigned int>             m_pseqDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_pseqComputeStats;
#endif

  std::string                        m_qseqDataInputFileName;
  std::string                        m_qseqDataInputFileType;
  unsigned int                       m_qseqSize;
  unsigned int                       m_qseqDisplayPeriod;
  bool                               m_qseqMeasureRunTimes;
  unsigned int                       m_qseqDataOutputPeriod;
  std::string                        m_qseqDataOutputFileName;
  std::string                        m_qseqDataOutputFileType;
  std::set<unsigned int>             m_qseqDataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  bool                               m_qseqComputeStats;
#endif

private:
#ifndef QUESO_DISABLE_BOOST_PROGRAM_OPTIONS
  ScopedPtr<BoostInputOptionsParser>::Type m_parser;
#endif  // QUESO_DISABLE_BOOST_PROGRAM_OPTIONS

  std::string                   m_option_help;
  std::string                   m_option_dataOutputFileName;
  std::string                   m_option_dataOutputAllowedSet;

  std::string                   m_option_pseq_dataOutputPeriod;
  std::string                   m_option_pseq_dataOutputFileName;
  std::string                   m_option_pseq_dataOutputFileType;
  std::string                   m_option_pseq_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_pseq_computeStats;
#endif

  std::string                   m_option_qseq_dataInputFileName;
  std::string                   m_option_qseq_dataInputFileType;
  std::string                   m_option_qseq_size;
  std::string                   m_option_qseq_displayPeriod;
  std::string                   m_option_qseq_measureRunTimes;
  std::string                   m_option_qseq_dataOutputPeriod;
  std::string                   m_option_qseq_dataOutputFileName;
  std::string                   m_option_qseq_dataOutputFileType;
  std::string                   m_option_qseq_dataOutputAllowedSet;
#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  std::string                   m_option_qseq_computeStats;
#endif

  //! Copies the option values from \c src to \c this.
  void copy(const McOptionsValues& src);

  void checkOptions();

  friend std::ostream & operator<<(std::ostream & os,
      const McOptionsValues & obj);

#ifdef QUESO_USES_SEQUENCE_STATISTICAL_OPTIONS
  SsOptionsValues             m_alternativePSsOptionsValues;
  SsOptionsValues             m_alternativeQSsOptionsValues;
#endif
};

}  // End namespace QUESO

#endif // UQ_MOC_SG_OPTIONS_H
