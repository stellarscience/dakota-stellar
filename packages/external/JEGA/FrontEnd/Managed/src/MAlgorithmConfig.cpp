/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Implementation of class MAlgorithmConfig.

    NOTES:

        See notes of MAlgorithmConfig.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Mon Feb 20 07:30:18 2006 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the MAlgorithmConfig class.
 */




/*
================================================================================
Includes
================================================================================
*/
#include <stdafx.h>
#include <MEvaluator.hpp>
#include <ManagedUtils.hpp>
#include <MAlgorithmConfig.hpp>
#include <MBasicParameterDatabaseImpl.hpp>

#pragma unmanaged
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/BasicParameterDatabaseImpl.hpp>
#pragma managed






/*
================================================================================
Namespace Using Directives
================================================================================
*/
using namespace JEGA::FrontEnd;
using namespace JEGA::Utilities;








/*
================================================================================
Begin Namespace
================================================================================
*/
namespace JEGA {
    namespace FrontEnd {
        namespace Managed {





/*
================================================================================
Static Member Data Definitions
================================================================================
*/








/*
================================================================================
Mutators
================================================================================
*/

bool
MAlgorithmConfig::SetAlgorithmName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetAlgorithmName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetPrintPopEachGen(
    bool doPrint
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetPrintPopEachGen(doPrint);
}

bool
MAlgorithmConfig::SetAlgorithmType(
    AlgType algType
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetAlgorithmType(
        static_cast<AlgorithmConfig::AlgType>(algType)
        );
}

bool
MAlgorithmConfig::SetLoggingFilename(
    System::String^ filename
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetLoggingFilename(ToStdStr(filename));
}

bool
MAlgorithmConfig::SetDefaultLoggingLevel(
    JEGA::Logging::LogLevel defLevel
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetDefaultLoggingLevel(defLevel);
}








/*
================================================================================
Accessors
================================================================================
*/

MParameterDatabase^
MAlgorithmConfig::GetParameterDB(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB;
}








/*
================================================================================
Public Methods
================================================================================
*/

const AlgorithmConfig&
MAlgorithmConfig::Manifest(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return *this->_theConfig;
}

bool
MAlgorithmConfig::SetMutatorName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->SetMutatorName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetConvergerName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetConvergerName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetCrosserName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetCrosserName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetNichePressureApplicatorName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetNichePressureApplicatorName(
        ToStdStr(name)
        );
}

bool
MAlgorithmConfig::SetFitnessAssessorName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetFitnessAssessorName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetSelectorName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetSelectorName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetInitializerName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetInitializerName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetPostProcessorName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetPostProcessorName(ToStdStr(name));
}

bool
MAlgorithmConfig::SetMainLoopName(
    System::String^ name
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->SetMainLoopName(ToStdStr(name));
}







/*
================================================================================
Subclass Visible Methods
================================================================================
*/








/*
================================================================================
Subclass Overridable Methods
================================================================================
*/
void
MAlgorithmConfig::DoDispose(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    delete this->_theConfig;
    this->_theConfig = 0x0;

    // The garbage collector will dispose it.  We need only
    // stop referencing it.
    this->_thePDB = nullptr;
}








/*
================================================================================
Private Methods
================================================================================
*/








/*
================================================================================
Structors
================================================================================
*/

MAlgorithmConfig::MAlgorithmConfig(
    MParameterDatabase^ thePDB,
    MEvaluator^ theEvaluator
    ) :
        _theConfig(0x0),
        _thePDB(thePDB)
{
    EDDY_FUNC_DEBUGSCOPE
    this->_theConfig = new AlgorithmConfig(
        theEvaluator->GetTheEvaluatorCreator(), this->_thePDB->Manifest()
        );
}

MAlgorithmConfig::MAlgorithmConfig(
    MEvaluator^ theEvaluator
    ) :
        _theConfig(0x0),
        _thePDB(gcnew MBasicParameterDatabaseImpl())
{
    EDDY_FUNC_DEBUGSCOPE
    this->_theConfig = new AlgorithmConfig(
        theEvaluator->GetTheEvaluatorCreator(), this->_thePDB->Manifest()
        );
}

MAlgorithmConfig::~MAlgorithmConfig(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->DoDispose();
}








/*
================================================================================
End Namespace
================================================================================
*/
        } // namespace Managed
    } // namespace FrontEnd
} // namespace JEGA


