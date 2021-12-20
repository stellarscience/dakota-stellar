/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Implementation of class MDriver.

    NOTES:

        See notes of MDriver.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Thu Feb 09 08:47:20 2006 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the MDriver class.
 */




/*
================================================================================
Includes
================================================================================
*/
#include <stdafx.h>
#include <MDriver.hpp>
#include <MSolution.hpp>
#include <ManagedUtils.hpp>
#include <MAlgorithmConfig.hpp>
#include <MGeneticAlgorithm.hpp>

#pragma unmanaged
#include <../FrontEnd/Core/include/Driver.hpp>
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/DesignMultiSet.hpp>
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








/*
================================================================================
Accessors
================================================================================
*/

System::UInt32
MDriver::GetRandomSeed(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return Driver::GetRandomSeed();
}

System::String^
MDriver::GetGlobalLogFilename(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ManagedUtils::ToSysString(Driver::GetGlobalLogFilename());
}



/*
================================================================================
Public Methods
================================================================================
*/
Driver&
MDriver::Manifest(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return *this->_theApp;
}

bool
MDriver::IsJEGAInitialized(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return Driver::IsJEGAInitialized();
}

System::String^
MDriver::GetXType(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ManagedUtils::ToSysString(Driver::GetXType());
}

System::String^
MDriver::GetGType(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ManagedUtils::ToSysString(Driver::GetGType());
}

System::String^
MDriver::GetFType(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ManagedUtils::ToSysString(Driver::GetFType());
}

bool
MDriver::InitializeJEGA(
    System::String^ globalLogFilename,
    JEGA::Logging::LogLevel globalLogDefLevel,
    System::UInt32 rSeed,
    MAlgorithmConfig::FatalBehavior fatalBehavior,
	bool registerSignalHandlers
    )
{
    EDDY_FUNC_DEBUGSCOPE
    try
    {
        return Driver::InitializeJEGA(
            ToStdStr(globalLogFilename), globalLogDefLevel, rSeed,
            static_cast<JEGA::Logging::Logger::FatalBehavior>(fatalBehavior),
			registerSignalHandlers
            );
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
}

System::UInt32
MDriver::ReSeed(
    System::UInt32 rSeed
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return Driver::ReSeed(rSeed);
}

void
MDriver::ResetGlobalLoggingLevel(
    JEGA::Logging::LogLevel globalLogDefLevel
    )
{
    EDDY_FUNC_DEBUGSCOPE
    try
    {
        Driver::ResetGlobalLoggingLevel(globalLogDefLevel);
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
}

void
MDriver::FlushGlobalLogStreams(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    try
    {
        Driver::FlushGlobalLogStreams();
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
}

SolutionVector^
MDriver::ExecuteAlgorithm(
    MAlgorithmConfig^ config
    )
{
    EDDY_FUNC_DEBUGSCOPE

    try
    {
        DesignOFSortSet solutions(
            this->_theApp->ExecuteAlgorithm(config->Manifest())
            );

        SolutionVector^ ret =
            gcnew SolutionVector(static_cast<int>(solutions.size()));

        for(DesignOFSortSet::const_iterator it(solutions.begin());
            it!=solutions.end(); ++it)
                ret->Add(gcnew MSolution(**it));

        // As per the Driver instructions, we must flush the solution
        // set now that we are done with it.
        solutions.flush();

        Driver::FlushGlobalLogStreams();

        return ret;
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
}

MGeneticAlgorithm^
MDriver::InitializeAlgorithm(
    MAlgorithmConfig^ config
    )
{
    EDDY_FUNC_DEBUGSCOPE
    try {
        return gcnew MGeneticAlgorithm(
            this->_theApp->InitializeAlgorithm(config->Manifest())
            );
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
}

bool
MDriver::PerformNextIteration(
    MGeneticAlgorithm^ theGA
    )
{
    EDDY_FUNC_DEBUGSCOPE
    try {
        return this->_theApp->PerformNextIteration(&theGA->Manifest());
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
}

SolutionVector^
MDriver::FinalizeAlgorithm(
    MGeneticAlgorithm^ theGA
    )
{
    EDDY_FUNC_DEBUGSCOPE
    try {

        DesignOFSortSet solutions(
            this->_theApp->FinalizeAlgorithm(&theGA->Manifest())
            );

        // finalize will delete the GA.  So we will dissociate from it.
        theGA->DissociateGA();

        SolutionVector^ ret =
            gcnew SolutionVector(static_cast<int>(solutions.size()));

        for(DesignOFSortSet::const_iterator it(solutions.begin());
            it!=solutions.end(); ++it)
                ret->Add(gcnew MSolution(**it));

        // As per the Driver instructions, we must flush the solution
        // set now that we are done with it.
        solutions.flush();

        return ret;
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
}

void
MDriver::DestroyAlgorithm(
    MGeneticAlgorithm^ theGA
    )
{
    EDDY_FUNC_DEBUGSCOPE
    try {
        this->_theApp->DestroyAlgorithm(&theGA->Manifest());
        theGA->DissociateGA();
    }
    catch(const std::exception& ex)
    {
        throw gcnew
            System::Exception(ManagedUtils::ToSysString(ex.what()));
    }
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
MDriver::DoDispose(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    delete this->_theApp;
    this->_theApp = 0x0;
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


MDriver::MDriver(
    MProblemConfig^ probConfig
    ) :
        _theApp(new Driver(probConfig->Manifest()))
{
    EDDY_FUNC_DEBUGSCOPE
}

MDriver::~MDriver(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    DoDispose();
}






/*
================================================================================
End Namespace
================================================================================
*/
        } // namespace Managed
    } // namespace FrontEnd
} // namespace JEGA


