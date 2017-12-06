/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Implementation of class MProblemConfig.

    NOTES:

        See notes of MProblemConfig.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Wed Feb 08 13:40:27 2006 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the MProblemConfig class.
 */




/*
================================================================================
Includes
================================================================================
*/
#include <stdafx.h>
#include <ManagedUtils.hpp>
#include <MProblemConfig.hpp>

#pragma unmanaged
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../FrontEnd/Core/include/ProblemConfig.hpp>
#pragma managed




/*
================================================================================
Namespace Using Directives
================================================================================
*/
using namespace System;
using namespace JEGA::FrontEnd;





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







/*
================================================================================
Public Methods
================================================================================
*/

const ProblemConfig&
MProblemConfig::Manifest(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return *_theConfig;
}

void
MProblemConfig::SetDiscardTracking(
    bool track
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->_theConfig->SetDiscardTracking(track);
}

bool
MProblemConfig::GetDiscardTracking(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_theConfig->GetDiscardTracking();
}

void
MProblemConfig::SetMaxGuffSize(
    System::UInt32 mgs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->_theConfig->SetMaxGuffSize(mgs);
}

System::UInt32
MProblemConfig::GetMaxGuffSize(
    )
{
    EDDY_FUNC_DEBUGSCOPE
	return static_cast<System::UInt32>(this->_theConfig->GetMaxGuffSize());
}

void
MProblemConfig::SetMaxDiscardCacheSize(
    System::UInt32 maxSize
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->_theConfig->SetMaxDiscardCacheSize(maxSize);
}

System::UInt32
MProblemConfig::GetMaxDiscardCacheSize(
    )
{
    EDDY_FUNC_DEBUGSCOPE
	return static_cast<System::UInt32>(
        this->_theConfig->GetMaxDiscardCacheSize()
        );
}

bool
MProblemConfig::AddContinuumRealVariable(
    System::String^ label,
    double lowerBound,
    double upperBound,
    int precision
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddContinuumRealVariable(
        ToStdStr(label), lowerBound, upperBound, precision
        );
}

bool
MProblemConfig::AddDiscreteRealVariable(
    System::String^ label,
    DoubleVector^ values
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddDiscreteRealVariable(
        ToStdStr(label), ToStdDoubleVector(values)
        );
}

bool
MProblemConfig::AddContinuumIntegerVariable(
    System::String^ label,
    int lowerBound,
    int upperBound
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddContinuumIntegerVariable(
        ToStdStr(label), lowerBound, upperBound
        );
}

bool
MProblemConfig::AddDiscreteIntegerVariable(
    System::String^ label,
    IntVector^ values
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddDiscreteIntegerVariable(
        ToStdStr(label), ToStdIntVector(values)
        );
}

bool
MProblemConfig::AddBooleanVariable(
    System::String^ label
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddBooleanVariable(ToStdStr(label));
}

bool
MProblemConfig::AddLinearMinimizeObjective(
    System::String^ label,
    DoubleVector^ coeffs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddLinearMinimizeObjective(
        ToStdStr(label), ToStdDoubleVector(coeffs)
        );
}

bool
MProblemConfig::AddLinearMaximizeObjective(
    System::String^ label,
    DoubleVector^ coeffs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddLinearMaximizeObjective(
        ToStdStr(label), ToStdDoubleVector(coeffs)
        );
}

bool
MProblemConfig::AddLinearSeekValueObjective(
    System::String^ label,
    double value,
    DoubleVector^ coeffs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddLinearSeekValueObjective(
        ToStdStr(label), value, ToStdDoubleVector(coeffs)
        );
}

bool
MProblemConfig::AddLinearSeekRangeObjective(
    System::String^ label,
    double lowerBound,
    double upperBound,
    DoubleVector^ coeffs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddLinearSeekRangeObjective(
        ToStdStr(label), lowerBound, upperBound,
        ToStdDoubleVector(coeffs)
        );
}

bool
MProblemConfig::AddNonlinearMinimizeObjective(
    System::String^ label
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddNonlinearMinimizeObjective(ToStdStr(label));
}

bool
MProblemConfig::AddNonlinearMaximizeObjective(
    System::String^ label
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddNonlinearMaximizeObjective(ToStdStr(label));
}

bool
MProblemConfig::AddNonlinearSeekValueObjective(
    System::String^ label,
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddNonlinearSeekValueObjective(ToStdStr(label), value);
}

bool
MProblemConfig::AddNonlinearSeekRangeObjective(
    System::String^ label,
    double lowerBound,
    double upperBound
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddNonlinearSeekRangeObjective(
        ToStdStr(label), lowerBound, upperBound
        );
}

bool
MProblemConfig::AddLinearInequalityConstraint(
    System::String^ label,
    con_val_t upperLimit,
    DoubleVector^ coeffs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddLinearInequalityConstraint(
        ToStdStr(label), upperLimit, ToStdDoubleVector(coeffs)
        );
}

bool
MProblemConfig::AddLinearEqualityConstraint(
    System::String^ label,
    con_val_t target,
    double allowedViol,
    DoubleVector^ coeffs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddLinearEqualityConstraint(
        ToStdStr(label), target, allowedViol, ToStdDoubleVector(coeffs)
        );
}

bool
MProblemConfig::AddLinearTwoSidedInequalityConstraint(
    System::String^ label,
    double lowerLimit,
    double upperLimit,
    DoubleVector^ coeffs
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddLinearTwoSidedInequalityConstraint(
        ToStdStr(label), lowerLimit, upperLimit, ToStdDoubleVector(coeffs)
        );
}

bool
MProblemConfig::AddNonlinearInequalityConstraint(
    System::String^ label,
    double upperLimit
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddNonlinearInequalityConstraint(
        ToStdStr(label), upperLimit
        );
}

bool
MProblemConfig::AddNonlinearEqualityConstraint(
    System::String^ label,
    double target,
    double allowedViol
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddNonlinearEqualityConstraint(
        ToStdStr(label), target, allowedViol
        );
}

bool
MProblemConfig::AddNonlinearTwoSidedInequalityConstraint(
    System::String^ label,
    double lowerLimit,
    double upperLimit
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return _theConfig->AddNonlinearTwoSidedInequalityConstraint(
        ToStdStr(label), lowerLimit, upperLimit
        );
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
MProblemConfig::DoDispose(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    delete _theConfig;
    _theConfig = 0x0;
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


MProblemConfig::MProblemConfig(
    ) :
        _theConfig(new ProblemConfig())
{
    EDDY_FUNC_DEBUGSCOPE
}

MProblemConfig::~MProblemConfig(
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


