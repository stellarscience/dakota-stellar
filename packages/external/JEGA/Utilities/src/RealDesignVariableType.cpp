/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class RealDesignVariableType.

    NOTES:

        See notes of RealDesignVariableType.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Jun 03 08:55:22 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the RealDesignVariableType class.
 */





/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <cfloat>
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/RealDesignVariableType.hpp>
#include <../Utilities/include/DesignVariableNatureBase.hpp>


/*
================================================================================
Namespace Using Directives
================================================================================
*/
using namespace std;







/*
================================================================================
Begin Namespace
================================================================================
*/
namespace JEGA {
    namespace Utilities {







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
string
RealDesignVariableType::ToString(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return string("Real");
}

bool
RealDesignVariableType::IsNatureLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return false;
}

DesignVariableTypeBase*
RealDesignVariableType::Clone(
    DesignVariableInfo& forDVI
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return new RealDesignVariableType(*this, forDVI);
}

double
RealDesignVariableType::GetValueOf(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetNature().GetValueOf(rep);
}

double
RealDesignVariableType::GetRandomValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetNature().GetRandomValue();
}

var_rep_t
RealDesignVariableType::GetRepOf(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetNature().GetRepOf(value);
}

var_rep_t
RealDesignVariableType::GetRandomRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetNature().GetRandomRep();
}

var_rep_t
RealDesignVariableType::GetRandomRep(
    const RegionOfSpace& within
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    // get the bounds of this variable from within.
    return this->GetNature().GetRandomRep(within);
}

double
RealDesignVariableType::GetNearestValidValue(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetNature().GetNearestValidValue(value);
}

var_rep_t
RealDesignVariableType::GetNearestValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetNature().GetNearestValidRep(rep);
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
RealDesignVariableType::RealDesignVariableType(
    DesignVariableInfo& info
    ) :
        DesignVariableTypeBase(info)
{
    EDDY_FUNC_DEBUGSCOPE
}

RealDesignVariableType::RealDesignVariableType(
    const RealDesignVariableType& copy,
    DesignVariableInfo& info
    ) :
        DesignVariableTypeBase(copy, info)
{
    EDDY_FUNC_DEBUGSCOPE
}







/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA
