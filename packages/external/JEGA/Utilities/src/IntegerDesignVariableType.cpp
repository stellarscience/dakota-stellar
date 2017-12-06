/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class IntegerDesignVariableType.

    NOTES:

        See notes of IntegerDesignVariableType.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Jun 03 08:55:28 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the IntegerDesignVariableType class.
 */



/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <cfloat>
#include <limits>
#include <utilities/include/Math.hpp>
#include <../Utilities/include/Logging.hpp>
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/DesignVariableInfo.hpp>
#include <../Utilities/include/DesignVariableNatureBase.hpp>
#include <../Utilities/include/IntegerDesignVariableType.hpp>



/*
================================================================================
Namespace Using Directives
================================================================================
*/
using namespace std;
using namespace JEGA::Logging;
using namespace eddy::utilities;







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
IntegerDesignVariableType::ToString(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return string("Integer");
}

bool
IntegerDesignVariableType::SetPrecision(
    eddy::utilities::int16_t prec
    )
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(prec <= 0);
    eddy::utilities::int16_t prevPrec = this->GetPrecision();
    if(!this->DesignVariableTypeBase::SetPrecision(prec)) return false;

    if(prec > 0)
    {
        JEGALOG_II_G(lquiet(), this,
            ostream_entry(lquiet(), "Precision for integer design variable "
                "type must be <= 0.  Supplied value of ") << prec << " for "
                << this->GetDesignVariableInfo().GetLabel() << " rejected."
            )

        this->DesignVariableTypeBase::SetPrecision(prevPrec);
        return false;
    }

    return true;
}

bool
IntegerDesignVariableType::IsNatureLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return false;
}

bool
IntegerDesignVariableType::IsValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return Math::IsWhole(rep) &&
        this->DesignVariableTypeBase::IsValidRep(rep);
}

DesignVariableTypeBase*
IntegerDesignVariableType::Clone(
    DesignVariableInfo& forDVI
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return new IntegerDesignVariableType(*this, forDVI);
}

double
IntegerDesignVariableType::GetValueOf(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return Math::IsWhole(rep) ?
		this->GetNature().GetValueOf(rep) :
		-std::numeric_limits<double>::max();
}

double
IntegerDesignVariableType::GetNearestValidValue(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(Math::IsWhole(GetMinValue()));
    EDDY_ASSERT(Math::IsWhole(GetMaxValue()));

    if(value == -std::numeric_limits<double>::max()) return value;

    double temp = this->GetNature().GetNearestValidValue(value);

    if(Math::IsWhole(temp) && this->IsValueInBounds(temp)) return temp;
    return this->GetNearestValidValue(Math::Round(temp));
}

var_rep_t
IntegerDesignVariableType::GetNearestValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

	if(rep == -std::numeric_limits<var_rep_t>::max())
		return -std::numeric_limits<var_rep_t>::max();

    return this->GetNature().GetNearestValidRep(
        static_cast<var_rep_t>(
            this->ubround(rep, this->GetMinRep(), this->GetMaxRep())
            )
        );
}

double
IntegerDesignVariableType::GetRandomValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    const double temp = this->GetNature().GetRandomValue();

    // If temp is valid, return it.  Otherwise continue.
    if(this->IsValidValue(temp)) return temp;

    return this->GetNearestValidValue(
        this->ubround(temp, this->GetMinValue(), this->GetMaxValue())
        );
}

var_rep_t
IntegerDesignVariableType::GetRepOf(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(Math::IsWhole(value));

    if(!Math::IsWhole(value)) return -std::numeric_limits<var_rep_t>::max();
    return this->GetNature().GetRepOf(value);
}

var_rep_t
IntegerDesignVariableType::GetRandomRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    var_rep_t temp = this->GetNature().GetRandomRep();
    return this->GetNearestValidRep(temp);
}

var_rep_t
IntegerDesignVariableType::GetRandomRep(
    const RegionOfSpace& within
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    var_rep_t temp = this->GetNature().GetRandomRep(within);
    return this->GetNearestValidRep(temp);
}

void
IntegerDesignVariableType::SetMinValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(Math::IsWhole(value));

    if(!Math::IsWhole(value))
    {
        JEGALOG_II_G(lquiet(), this,
            ostream_entry(lquiet(), "Integral lower bound value ")
                << value << " being rounded off to" << Math::Round(value)
                << "."
            )
        value = Math::Round(value);
    }
    this->DesignVariableTypeBase::SetMinValue(value);
}


void
IntegerDesignVariableType::SetMaxValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(Math::IsWhole(value));

    if(!Math::IsWhole(value))
    {
        JEGALOG_II_G(lquiet(), this,
            ostream_entry(lquiet(), "Integral upper bound value ")
                << value << " being rounded off to" << Math::Round(value)
                << "."
            )
        value = Math::Round(value);
    }
    this->DesignVariableTypeBase::SetMaxValue(value);
}




/*
================================================================================
Private Methods
================================================================================
*/

double
IntegerDesignVariableType::ubround(
    const double& value,
    const double& min,
    const double& max
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    const double pct = (value - min) / (max - min);
    const double adjMinVal = min - 0.5;

    // There is a slight chance that the provided value is ub and so will adjust
    // wind up rounding to 1 greater than ub after adjustment.  All calling
    // sites to this function will correct for that.  Alternatives would be
    // to avoid that by shortening the adjustment to as close to 0.5 as we can
    // such that adjMinVal != min-0.5 and adjMaxVal != max + 0.5.
    return Math::Round(
        (pct * ((max + 0.5) - adjMinVal)) + adjMinVal
        );
}





/*
================================================================================
Structors
================================================================================
*/

IntegerDesignVariableType::IntegerDesignVariableType(
    DesignVariableInfo& info
    ) :
        DesignVariableTypeBase(info)
{
    EDDY_FUNC_DEBUGSCOPE
    this->GetNature().SetPrecision(0);
}

IntegerDesignVariableType::IntegerDesignVariableType(
    const IntegerDesignVariableType& copy,
    DesignVariableInfo& info
    ) :
        DesignVariableTypeBase(copy, info)
{
    EDDY_FUNC_DEBUGSCOPE
    this->GetNature().SetPrecision(0);
}






/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA
