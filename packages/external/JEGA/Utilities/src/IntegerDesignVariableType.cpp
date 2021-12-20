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

#include <cmath>
#include <cfloat>
#include <limits>
#include <type_traits>
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
struct SignedRounder
{
    static
    double
    Round(
        const double& value,
        const double& min,
        const double& max
        )
    {
        EDDY_FUNC_DEBUGSCOPE

        const double pct = (value - min) / (max - min);
        const double adjMinVal = min - 0.5;

        // There is a slight chance that the provided value is ub and so will
        // wind up rounding to 1 greater than ub after adjustment.  All calling
        // sites to this function will correct for that.  Alternatives would be
        // to avoid that by shortening the adjustment to as close to 0.5 as we
        // can such that adjMinVal != min-0.5 and adjMaxVal != max + 0.5.
        return Math::Round(
            (pct * ((max + 0.5) - adjMinVal)) + adjMinVal
            );
    }
};

struct IntegralRounder
{
    static
    double
    Round(
        const double& value,
        const double& min,
        const double& max
        )
    {
        EDDY_FUNC_DEBUGSCOPE
        return value;
    }
};


struct UnsignedRounder
{
    static
    double
    Round(
        const double& value,
        const double& min,
        const double& max
        )
    {
        EDDY_FUNC_DEBUGSCOPE

        // If the range is such that subtracting 0.5 from the min will still be
        // an unsigned value, then use the signed rounder.
        if(min >= 0.5) return SignedRounder::Round(value, min, max);

        // Otherwise, use a different technique.  NOTE: Assume there are no
        // unsigned real valued types.

        // There is a slight chance that the provided value is ub and so will
        // wind up as 1 greater than ub after adjustment.  All calling
        // sites to this function will correct for that.
        const double pct = (value - min) / (max - min);
        return std::trunc(pct * (max + 1.0 - min) + min);
    }
};





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

    typedef std::conditional<
        std::is_integral<var_rep_t>::value, IntegralRounder, UnsignedRounder
        >::type rndr_t;

	if(rep == -std::numeric_limits<var_rep_t>::max())
		return -std::numeric_limits<var_rep_t>::max();

    return this->GetNature().GetNearestValidRep(
        static_cast<var_rep_t>(
            rndr_t::Round(rep, this->GetMinRep(), this->GetMaxRep())
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
        SignedRounder::Round(temp, this->GetMinValue(), this->GetMaxValue())
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
