/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class ContinuumDesignVariableNature.

    NOTES:

        See notes of ContinuumDesignVariableNature.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Jun 03 08:55:45 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the ContinuumDesignVariableNature
 *        class.
 */




/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <cfloat>
#include <type_traits>
#include <utilities/include/Math.hpp>
#include <../Utilities/include/Logging.hpp>
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/DesignVariableInfo.hpp>
#include <utilities/include/RandomNumberGenerator.hpp>
#include <../Utilities/include/ContinuumDesignVariableNature.hpp>


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
struct YesIntCaster
{
    inline static long long Cast(var_rep_t val)
    {
        return static_cast<long long>(val);
    }
};

struct NoIntCaster
{
    inline static var_rep_t Cast(var_rep_t val)
    {
        return val;
    }
};

struct IntRepTRNG
{
    inline static var_rep_t RNG(var_rep_t lb, var_rep_t ub)
    {
        typedef std::conditional<
            std::is_integral<var_rep_t>::value, NoIntCaster, YesIntCaster
            >::type caster;

        return static_cast<var_rep_t>(RandomNumberGenerator::UniformInt(
            caster::Cast(lb), caster::Cast(ub)
            ));
    }
};


struct RealRepTRNG
{
    inline static var_rep_t RNG(var_rep_t lb, var_rep_t ub)
    {
        return static_cast<var_rep_t>(
            RandomNumberGenerator::UniformReal(lb, ub)
            );
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
ContinuumDesignVariableNature::ToString(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return string("Continuum");
}

DesignVariableNatureBase*
ContinuumDesignVariableNature::Clone(
    DesignVariableTypeBase& forType
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return new ContinuumDesignVariableNature(*this, forType);
}

var_rep_t
ContinuumDesignVariableNature::GetMaxRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
	return this->GetRepOf(this->_maxVal);
}

var_rep_t
ContinuumDesignVariableNature::GetMinRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
	return this->GetRepOf(this->_minVal);
}

var_rep_t
ContinuumDesignVariableNature::GetRandomRep(
    var_rep_t lb,
    var_rep_t ub
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    typedef std::conditional<
        std::is_integral<var_rep_t>::value, IntRepTRNG, RealRepTRNG
        >::type rng_t;

    return static_cast<var_rep_t>(rng_t::RNG(lb, ub));
}

var_rep_t
ContinuumDesignVariableNature::GetRepOf(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return static_cast<var_rep_t>(value);
}

double
ContinuumDesignVariableNature::GetRandomValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return RandomNumberGenerator::UniformReal(
        this->GetMinValue(), this->GetMaxValue()
        );
}

double
ContinuumDesignVariableNature::GetMaxValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_maxVal;
}

double
ContinuumDesignVariableNature::GetMinValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_minVal;
}

double
ContinuumDesignVariableNature::GetValueOf(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return static_cast<double>(rep);
}

var_rep_t
ContinuumDesignVariableNature::GetDistanceBetweenReps(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return static_cast<var_rep_t>(
		Math::Pow(10, static_cast<double>(-this->GetPrecision()))
		);
}

bool
ContinuumDesignVariableNature::AddDiscreteValue(
    double
    )
{
    EDDY_FUNC_DEBUGSCOPE

    JEGALOG_II_G_F(this,
        text_entry(lfatal(), this->GetDesignVariableInfo().GetLabel() +
            ": Continuum natured variable cannot accept discrete values.")
        )

    return false;
}

void
ContinuumDesignVariableNature::ClearDiscreteValues(
    )
{
    EDDY_FUNC_DEBUGSCOPE

    JEGALOG_II_G_F(this,
        text_entry(lfatal(), this->GetDesignVariableInfo().GetLabel() +
            ": Cannot clear discrete values for continuum natured variable.")
        )
}

bool
ContinuumDesignVariableNature::RemoveDiscreteValue(
    double
    )
{
    EDDY_FUNC_DEBUGSCOPE

    JEGALOG_II_G_F(this,
        text_entry(lfatal(), this->GetDesignVariableInfo().GetLabel() +
            ": Continuum natured variable has no discrete values to remove.")
        )

    return false;
}

void
ContinuumDesignVariableNature::SetMaxValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(value != -std::numeric_limits<double>::max());
    this->_maxVal = value;
}

void
ContinuumDesignVariableNature::SetMinValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(value != -std::numeric_limits<double>::max());
    this->_minVal = value;
}

bool
ContinuumDesignVariableNature::IsDiscreteValueLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return true;
}

bool
ContinuumDesignVariableNature::IsValueInBounds(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(this->GetMinValue() <= this->GetMaxValue());
    return (value >= this->GetMinValue()) && (value <= this->GetMaxValue());
}

bool
ContinuumDesignVariableNature::IsRepInBounds(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(this->GetMinRep() <= this->GetMaxRep());
    return (rep >= this->GetMinRep()) && (rep <= this->GetMaxRep());
}

bool
ContinuumDesignVariableNature::IsOutOfBoundsDefined(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return true;
}

bool
ContinuumDesignVariableNature::IsPrecisionLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return false;
}

bool
ContinuumDesignVariableNature::IsValidValue(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->IsValueInBounds(value);
}

bool
ContinuumDesignVariableNature::IsValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->IsRepInBounds(rep);
}

bool
ContinuumDesignVariableNature::IsDiscrete(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return false;
}

bool
ContinuumDesignVariableNature::IsContinuum(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return true;
}

double
ContinuumDesignVariableNature::GetNearestValidValue(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    // The best guess that this nature can make is
    // to return something in bounds
    return Math::Max(
        Math::Min(value, this->GetMaxValue()), this->GetMinValue()
        );
}

var_rep_t
ContinuumDesignVariableNature::GetNearestValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    // The best guess that this nature can make is
    // to return something in bounds
    return Math::Max(
        Math::Min(rep, this->GetMaxRep()),
        this->GetMinRep()
        );
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
ContinuumDesignVariableNature::ContinuumDesignVariableNature(
    DesignVariableTypeBase& type
    ) :
        DesignVariableNatureBase(type),
        _maxVal(0.0),
        _minVal(0.0)
{
    EDDY_FUNC_DEBUGSCOPE
}

ContinuumDesignVariableNature::ContinuumDesignVariableNature(
    const ContinuumDesignVariableNature& copy,
    DesignVariableTypeBase& type
    ) :
        DesignVariableNatureBase(copy, type),
        _maxVal(0.0),
        _minVal(0.0)
{
    EDDY_FUNC_DEBUGSCOPE
}

ContinuumDesignVariableNature::~ContinuumDesignVariableNature(
    )
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
