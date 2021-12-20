/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Inline methods of class DesignVariableInfo.

    NOTES:

        See notes of DesignVariableInfo.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Jun 03 07:34:36 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the inline methods of the DesignVariableInfo class.
 */



/*
================================================================================
Includes
================================================================================
*/








/*
================================================================================
Begin Namespace
================================================================================
*/
namespace JEGA {
    namespace Utilities {







/*
================================================================================
Inline Mutators
================================================================================
*/








/*
================================================================================
Inline Accessors
================================================================================
*/
inline
DesignVariableTypeBase&
DesignVariableInfo::GetType(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(this->_type != 0x0)
    return *this->_type;
}







/*
================================================================================
Inline Public Methods
================================================================================
*/

inline
double
DesignVariableInfo::GetRepBoundViolation(
    const Design& des
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetRepBoundViolation(this->WhichRep(des));
}

inline
double
DesignVariableInfo::GetValueBoundViolation(
    const Design& des
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetValueBoundViolation(this->WhichValue(des));
}

inline
double
DesignVariableInfo::AssertPrecision(
    double val
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().AssertPrecision(val);
}

inline
double
DesignVariableInfo::GetDefaultValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetDefaultValue();
}

inline
double
DesignVariableInfo::GetMaxValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetMaxValue();
}

inline
double
DesignVariableInfo::GetMinValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetMinValue();
}

inline
double
DesignVariableInfo::GetValueRange(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetMaxValue() - this->GetMinValue();
}

inline
double
DesignVariableInfo::GetValueOf(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetValueOf(rep);
}

inline
double
DesignVariableInfo::GetRandomValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetRandomValue();
}

inline
var_rep_t
DesignVariableInfo::GetDefaultRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetDefaultRep();
}

inline
var_rep_t
DesignVariableInfo::GetMaxRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetMaxRep();
}

inline
var_rep_t
DesignVariableInfo::GetMinRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetMinRep();
}

inline
var_rep_t
DesignVariableInfo::GetRepRange(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetMaxRep() - this->GetMinRep();
}

inline
var_rep_t
DesignVariableInfo::GetRepOf(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetRepOf(value);
}

inline
var_rep_t
DesignVariableInfo::GetRandomRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetRandomRep();
}

inline
var_rep_t
DesignVariableInfo::GetRandomRep(
    const RegionOfSpace& within
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetRandomRep(within);
}

inline
var_rep_t
DesignVariableInfo::GetDistanceBetweenReps(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetDistanceBetweenReps();
}

inline
bool
DesignVariableInfo::IsDiscreteValueLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsDiscreteValueLocked();
}

inline
bool
DesignVariableInfo::IsValueInBounds(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsValueInBounds(value);
}

inline
bool
DesignVariableInfo::IsRepInBounds(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsRepInBounds(rep);
}

inline
bool
DesignVariableInfo::IsNatureLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsNatureLocked();
}

inline
bool
DesignVariableInfo::IsOutOfBoundsDefined(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsOutOfBoundsDefined();
}

inline
bool
DesignVariableInfo::IsPrecisionLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsPrecisionLocked();
}

inline
bool
DesignVariableInfo::IsValidValue(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsValidValue(value);
}

inline
bool
DesignVariableInfo::IsValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsValidRep(rep);
}

inline
std::string
DesignVariableInfo::GetTypeString(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().ToString();
}

inline
std::string
DesignVariableInfo::GetNatureString(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetNatureString();
}

inline
bool
DesignVariableInfo::AddDiscreteValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().AddDiscreteValue(value);
}

inline
void
DesignVariableInfo::ClearDiscreteValues(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->GetType().ClearDiscreteValues();
}

inline
bool
DesignVariableInfo::RemoveDiscreteValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().RemoveDiscreteValue(value);
}

inline
void
DesignVariableInfo::SetMinValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->GetType().SetMinValue(value);
}

inline
void
DesignVariableInfo::SetMaxValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->GetType().SetMaxValue(value);
}

inline
bool
DesignVariableInfo::SetPrecision(
    eddy::utilities::int16_t prec
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().SetPrecision(prec);
}

inline
eddy::utilities::int16_t
DesignVariableInfo::GetPrecision(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetPrecision();
}

inline
bool
DesignVariableInfo::SetNature(
    DesignVariableNatureBase* nature
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().SetNature(nature);
}

inline
DesignVariableNatureBase&
DesignVariableInfo::GetNature(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().GetNature();
}

inline
bool
DesignVariableInfo::IsDiscrete(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsDiscrete();
}

inline
bool
DesignVariableInfo::IsContinuum(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->GetType().IsContinuum();
}






/*
================================================================================
Inline Subclass Visible Methods
================================================================================
*/








/*
================================================================================
Inline Private Methods
================================================================================
*/








/*
================================================================================
Inline Structors
================================================================================
*/








/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA
