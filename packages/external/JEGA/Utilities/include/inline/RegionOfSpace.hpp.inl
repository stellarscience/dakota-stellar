/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Inline methods of class RegionOfSpace.

    NOTES:

        See notes of RegionOfSpace.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        2.0.0

    CHANGES:

        Thu Apr 13 07:42:48 2006 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the inline methods of the RegionOfSpace class.
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

inline
void
RegionOfSpace::SetLowerLimit(
    eddy::utilities::extremes<var_rep_t>::size_type dim,
    var_rep_t value
    )
{
    this->_limits.set_min(dim, value);
}

inline
void
RegionOfSpace::SetUpperLimit(
    eddy::utilities::extremes<var_rep_t>::size_type dim,
    var_rep_t value
    )
{
    this->_limits.set_max(dim, value);
}

inline
void
RegionOfSpace::SetLimits(
    eddy::utilities::extremes<var_rep_t>::size_type dim,
    var_rep_t lowerLimit,
    var_rep_t upperLimit
    )
{
    this->_limits.set_min(dim, lowerLimit);
    this->_limits.set_max(dim, upperLimit);
}







/*
================================================================================
Inline Accessors
================================================================================
*/








/*
================================================================================
Inline Public Methods
================================================================================
*/
inline
eddy::utilities::extremes<var_rep_t>::size_type
RegionOfSpace::Dimensionality(
    ) const
{
    return this->_limits.size();
}

inline
std::vector<var_rep_t>
RegionOfSpace::GetLowerLimits(
    ) const
{
    return this->_limits.get_mins();
}

inline
std::vector<var_rep_t>
RegionOfSpace::GetUpperLimits(
    ) const
{
    return this->_limits.get_maxs();
}

inline
var_rep_t
RegionOfSpace::GetLowerLimit(
    eddy::utilities::extremes<var_rep_t>::size_type dim
    ) const
{
    return this->_limits.get_min(dim);
}

inline
var_rep_t
RegionOfSpace::GetUpperLimit(
    eddy::utilities::extremes<var_rep_t>::size_type dim
    ) const
{
    return this->_limits.get_max(dim);
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

