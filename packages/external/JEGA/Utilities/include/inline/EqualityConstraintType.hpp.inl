/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Inline methods of class EqualityConstraintType.

    NOTES:

        See notes of EqualityConstraintType.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Jun 10 08:43:33 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the inline methods of the EqualityConstraintType class.
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
EqualityConstraintType::SetTargetValue(
    con_val_t val
    )
{
    this->_value = val;
}

inline
void
EqualityConstraintType::SetAllowableViolation(
    double viol
    )
{
    this->_viol = viol;
}





/*
================================================================================
Inline Accessors
================================================================================
*/

inline
con_val_t
EqualityConstraintType::GetTargetValue(
    ) const
{
    return this->_value;
}

inline
double
EqualityConstraintType::GetAllowableViolation(
    ) const
{
    return this->_viol;
}






/*
================================================================================
Inline Public Methods
================================================================================
*/








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
