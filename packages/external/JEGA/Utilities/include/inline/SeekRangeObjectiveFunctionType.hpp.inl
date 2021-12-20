/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Inline methods of class SeekRangeObjectiveFunctionType.

    NOTES:

        See notes of SeekRangeObjectiveFunctionType.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Sun Oct 12 17:36:14 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the inline methods of the SeekRangeObjectiveFunctionType
 *        class.
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
SeekRangeObjectiveFunctionType::SetLowerBound(
    obj_val_t value
    )
{
    this->_lowerBound = value;
}

inline
void
SeekRangeObjectiveFunctionType::SetUpperBound(
    obj_val_t value
    )
{
    this->_upperBound = value;
}





/*
================================================================================
Inline Accessors
================================================================================
*/

inline
obj_val_t
SeekRangeObjectiveFunctionType::GetLowerBound(
    ) const
{
    return this->_lowerBound;
}

inline
obj_val_t
SeekRangeObjectiveFunctionType::GetUpperBound(
    ) const
{
    return this->_upperBound;
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

