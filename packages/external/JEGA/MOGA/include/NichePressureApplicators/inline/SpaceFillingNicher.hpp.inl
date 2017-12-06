/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Inline methods of class SpaceFillingNicher.

    NOTES:

        See notes of SpaceFillingNicher.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public
        License as published by the Free Software Foundation; either
        version 2.1 of the License, or (at your option) any later version.
        
        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Lesser General Public License for more details.
        
        You should have received a copy of the GNU Lesser General Public
        License along with this library; if not, write to the Free Software
        Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
        USA

    VERSION:

        2.7.0

    CHANGES:

        Thu Aug 20 08:50:56 2015 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the inline methods of the SpaceFillingNicher class.
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
    namespace Algorithms {





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
const std::size_t&
SpaceFillingNicher::GetNumDesigns2Keep(
    ) const
{
    return this->_nDes2Keep;
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

double
SpaceFillingNicher::NormalizedObjVal(
    const JEGA::Utilities::Design& des,
    const size_t& of,
	const eddy::utilities::extremes<obj_val_t>& popExtremes
    )
{
    EDDY_FUNC_DEBUGSCOPE
	return (des.GetObjective(of) - popExtremes.get_min(of)) /
        popExtremes.get_range(of);
}



/*
================================================================================
Inline Private Methods
================================================================================
*/

inline
std::size_t
SpaceFillingNicher::hash(
	const JEGA::Utilities::Design& des1, 
	const JEGA::Utilities::Design& des2
	)
{
	const std::size_t d1h = des1.GetID();
	return d1h ^ (des2.GetID() + 0x9e3779b9 + (d1h<<6) + (d1h>>2));
}






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
    } // namespace Algorithms
} // namespace JEGA

