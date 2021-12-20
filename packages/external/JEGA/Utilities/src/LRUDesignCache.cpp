/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class LRUDesignCache.

    NOTES:

        See notes of LRUDesignCache.hpp.

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

        Mon Apr 20 08:37:40 2015 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the LRUDesignCache class.
 */




/*
================================================================================
Includes
================================================================================
*/
#include <limits>
#include "../include/LRUDesignCache.hpp"







/*
================================================================================
Namespace Using Directives
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
void
LRUDesignCache::indexed_list::remove(
    const key_type& key
    )
{
    EDDY_FUNC_DEBUGSCOPE
    set_iter_t it(this->_indices.find(key));
    if(it == this->_indices.end()) return;
    this->_data.erase(it->second);
    this->_indices.erase(it);
}


DesignDVSortSet::size_type
LRUDesignCache::erase_exacts(
    const key_type key
    )
{
    EDDY_FUNC_DEBUGSCOPE

    // store the initial size so we can easily return the number erased.
    const size_type isize = this->_data.size();

    // start by bounding the search region.
    DesignDVSortSet::iterator_pair bounds(this->_data.equal_range(key));

    // now, look at every member in that range.
    // when any exact matches to key are found, erase them.
    for(; bounds.first!=bounds.second;)
    {
        if(*bounds.first == key) this->erase(bounds.first++);
        else ++bounds.first;
    }

    // return the difference in size from the beginning.
    return isize - this->_data.size();
}









/*
================================================================================
Subclass Visible Methods
================================================================================
*/

void
LRUDesignCache::manage_size(
    )
{
    EDDY_FUNC_DEBUGSCOPE

    if(this->_doCache)
	{
		while(this->_data.size() > this->_maxSize)
		{
			this->_data.erase(this->_lruList.front());
			this->_lruList.remove_first();
		}
	}
}

/*
================================================================================
Subclass Overridable Methods
================================================================================
*/








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

LRUDesignCache::LRUDesignCache(
    const size_type& maxSize
    ) :
        _lruList(maxSize == std::numeric_limits<size_type>::max() ? 0 : 100),
        _data(),
        _maxSize(maxSize),
        _doCache(maxSize < std::numeric_limits<size_type>::max())
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

