/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Inline methods of class LRUDesignCache.

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
 * \brief Contains the inline methods of the LRUDesignCache class.
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








/*
================================================================================
Inline Public Methods
================================================================================
*/

inline
void
LRUDesignCache::indexed_list::on_accessed(
    const key_type& key
    )
{
    this->remove(key);
    this->add(key);
}

inline
void
LRUDesignCache::indexed_list::add(
    const key_type& key
    )
{
    this->_indices[key] = this->_data.insert(this->_data.end(), key);
}


inline
void
LRUDesignCache::indexed_list::remove_first(
    )
{
    this->_indices.erase(*this->_data.begin());
    this->_data.pop_front();
}

inline
std::size_t
LRUDesignCache::indexed_list::size(
    )
{
    return this->_data.size();
}


inline
void
LRUDesignCache::indexed_list::clear(
    )
{
    this->_indices.clear();
    this->_data.clear();
}


inline
LRUDesignCache::const_reference
LRUDesignCache::indexed_list::front(
    )
{
    return this->_data.front();
}


inline
LRUDesignCache::const_reference
LRUDesignCache::indexed_list::back(
    )
{
    return this->_data.back();
}

inline
const DesignDVSortSet&
LRUDesignCache::DVSortSet(
    ) const
{
    return this->_data;
}

inline
const LRUDesignCache::size_type&
LRUDesignCache::max_size(
    ) const
{
    return this->_maxSize;
}
    
inline
void
LRUDesignCache::max_size(
    const size_type& newSize
    )
{
    this->_maxSize = newSize;
    this->manage_size();
}

inline
LRUDesignCache::size_type
LRUDesignCache::size(
    ) const
{
    return this->_data.size();
}

inline
bool
LRUDesignCache::empty(
    ) const
{
    return this->_data.empty();
}

inline
void
LRUDesignCache::insert(
    const key_type& key
    )
{
    this->_data.insert(key);
    this->_lruList.add(key);
    this->manage_size();
}

inline
void
LRUDesignCache::erase(
    const_iterator loc
    )
{
    if(this->_doCache) this->_lruList.remove(*loc);
    this->_data.erase(loc);
}

inline
LRUDesignCache::size_type
LRUDesignCache::erase(
    const key_type& key
    )
{
    if(this->_doCache) this->_lruList.remove(key);
    return this->_data.erase(key);
}

inline
void
LRUDesignCache::erase(
    const_iterator b,
    const_iterator e
    )
{
    if(this->_doCache) for(const_iterator it(b); it!=e; ++it)
        this->_lruList.remove(*it);

    this->_data.erase(b, e);
}

inline
DesignDVSortSet::const_iterator
LRUDesignCache::find_exact(
    const key_type& key
    ) const
{
    return this->on_maybe_accessed(this->_data.find_exact(key));
}

inline
DesignDVSortSet::const_iterator
LRUDesignCache::find_not_exact(
    const key_type& key
    ) const
{
    return this->on_maybe_accessed(this->_data.find_not_exact(key));
}

inline
DesignDVSortSet::iterator
LRUDesignCache::find_exact(
    const key_type& key
    )
{
    return this->on_maybe_accessed(this->_data.find_exact(key));
}

inline
DesignDVSortSet::iterator
LRUDesignCache::find_not_exact(
    const key_type& key
    )
{
    return this->on_maybe_accessed(this->_data.find_not_exact(key));
}

inline
DesignDVSortSet::size_type
LRUDesignCache::count_non_unique(
    ) const
{
    return this->_data.count_non_unique();
}

inline
std::ostream&
LRUDesignCache::stream_out(
    std::ostream& stream
    ) const
{
    return this->_data.stream_out(stream);
}

inline
std::ostream&
LRUDesignCache::stream_out(
    const key_type des,
    std::ostream& stream
    )
{
    return DesignDVSortSet::stream_out(des, stream);
}

inline
void
LRUDesignCache::flush(
    )
{
    this->_data.flush();
    this->_lruList.clear(); // doCache check not needed.
}

inline
DesignDVSortSet::const_iterator
LRUDesignCache::test_for_clone(
    const key_type key
    ) const
{
    return this->on_maybe_accessed(this->_data.test_for_clone(key));
}

inline
DesignDVSortSet::size_type
LRUDesignCache::test_within_list_for_clones(
    ) const
{
    return this->_data.test_within_list_for_clones();
}

inline
LRUDesignCache::const_iterator
LRUDesignCache::begin(
    ) const
{
    return this->_data.begin();
}

inline
LRUDesignCache::iterator
LRUDesignCache::begin(
    )
{
    return this->_data.begin();
}

inline
LRUDesignCache::const_iterator
LRUDesignCache::end(
    ) const
{
    return this->_data.end();
}

inline
LRUDesignCache::iterator
LRUDesignCache::end(
    )
{
    return this->_data.end();
}




//
//const LRUDesignCache::key_type&
//LRUDesignCache::operator [](
//    const key_type& key
//    ) const
//{
//    return this->retrieve(key);
//}
//
//
//void
//LRUDesignCache::insert(
//    const key_type& key
//    )
//{
//    this->_data.insert(key);
//    this->_lruList.on_accessed(key);
//    this->manage_size();
//}
//
//
//const LRUDesignCache::key_type&
//LRUDesignCache::retrieve(
//    const key_type& key
//    ) const
//{
//    const_iterator it(this->_data.find(key));
//    
//    if(it == this->_data.end()) throw std::exception(
//        "Key not found in LRU cache."
//        );
//
//    this->_lruList.on_accessed(key);
//    
//    return *it;
//}

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
template <typename ItT>
const ItT&
LRUDesignCache::on_accessed(
    const ItT& it
    ) const
{
    if(this->_doCache) this->_lruList.on_accessed(*it);
    return it;
}

template <typename ItT>
const ItT&
LRUDesignCache::on_maybe_accessed(
    const ItT& it
    ) const
{
    if(this->_doCache && (it != this->_data.end()))
        this->_lruList.on_accessed(*it);
    return it;
}








/*
================================================================================
Inline Structors
================================================================================
*/

inline
LRUDesignCache::indexed_list::indexed_list(
    size_type initCap 
    ) :
        _data(),
#ifdef JEGA_HAVE_BOOST
        _indices(initCap)
#else
        _indices()
#endif
{}














/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA

