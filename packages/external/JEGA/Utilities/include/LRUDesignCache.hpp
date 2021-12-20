/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Definition of class LRUDesignCache.

    NOTES:

        See notes under "Document this File" section of this file.

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
 * \brief Contains the definition of the LRUDesignCache class.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_UTILITIES_LRUDESIGNCACHE_HPP
#define JEGA_UTILITIES_LRUDESIGNCACHE_HPP

#pragma once





/*
================================================================================
Includes
================================================================================
*/
#include <../Utilities/include/JEGAConfig.hpp>
#include <../Utilities/include/DesignMultiSet.hpp>

#ifdef JEGA_HAVE_BOOST
#   include <boost/unordered_map.hpp>
#else
#   include <map>
#endif

#include <list>
#include <cstddef>




/*
================================================================================
Pre-Namespace Forward Declares
================================================================================
*/
#ifdef JEGA_HAVE_BOOST
#define MAP_TYPE boost::unordered_map
#else
#define MAP_TYPE std::map
#endif







/*
================================================================================
Namespace Aliases
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
In-Namespace Forward Declares
================================================================================
*/







/*
================================================================================
In-Namespace File Scope Typedefs
================================================================================
*/








/*
================================================================================
Class Definition
================================================================================
*/
/**
 * \brief
 *
 *
 */
class LRUDesignCache
{
    /*
    ============================================================================
    Class Scope Typedefs
    ============================================================================
    */
    public:

        typedef DesignDVSortSet::key_type key_type;
        typedef DesignDVSortSet::size_type size_type;

        typedef DesignDVSortSet::iterator iterator;
        typedef DesignDVSortSet::const_iterator const_iterator;
        
        typedef DesignDVSortSet::reference reference;
        typedef DesignDVSortSet::const_reference const_reference;



    protected:


    private:

        class indexed_list
        {
            private:

                typedef std::list<key_type> list_t;
                typedef list_t::iterator list_iter_t;
                typedef MAP_TYPE<key_type, list_iter_t> map_t;
                typedef map_t::iterator set_iter_t;
                typedef map_t::const_iterator set_citer_t;

                list_t _data;
                map_t _indices;

            public:

                inline
                void
                on_accessed(
                    const key_type& key
                    );

                inline
                void
                add(
                    const key_type& key
                    );

                inline
                void
                remove_first(
                    );

                void
                remove(
                    const key_type& key
                    );

                inline
                std::size_t
                size(
                    );

                inline
                void
                clear(
                    );

                inline
                const key_type&
                front(
                    );

                inline
                const key_type&
                back(
                    );

                inline
                indexed_list(
                    map_t::size_type initCap 
                    );
        };


    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:

        mutable indexed_list _lruList;

        DesignDVSortSet _data;

        std::size_t _maxSize;

        bool _doCache;

    /*
    ============================================================================
    Mutators
    ============================================================================
    */
    public:


    protected:


    private:


    /*
    ============================================================================
    Accessors
    ============================================================================
    */
    public:


    protected:


    private:


    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        ///// Returns the last element in the container (constant)
        ///**
        // * Calling this on an empty container will result in a crash or at
        // * least undefined behavior.
        // *
        // * \return An immutable reference to the last element.
        // */
        //inline
        //DesignDVSortSet::const_reference
        //back(
        //    ) const;

        ///// Returns the first element in the container (constant)
        ///**
        // * Calling this on an empty container will result in a crash or at
        // * least undefined behavior.
        // *
        // * \return An immutable reference to the first element.
        // */
        //inline
        //DesignDVSortSet::const_reference
        //front(
        //    ) const;

        ///// Returns the last element in the container (constant)
        ///**
        // * Calling this on an empty container will result in a crash or at
        // * least undefined behavior.
        // *
        // * \return A reference to the last element.
        // */
        //inline
        //DesignDVSortSet::reference
        //back(
        //    );

        ///// Returns the first element in the container (constant)
        ///**
        // * Calling this on an empty container will result in a crash or at
        // * least undefined behavior.
        // *
        // * \return A reference to the first element.
        // */
        //inline
        //DesignDVSortSet::reference
        //front(
        //    );

        inline
        const DesignDVSortSet&
        DVSortSet(
            ) const;

        inline
        const_iterator
        begin(
            ) const;
            
        inline
        iterator
        begin(
            );
            
        inline
        const_iterator
        end(
            ) const;
            
        inline
        iterator
        end(
            );
            
        inline
        const size_type&
        max_size(
            ) const;
            
        inline
        void
        max_size(
            const size_type& newSize
            );

        inline
        size_type
        size(
            ) const;

        inline
        bool
        empty(
            ) const;

        inline
	    void
        erase(
            const_iterator loc
            );

        inline
	    void
        insert(
            const key_type& key
            );

        inline
	    size_type
        erase(
            const key_type& key
            );

        inline
	    void
        erase(
            const_iterator b,
            const_iterator e
            );

        /// Finds the actual passed in design without regard to the predicate.
        /**
         * If the passed in pointer is not kept in here, the return is end.
         * Otherwise it is an iterator to the location of the Design.
         *
         * The regular find method will find an equivalent Design which may
         * or may not be the actual passed in design.
         *
         * \param key The Design for which an exact match is sought.
         * \return An immutable iterator to the location of "key" or end()
         *         if not found.
         */
        inline
        DesignDVSortSet::const_iterator
        find_exact(
            const key_type& key
            ) const;

        /**
         * \brief Finds a logical equivalent to the passed in Design but
         *        disregards the Design itself.
         *
         * The pointer value of the returned iterator (if not end()) will not
         * be equal to key. If no equivalent matches other than an exact match
         * of key are present, the return is end().  This is just like the
         * regular find method except that a return of the exact same object as
         * is pointed to by "key" is not allowed.
         *
         * The regular find method will find an equivalent Design which may
         * or may not be the actual passed in design and find_exact will find
         * the actual design.
         *
         * \param key The Design for which a duplicate is sought.
         * \return An immutable iterator to the location of a duplicate or
         *         end() if not found.
         */
        inline
        DesignDVSortSet::const_iterator
        find_not_exact(
            const key_type& key
            ) const;

        /// Finds the actual passed in design without regard to the predicate.
        /**
         * If the passed in pointer is not kept in here, the return is end.
         * Otherwise it is an iterator to the location of the Design.
         *
         * The regular find method will find an equivalent Design which may
         * or may not be the actual passed in design and the find_not_exact
         * method will find one that is definitely not the same object as key.
         *
         * \param key The Design for which an exact match is sought.
         * \return An iterator to the location of "key" or end()
         *         if not found.
         */
        inline
        DesignDVSortSet::iterator
        find_exact(
            const key_type& key
            );

        /// Finds the actual passed in design without regard to the predicate.
        /**
         * The pointer value of the returned iterator (if not end()) will not
         * be equal to key. If no equivalent matches other than an exact match
         * of key are present, the return is end().  This is just like the
         * regular find method except that a return of the exact same object as
         * is pointed to by "key" is not allowed.
         *
         * The regular find method will find an equivalent Design which may
         * or may not be the actual passed in design and find_exact will find
         * the actual design.
         *
         * \param key The Design for which a duplicate is sought.
         * \return An iterator to the location of a duplicate or end()
         *         if not found.
         */
        inline
        DesignDVSortSet::iterator
        find_not_exact(
            const key_type& key
            );

        /**
         * \brief Counts the number of duplicate Design configurations in this
         *        container.
         *
         * For each set of duplicate Design configurations, 1 is considered
         * unique and the rest are considered non-unique.  No tagging of clones
         * takes place. Comparisons are made according to the predicate.
         *
         * \return The number of non-unique Designs found.
         */
        inline
        DesignDVSortSet::size_type
        count_non_unique(
            ) const;

        /**
         * \brief Prints each Design class object in this container to "stream"
         * using the overload.
         *
         * \param stream The stream into which to write the designs of this
         *               container.
         * \return The supplied \a stream is returned for convenience.
         *         It allows for chaining.
         */
        inline
        std::ostream&
        stream_out(
            std::ostream& stream
            ) const;

        /// Writes the Design "des" to "stream".
        /**
         * The design is written in tab-delimited flat file format.  No matter
         * what, all design variables are written.  After that, objective and
         * constraint values are written iff the design has been evaluated and
         * is not ill-conditioned.
         *
         * \param des The design to write into the supplied "stream".
         * \param stream The stream into which to write the designs of this
         *               container.
         * \return The supplied \a stream is returned for convenience.
         *         It allows for chaining.
         */
        inline static
        std::ostream&
        stream_out(
            const key_type des,
            std::ostream& stream
            );

        /**
         * \brief Deletes each Design pointer stored in this container then
         *        clears the set.
         */
        void
        flush(
            );

        /// Erases all occurrences of the exact key specified.
        /**
         * This does not mean logically equivalent according to the predicate.
         * This erases all occurrences of the design pointer specified.  Clones
         * of key may remain in the set.
         *
         * \param key The Design for which instances should be removed.
         * \return The number of occurrences found and marked.
         */
        DesignDVSortSet::size_type
        erase_exacts(
            const key_type key
            );

        /// Copies all elements of "other" into this.  Leaves "other" in tact.
        /**
         * \param other The container from which to copy Designs into this.
         */
        template <typename DesignContainer>
        void
        copy_in(
            const DesignContainer& other
            )
        {
            EDDY_FUNC_DEBUGSCOPE

            typedef typename DesignContainer::const_iterator it_t;

            // look at each member of other in turn.
            const it_t oe(other.end());

            for(it_t it(other.begin()); it!=oe; ++it)
            {
                this->_data.insert(*it);
                if(this->_doCache) this->_lruList.add(*it);
            }
        }

        /// Looks for a duplicate of "key" existing within this container.
        /**
         * Comparisons are made according to the predicate.  If found, this
         * method tags the Designs as clones.
         *
         * "key" should not be a member of this set.
         *
         * \param key The Design for which a clone is sought.
         * \return An immutable iterator to the found clone or end() if not
         *         one.
         */
        DesignDVSortSet::const_iterator
        test_for_clone(
            const key_type key
            ) const;

        ///// Detects and sets dead each non-unique Design configuration.
        ///**
        // * For each set of duplicate Design configurations, 1 is considered
        // * unique and the rest are considered non-unique.  No tagging of clones
        // * takes place. Comparisons are made according to the predicate.
        // *
        // * \return The number of Designs marked.
        // */
        //DesignDVSortSet::size_type
        //mark_non_unique(
        //    std::size_t mark = MARK
        //    ) const;

        ///**
        // * \brief Goes through this container and marks all designs for which
        // *        the predicate returns true.
        // *
        // * \param pred The function object to call on each design as it is
        // *             considered.
        // * \return The number of Designs marked.
        // */
        //template <typename MarkPred>
        //typename DesignDVSortSet::size_type
        //mark_if(
        //    MarkPred pred,
        //    std::size_t mark = MARK
        //    )
        //{
        //    EDDY_FUNC_DEBUGSCOPE

        //    // store the initial size so we can return the removal count.
        //    size_type nmarked = 0;
        //    const const_iterator te(this->end());

        //    // iterate the set and remove all designs for which func evaluates
        //    // to evaluatesTo.
        //    for(iterator it(this->begin()); it!=te; ++it)
        //    {
        //        bool marked = pred(*it);
        //        (*it)->Design::ModifyAttribute(mark, marked);
        //        if(marked) ++nmarked;
        //    }

        //    // return the difference between the current size and the old size.
        //    return nmarked;

        //} // DesignMultiSet::mark_if

        /// Checks all the members of "other" for duplicates in this list.
        /**
         * Duplicates found are tagged as clones.
         *
         * Returns the number of times duplicates were detected which is not
         * the same as the number of newly found clones.  It is possible that
         * some clones detected were already clones of other designs and
         * therefore were not "newly" found.
         *
         * Comparisons are made according to the predicate.
         *
         * The design container must be a forward iteratable STL compliant
         * container of Design* 's with begin, end, and empty methods and
         * define the const_iterator type.
         *
         * \param other The other container from which to seek duplicates in
         *              this.
         * \return The number of times duplicates were detected.
         */
        template <typename DesignContainer>
        size_type
        test_for_clones(
            const DesignContainer& other
            ) const
        {
            EDDY_FUNC_DEBUGSCOPE

            // check for the trivial abort conditions
            if(other.empty() || this->empty()) return 0;

            // prepare to store the number of successful clone tests.
            size_type clonecount = 0;
            const typename DesignContainer::const_iterator oe(other.end());
            const const_iterator te(this->end());

            // look at each member of other in turn.
            for(typename DesignContainer::const_iterator it(other.begin());
                it!=oe; ++it)
                    clonecount += (this->test_for_clone(*it)==te) ? 0 : 1;

            // return the number of successful clone tests.
            return clonecount;

        } // test_for_clones

        /**
         * \brief A specialization of the test_for_clones method for containers
         *        of this type.
         *
         * The knowledge that the sorts are the same can be used to expedite
         * this search.  This must be re-implemented in any derived classes and
         * called back on in order for the specialization to work properly.
         *
         * \param other The other container from which to seek duplicates in
         *              this.
         * \return The number of times duplicates were detected.
         */
        DesignDVSortSet::size_type
        test_for_clones(
            const DesignDVSortSet& other
            ) const
        {
            EDDY_FUNC_DEBUGSCOPE

            // check for the trivial abort conditions
            if(other.empty() || this->empty()) return 0;

            // Check to be sure that the passed in deque is not this.
            // If it is, this method is not appropriate.
            if(&(this->_data) == &other)
                return this->test_within_list_for_clones();

            // remember that the complexity of searches in this container is
            // bounded as O(log10(size)).  So it is probably safe to assume
            // that searching in the larger of the lists is nearly as cheap
            // as searching in the smaller.

            // We are going to iterate the smaller list and look at each member
            // but we are going to start with the first possible match and
            // finish with the last.  We can determine those matches using the
            // lower_bound and upper_bound methods.

            // figure out which is the smaller and which is the larger.
            const DesignDVSortSet& smaller =
                (this->size() <= other.size()) ? this->_data : other;
            const DesignDVSortSet& larger =
                (&smaller == &this->_data) ? other : this->_data;

            // The lower_bound of the first element of the larger in the
            // smaller will be where we start looking and the upper_bound of
            // the last of larger in the smaller will be where we stop.
            const_iterator s(smaller.lower_bound(*(larger.begin())));
            const const_iterator e(smaller.upper_bound(*(larger.rbegin())));
            const const_iterator le(larger.end());

            // now do our thing with these elements.
            // prepare to store the number of successful clone tests.
            size_type clonecount = 0;

            // look at each member of other in turn.
            for(; s!=e; ++s)
            {
                const_iterator it(larger.test_for_clone(*s));
                if(it != le)
                {
                    this->on_accessed(it);
                    ++clonecount;
                }
            }

            // return the number of successful clone tests.
            return clonecount;

        } // test_for_clones

        /// Performs clone testing within this list.
        /**
         * This method should be used instead of comparing this container
         * to itself using test_for_clones.
         *
         * Found clones are tagged as such using Design::TagAsClones.
         * The return is the number of newly found clones (those that were
         * not previously marked as a clone of any other Design).
         *
         * Comparisons are made according to the predicate.
         *
         * \return The number of newly found clones.
         */
        inline
        DesignDVSortSet::size_type
        test_within_list_for_clones(
            ) const;

        ///// Marks all occurrences of logical equivalents to the key specified.
        ///**
        // * This includes the key itself.
        // *
        // * \param key The Design for which instances should be removed.
        // * \return The number of occurrences found and marked.
        // */
        //DesignDVSortSet::size_type
        //mark(
        //    const key_type key,
        //    std::size_t mark = MARK
        //    );

        ///**
        // * \brief Marks all occurrences of logical equivalents to the key
        // *        with the exception of the key itself.
        // *
        // * \param key The Design for which instances should be removed.
        // * \return The number of occurrences found and marked.
        // */
        //DesignDVSortSet::size_type
        //mark_not_exact(
        //    const key_type key,
        //    std::size_t mark = MARK
        //    );







































        //const key_type&
        //operator [](
        //    const key_type& key
        //    ) const;

        //void
        //insert(
        //    const key_type& key
        //    );

        //const key_type&
        //retrieve(
        //    const key_type& key
        //    ) const;

    /*
    ============================================================================
    Subclass Visible Methods
    ============================================================================
    */
    protected:

        void
        manage_size(
            );



    /*
    ============================================================================
    Subclass Overridable Methods
    ============================================================================
    */
    public:


    protected:


    private:


    /*
    ============================================================================
    Private Methods
    ============================================================================
    */
    private:

        template <typename ItT>
        const ItT&
        on_accessed(
            const ItT& it
            ) const;

        template <typename ItT>
        const ItT&
        on_maybe_accessed(
            const ItT& it
            ) const;



    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        /// Constructs an LRUDesignCache limited to the supplied size.
        /**
         * \param pred The maximum allowable size for the cache.
         */
        LRUDesignCache(
            const size_type& maxSize
            );

}; // class LRUDesignCache



/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA







/*
================================================================================
Include Inlined Functions File
================================================================================
*/
#include "inline/LRUDesignCache.hpp.inl"



/*
================================================================================
End of Multiple Inclusion Check
================================================================================
*/
#endif // JEGA_UTILITIES_LRUDESIGNCACHE_HPP
