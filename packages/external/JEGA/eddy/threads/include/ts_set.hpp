/*
================================================================================
    PROJECT:

        Eddy C++ Thread Safety Project

    CONTENTS:

        Definition of class ts_set.

    NOTES:

        See notes under Class Definition section of this file.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Thu Apr 08 16:07:28 2004 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of the ts_set class.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef EDDY_THREADS_TS_SET_HPP
#define EDDY_THREADS_TS_SET_HPP







/*
================================================================================
Includes
================================================================================
*/
// config.hpp should be the first include.
#include "../include/config.hpp"

#include <set>
#include "mutex_lock.hpp"
#include "scoped_lock_concept.hpp"







/*
================================================================================
Pre-Namespace Forward Declares
================================================================================
*/








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
namespace eddy {
    namespace threads {





/*
================================================================================
In-Namespace Forward Declares
================================================================================
*/
template <
            typename Key,
            typename Compare,
            typename Allocator
          >
class ts_set;







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
/// A set class implementation with mutex protection.
/**
 * All methods lock for their duration.  To lock over multiple method
 * calls, use the lock and unlock methods or the scoped_lock.
 *
 * A note on iterating this container.  To do so, you should lock and unlock
 * or use the scoped_lock.  There is no iterator level mutex protection.
 * The begin and end methods will lock for their duration but to be sure
 * that it is not modified by another thread for the entire time you are
 * iterating it, lock, and when you are done, unlock.
 *
 * Certain operations, such as assignment, are defined for both a
 * right-hand-side of this type as well as for one of type
 * std::set<Key, Compare, Allocator>
 *
 * \param Key The type of the items stored in this set.
 * \param Compare the comparator used to order the items in this set.
 * \param Allocator the object used to allocate space for new items in this
 *        set.
 */
template <
            typename Key,
            typename Compare = std::less<Key>,
            typename Allocator = std::allocator<Key>
          >
class EDDY_SL_IEDECL ts_set
{
    /*
    ============================================================================
    Typedefs
    ============================================================================
    */
    public:

        /// A shorthand for this type.
        typedef
        ts_set<Key, Compare, Allocator>
        my_type;

        /// The type of the scoped lock that can be used with this container.
        typedef
        scoped_lock_concept<my_type>
        scoped_lock;

        /// A shorthand for the underlying set type.
        typedef
        std::set<Key, Compare, Allocator>
        container_type;

        /// An type that represents the allocator class for this set.
        typedef
        typename container_type::allocator_type
        allocator_type;

        /// An iterator type for this set.
        typedef
        typename container_type::iterator
        iterator;

        /// A constant iterator type for this set.
        typedef
        typename container_type::const_iterator
        const_iterator;

        /// A type for keeping the size of this set.
        typedef
        typename container_type::size_type
        size_type;

        /// The type of the allocator this set is using.
        typedef
        typename container_type::allocator_type
        allocator_type;

        /**
         * \brief A type to represent the difference between the address of two
         *        elements.
         */
        typedef
        typename container_type::difference_type
        difference_type;

        /// A type to serve as a non-const pointer to data.
        typedef
        typename container_type::pointer
        pointer;

        /// A type to serve as a constant pointer to an element in the set.
        typedef
        typename container_type::const_pointer
        const_pointer;

        /// A type to serve as a reference to an element in the set.
        typedef
        typename container_type::reference
        reference;

        /// A type to serve as a constant-iterator for the set.
        typedef
        typename container_type::const_reference
        const_reference;

        /// A type to serve as a non-constant reverse iterator for the set.
        typedef
        typename container_type::reverse_iterator
        reverse_iterator;

        /// A type to serve as a constant reverse iterator for the set.
        typedef
        typename container_type::const_reverse_iterator
        const_reverse_iterator;

        /// A type to serve as a key and value type for this set.
        typedef
        typename container_type::key_type
        key_type;

        /**
         * \brief A type to serve as a comparison function for values of type
         *        key_type.
         */
        typedef
        typename container_type::key_compare
        key_compare;

        /// A type to serve as a value stored in this set.
        typedef
        typename container_type::value_type
        value_type;

        /**
         * \brief A type to serve as a comparison function for values of type
         *        value_type.
         */
        typedef
        typename container_type::value_compare
        value_compare;

        /// A pair of const_iterators.
        /**
         * This is the return type for the const version of the equal_range
         * method.
         */
        typedef
        std::pair<const_iterator, const_iterator>
        const_iterator_pair;

        /// A pair of iterators.
        /**
         * This is the return type for the non-const version of the equal_range
         * method.
         */
        typedef
        std::pair<iterator, iterator>
        iterator_pair;

    protected:

    private:



    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:

        /// The underlying un-protected set for this set.
        container_type _container;

        /// The mutex used to protect this set.
        mutable mutex _mutex;



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

        /// Locks this set until an unlock is issued.
        /**
         * All methods lock the set for their duration.
         * Use this method if you need to lock the set over multiple actions.
         *
         * The set may be locked recursively.
         */
        inline
        void
        lock(
            ) const;

        /// Unlocks this set.
        /**
         * Only use this method to counter a lock call.
         */
        inline
        void
        unlock(
            ) const;

        /**
         * \brief Returns a const_iterator to the first element in this set
         *        or end() if empty().
         *
         * \return Iterator to the first element in the container.
         */
        inline
        const_iterator
        begin(
            ) const;

        /**
         * \brief Returns an iterator to the first element in this set
         *        or end() if empty().
         *
         * \return Iterator to the first element in the container.
         */
        inline
        iterator
        begin(
            );

        /// Erases all elements of this set.
        inline
        void
        clear(
            );

        /// Counts the number of elements whose key matches \a key
        /**
         * Given that this is not a multiset implementation, the return
         * should always be 1.
         *
         * \param key The key of interest.
         * \return The count of items in this set that match \a key.
         */
        inline
        size_type
        count(
            const key_type& key
            ) const;

        /// Returns true if the size of this set is 0 and false otherwise.
        /**
         * \return true if the set has no elements and false otherwise.
         */
        bool
        empty(
            ) const;

        /// Returns a const_iterator to one past the last element in this set.
        /**
         * \return Iterator to 1 past the last element in this set.
         */
        inline
        const_iterator
        end(
            ) const;

        /// Returns an iterator to one past the last element in this set.
        /**
         * \return Iterator to 1 past the last element in this set.
         */
        inline
        iterator
        end(
            );

        /**
         * \brief Returns a pair of iterators where first is
         *        lower_bound(\a key) and second is upper_bound(\a key)
         *
         * \param key The key to find the upper and lower bounds for.
         * \return A pair of iterators comprised of the lower and upper bounds
         *         of \a key in this set.
         */
        inline
        const_iterator_pair
        equal_range(
            const key_type& key
            ) const;

        /**
         * \brief Returns a pair of iterators where first is
         *        lower_bound(\a key) and second is upper_bound(\a key)
         *
         * \param key The key to find the upper and lower bounds for.
         * \return A pair of iterators comprised of the lower and upper bounds
         *         of \a key in this set.
         */
        inline
        iterator_pair
        equal_range(
            const key_type& key
            );

        /// Erases the element pointed to by \a where from the set.
        /**
         * \param where The location of the element to remove from this set.
         * \return An iterator to the next element after the one removed
         *         or end() if such an element does not exist.
         */
        inline
        iterator
        erase(
            iterator where
            );

        /// Erases the elements in the range [first, last).
        /**
         *
         * \param first The location of the first element to remove from this
         *              set.
         * \param last 1 past the last location of an element to remove from
         *             this set.
         * \return An iterator to the next element after those removed or
         *         end() if such an element does not exist.
         */
        inline
        iterator
        erase(
            iterator first,
            iterator last
            );

        /// Erases all elements for that compare equal to \a key.
        /**
         * \param key The key for which all matching elements should be
         *            removed.
         * \return The number of elements removed.
         */
        inline
        size_type
        erase(
            const key_type& key
            );

        /**
         * \brief Returns an iterator to the first element found whose key
         *        matches \a key or end() if no such element can be found.
         *
         * \param key The key sought in this set.
         * \return An iterator to the first occurrence of a match to \a key
         *         or end if not found.
         */
        inline
        iterator
        find(
            const key_type& key
            );

        /**
         * \brief Returns a const_iterator to the first element found whose key
         *        matches \a key or end() if no such element can be found.
         *
         * \param key The key sought in this set.
         * \return An iterator to the first occurrence of a match to \a key
         *         or end if not found.
         */
        inline
        const_iterator
        find(
            const key_type& key
            ) const;

        /// Returns a duplicate of the allocator being used by this set.
        /**
         * \return A copy of the allocator being used by this set.
         */
        inline
        allocator_type
        get_allocator(
            ) const;

        /**
         * \brief Inserts an element with key \a key into the proper location
         *        in this set.
         *
         * \param value The item to insert into the set.
         * \return An iterator to the inserted value.
         */
        inline
        std::pair <iterator, bool>
        insert(
            const value_type& value
            );

        /**
         * \brief Inserts the supplied value after the location indicated by
         *        \a where.
         *
         * \param where The iterator to the location after which the new items
         *              location should be found.
         * \param value The item to insert into the set.
         * \return An iterator to the inserted value.
         */
        inline
        iterator
        insert(
            iterator where,
            const value_type& value
            );

        /// Inserts the values in the supplied range to this set.
        /**
         * \a last must be reachable by forward iteration from \a first.
         *
         * \param first The iterator to the first element in the range to be
         *              inserted into this set.
         * \param last The iterator to 1 past the last element in the range to
         *             be inserted into this set.
         */
        template<typename InputIterator>
        inline
        void
        insert(
            InputIterator first,
            InputIterator last
            )
        {
            mutex_lock lock(this->_mutex);
            this->_container.insert(first, last);
        }

        /**
         * \brief Returns a duplicate of the comparison object used to order
         *        this set.
         *
         * \return The comparison predicate used by this set (by value).
         */
        inline
        key_compare
        key_comp(
            ) const;

        /**
         * \brief Returns an iterator to the first element in this set with a
         *        key that is equal to or greater than a specified \a key.
         *
         * \param key The value for which the lower bound is sought.
         * \return The location in this set that contains the first element
         *         with a key that is equal to or greater than the supplied
         *         key.
         */
        inline
        const_iterator
        lower_bound(
            const key_type& key
            ) const;

        /**
         * \brief Returns an iterator to the first element in this set with a
         *        key that is equal to or greater than the specified \a key.
         *
         * \param key The value for which the lower bound is sought.
         * \return The location in this set that contains the first element
         *         with a key that is equal to or greater than the supplied
         *         key.
         */
        inline
        iterator
        lower_bound(
            const key_type& key
            );

        /// Returns the maximum possible size of this set.
        /**
         * \return The maximum allowable size of this set.
         */
        inline
        size_type
        max_size(
            ) const;

        /**
         * \brief Returns a const_reverse_iterator to the first element in this
         *        set reversed or rend() if empty().
         *
         * \return Iterator to the first element in the reverse sequence of
         *         this container.
         */
        inline
        const_reverse_iterator
        rbegin(
            ) const;

        /**
         * \brief Returns a reverse_iterator to the first element in this set
         *        reversed or rend() if empty().
         *
         * \return Iterator to the first element in the reverse sequence of
         *         this container.
         */
        inline
        reverse_iterator
        rbegin(
            );

        /**
         * \brief Returns a const_reverse_iterator to one past the last element
         *        in this set reversed.
         *
         * \return Iterator to 1 past the last element in the reverse sequence
         *         of this container.
         */
        inline
        const_reverse_iterator
        rend(
            ) const;

        /**
         * \brief Returns a reverse_iterator to one past the last element in
         *        this set reversed.
         *
         * \return Iterator to 1 past the last element in the reverse sequence
         *         of this container.
         */
        inline
        reverse_iterator
        rend(
            );

        /// Returns the current number of elements in this set.
        /**
         * \return The number of elements currently in this set.
         */
        inline
        size_type
        size(
            ) const;

        /// Exchanges the elements of this with those of right.
        /**
         * \param right The existing container type with which to swap this.
         */
        inline
        void
        swap(
            container_type& right
            );

        /// Exchanges the elements of this with those of right.
        /**
         * \param right The existing set with which to swap this.
         */
        inline
        void
        swap(
            my_type& right
            );

        /**
         * \brief Returns an iterator to the last element in this set with a
         *        key that is equal to the specified \a key or the first one
         *        with a key greater than the specified \a key if no equals are
         *        found.
         *
         * \param key The value for which the upper bound is sought.
         * \return The location in this set that contains the last element
         *         with a key that is equal to the specified \a key or the
         *         first one with a key greater than the specified \a key if no
         *         equals are found.
         */
        inline
        const_iterator
        upper_bound(
            const key_type& key
            ) const;

        /**
         * \brief Returns an iterator to the last element in this set with a
         *        key that is equal to the specified \a key or the first one
         *        with a key greater than the specified \a key if no equals are
         *        found.
         *
         * \param key The value for which the upper bound is sought.
         * \return The location in this set that contains the last element
         *         with a key that is equal to the specified \a key or the
         *         first one with a key greater than the specified \a key if no
         *         equals are found.
         */
        inline
        iterator
        upper_bound(
            const key_type& key
            );

        /// Returns a copy of the object used to order elements in this set.
        /**
         * \return The comparator used to order the values of this set.
         */
        inline
        value_compare
        value_comp(
            ) const;

        /// Assigns the contents of this set to those of \a rhs
        /**
         * \param rhs An existing set from which to copy elements into this.
         * \return This set after assignment.
         */
        inline
        const my_type&
        operator =(
            const my_type& rhs
            );

        /// Assigns the contents of this set to those of \a rhs
        /**
         * \param rhs An existing container from which to copy elements into
         *            this.
         * \return This set after assignment.
         */
        inline
        const my_type&
        operator =(
            const container_type& rhs
            );

    /*
    ============================================================================
    Subclass Visible Methods
    ============================================================================
    */
    protected:





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





    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:


        /// Default constructs a ts_set which is initially empty.
        inline
        ts_set(
            );

        /// Constructs a ts_set which used \a comp to order elements.
        explicit
        inline
        ts_set(
            const key_compare& comp
            );

        /**
         * \brief Constructs a ts_set which used \a comp to order elements
         *        and \a al for allocation.
         */
        explicit
        inline
        ts_set(
            const key_compare& comp,
            const allocator_type& al
            );

        /// Copy constructs a ts_set from another.
        inline
        ts_set(
            const my_type& copy
            );

        /// Copy constructs a ts_set from a non-thread safe set.
        inline
        ts_set(
            const container_type& copy
            );

        /**
         * \brief Constructs a ts_set composed initially of the elements
         *        in the range [first, last)
         */
        template<typename InputIterator>
        inline
        ts_set(
            InputIterator first,
            InputIterator last
            ) :
                _container(first, last),
                _mutex(PTHREAD_MUTEX_RECURSIVE)
        {
        }

        /**
         * \brief Constructs a ts_set composed initially of the elements
         *        in the range [first, last) which uses \a comp to order
         *        elements.
         */
        template<typename InputIterator>
        inline
        ts_set(
            InputIterator first,
            InputIterator last,
            const key_compare& comp
            ) :
                _container(first, last, comp),
                _mutex(PTHREAD_MUTEX_RECURSIVE)
        {
        }

        /**
         * \brief Constructs a ts_set composed initially of the elements
         *        in the range [first, last) which uses \a comp to order
         *        elements and \a al to do allocation.
         */
        template<typename InputIterator>
        inline
        ts_set(
            InputIterator first,
            InputIterator last,
            const key_compare& comp,
            const allocator_type& al
            ) :
                _container(first, last, comp, al),
                _mutex(PTHREAD_MUTEX_RECURSIVE)
        {
        }


}; // class ts_set



/*
================================================================================
End Namespace
================================================================================
*/
    } //  namespace threads
} // namespace eddy







/*
================================================================================
Include Inlined Functions File
================================================================================
*/
#include "./inline/ts_set.hpp.inl"



/*
================================================================================
End of Multiple Inclusion Check
================================================================================
*/
#endif // EDDY_THREADS_TS_SET_HPP
