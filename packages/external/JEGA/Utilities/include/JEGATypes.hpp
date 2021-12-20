/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        All data types used in the JEGA project.

    NOTES:



    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Wed Dec 14 09:32:21 2005 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the data types used in the JEGA project.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_JEGATYPES_HPP
#define JEGA_JEGATYPES_HPP







/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <set>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <utility>
#include <cstddef>
#include <utilities/include/int_types.hpp>






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
namespace JEGA {





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
/// Pair of strings for general use.
typedef
std::pair<std::string, std::string>
StringPair;

/// A vector of StringPairs for general use.
typedef
std::map<std::string, std::string>
StringMap;

/// A set of strings for a general purpose sorted, unique string collection.
typedef
std::set<std::string>
StringSet;

/// A vector of strings for general use.
typedef
std::vector<std::string>
StringVector;

/// A linked list of strings for general use.
typedef
std::list<std::string>
StringList;

/// A vector of doubles for general use
typedef
std::vector<double>
DoubleVector;

/// A linked list of doubles for general use.
typedef
std::list<double>
DoubleList;

/// A vector of size_t's for general use
typedef
std::vector<std::size_t>
SizeTVector;

/// A vector of ints for general use
typedef
std::vector<int>
IntVector;

/// A vector of 8 bit ints for general use
typedef
std::vector<eddy::utilities::int8_t>
Int8Vector;

/// A vector of 16 bit ints for general use
typedef
std::vector<eddy::utilities::int16_t>
Int16Vector;

/// A vector of 32 bit ints for general use
typedef
std::vector<eddy::utilities::int32_t>
Int32Vector;

/// A vector of 64 bit ints for general use
typedef
std::vector<eddy::utilities::int64_t>
Int64Vector;

/// A linked list of ints for general use
typedef
std::list<int>
IntList;

/// A matrix of doubles for general use
typedef
std::vector<DoubleVector>
DoubleMatrix;

/// A matrix of ints for general use
typedef
std::vector<IntVector>
IntMatrix;


#ifndef X_TYPE
#	define X_TYPE double
#endif

#ifndef F_TYPE
#	define F_TYPE double
#endif

#ifndef G_TYPE
#	define G_TYPE double
#endif

typedef X_TYPE var_rep_t;
typedef F_TYPE obj_val_t;
typedef G_TYPE con_val_t;


template<
    typename KeyT,
    typename ValT,
	typename PredT = std::less<KeyT>,
	typename AllocT = std::allocator<std::pair<const KeyT, ValT> >
    >
class icmap : public std::map<KeyT, ValT, PredT, AllocT>
{
    private:

        typedef icmap<KeyT, ValT, PredT, AllocT> my_type;
        typedef std::map<KeyT, ValT, PredT, AllocT> base_type;

    public:
	    icmap() : base_type() {}
        icmap(const typename icmap::size_type&) : base_type() {};
};

template<
    typename ValT,
	typename PredT = std::less<ValT>,
	typename AllocT = std::allocator<ValT>
    >
class icset : public std::set<ValT, PredT, AllocT>
{
    private:

        typedef icset<ValT, PredT, AllocT> my_type;
        typedef std::set<ValT, PredT, AllocT> base_type;
            
    public:
	    icset() : base_type() {}
        icset(const typename icset::size_type&) : base_type() {};
};

/*
================================================================================
End Namespace
================================================================================
*/
} // namespace JEGA







/*
================================================================================
Include Inlined Functions File
================================================================================
*/
// Not using an Inlined Functions File.



/*
================================================================================
End of Multiple Inclusion Check
================================================================================
*/
#endif // JEGA_JEGATYPES_HPP
