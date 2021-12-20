/*
================================================================================
    PROJECT:

        Eddy C++ Utilities Project

    CONTENTS:

        Inline methods of class Math.

    NOTES:

        See notes of Math.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Jun 03 10:30:43 2003 - Original Version (JE)

================================================================================
*/



/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the inline methods of the Math class.
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
namespace eddy {
    namespace utilities {







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
double
Math::Pi(
    )
{
    return 3.1415926535897932384;

} // Math::Pi

template <typename T>
inline
bool
Math::IsWhole(
    T val
    )
{
	return true;
}

inline
bool
Math::IsWhole(
    double val
    )
{
    return Round(val) == val;
}

inline
bool
Math::IsWhole(
    float val
    )
{
    return Round(val) == val;
}

template<typename T>
inline
const T&
Math::Max(
    const T& arg1,
    const T& arg2
    )
{
    return (arg1 < arg2) ? arg2 : arg1;

} // Math::Max

template<typename T>
inline
const T&
Math::Min(
    const T& arg1,
    const T& arg2
    )
{
    return (arg1 < arg2) ? arg1 : arg2;

} // Math::Min


template <typename T>
T
Math::Round(
    T val,
    int prec
    )
{
	return static_cast<T>(Round(static_cast<double>(val), prec));
}

inline
unsigned int
Math::Abs(
    unsigned int val
    )
{
    return val;
}

inline
unsigned long
Math::Abs(
    unsigned long val
    )
{
    return val;
}

inline
unsigned short
Math::Abs(
    unsigned short val
    )
{
    return val;
}

inline
unsigned long long
Math::Abs(
    unsigned long long val
    )
{
    return val;
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
    } // namespace utilities
} // namespace eddy
