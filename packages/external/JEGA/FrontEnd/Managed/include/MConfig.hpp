/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        General configuration for the managed project.

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

        Mon Sep 25 10:51:11 2006 - Original Version (JE)

================================================================================
*/
#pragma once


/// An empty argument for macros.
#define EMPTY_ARG()

#pragma unmanaged
#include <../Utilities/include/Logging.hpp>
#include <cstddef>
#include <utilities/include/int_types.hpp>
#pragma managed


/*
================================================================================
Begin Namespace
================================================================================
*/
namespace JEGA {
    namespace FrontEnd {
        namespace Managed {



/*
================================================================================
In-Namespace Forward Declares
================================================================================
*/
// Forward declare the MDesign for creation of the DesignVector below.
ref class MDesign;

// Do the same for the MSolution
ref class MSolution;

#define VECTOR_CLASS(dtype, name) \
    public ref class name : \
        public System::Collections::Generic::List< dtype > \
    { \
        public: \
            name() : System::Collections::Generic::List< dtype >() {} \
            name(System::Int32 initCap) : \
                System::Collections::Generic::List< dtype >(initCap) {} \
            name(System::Collections::Generic::IEnumerable< dtype >^ copy) : \
                System::Collections::Generic::List< dtype >(copy) {} \
    };

VECTOR_CLASS(System::Double, DoubleVector)
VECTOR_CLASS(DoubleVector^, DoubleMatrix)
VECTOR_CLASS(System::Int32, IntVector)
VECTOR_CLASS(System::String^, StringVector)
VECTOR_CLASS(MDesign^, DesignVector)
VECTOR_CLASS(MSolution^, SolutionVector)


/*
================================================================================
End Namespace
================================================================================
*/
        } // namespace Managed
    } // namespace FrontEnd
} // namespace JEGA
