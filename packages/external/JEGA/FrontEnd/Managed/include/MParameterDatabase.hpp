/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Definition of class MParameterDatabase.

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

        Wed Feb 08 16:29:50 2006 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of the MParameterDatabase class.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_FRONTEND_MANAGED_MPARAMETERDATABASE_HPP
#define JEGA_FRONTEND_MANAGED_MPARAMETERDATABASE_HPP

#pragma once


/*
================================================================================
Includes
================================================================================
*/
#include <MConfig.hpp>







/*
================================================================================
Pre-Namespace Forward Declares
================================================================================
*/
namespace JEGA
{
    namespace Utilities
    {
        class ParameterDatabase;
    }
}






/*
================================================================================
Namespace Aliases
================================================================================
*/
#pragma managed
#using <mscorlib.dll>








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
 * \brief A managed version of the ParameterDatabase interface for use in the
 *        managed front end.
 *
 * This interface only requires the "Add" methods but also requires that a
 * manifestation to a JEGA::Utilities::ParameterDatabase be possible for
 * extraction of parameters from within JEGA core.
 */
public interface class MParameterDatabase :
    public System::IDisposable
{
    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        /**
         * \brief Override to manifest this managed parameter database as a
         *        JEGA::Utilities::ParameterDatabase object.
         *
         * This is necessary so that the data pushed into this keep can be
         * accessed by JEGA core. This is typically simple because typically,
         * managed databases will internally use
         * JEGA::Utilities::ParameterDatabase objects.
         *
         * \return The JEGA::Utilities::ParameterDatabase representation of
         *         this managed database.
         */
        virtual
        JEGA::Utilities::ParameterDatabase&
        Manifest(
            ) = 0;

        /// Tests the integer param map for a value mapped to the supplied tag.
        /**
         * \param tag The key to search for in the integer param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasIntegralParam(
            System::String^ tag
            ) = 0;

        /// Tests the short param map for a value mapped to the supplied tag.
        /**
         * \param tag The key to search for in the short param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasShortParam(
            System::String^ tag
            ) = 0;

        /// Tests the double param map for a value mapped to the supplied tag.
        /**
         * \param tag The key to search for in the double param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasDoubleParam(
            System::String^ tag
            ) = 0;

        /// Tests the size_t param map for a value mapped to the supplied tag.
        /**
         * \param tag The key to search for in the size_t param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasSizeTypeParam(
            System::String^ tag
            ) = 0;

        /// Tests the Boolean param map for a value mapped to the supplied tag.
        /**
         * \param tag The key to search for in the Boolean param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasBooleanParam(
            System::String^ tag
            ) = 0;

        /// Tests the string param map for a value mapped to the supplied tag.
        /**
         * \param tag The key to search for in the string param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasStringParam(
            System::String^ tag
            ) = 0;

        /**
         * \brief Tests the vector of doubles param map for a value mapped to
         *        the supplied tag.
         *
         * \param tag The key to search for in the vector of doubles param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasDoubleVectorParam(
            System::String^ tag
            ) = 0;

        /**
         * \brief Tests the vector of ints param map for a value mapped to the
         *        supplied tag.
         *
         * \param tag The key to search for in the vector of ints param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasIntVectorParam(
            System::String^ tag
            ) = 0;

        /**
         * \brief Tests the matrix of doubles param map for a value mapped to
         *        the supplied tag.
         *
         * \param tag The key to search for in the matrix of doubles param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasDoubleMatrixParam(
            System::String^ tag
            ) = 0;

        /**
         * \brief Tests the vector of strings param map for a value mapped to
         *        the supplied tag.
         *
         * \param tag The key to search for in the vector of strings param map.
         * \return true if \a tag is found as a key and false otherwise.
         */
        virtual
        bool
        HasStringVectorParam(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as an integer from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The int value associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with an
         *                                 int value.
         */
        virtual
        int
        GetIntegral(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a short from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The short value associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with an
         *                                 short value.
         */
        virtual
        short
        GetShort(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a double from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The double value associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with a
         *                                 double value.
         */
        virtual
        double
        GetDouble(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a size_t from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The size_t value associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with a
         *                                 size_t value.
         */
        virtual
        eddy::utilities::uint64_t
        GetSizeType(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a bool from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The bool value associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with a
         *                                 bool value.
         */
        virtual
        bool
        GetBoolean(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a string from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The std::string value associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with a
         *                                 string value.
         */
        virtual
        System::String^
        GetString(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a DoubleVector from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The vector of doubles associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with a
         *                                 DoubleVector value.
         */
        virtual
        JEGA::FrontEnd::Managed::DoubleVector^
        GetDoubleVector(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as an IntVector from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The vector of ints associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with an
         *                                 IntVector value.
         */
        virtual
        JEGA::FrontEnd::Managed::IntVector^
        GetIntVector(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a DoubleMatrix from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The matrix of doubles associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with a
         *                                 DoubleMatrix value.
         */
        virtual
        JEGA::FrontEnd::Managed::DoubleMatrix^
        GetDoubleMatrix(
            System::String^ tag
            ) = 0;

        /// Supplies the requested parameter as a StringVector from this DB.
        /**
         * \param tag The key by which the requested value is to be retrieved.
         * \return The vector of std::strings associated with the supplied tag.
         * \throws no_such_parameter_error If tag is not associated with a
         *                                 StringVector value.
         */
        virtual
        JEGA::FrontEnd::Managed::StringVector^
        GetStringVector(
            System::String^ tag
            ) = 0;

        /// Maps the supplied integer value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other integral value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddIntegralParam(
            System::String^ tag,
            int value
            ) = 0;

        /// Maps the supplied double value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other double value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddDoubleParam(
            System::String^ tag,
            double value
            ) = 0;

        /// Maps the supplied size_t value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other size_t value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddSizeTypeParam(
            System::String^ tag,
            eddy::utilities::uint64_t value
            ) = 0;

        /// Maps the supplied Boolean value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other Boolean value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddBooleanParam(
            System::String^ tag,
            bool value
            ) = 0;

        /// Maps the supplied string value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other string value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddStringParam(
            System::String^ tag,
            System::String^ value
            ) = 0;

        /// Maps the supplied vector of doubles value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other vector of doubles value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddDoubleVectorParam(
            System::String^ tag,
            DoubleVector^ value
            ) = 0;

        /// Maps the supplied vector of ints value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other vector of ints value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddIntVectorParam(
            System::String^ tag,
            IntVector^ value
            ) = 0;

        /// Maps the supplied matrix of doubles value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other matrix of doubles value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddDoubleMatrixParam(
            System::String^ tag,
            DoubleMatrix^ value
            ) = 0;

        /// Maps the supplied vector of strings value to the supplied tag.
        /**
         * The mapping will only succeed if \a tag is not already the key for
         * some other vector of strings value.
         *
         * \param tag The key to which the requested value is to be mapped.
         * \param value The value that is to be mapped by key.
         * \return True of the mapping succeeds and false otherwise.
         */
        virtual
        bool
        AddStringVectorParam(
            System::String^ tag,
            StringVector^ value
            ) = 0;

}; // class MParameterDatabase



/*
================================================================================
End Namespace
================================================================================
*/
        } // namespace Managed
    } // namespace FrontEnd
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
#endif // JEGA_FRONTEND_MANAGED_MPARAMETERDATABASE_HPP
