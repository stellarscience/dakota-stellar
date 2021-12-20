/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Definition of class MStreamHijacker.

    NOTES:

        See notes under "Document this File" section of this file.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Dec 13 11:51:21 2007 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of managed stream hijacking classes.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_FRONTEND_MANAGED_MSTREAMHIJACKER_HPP
#define JEGA_FRONTEND_MANAGED_MSTREAMHIJACKER_HPP







/*
================================================================================
Includes
================================================================================
*/
#include <MConfig.hpp>
#include <vcclr.h>

#pragma unmanaged
#include <iostream>





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
class StreamHijacker;
ref class MStreamHijacker;
ref class MCoutHijacker;
ref class MCerrHijacker;







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
class StreamHijacker : public std::streambuf
{

    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:

        std::ostream& _hijacked;
        std::streambuf* _prevBuf;
        gcroot<System::IO::TextWriter^> _out;

    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        inline
        void
        detach(
            )
        {
            this->_hijacked.rdbuf(this->_prevBuf);
        }

        inline
        void
        reattach(
            )
        {
            this->_prevBuf = this->_hijacked.rdbuf(this);
        }

        inline
        System::IO::TextWriter^
        text_writer(
            ) const
        {
            return this->_out;
        }

        inline
        std::ostream&
        hijacked_stream(
            ) const
        {
            return this->_hijacked;
        }

    /*
    ============================================================================
    Subclass Overridable Methods
    ============================================================================
    */
    public:


    protected:

        virtual int
        overflow(
            int c = EOF
            );

        virtual
        int
        sync(
            );

    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        StreamHijacker(
            std::ostream& toHijack,
            gcroot<System::IO::TextWriter^> hijackStr
            );

        ~StreamHijacker(
            );
};

/**
 * \brief
 *
 *
 */
public ref class MStreamHijacker abstract
{

    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:

        StreamHijacker* _guts;

    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        void
        detach(
            );

        void
        reattach(
            );

        System::IO::TextWriter^
        text_writer(
            );

        virtual
        void
        on_hijacked_sync(
            );

    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        MStreamHijacker(
            std::ostream& toHijack,
            System::IO::TextWriter^ hijackStr
            );

        ~MStreamHijacker(
            );

};

/**
 * \brief
 *
 *
 */
public ref class MCoutHijacker : public MStreamHijacker
{
    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        MCoutHijacker(
            System::IO::TextWriter^ hijackStr
            );
};

/**
 * \brief
 *
 *
 */
public ref class MCerrHijacker : public MStreamHijacker
{
    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        MCerrHijacker(
            System::IO::TextWriter^ hijackStr
            );
};


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
#endif // JEGA_FRONTEND_MANAGED_MSTREAMHIJACKER_HPP
