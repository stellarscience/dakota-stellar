/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Definition of class MMessageBoard.

    NOTES:

        See notes under "Document this File" section of this file.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        2.7.0

    CHANGES:

        Wed Feb 01 15:12:07 2012 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of the MMessageBoard class.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_FRONTEND_MANAGED_MMESSAGEBOARD_HPP
#define JEGA_FRONTEND_MANAGED_MMESSAGEBOARD_HPP

#pragma once





/*
================================================================================
Includes
================================================================================
*/
#pragma unmanaged
#include <../Utilities/include/MessageBoard.hpp>

#ifdef JEGA_MESSAGE_BOARD

#pragma managed
#include <vcclr.h>
#include <MConfig.hpp>
#include <ManagedUtils.hpp>







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
ref class MMessageBoard;







/*
================================================================================
In-Namespace File Scope Typedefs
================================================================================
*/

namespace detail
{
class PredicateTie
{
    private:

        gcroot<System::Delegate^> _mPred;

    public:

        inline
        bool
        operator()(
            const JEGA::Utilities::MessageBoard::MessageIdentifier& msgId
            ) const;

        PredicateTie(
            gcroot<System::Delegate^> mPred
            ) :
                _mPred(mPred)
        {}

        PredicateTie(
            const PredicateTie& copy
            ) :
                _mPred(copy._mPred)
        {}
};

class CallbackTie
{
    private:

        gcroot<System::Delegate^> _mCallback;

    public:

        inline
        void
        operator()(
            const JEGA::Utilities::MessageBoard::Message& msg
            ) const;

        CallbackTie(
            gcroot<System::Delegate^> mCallback
            ) :
                _mCallback(mCallback)
        {}

        CallbackTie(
            const CallbackTie& copy
            ) :
                _mCallback(copy._mCallback)
        {}
};
}





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
public ref class MMessageBoard
{
    /*
    ============================================================================
    Class Scope Typedefs
    ============================================================================
    */
    public:

        ref class MessageIdentifier
        {
            private:

                const JEGA::Utilities::MessageBoard::MessageIdentifier& _msgId;

            public:

                System::String^
                Tag(
                    )
                {
                    return ManagedUtils::ToSysString(this->_msgId.Tag());
                }

                System::String^
                SenderType(
                    )
                {
                    return ManagedUtils::ToSysString(this->_msgId.SenderType());
                }

                std::size_t
                Sender(
                    )
                {
                    return reinterpret_cast<std::size_t>(this->_msgId.Sender());
                }

                const JEGA::Utilities::MessageBoard::MessageIdentifier&
                Manifest(
                    )
                {
                    return this->_msgId;
                }

            public:

                MessageIdentifier(
                    const JEGA::Utilities::MessageBoard::MessageIdentifier& msgId
                    ) :
                        _msgId(msgId)
                {}

        };

        ref class Message
        {
            private:

                const JEGA::Utilities::MessageBoard::Message& _msg;
                MessageIdentifier^ _msgId;

            public:

                System::String^
                Text(
                    )
                {
                    return ManagedUtils::ToSysString(this->_msg.Text());
                }

                MessageIdentifier^
                Identifier(
                    )
                {
                    return this->_msgId;
                }

                const JEGA::Utilities::MessageBoard::Message&
                Manifest(
                    )
                {
                    return this->_msg;
                }

            public:

                Message(
                    const JEGA::Utilities::MessageBoard::Message& msg
                    ) :
                        _msg(msg),
                        _msgId(gcnew MessageIdentifier(msg.Identifier()))
                {}

        };

        delegate
        bool
        Predicate(
            MessageIdentifier^ msgId
            );

        delegate
        void
        Callback(
            Message^ msg
            );

        ref class Subscription
        {
            private:

                const JEGA::Utilities::MessageBoard::Subscription* _subsc;

            public:

                bool
                Check(
                    Message^ msg
                    )
                {
                    return this->_subsc->Check(msg->Manifest());
                }

                void
                Service(
                    Message^ msg
                    )
                {
                    this->_subsc->Service(msg->Manifest());
                }

                const JEGA::Utilities::MessageBoard::Subscription&
                Manifest(
                    )
                {
                    return *this->_subsc;
                }

                bool
                operator ==(
                    Subscription^ other
                    )
                {
                    return *this->_subsc == *other->_subsc;
                }

            public:

                Subscription(
                    MMessageBoard::Predicate^ pred,
                    MMessageBoard::Callback^ cb
                    ) :
                        _subsc(0x0)
                {
                    detail::PredicateTie pt(pred);
                    detail::CallbackTie ct(cb);

                    this->_subsc =
                        new JEGA::Utilities::MessageBoard::Subscription(
                            JEGA::Utilities::MessageBoard::Subscribe(pt, ct)
                            );
                }


                Subscription(
                    JEGA::Utilities::MessageBoard::Subscription* subsc
                    ) :
                        _subsc(subsc)
                {}

                ~Subscription(
                    )
                {
                    delete _subsc;
                }

        };

    protected:


    private:

    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:





    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:


        static
        Subscription^
        Subscribe(
            Predicate^ pred,
            Callback^ cb
            )
        {
            return gcnew Subscription(pred, cb);
        }

        static
        bool
        Unsubscribe(
            Subscription^ subscr
            )
        {
            return JEGA::Utilities::MessageBoard::Unsubscribe(
                subscr->Manifest()
                );
        }


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
    private:

        MMessageBoard(
            )
        {}

}; // class MMessageBoard

namespace detail
{
inline
bool
PredicateTie::operator()(
    const JEGA::Utilities::MessageBoard::MessageIdentifier& msgId
    ) const
{
    return static_cast<bool>(this->_mPred->DynamicInvoke(
        gcnew MMessageBoard::MessageIdentifier(msgId)
        ));
}

inline
void
CallbackTie::operator()(
    const JEGA::Utilities::MessageBoard::Message& msg
    ) const
{
    this->_mCallback->DynamicInvoke(
        gcnew MMessageBoard::Message(msg)
        );
}

}


/*
================================================================================
End Namespace
================================================================================
*/
        } // namespace Managed
    } // namespace FrontEnd
} // namespace JEGA




#endif


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
#endif // JEGA_FRONTEND_MANAGED_MMESSAGEBOARD_HPP
