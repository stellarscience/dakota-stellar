/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Implementation of class MBasicParameterDatabaseImpl.

    NOTES:

        See notes of MBasicParameterDatabaseImpl.hpp.

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
 * \brief Contains the implementation of the MBasicParameterDatabaseImpl class.
 */




/*
================================================================================
Includes
================================================================================
*/
#include <stdafx.h>
#include <ManagedUtils.hpp>
#include <MBasicParameterDatabaseImpl.hpp>

#pragma unmanaged
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/BasicParameterDatabaseImpl.hpp>
#pragma managed







/*
================================================================================
Namespace Using Directives
================================================================================
*/
using namespace std;
using namespace System;
using namespace JEGA::Utilities;
using namespace System::Collections;







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
ParameterDatabase&
MBasicParameterDatabaseImpl::Manifest(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return *this->_thePDB;
}

bool
MBasicParameterDatabaseImpl::HasIntegralParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasIntegralParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasShortParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasShortParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasDoubleParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasDoubleParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasSizeTypeParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasSizeTypeParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasBooleanParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasBooleanParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasStringParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasStringParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasDoubleVectorParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasDoubleVectorParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasIntVectorParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasIntVectorParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasDoubleMatrixParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasDoubleMatrixParam(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::HasStringVectorParam(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->HasStringVectorParam(ToStdStr(tag));
}

int
MBasicParameterDatabaseImpl::GetIntegral(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->GetIntegral(ToStdStr(tag));
}

short
MBasicParameterDatabaseImpl::GetShort(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->GetShort(ToStdStr(tag));
}

double
MBasicParameterDatabaseImpl::GetDouble(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->GetDouble(ToStdStr(tag));
}

eddy::utilities::uint64_t
MBasicParameterDatabaseImpl::GetSizeType(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->GetSizeType(ToStdStr(tag));
}

bool
MBasicParameterDatabaseImpl::GetBoolean(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->GetBoolean(ToStdStr(tag));
}

System::String^
MBasicParameterDatabaseImpl::GetString(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ManagedUtils::ToSysString(this->_thePDB->GetString(ToStdStr(tag)));
}

JEGA::FrontEnd::Managed::DoubleVector^
MBasicParameterDatabaseImpl::GetDoubleVector(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ToSysDoubleVector(this->_thePDB->GetDoubleVector(ToStdStr(tag)));
}

JEGA::FrontEnd::Managed::IntVector^
MBasicParameterDatabaseImpl::GetIntVector(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ToSysIntVector(this->_thePDB->GetIntVector(ToStdStr(tag)));
}

JEGA::FrontEnd::Managed::DoubleMatrix^
MBasicParameterDatabaseImpl::GetDoubleMatrix(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return nullptr;
}

JEGA::FrontEnd::Managed::StringVector^
MBasicParameterDatabaseImpl::GetStringVector(
    System::String^ tag
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return ToSysStringVector(this->_thePDB->GetStringVector(ToStdStr(tag)));
}

bool
MBasicParameterDatabaseImpl::AddIntegralParam(
    String^ tag,
    int value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddIntegralParam(ToStdStr(tag), value);
}

bool
MBasicParameterDatabaseImpl::AddDoubleParam(
    String^ tag,
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddDoubleParam(ToStdStr(tag), value);
}

bool
MBasicParameterDatabaseImpl::AddSizeTypeParam(
    String^ tag,
    eddy::utilities::uint64_t value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddSizeTypeParam(
        ToStdStr(tag), static_cast<std::size_t>(value)
        );
}

bool
MBasicParameterDatabaseImpl::AddBooleanParam(
    String^ tag,
    bool value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddBooleanParam(ToStdStr(tag), value);
}

bool
MBasicParameterDatabaseImpl::AddStringParam(
    String^ tag,
    String^ value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddStringParam(ToStdStr(tag), ToStdStr(value));
}

bool
MBasicParameterDatabaseImpl::AddDoubleVectorParam(
    String^ tag,
    DoubleVector^ value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddDoubleVectorParam(
        ToStdStr(tag), ToStdDoubleVector(value)
        );
}

bool
MBasicParameterDatabaseImpl::AddIntVectorParam(
    String^ tag,
    IntVector^ value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddIntVectorParam(
        ToStdStr(tag), ToStdIntVector(value)
        );
}

bool
MBasicParameterDatabaseImpl::AddDoubleMatrixParam(
    String^ tag,
    DoubleMatrix^ value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_thePDB->AddDoubleMatrixParam(
        ToStdStr(tag), ToStdDoubleMatrix(value)
        );
}

bool
MBasicParameterDatabaseImpl::AddStringVectorParam(
    String^ tag,
    StringVector^ value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    JEGA::StringVector svec;
    svec.reserve(value->Count);

    System::Collections::Generic::IEnumerator<String^>^ oe =
        value->GetEnumerator();

    while(oe->MoveNext())
        svec.push_back(ToStdStr(safe_cast<String^>(oe->Current)));

    return this->_thePDB->AddStringVectorParam(ToStdStr(tag), svec);
}






/*
================================================================================
Subclass Visible Methods
================================================================================
*/








/*
================================================================================
Subclass Overridable Methods
================================================================================
*/
void
MBasicParameterDatabaseImpl::DoDispose(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    delete this->_thePDB;
    this->_thePDB = 0x0;
}








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


MBasicParameterDatabaseImpl::MBasicParameterDatabaseImpl(
    ) :
        _thePDB(new BasicParameterDatabaseImpl())
{
    EDDY_FUNC_DEBUGSCOPE
}

MBasicParameterDatabaseImpl::~MBasicParameterDatabaseImpl(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    DoDispose();
}

/*
================================================================================
End Namespace
================================================================================
*/
        } // namespace Managed
    } // namespace FrontEnd
} // namespace JEGA


