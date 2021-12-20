' ==============================================================================
#Region "File Preamble"
' ==============================================================================
'     PROJECT:
'
'         John Eddy's Genetic Algorithms (JEGA) Visual Basic Front End
'
'     CONTENTS:
'
'         Declaration of Interface VBJEGAAlgorithmConfigLoader.
'
'     NOTES:
'
'
'
'     PROGRAMMERS:
'
'         John Eddy (jpeddy@sandia.gov) (JE)
'
'     ORGANIZATION:
'
'         Sandia National Laboratories
'
'     COPYRIGHT:
'
'         This library is free software; you can redistribute it and/or
'         modify it under the terms of the GNU Lesser General Public
'         License as published by the Free Software Foundation; either
'         version 2.1 of the License, or (at your option) any later version.
'
'         This library is distributed in the hope that it will be useful,
'         but WITHOUT ANY WARRANTY; without even the implied warranty of
'         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
'         Lesser General Public License for more details.
'
'         You should have received a copy of the GNU Lesser General Public
'         License along with this library; if not, write to the Free Software
'         Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
'         USA
'
'     VERSION:
'
'         1.0.0
'
'     CHANGES:
'
'         Fri Jan 06 10:06:49 2006 - Original Version (JE)
'
' ==============================================================================
#End Region




'===============================================================================
#Region "Options"

Option Explicit On
Option Strict On

#End Region ' Options




'===============================================================================
#Region "Imports"

Imports JEGA.FrontEnd.Managed

#End Region ' Imports




'===============================================================================
' Namespace Declaration
'==============================================================================
Namespace JEGA.FrontEnd.VisualBasic



    '===========================================================================
    ' Interface Declaration
    '===========================================================================
    ''' <summary>
    ''' An interface for a class that wishes to be used as an Algorithm
    ''' configuration object loader for the JEGA VB front end.
    ''' </summary>
    ''' 
    ''' <remarks>
    ''' Use of this interface is not required as there is a version of
    ''' VBJEGADriver.Run that accepts an MAlgorithmConfig directly.  This is
    ''' useful when you would like to neatly encapsulate that functionality
    ''' within a separate class.
    ''' </remarks>
    Public Interface VBJEGAAlgorithmConfigLoader

        '=======================================================================
#Region "Delegates"

#End Region ' Delegates




        '=======================================================================
#Region "Events"

#End Region ' Events




        '=======================================================================
#Region "Methods"

        ''' <summary>
        ''' This is the required method of this interface in which a
        ''' configuration object is to be loaded.
        ''' </summary>
        ''' 
        ''' <remarks>
        ''' The parameter database <paramref>pdb</paramref> can be acquired by
        ''' a call to config.GetParameterDB() but is provided as a separate
        ''' argument for convenience b/c it will pretty much always be needed.
        ''' </remarks>
        ''' 
        ''' <param name="config">
        ''' The configuration object to load up.
        ''' </param>
        ''' <param name="pdb">
        ''' The problem database into which to load algorithm parameters.
        ''' </param>
        ''' <returns>
        ''' True if the loading completes successfully and false otherwise.
        ''' </returns>
        Function _
        LoadConfiguration( _
            ByVal config As MAlgorithmConfig, _
            ByVal pdb As MParameterDatabase _
            ) As Boolean

#End Region ' Methods


    End Interface ' VBJEGAAlgorithmConfigLoader

End Namespace ' JEGA.FrontEnd.VisualBasic