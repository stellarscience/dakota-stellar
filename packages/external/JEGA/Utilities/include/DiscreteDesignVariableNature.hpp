/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class DiscreteDesignVariableNature.

    NOTES:

        See notes under section "Class Definition" of this file.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Tue Jun 03 08:55:37 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of the DiscreteDesignVariableNature class.
 */



/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_UTILITIES_DISCRETEDESIGNVARIABLENATURE_HPP
#define JEGA_UTILITIES_DISCRETEDESIGNVARIABLENATURE_HPP







/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <../Utilities/include/JEGATypes.hpp>

#include <../Utilities/include/DesignVariableNatureBase.hpp>






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
    namespace Utilities {









/*
================================================================================
In Namespace File Scope Typedefs
================================================================================
*/






/*
================================================================================
In-Namespace Forward Declares
================================================================================
*/
class DiscreteDesignVariableNature;







/*
================================================================================
Class Definition
================================================================================
*/

/// A nature for discrete design variables.
/**
 * This nature is used to store discrete values for any type variable type.
 * By forwarding work onto this class, DesignVariableInfo objects can get
 * valid random values, boundary values, etc.
 *
 * To be a discrete variable for the purposes of this project,
 * it must be the case that there is a finite number of values for
 * which there is no regular and patterned progression from one to
 * the next.  For example, it is not intended that you use this nature
 * to represent all the integers from some lower bound to some upper
 * bound.  Use the continuum nature with the integer type for that.
 */
class JEGA_SL_IEDECL DiscreteDesignVariableNature :
    public DesignVariableNatureBase
{
    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:

        /**
         * \brief This vector holds all the allowed discrete values for this
         *        variable.
         *
         * They are stored as doubles and their double representation is
         * the index of their location in the vector.
         */
        JEGA::DoubleVector _disVals;




    /*
    ============================================================================
    Mutators
    ============================================================================
    */
    public:





    /*
    ============================================================================
    Accessors
    ============================================================================
    */
    public:

        /// Allows immutable access to the discrete values for this variable.
        /**
         * \return The list of discrete values for this variable.
         */
        inline
        const JEGA::DoubleVector&
        GetDiscreteValues(
            ) const;




    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:



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

        /// Returns the string name of this nature.
        /**
         * \return The string "Discrete".
         */
        virtual
        std::string
        ToString(
            ) const;

        /// Returns an exact duplicate of this nature object.
        /**
         * \param forType The type base with which this nature is being used.
         * \return An exact duplicate of this nature created for use with the
         *         supplied type.
         */
        virtual
        DesignVariableNatureBase*
        Clone(
            DesignVariableTypeBase& forType
            ) const;

        /**
         * \brief Returns the representation of the max value for this nature
         *        as a double.
         *
         * A return of -limits::max indicates failure.
         *
         * \return The representation of the maximum value for this variable or
         *         -limits::max if error.
         */
        virtual
        var_rep_t
        GetMaxRep(
            ) const;

        /**
         * \brief Returns the representation of the min value for this nature
         *        as a double.
         *
         * A return of -limits::max indicates failure.
         *
         * \return The representation of the minimum value for this variable or
         *         -limits::max if error.
         */
        virtual
        var_rep_t
        GetMinRep(
            ) const;

        /**
         * \brief Returns a random representation existing within
         *        the supplied range.
         *
         * A return of -limits::max indicates failure.
         *
         * \param lb The lower bound on the desired random value.
         * \param ub The upper bound on the desired value.
         * \return The representation of a random value inside \a within.
         */
        virtual
        var_rep_t
        GetRandomRep(
            var_rep_t lb,
            var_rep_t ub
            ) const;

        /// Returns the proper representation of \a value as a double.
        /**
         * A return of -limits::max indicates failure.
         *
         * \param value The value to retrieve the representation of.
         * \return The representation of the value \a value.
         */
        virtual
		var_rep_t
        GetRepOf(
            double value
            ) const;

        /// Returns a random valid value for this type as a double.
        /**
         * A return of -limits::max indicates failure.
         *
         * \return A random value for this variable for which IsValidValue will
         *         return true;
         */
        virtual
        double
        GetRandomValue(
            ) const;

        /// Returns the maximum value this nature may have as a double.
        /**
         * A return of -limits::max indicates failure.
         *
         * \return The largest legitimate value for this variable.
         */
        virtual
        double
        GetMaxValue(
            ) const;

        /// Returns the minimum value this nature may have as a double.
        /**
         * A return of -limits::max indicates failure.
         *
         * \return The smallest legitimate value for this variable.
         */
        virtual
        double
        GetMinValue(
            ) const;

        /// Returns the value represented by \a rep as a double.
        /**
         * A return of -limits::max indicates failure.
         *
         * \param rep The representation to convert to a value.
         * \return The value associated with or represented by \a rep.
         */
        virtual
        double
        GetValueOf(
            var_rep_t rep
            ) const;

        /// Returns the nearest valid value to \a value.
        /**
         * A return of -limits::max indicates failure.
         *
         * \param value The value to correct to a valid value.
         * \return The nearest value to \a value for which IsValidValue will
         *         return true;
         */
        virtual
        double
        GetNearestValidValue(
            double value
            ) const;

        /// Returns the nearest valid double rep to \a rep.
        /**
         * A return of -limits::max indicates failure.
         *
         * \param rep The representation to correct to a valid representation.
         * \return The nearest representation to \a rep for which
         *         IsValidRep will return true;
         */
        virtual
        var_rep_t
        GetNearestValidRep(
            var_rep_t rep
            ) const;

        /// Returns the distance between valid representations as a double.
        /**
         * A return of -limits::max indicates failure.
         *
         * \return The increment that exists between consecutive
         *         representations according to the decimal precision.
         */
        virtual
        var_rep_t
        GetDistanceBetweenReps(
            ) const;

        /**
         * \brief The mechanism by which discrete values can be added to this
         *        object.
         *
         * \a value is checked to see that it is unique in the list of discrete
         * values.  Returns true if \a value was successfully added and false
         * otherwise.
         *
         * \param value The value that would be added to the list of discrete
         *              values.
         * \return true if the value is added and false otherwise.
         */
        virtual
        bool
        AddDiscreteValue(
            double value
            );

        /// This method empties the list of discrete values.
        virtual
        void
        ClearDiscreteValues(
            );

        /**
         * \brief This method allows the removal of the specified discrete
         *        value \a value.
         *
         * Returns true if \a value was removed and false otherwise (which
         * usually means that \a value was not found).
         *
         * \param value The value that to be removed from the list of
         *              discrete values.
         * \return true if the value is found and removed and false otherwise.
         */
        virtual
        bool
        RemoveDiscreteValue(
            double value
            );

        /**
         * \brief Sets the upper bound or largest value for this variable to
         *        \a value.
         *
         * This method causes the removal of any discrete values larger than
         * value.
         *
         * \param value The new maximum value for this variable.
         */
        virtual
        void
        SetMaxValue(
            double value
            );

        /**
         * \brief Sets the lower bound or smallest value for this variable to
         *        \a value.
         *
         * This method causes the removal of any discrete values smaller than
         * value.
         *
         * \param value The new minimum value for this variable.
         */
        virtual
        void
        SetMinValue(
            double value
            );

        /**
         * \brief Always returns false b/c discrete values can be handled
         *        by this nature.
         *
         * \return false, always.
         */
        virtual
        bool
        IsDiscreteValueLocked(
            ) const;

        /**
         * \brief Returns true if \a value is one of the discrete values known
         *        to this nature.
         *
         * \param value The value to check to see if it is in bounds or not.
         * \return True if value is in bounds and false otherwise.
         */
        virtual
        bool
        IsValueInBounds(
            double value
            ) const;

        /**
         * \brief Returns true if \a rep lies within the upper and lower bounds
         *        of this variable in terms of representations ( not values ).
         *
         * In order for a rep to be in bounds for a discrete value, it must
         * not only be within the upper and lower bound but must also be a
         * whole number.
         *
         * \param rep The representation to check to see if it is in bounds or
         *            not.
         * \return True if rep is in bounds and false otherwise.
         */
        virtual
        bool
        IsRepInBounds(
            var_rep_t rep
            ) const;

        /**
         * \brief Returns true if this nature can take on values outside of the
         *        bounds and still be evaluated (even though it will no doubt
         *        be infeasible).
         *
         * This is useful to be sure that an algorithm does not generate values
         * that cannot be used.  In this case of a discrete variable, it
         * makes no sense to attempt to use the 10th variable if only 9 exist.
         *
         * \return false, always.
         */
        virtual
        bool
        IsOutOfBoundsDefined(
            ) const;

        /**
         * \brief Returns true b/c the precision refers to representations
         *        which are always integral for this nature.
         *
         * \return true, always.
         */
        virtual
        bool
        IsPrecisionLocked(
            ) const;

        /// Returns true if \a value is a valid value for this variable.
        /**
         * Valid values are those that may be returned by GetRandomValue.
         * This method considers a value to be valid if it is in bounds.
         *
         * \param value The value to check for validity with this variable
         *              nature.
         * \return true if \a value is a valid value for this variable and
         *         false otherwise.
         */
        virtual
        bool
        IsValidValue(
            double value
            ) const;

        /**
         * \brief Returns true if \a rep is the representation of a valid value
         *        for this variable.
         *
         * Valid representations are those that may be returned by
         * GetRandomRep.  This method considers a representation to be
         * valid if it is in bounds.
         *
         * \param rep The representation to check for validity with this
         *            variable nature.
         * \return true if \a rep is a valid representation for this variable
         *         and false otherwise.
         */
        virtual
        bool
        IsValidRep(
            var_rep_t rep
            ) const;

        /// Returns true b/c this is the discrete nature.
        /**
         * \return true, always.
         */
        virtual
        bool
        IsDiscrete(
            ) const;


        /// Returns false b/c this is the discrete nature.
        /**
         * \return false, always.
         */
        virtual
        bool
        IsContinuum(
            ) const;

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

        /// Constructs a DiscreteDesignVariableNature known by \a type.
        /**
         * \param type The type along with which this nature will be used to
         *             describe/define the behavior of a design variable.
         */
        DiscreteDesignVariableNature(
            DesignVariableTypeBase& type
            );

        /// Copy constructs a DiscreteDesignVariableNature known by \a type.
        /**
         * \param copy The existing DesignVariableNatureBase from which to copy
         *             properties into this.
         * \param type The type along with which this nature will be used to
         *             describe/define the behavior of a design variable.
         */
        DiscreteDesignVariableNature(
            const DiscreteDesignVariableNature& copy,
            DesignVariableTypeBase& type
            );

        /// Destructs a DiscreteDesignVariableNature object.
        virtual
        ~DiscreteDesignVariableNature(
            );




}; // class DiscreteDesignVariableNature


/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA







/*
================================================================================
Include Inlined Methods File
================================================================================
*/
#include "./inline/DiscreteDesignVariableNature.hpp.inl"



/*
================================================================================
End of Multiple Inclusion Check
================================================================================
*/
#endif // JEGA_UTILITIES_DISCRETEDESIGNVARIABLENATURE_HPP
