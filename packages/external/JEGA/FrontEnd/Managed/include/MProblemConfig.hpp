/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA) Managed Front End

    CONTENTS:

        Definition of class MProblemConfig.

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

        Wed Feb 08 13:40:27 2006 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of the MProblemConfig class.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_FRONTEND_MANAGED_MPROBLEMCONFIG_HPP
#define JEGA_FRONTEND_MANAGED_MPROBLEMCONFIG_HPP

#pragma once


/*
================================================================================
Includes
================================================================================
*/
#include <MConfig.hpp>

#pragma unmanaged
#include <../Utilities/include/JEGATypes.hpp>





/*
================================================================================
Pre-Namespace Forward Declares
================================================================================
*/
namespace JEGA
{
    namespace FrontEnd
    {
        class ProblemConfig;
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
 * \brief A managed wrapper for the JEGA::FrontEnd::ProblemConfig object.
 *
 * This class uses a JEGA::FrontEnd::ProblemConfig instance internally and
 * exposes the required methods in a managed compliant manner for use from
 * .NET languages.
 */
public ref class MProblemConfig :
    public System::IDisposable
{
    /*
    ============================================================================
    Class Scope Enumerations
    ============================================================================
    */
    public:


    protected:


    private:


    /*
    ============================================================================
    Class Scope Typedefs
    ============================================================================
    */
    public:


    protected:


    private:


    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:

        /// A pointer to the actual underlying problem configuration object.
        JEGA::FrontEnd::ProblemConfig* _theConfig;
        

    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        /**
         * \brief A method that can be used to manifest this managed wrapper
         *        into an actual JEGA::FrontEnd::ProblemConfig.
         *
         * The object returned by this method can be used directly in a call
         * to JEGA::FrontEnd::Driver::Driver.
         *
         * \return The JEGA::FrontEnd::ProblemConfig that is actually used by
         *         this wrapper.
         */
        const JEGA::FrontEnd::ProblemConfig&
        Manifest(
            );

        void
        SetDiscardTracking(
            bool track
            );

        bool
        GetDiscardTracking(
            );

        void
        SetMaxGuffSize(
            System::UInt32 mgs
            );

		System::UInt32
        GetMaxGuffSize(
            );

        void
        SetMaxDiscardCacheSize(
            System::UInt32 maxSize
            );

        System::UInt32
        GetMaxDiscardCacheSize(
            );

        /**
         * \brief Adds a continuum, real valued variable to the target
         *        using the supplied information.
         *
         * \param label The label for the new design variable description.
         * \param lowerBound The lowest legal value for the new design
         *                   variable.
         * \param upperBound The highest legal value for the new design
         *                   variable.
         * \param precision The decimal precision to be used by this new
         *                  variable (number of decimal places of interest).
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddContinuumRealVariable method.
         */
        bool
        AddContinuumRealVariable(
            System::String^ label,
            double lowerBound,
            double upperBound,
            int precision
            );

        /**
         * \brief Adds a discrete, real valued variable to the target
         *        using the supplied information.
         *
         * \param label The label for the new design variable description.
         * \param values The initial set of discrete values for this variable.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddDiscreteRealVariable method.
         */
        bool
        AddDiscreteRealVariable(
            System::String^ label,
            DoubleVector^ values
            );

        /**
         * \brief Adds a continuum, integer valued variable to the
         *        target using the supplied information.
         *
         * \param label The label for the new design variable description.
         * \param lowerBound The lowest legal value for the new design
         *                   variable.
         * \param upperBound The highest legal value for the new design
         *                   variable.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddContinuumIntegerVariable method.
         */
        bool
        AddContinuumIntegerVariable(
            System::String^ label,
            int lowerBound,
            int upperBound
            );

        /**
         * \brief Adds a discrete, integer valued variable to the
         *        target using the supplied information.
         *
         * \param label The label for the new design variable description.
         * \param values The initial set of discrete values for this variable.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddDiscreteIntegerVariable method.
         */
        bool
        AddDiscreteIntegerVariable(
            System::String^ label,
            IntVector^ values
            );

        /**
         * \brief Adds a Boolean valued variable to the target using
         *        the supplied information.
         *
         * \param label The label for the new design variable description.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddDesignVariableInfo method.
         */
        bool
        AddBooleanVariable(
            System::String^ label
            );

        /**
         * \brief Adds a linear, minimization objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \param coeffs The coefficients for the linear equation of this
         *               objective.  If not supplied, the created info object
         *               will be unable to evaluate itself.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddLinearMinimizeObjective method.
         */
        bool
        AddLinearMinimizeObjective(
            System::String^ label,
            DoubleVector^ coeffs
            );

        /**
         * \brief Adds a linear, maximization objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \param coeffs The coefficients for the linear equation of this
         *               objective.  If not supplied, the created info object
         *               will be unable to evaluate itself.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddLinearMaximizeObjective method.
         */
        bool
        AddLinearMaximizeObjective(
            System::String^ label,
            DoubleVector^ coeffs
            );

        /**
         * \brief Adds a linear, seek value objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \param value The sought value for this objective.
         * \param coeffs The coefficients for the linear equation of this
         *               objective.  If not supplied, the created info object
         *               will be unable to evaluate itself.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddLinearSeekValueObjective method.
         */
        bool
        AddLinearSeekValueObjective(
            System::String^ label,
            double value,
            DoubleVector^ coeffs
            );

        /**
         * \brief Adds a linear, seek range objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \param lowerBound The lower bound of the range of sought values.
         * \param upperBound The upper bound of the range of sought values.
         * \param coeffs The coefficients for the linear equation of this
         *               objective.  If not supplied, the created info object
         *               will be unable to evaluate itself.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddLinearSeekRangeObjective method.
         */
        bool
        AddLinearSeekRangeObjective(
            System::String^ label,
            double lowerBound,
            double upperBound,
            DoubleVector^ coeffs
            );

        /**
         * \brief Adds a non-linear, minimization objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddNonlinearMinimizeObjective method.
         */
        bool
        AddNonlinearMinimizeObjective(
            System::String^ label
            );

        /**
         * \brief Adds a non-linear, maximization objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddNonlinearMaximizeObjective method.
         */
        bool
        AddNonlinearMaximizeObjective(
            System::String^ label
            );

        /**
         * \brief Adds a non-linear, seek value objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \param value The sought value for this objective.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddNonlinearSeekValueObjective method.
         */
        bool
        AddNonlinearSeekValueObjective(
            System::String^ label,
            double value
            );

        /**
         * \brief Adds a non-linear, seek range objective to the
         *        target using the supplied information.
         *
         * \param label The label for the new objective function description.
         * \param lowerBound The lower bound of the range of sought values.
         * \param upperBound The upper bound of the range of sought values.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddNonlinearSeekRangeObjective method.
         */
        bool
        AddNonlinearSeekRangeObjective(
            System::String^ label,
            double lowerBound,
            double upperBound
            );

        /**
         * \brief Adds a linear, single-sided inequality constraint to the
         *        target using the supplied information.
         *
         * \param label The label for the new constraint description.
         * \param upperLimit The value above which the constraint is considered
         *                   to be violated.
         * \param coeffs The coefficients for the linear equation of this
         *               constraint.  If not supplied, the created info object
         *               will be unable to evaluate itself.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddLinearInequalityConstraint method.
         */
        bool
        AddLinearInequalityConstraint(
            System::String^ label,
            con_val_t upperLimit,
            DoubleVector^ coeffs
            );

        /**
         * \brief Adds a linear, equality constraint to the target
         *        using the supplied information.
         *
         * \param label The label for the new constraint description.
         * \param target The value away from which the constraint is
         *               considered to be violated.
         * \param allowedViol The amount by which strict equality can be
         *                    violated without infeasibility.
         * \param coeffs The coefficients for the linear equation of this
         *               constraint.  If not supplied, the created info object
         *               will be unable to evaluate itself.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddLinearEqualityConstraint method.
         */
        bool
        AddLinearEqualityConstraint(
            System::String^ label,
            con_val_t target,
            double allowedViol,
            DoubleVector^ coeffs
            );

        /**
         * \brief Adds a linear, two-sided inequality constraint to the
         *        target using the supplied information.
         *
         * \param label The label for the new constraint description.
         * \param lowerLimit The value below which the constraint is considered
         *                   to be violated.
         * \param upperLimit The value above which the constraint is considered
         *                   to be violated.
         * \param coeffs The coefficients for the linear equation of this
         *               constraint.  If not supplied, the created info object
         *               will be unable to evaluate itself.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddLinearTwoSidedInequalityConstraint method.
         */
        bool
        AddLinearTwoSidedInequalityConstraint(
            System::String^ label,
            double lowerLimit,
            double upperLimit,
            DoubleVector^ coeffs
            );

        /**
        * \brief Adds a linear, not-equality constraint to the target
        *        using the supplied information.
        *
        * \param label The label for the new constraint description.
        * \param target The value that this constraint cannot take on. Any other
        *               value is acceptable.
        * \param coeffs The coefficients for the linear equation of this
        *               constraint.  If not supplied, the created info object
        *               will be unable to evaluate itself.
        * \return True if the new description is added successfully and false
        *         if not.  This comes directly from the
        *         ConfigHelper::AddLinearNotEqualityConstraint method.
        */
        bool
        AddLinearNotEqualityConstraint(
            System::String^ label,
            double target,
            DoubleVector^ coeffs
            );

        /**
         * \brief Adds a non-linear, single-sided inequality constraint to the
         *        target using the supplied information.
         *
         * \param label The label for the new constraint description.
         * \param upperLimit The value above which the constraint is considered
         *                   to be violated.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddNonlinearInequalityConstraint method.
         */
        bool
        AddNonlinearInequalityConstraint(
            System::String^ label,
            double upperLimit
            );

        /**
         * \brief Adds a non-linear, equality constraint to the target
         *        using the supplied information.
         *
         * \param label The label for the new constraint description.
         * \param target The value away from which the constraint is
         *               considered to be violated.
         * \param allowedViol The amount by which strict equality can be
         *                    violated without infeasibility.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddNonlinearEqualityConstraint method.
         */
        bool
        AddNonlinearEqualityConstraint(
            System::String^ label,
            double target,
            double allowedViol
            );

        /**
        * \brief Adds a non-linear, not-equality constraint to the target
        *        using the supplied information.
        *
        * \param label The label for the new constraint description.
        * \param target The value that this constraint cannot take on. Any other
        *               value is acceptable
        * \return True if the new description is added successfully and false
        *         if not.  This comes directly from the
        *         ConfigHelper::AddNonlinearNotEqualityConstraint method.
        */
        bool
        AddNonlinearNotEqualityConstraint(
            System::String^ label,
            double target
            );

        /**
         * \brief Adds a non-linear, two-sided inequality constraint to the
         *        target using the supplied information.
         *
         * \param label The label for the new constraint description.
         * \param lowerLimit The value below which the constraint is considered
         *                   to be violated.
         * \param upperLimit The value above which the constraint is considered
         *                   to be violated.
         * \return True if the new description is added successfully and false
         *         if not.  This comes directly from the
         *         ConfigHelper::AddNonlinearTwoSidedInequalityConstraint
         *         method.
         */
        bool
        AddNonlinearTwoSidedInequalityConstraint(
            System::String^ label,
            double lowerLimit,
            double upperLimit
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

        /// Implementation of the Dispose method for the IDisposable interface.
        /**
         * This method destroys the underlying JEGA problem configuration.
         * It is safe to call this more than once but it is not safe to use
         * the other methods of this class once this has been called.
         *
         * The destructor should call this method.
         */
        virtual
        void
        DoDispose(
            );

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

        /// Default constructs an MProblemConfig.
        MProblemConfig(
            );

        /// Destructs an MProblemConfig.
        /**
         * The body of this method calls the Dispose method to destroy the
         * JEGA core problem configuration object.
         */
        ~MProblemConfig(
            );


}; // class MProblemConfig



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
#endif // JEGA_FRONTEND_MANAGED_MPROBLEMCONFIG_HPP
