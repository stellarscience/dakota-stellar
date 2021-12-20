/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Definition of class NullOperatorGroup.

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

        Thu May 29 10:25:04 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of the NullOperatorGroup class.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_ALGORITHMS_NULLOPERATORGROUP_HPP
#define JEGA_ALGORITHMS_NULLOPERATORGROUP_HPP







/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <GeneticAlgorithmOperatorGroup.hpp>








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
    namespace Algorithms {







/*
================================================================================
In-Namespace Forward Declares
================================================================================
*/
class NullOperatorGroup;







/*
================================================================================
Class Definition
================================================================================
*/

/// Group containing all the null operators.
/**
 * This group is intended to contain all the null operators.
 */
class NullOperatorGroup :
    public GeneticAlgorithmOperatorGroup
{
    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    private:

        /**
         * \brief Flag indicating whether or not AcquireOperators has been
         *        called for this class yet or not.
         */
        static bool _acquired_operators;




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





    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        /// Returns the one instance of this group as a base class object.
        /**
         * This method is useful for registering groups in operator group
         * registries.  Otherwise, it is generally more appropriate to use the
         * FullInstance() method.
         *
         * \return The one instance of this operator group as a
         *         GeneticAlgorithmOperatorGroup instance.
         */
        inline static
        const GeneticAlgorithmOperatorGroup&
        Instance(
            );

        /// Returns the one instance of this group as its full type.
        /**
         * \return The one instance of this operator group.
         */
        inline static
        const NullOperatorGroup&
        FullInstance(
            );

        /// Returns the registry of mutators valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the mutators allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        MutatorRegistry(
            );

        /// Returns the registry of convergers valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the convergers allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        ConvergerRegistry(
            );

        /// Returns the registry of crossers valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the crossers allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        CrosserRegistry(
            );

        /// Returns the registry of fitness assessors valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the fitness assessors allowed in this
         *         group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        FitnessAssessorRegistry(
            );

        /// Returns the registry of selectors valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the selectors allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        SelectorRegistry(
            );

        /**
         * \brief Returns the registry of niche pressure applicators valid for
         *        this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the niche pressure applicators allowed in
         *         this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        NichePressureApplicatorRegistry(
            );

        /// Returns the registry of initializers valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the initializers allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        InitializerRegistry(
            );

        /// Returns the registry of evaluators valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the evaluators allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        EvaluatorRegistry(
            );

        /// Returns the registry of main loops valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the main loops allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        MainLoopRegistry(
            );

        /// Returns the registry of post processors valid for this group.
        /**
         * \return The operator registry holding the names and creation
         *         functions of all the post processors allowed in this group.
         */
        inline static
        GeneticAlgorithmOperatorRegistry&
        PostProcessorRegistry(
            );

        /// Returns the unique name for this operator group.
        /**
         * \return The string "Null Operator Group".
         */
        static
        const std::string&
        Name(
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

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmMutators for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the mutators allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetMutatorRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmConvergers for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the convergers allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetConvergerRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmCrossers for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the crossers allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetCrosserRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmFitnessAssessors for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the fitness assessors allowed in this
         *         group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetFitnessAssessorRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmSelectors for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the selectors allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetSelectorRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmInitializers for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the initializers allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetInitializerRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmEvaluators for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the evaluators allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetEvaluatorRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmNichePressureApplicators for this
         *        group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the niche pressure applicators allowed in
         *         this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetNichePressureApplicatorRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmMainLoops for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the main loops allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetMainLoopRegistry(
            ) const;

        /**
         * \brief This method returns the operator registry containing all
         *        allowable GeneticAlgorithmPostProcessors for this group.
         *
         * \return The operator registry holding the names and creation
         *         functions of all the post processors allowed in this group.
         */
        virtual
        GeneticAlgorithmOperatorRegistry&
        GetPostProcessorRegistry(
            ) const;

        /// This method returns the unique name of this operator group.
        /**
         * \return See Name().
         */
        virtual
        std::string
        GetName(
            ) const;

    protected:


    private:





    /*
    ============================================================================
    Private Methods
    ============================================================================
    */
    private:
        
        /**
         * \brief This method performs the registration of the operators
         *        that are supposed to be part of this group.
         *
         * It is called from the constructor and only acts on the first
         * instantiation of an object of this type.
         */
        static
        void
        AcquireOperators(
            );

    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

    private:

        /// Constructs a NullOperatorGroup
        /**
         * Construction causes acquisition of operators and absorption of
         * various operator groups.
         */
        NullOperatorGroup(
            );

}; // class NullOperatorGroup


/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Algorithms
} // namespace JEGA







/*
================================================================================
Include Inlined Methods File
================================================================================
*/
#include "./inline/NullOperatorGroup.hpp.inl"



/*
================================================================================
End of Multiple Inclusion Check
================================================================================
*/
#endif // JEGA_ALGORITHMS_NULLOPERATORGROUP_HPP
