/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Definition of class SpaceFillingNicher.

    NOTES:

        See notes under "Document this File" section of this file.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public
        License as published by the Free Software Foundation; either
        version 2.1 of the License, or (at your option) any later version.
        
        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Lesser General Public License for more details.
        
        You should have received a copy of the GNU Lesser General Public
        License along with this library; if not, write to the Free Software
        Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
        USA

    VERSION:

        2.7.0

    CHANGES:

        Thu Aug 20 08:50:56 2015 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the definition of the SpaceFillingNicher class.
 */




/*
================================================================================
Prevent Multiple Inclusions
================================================================================
*/
#ifndef JEGA_ALGORITHMS_SPACEFILLINGNICHER_HPP
#define JEGA_ALGORITHMS_SPACEFILLINGNICHER_HPP






/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <../Utilities/include/JEGATypes.hpp>
#include <../Utilities/include/DesignValueMap.hpp>
#include <GeneticAlgorithmNichePressureApplicator.hpp>

#ifdef JEGA_HAVE_BOOST
#   include <boost/unordered_map.hpp>
#else
#   include <map>
#endif






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
class SpaceFillingNicher;







/*
================================================================================
In-Namespace File Scope Typedefs
================================================================================
*/
#ifdef JEGA_HAVE_BOOST
#   define MAP_BASE boost::unordered_map
#else
#   define MAP_BASE std::map
#endif







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
class SpaceFillingNicher :
    public GeneticAlgorithmNichePressureApplicator
{
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
    public:
        
        static const std::size_t DEFAULT_NUM_2_KEEP;

    private:

        std::size_t _nDes2Keep;

		mutable MAP_BASE<std::size_t, double> _distCache;

    /*
    ============================================================================
    Mutators
    ============================================================================
    */
    public:

        void
        SetNumDesigns2Keep(
            std::size_t numDesigns
            );


    protected:


    private:


    /*
    ============================================================================
    Accessors
    ============================================================================
    */
    public:

        inline
        const std::size_t&
        GetNumDesigns2Keep(
            ) const;


    protected:


    private:


    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        /// Returns the proper name of this operator.
        /**
         * \return The string "space_filling".
         */
        static
        const std::string&
        Name(
            );

        static
        const std::string&
        Description(
            );

        static
        GeneticAlgorithmOperator*
        Create(
            GeneticAlgorithm& algorithm
            );


    /*
    ============================================================================
    Subclass Visible Methods
    ============================================================================
    */
    protected:
		
        double
        ComputeDistance(
            const JEGA::Utilities::Design& des1,
            const JEGA::Utilities::Design& des2,
			const eddy::utilities::extremes<obj_val_t>& popExtremes
            ) const;

		static inline
        double
        NormalizedObjVal(
            const JEGA::Utilities::Design& des,
            const std::size_t& of,
			const eddy::utilities::extremes<obj_val_t>& popExtremes
            );



    /*
    ============================================================================
    Subclass Overridable Methods
    ============================================================================
    */
    public:

        /**
         * \brief Called prior to the selection operators to prepare for the
         *        selection operations.
         *
         * The selection operators include fitness assessment, selection, and
         * niche pressure,  This method is called prior to them all.  This
         * implementation of it re-assimilates any buffered designs back into
         * the population.
         *
         * \param population The current population prior to selection.  It is
         *                   into this group that buffered designs will be
         *                   placed.
         */
        virtual
        void
        PreSelection(
            JEGA::Utilities::DesignGroup& population
            );

        /**
         * \brief Overridden to carry out the specific niche pressure algorithm.
         *
         * Applies niche pressure to the supplied group according to the
         * assigned fitnesses as can be found in the supplied fitness record.
         *
         * \param population The group (presumably the population of the GA)
         *                   to which to apply niche pressure.
         * \param fitnesses A record of the fitness values assigned to each of
         *                  the members of the supplied group.
         */
        virtual
        void
        ApplyNichePressure(
            JEGA::Utilities::DesignGroup& population,
            const FitnessRecord& fitnesses
            );

        /// Returns the proper name of this operator.
        /**
         * \return See Name().
         */
        virtual
        std::string
        GetName(
            ) const;

        /// Returns a full description of what this operator does and how.
        /**
         * \return See Description().
         */
        virtual
        std::string
        GetDescription(
            ) const;

        /**
         * \brief Creates and returns a pointer to an exact duplicate of this
         *        operator.
         *
         * \param algorithm The GA for which the clone is being created.
         * \return A clone of this operator.
         */
        virtual
        GeneticAlgorithmOperator*
        Clone(
            GeneticAlgorithm& algorithm
            ) const;

        /// Retrieves specific parameters using Get...FromDB methods.
        /**
         * This method is used to extract needed information for this operator.
         * It does so using the "Get...FromDB" class of methods from the
         * GeneticAlgorithmOperator base class.  The return value is indicative
         * of the success of the method.
         *
         * This implementation retrieves the distance percentages for each
         * objective.
         *
         * \param db The database of parameter values from which to do
         *           retrieval.
         * \return true if polling completed successfully and false otherwise.
         */
        virtual
        bool
        PollForParameters(
            const JEGA::Utilities::ParameterDatabase& db
            );


    protected:


    private:


    /*
    ============================================================================
    Private Methods
    ============================================================================
    */
    private:

		static inline
		std::size_t
		hash(
			const JEGA::Utilities::Design& des1, 
			const JEGA::Utilities::Design& des2
			);

    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        /**
         * \brief Constructs an SpaceFillingNicher for use by
         *        \a algorithm.
         *
         * \param algorithm The GA for which this niche pressure applicator is
         *                  being constructed.
         */
        SpaceFillingNicher(
            GeneticAlgorithm& algorithm
            );

        /// Copy constructs an SpaceFillingNicher.
        /**
         * \param copy The instance from which properties should be copied into
         *             this.
         */
        SpaceFillingNicher(
            const SpaceFillingNicher& copy
            );

        /**
         * \brief Copy constructs an SpaceFillingNicher for use by
         *        \a algorithm.
         *
         * \param copy The instance from which properties should be copied into
         *             this.
         * \param algorithm The GA for which this niche pressure  applicator is
         *                  being constructed.
         */
        SpaceFillingNicher(
            const SpaceFillingNicher& copy,
            GeneticAlgorithm& algorithm
            );



}; // class SpaceFillingNicher



/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Algorithms
} // namespace JEGA







/*
================================================================================
Include Inlined Functions File
================================================================================
*/
#include "inline/SpaceFillingNicher.hpp.inl"



/*
================================================================================
End of Multiple Inclusion Check
================================================================================
*/
#endif // JEGA_ALGORITHMS_SPACEFILLINGNICHER_HPP
