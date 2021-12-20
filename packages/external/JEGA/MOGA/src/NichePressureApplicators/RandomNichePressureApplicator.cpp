/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class RandomNichePressureApplicator.

    NOTES:

        See notes of RandomNichePressureApplicator.hpp.

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

        Thu May 10 11:06:08 2018 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the RandomNichePressureApplicator class.
 */




/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <../Utilities/include/Logging.hpp>
#include <../Utilities/include/DesignGroup.hpp>
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/ParameterExtractor.hpp>
#include <utilities/include/RandomNumberGenerator.hpp>
#include <../Utilities/include/MultiObjectiveStatistician.hpp>
#include "../../include/NichePressureApplicators/RandomNichePressureApplicator.hpp"







/*
================================================================================
Namespace Using Directives
================================================================================
*/
using namespace std;
using namespace JEGA::Logging;
using namespace JEGA::Utilities;
using namespace eddy::utilities;






/*
================================================================================
Begin Namespace
================================================================================
*/
namespace JEGA {
    namespace Algorithms {





/*
================================================================================
Static Member Data Definitions
================================================================================
*/
const size_t RandomNichePressureApplicator::DEFAULT_MAX_DESIGNS(100);







/*
================================================================================
Mutators
================================================================================
*/

void
RandomNichePressureApplicator::SetMaximumDesigns(
    std::size_t maxDesigns
    )
{
    EDDY_FUNC_DEBUGSCOPE

    this->_maxDesigns = maxDesigns;

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(),
            this->GetName() + ": Maximum designs now = "
            ) << this->_maxDesigns
        )
}



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

const string&
RandomNichePressureApplicator::Name(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    static const string ret("random");
    return ret;
}

const string&
RandomNichePressureApplicator::Description(
    )
{
    EDDY_FUNC_DEBUGSCOPE

    static const string ret(
        "This niche pressure applicator is designed to choose a limited "
        "number of solutions to remain in the population at random."
        );
    return ret;
}

GeneticAlgorithmOperator*
RandomNichePressureApplicator::Create(
    GeneticAlgorithm& algorithm
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return new RandomNichePressureApplicator(algorithm);
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

string
RandomNichePressureApplicator::GetName(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return RandomNichePressureApplicator::Name();
}

string
RandomNichePressureApplicator::GetDescription(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return RandomNichePressureApplicator::Description();
}

GeneticAlgorithmOperator*
RandomNichePressureApplicator::Clone(
    GeneticAlgorithm& algorithm
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return new RandomNichePressureApplicator(*this, algorithm);
}

bool
RandomNichePressureApplicator::PollForParameters(
    const ParameterDatabase& db
    )
{
    EDDY_FUNC_DEBUGSCOPE
        
    bool success = ParameterExtractor::GetSizeTypeFromDB(
        db, "method.population_size", this->_maxDesigns
        );

    // If we did not find the population size either, warn about it and move
    // on to using the default size.  Note that if !success, then
    // _maxDesigns has not been altered.
    JEGAIFLOG_CF_II(!success, this->GetLogger(), lverbose(), this,
        text_entry(lverbose(), this->GetName() + ": The population size "
            "was not found in the parameter database either.  Using the "
            "default value for the maximum post niching design count.")
        )

    this->SetMaximumDesigns(this->_maxDesigns);

    return this->GeneticAlgorithmNichePressureApplicator::PollForParameters(db);
}

void
RandomNichePressureApplicator::PreSelection(
    DesignGroup& population
    )
{
    EDDY_FUNC_DEBUGSCOPE

    // if we are not caching designs, we needn't do anything here.
    if(!this->GetCacheDesigns()) return;

    // Synchronize the lists just in case.
    population.SynchronizeOFAndDVContainers();

    JEGA_LOGGING_IF_ON(
        const DesignOFSortSet::size_type initPSize = population.SizeOF();
        )

    // Re-assimilate the buffered designs into the population so that they
    // can be considered when making the initial selection.  We will cull
    // them out again when ApplyNichePressure is called later if appropriate.
    this->ReAssimilateBufferedDesigns(population);

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(), this->GetName() + ": Returned ")
            << (population.SizeOF() - initPSize) << " designs during "
               "pre-selection phase of niche pressure application."
        )
}


void
RandomNichePressureApplicator::ApplyNichePressure(
    DesignGroup& population,
    const FitnessRecord& fitnesses
    )
{
    EDDY_FUNC_DEBUGSCOPE

    JEGALOG_II(GetLogger(), ldebug(), this, text_entry(ldebug(),
        "max designs nicher: in use."))

    // If the population is empty, we needn't go any further.
    if(population.IsEmpty()) return;

    // Make sure that the Taboo mark is clear on all designs.
    for(DesignDVSortSet::const_iterator it(population.BeginDV());
        it!=population.EndDV(); ++it) (*it)->ModifyAttribute(TABOO_MARK, false);

    // in case we are not caching, we will need the target below.
    DesignTarget& target = this->GetDesignTarget();

    // we will need the number of objectives for a few things here.
    const size_t nof = target.GetNOF();

    // Synchronize the lists just in case.
    population.SynchronizeOFAndDVContainers();

    // The number of designs to keep for repeated use below.
    const size_t n2Keep = this->GetMaximumDesigns();

    // See if there are fewer solutions in the population than we are to
    // niche to.  If so, we keep them all.
    if(population.SizeOF() < n2Keep) return;

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(), this->GetName() + ": Population size "
            "before niching is ") << population.GetSize() << "."
        )

    const DesignOFSortSet& popByOf = population.GetOFSortContainer();

    JEGA_LOGGING_IF_ON(std::size_t prevPopSize = popByOf.size();)

    vector<const Design*> allDesVec;
    allDesVec.reserve(popByOf.size());

    for(DesignOFSortSet::const_iterator curr(popByOf.begin());
        curr!=popByOf.end(); ++curr) allDesVec.push_back(*curr);
    
    // Start at the end of the list and remove the required number all the while
    // skipping any taboo.
    size_t n2Remove = allDesVec.size() - n2Keep;

    while(n2Remove > 0)
    {
        const int index = RandomNumberGenerator::UniformInt<int>(
			0, static_cast<int>(allDesVec.size() - 1)
			);
        const Design* des = allDesVec[index];
        if(des->HasAttribute(TABOO_MARK)) continue;

        const bool buffered = this->BufferDesign(des);
        population.Erase(des);
        if(!buffered) target.TakeDesign(const_cast<Design*>(des));
        --n2Remove;
    }

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(), this->GetName() + ": Final population size "
            "after niching is ") << population.GetSize() << "."
        )

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
RandomNichePressureApplicator::RandomNichePressureApplicator(
    GeneticAlgorithm& algorithm
    ) :
        GeneticAlgorithmNichePressureApplicator(algorithm),
        _maxDesigns(DEFAULT_MAX_DESIGNS)
{
    EDDY_FUNC_DEBUGSCOPE
}

RandomNichePressureApplicator::RandomNichePressureApplicator(
    const RandomNichePressureApplicator& copy
    ) :
        GeneticAlgorithmNichePressureApplicator(copy),
        _maxDesigns(copy._maxDesigns)
{
    EDDY_FUNC_DEBUGSCOPE
}

RandomNichePressureApplicator::RandomNichePressureApplicator(
    const RandomNichePressureApplicator& copy,
    GeneticAlgorithm& algorithm
    ) :
        GeneticAlgorithmNichePressureApplicator(copy, algorithm),
        _maxDesigns(copy._maxDesigns)
{
    EDDY_FUNC_DEBUGSCOPE
}


/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Algorithms
} // namespace JEGA

