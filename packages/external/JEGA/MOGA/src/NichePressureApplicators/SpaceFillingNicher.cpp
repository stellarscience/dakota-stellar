/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class SpaceFillingNicher.

    NOTES:

        See notes of SpaceFillingNicher.hpp.

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
 * \brief Contains the implementation of the SpaceFillingNicher class.
 */




/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <FitnessRecord.hpp>
#include <GeneticAlgorithmSelector.hpp>
#include <../Utilities/include/Logging.hpp>
#include <../Utilities/include/DesignGroup.hpp>
#include <utilities/include/EDDY_DebugScope.hpp>
#include <../Utilities/include/ParameterExtractor.hpp>
#include <utilities/include/RandomNumberGenerator.hpp>
#include <../Utilities/include/MultiObjectiveStatistician.hpp>
#include "../../include/NichePressureApplicators/SpaceFillingNicher.hpp"







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
const size_t SpaceFillingNicher::DEFAULT_NUM_2_KEEP(100);








/*
================================================================================
Mutators
================================================================================
*/

void
SpaceFillingNicher::SetNumDesigns2Keep(
    size_t numDesigns
    )
{
    EDDY_FUNC_DEBUGSCOPE

    this->_nDes2Keep = numDesigns;

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(),
            this->GetName() + ": Number of designs now = "
            ) << this->_nDes2Keep
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
SpaceFillingNicher::Name(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    static const string ret("space_filling");
    return ret;
}

const string&
SpaceFillingNicher::Description(
    )
{
    EDDY_FUNC_DEBUGSCOPE

    static const string ret(
        "This niche pressure applicator is designed to ..."
        );
    return ret;
}

GeneticAlgorithmOperator*
SpaceFillingNicher::Create(
    GeneticAlgorithm& algorithm
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return new SpaceFillingNicher(algorithm);
}









/*
================================================================================
Subclass Visible Methods
================================================================================
*/
double
SpaceFillingNicher::ComputeDistance(
    const Design& des1,
    const Design& des2,
	const eddy::utilities::extremes<obj_val_t>& popExtremes
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

	typedef MAP_BASE<size_t, double>::iterator dmIt_T;

    // Start by testing the des1 then des2 hash.  If found, great.
    // If it is not, then test the des2 then des1 hash.

    size_t hVal = hash(des1, des2);

	std::pair<dmIt_T, bool> insDat =
        this->_distCache.insert(std::make_pair(hVal, 0.0));

	if(!insDat.second) return insDat.first->second;

    // The item was not found with des1 before des2.  So we must test des2
    // before des1.  Don't do an insert though and be sure to update the
    // 1 then 2 entry that was just added.
    hVal = hash(des2, des1);
    dmIt_T d2t1 = this->_distCache.find(hVal);
    if(d2t1 != this->_distCache.end())
    {
        insDat.first->second = d2t1->second;
        return d2t1->second;
    }

    const size_t nof = this->GetDesignTarget().GetNOF();

    double dist = 0;
    
    for(size_t of=0; of<nof; ++of)
    {
        const double delta =
			NormalizedObjVal(des1, of, popExtremes) -
			NormalizedObjVal(des2, of, popExtremes);

        dist += delta * delta;
    }

	const double ret = sqrt(dist);
	insDat.first->second = ret;
    return ret;
}






/*
================================================================================
Subclass Overridable Methods
================================================================================
*/

string
SpaceFillingNicher::GetName(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return SpaceFillingNicher::Name();
}

string
SpaceFillingNicher::GetDescription(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return SpaceFillingNicher::Description();
}

GeneticAlgorithmOperator*
SpaceFillingNicher::Clone(
    GeneticAlgorithm& algorithm
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return new SpaceFillingNicher(*this, algorithm);
}

bool
SpaceFillingNicher::PollForParameters(
    const ParameterDatabase& db
    )
{
    EDDY_FUNC_DEBUGSCOPE

    bool success = ParameterExtractor::GetSizeTypeFromDB(
        db, "method.jega.max_designs", this->_nDes2Keep
        );

    // If we did not find the maximum number of designs, warn about it and move
    // on to trying the population size.  Note that if !success, then
    // _nDes2Keep has not been altered.
    JEGAIFLOG_CF_II(!success, this->GetLogger(), lverbose(), this,
        text_entry(lverbose(), this->GetName() + ": The maximum post niching "
            "design count was not found in the parameter database.  Attempting "
            "to find the population size to use it.")
        )

    if(!success) success = ParameterExtractor::GetSizeTypeFromDB(
        db, "method.population_size", this->_nDes2Keep
        );

    // If we did not find the population size either, warn about it and move
    // on to using the default size.  Note that if !success, then
    // _nDes2Keep has not been altered.
    JEGAIFLOG_CF_II(!success, this->GetLogger(), lverbose(), this,
        text_entry(lverbose(), this->GetName() + ": The population size "
            "was not found in the parameter database either.  Using the "
            "default value for the maximum post niching design count.")
        )

    this->SetNumDesigns2Keep(this->_nDes2Keep);

    return this->GeneticAlgorithmNichePressureApplicator::PollForParameters(db);
}

void
SpaceFillingNicher::PreSelection(
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
SpaceFillingNicher::ApplyNichePressure(
    DesignGroup& population,
    const FitnessRecord& fitnesses
    )
{
    EDDY_FUNC_DEBUGSCOPE

    JEGALOG_II(GetLogger(), ldebug(), this, text_entry(ldebug(),
        "Space filling nicher: in use."))

    // If the population is empty, we needn't go any further.
    if(population.IsEmpty()) return;

    // in case we are not caching, we will need the target below.
    DesignTarget& target = this->GetDesignTarget();

    //// For now, since we will only pass Pareto optimal points, we can go ahead
    //// and strip out all dominated.  This will have to go if we modify the
    //// algorithm to take fitness into account.
    //for(DesignDVSortSet::iterator it(population.BeginDV());
    //    it!=population.EndDV();)
    //{
    //    double currFit = fitnesses.GetFitness(**it);
    //    if(currFit < fitnesses.GetMaxFitness())
    //    {
    //        Design* des = *it;

    //        it = population.EraseRetDV(it);

    //        // if we are caching, put design in our buffer.  If not,
    //        // give it back to the target.
    //        if(!this->BufferDesign(des)) target.TakeDesign(des);
    //    }
    //    else ++it;
    //}

    for(DesignDVSortSet::const_iterator it(population.BeginDV());
        it!=population.EndDV(); ++it) (*it)->ModifyAttribute(TABOO_MARK, false);

    // we will need the number of objectives for a few things here.
    const size_t nof = target.GetNOF();

    // Synchronize the lists just in case.
    population.SynchronizeOFAndDVContainers();

    // The the number of designs to keep for repeated use below.
    const size_t n2Keep = this->GetNumDesigns2Keep();

    // See if there are fewer solutions in the population than we are to
    // niche to.  If so, we keep them all.
    if(population.SizeOF() < n2Keep) return;

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(), this->GetName() + ": Population size "
            "before niching is ") << population.GetSize() << "."
        )

    const DesignOFSortSet& popByOf = population.GetOFSortContainer();

    JEGA_LOGGING_IF_ON(size_t prevPopSize = popByOf.size();)

    // Now continue by extracting the populations objective function extremes
    eddy::utilities::extremes<obj_val_t> popExtremes(
        DesignStatistician::GetObjectiveFunctionExtremes(popByOf)
        );

    // we are only going to consider the "best" (which should be the
    // non-dominated) designs as defined by the fitnesses.  We will call them
    // the Pareto even though they may not be.

    // In the future, this may need to be a call to "GetBest".
    const DesignOFSortSet& pareto = population.GetOFSortContainer();

    if(pareto.size() < n2Keep)
    {
        JEGALOG_II(this->GetLogger(), lverbose(), this,
            ostream_entry(lverbose(), this->GetName() + ": Final population "
                "size after niching is ") << population.GetSize() << "."
            )

        this->_distCache.clear();

        return;
    }

    // Now continue by extracting the Pareto extremes
    size_t nTagged = this->TagTabooNicheDesigns(pareto);

    while(nTagged < n2Keep)
    {
        Design* ldDes = 0x0;
        double largeDist = 0.0;

        for(DesignOFSortSet::const_iterator it(pareto.begin());
            it!=pareto.end(); ++it)
        {
            Design& des = **it;
            if(des.HasAttribute(TABOO_MARK)) continue;

            // At this point, we know that des does not have the taboo mark.
            // Therefore, it is a viable candidate to receive it.

            // Find the smallest distance to a marked point.
            double smallDist = std::numeric_limits<double>::max();

            for(DesignOFSortSet::const_iterator jt(pareto.begin());
                jt!=pareto.end(); ++jt)
            {
                // We are looking for the smallest distance to a tagged design.
                // So if jt is not tagged, then we are not interested.
                if(!(*jt)->HasAttribute(TABOO_MARK)) continue;

                // We don't need to test for it!=jt b/c we know "it" is not
                // tagged and "jt" is so they cannot be the same.
                double dist = this->ComputeDistance(des, **jt, popExtremes);
                if(dist < smallDist) smallDist = dist;
            }

            if(smallDist > largeDist)
            {
                largeDist = smallDist;
                ldDes = *it;
            }
        }

        if(ldDes != 0x0)
        {
            ldDes->ModifyAttribute(TABOO_MARK, true);
            ++nTagged;
        }
    }

    for(DesignDVSortSet::iterator it(population.BeginDV());
        it!=population.EndDV();)
    {
        Design* ldDes = *it;

        if(!ldDes->HasAttribute(TABOO_MARK))
        {
            it = population.EraseRetDV(it);

            // if we are caching, put design in our buffer.  If not,
            // give it back to the target.
            if(!this->BufferDesign(ldDes)) target.TakeDesign(ldDes);
        }
        else
        {
            ldDes->ModifyAttribute(TABOO_MARK, false);
            ++it;
        }
    }

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(), this->GetName() + ": Final population size "
            "after niching is ") << population.GetSize() << "."
        )

	this->_distCache.clear();
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



SpaceFillingNicher::SpaceFillingNicher(
    GeneticAlgorithm& algorithm
    ) :
        GeneticAlgorithmNichePressureApplicator(algorithm),
        _nDes2Keep(DEFAULT_NUM_2_KEEP)
{
    EDDY_FUNC_DEBUGSCOPE
}

SpaceFillingNicher::SpaceFillingNicher(
    const SpaceFillingNicher& copy
    ) :
        GeneticAlgorithmNichePressureApplicator(copy),
        _nDes2Keep(copy._nDes2Keep)
{
    EDDY_FUNC_DEBUGSCOPE
}

SpaceFillingNicher::SpaceFillingNicher(
    const SpaceFillingNicher& copy,
    GeneticAlgorithm& algorithm
    ) :
        GeneticAlgorithmNichePressureApplicator(copy, algorithm),
        _nDes2Keep(copy._nDes2Keep)
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

