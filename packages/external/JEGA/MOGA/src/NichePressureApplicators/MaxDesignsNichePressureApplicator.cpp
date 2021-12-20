/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class MaxDesignsNichePressureApplicator.

    NOTES:

        See notes of MaxDesignsNichePressureApplicator.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        2.7.0

    CHANGES:

        Wed Dec 21 16:25:44 2011 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the implementation of the MaxDesignsNichePressureApplicator
 *        class.
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
#include <../MOGA/include/NichePressureApplicators/MaxDesignsNichePressureApplicator.hpp>







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
const double MaxDesignsNichePressureApplicator::DEFAULT_DIST_PCT(0.01);
const size_t MaxDesignsNichePressureApplicator::DEFAULT_MAX_DESIGNS(100);







/*
================================================================================
Mutators
================================================================================
*/

void
MaxDesignsNichePressureApplicator::SetDistancePercentages(
    const JEGA::DoubleVector& pcts
    )
{
    EDDY_FUNC_DEBUGSCOPE

    const std::size_t nof = this->GetDesignTarget().GetNOF();

    JEGAIFLOG_CF_II(nof < pcts.size(), this->GetLogger(), lquiet(), this,
        text_entry(lquiet(),
            this->GetName() + ": Received more percentages than there are "
            "objective functions.  Extras will be ignored.")
        )

    JEGAIFLOG_CF_II(nof > pcts.size() && pcts.size() > 1, this->GetLogger(),
        lquiet(), this, ostream_entry(lquiet(),
            this->GetName() + ": Received fewer percentages (") << pcts.size()
            << ") than there are objective functions (" << nof << ").  "
            "Using default value of " << DEFAULT_DIST_PCT << " to fill in."
        )

    JEGAIFLOG_CF_II(nof > pcts.size() && pcts.size() == 1, this->GetLogger(),
        lquiet(), this, ostream_entry(lquiet(),
            this->GetName() + ": Received a single distance percentage for a ")
            << nof << " objective function problem.  Using the supplied value "
            "of " << pcts[0] << " for all objectives."
        )

    this->_distPcts = pcts;

    const double fill_val =
        (this->_distPcts.size() == 1) ? this->_distPcts[0] : DEFAULT_DIST_PCT;

    if(nof > this->_distPcts.size())
        this->_distPcts.resize(
            static_cast<JEGA::DoubleVector::size_type>(nof), fill_val
            );

    // now go through and set each one individually so that they can be checked
    // for legitimacy.
    for(JEGA::DoubleVector::size_type i=0; i<nof; ++i)
        this->SetDistancePercentage(i, this->_distPcts[i]);
}

void
MaxDesignsNichePressureApplicator::SetDistancePercentages(
    double pct
    )
{
    EDDY_FUNC_DEBUGSCOPE

    JEGA_LOGGING_IF_ON(
        static const double minPct = std::numeric_limits<double>::min();
        )

    const std::size_t nof = this->GetDesignTarget().GetNOF();

    JEGAIFLOG_CF_II(pct < 0.0, this->GetLogger(), lquiet(), this,
        ostream_entry(lquiet(),
            this->GetName() + ": Distance percentages must be at least ")
            << minPct << " Supplied value of " << pct
            << " will be replaced by the minimum."
        )

    JEGAIFLOG_CF_II(pct > 1.0, this->GetLogger(), lquiet(), this,
        ostream_entry(lquiet(),
            this->GetName() + ": Distance percentages cannot exceed 100%.  "
            "Supplied value of ") << pct << " will be replaced by 100%."
        )

    pct = Math::Max(0.0, Math::Min(pct, 1.0));

    this->_distPcts.assign(
        static_cast<JEGA::DoubleVector::size_type>(nof), pct
        );

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(),
            this->GetName() + ": All distance percentages now = ") << pct
        )
}

void
MaxDesignsNichePressureApplicator::SetDistancePercentage(
    std::size_t of,
    double pct
    )
{
    EDDY_FUNC_DEBUGSCOPE

    JEGA_LOGGING_IF_ON(
        static const double minPct = std::numeric_limits<double>::min();
        )

    const DesignTarget& target = GetDesignTarget();
    const std::size_t nof = target.GetNOF();
    const JEGA::DoubleVector::size_type dvsof =
        static_cast<JEGA::DoubleVector::size_type>(of);

    // make sure we have enough locations in the percentages vector.
    this->_distPcts.resize(
        static_cast<JEGA::DoubleVector::size_type>(nof), DEFAULT_DIST_PCT
        );

    // now verify the supplied objective function index.
    JEGAIFLOG_CF_II_F(of >= nof, this->GetLogger(), this,
        ostream_entry(lfatal(),
            this->GetName() + ": Request to change objective with index #")
            << of << ".  Valid indices are 0 through " << (nof-1) << "."
        )

    // now verify the supplied value.
    JEGAIFLOG_CF_II(pct < 0.0, this->GetLogger(), lquiet(), this,
        ostream_entry(lquiet(),
            this->GetName() + ": Distance percentages must be at least ")
            << minPct << " Supplied value of " << pct << " for objective \""
            << target.GetObjectiveFunctionInfos()[dvsof]->GetLabel()
            << "\" will be replaced by the minimum."
        )

    JEGAIFLOG_CF_II(pct < 0.0, this->GetLogger(), lquiet(), this,
        ostream_entry(lquiet(),
            this->GetName() + ": Distance percentages cannot exceed 100%.  "
            "Supplied value of ") << pct << " for objective \""
            << target.GetObjectiveFunctionInfos()[dvsof]->GetLabel()
            << "\" will be replaced by 100%."
        )

    pct = Math::Max(0.0, Math::Min(pct, 1.0));

    this->_distPcts[dvsof] = pct;

    JEGALOG_II(this->GetLogger(), lverbose(), this,
        ostream_entry(lverbose(),
            this->GetName() + ": Distance for objective \"")
            << target.GetObjectiveFunctionInfos()[dvsof]->GetLabel()
            << "\" now = " << pct << "."
        )
}

void
MaxDesignsNichePressureApplicator::SetMaximumDesigns(
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
MaxDesignsNichePressureApplicator::Name(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    static const string ret("max_designs");
    return ret;
}

const string&
MaxDesignsNichePressureApplicator::Description(
    )
{
    EDDY_FUNC_DEBUGSCOPE

    static const string ret(
        "This niche pressure applicator is designed to choose a limited "
        "number of solutions to remain in the population.  It does so in "
        "order to balance the tendency for populations to grow very large "
        "and thus consuming too many computer resources.  It operates by "
        "ranking designs according to their fitness standing and a "
        "computed count of how many other designs are too close to them.  Too "
        "close is a function of the supplied niche_vector.  Once the designs "
        "are all ranked, the top max_designs designs are kept in the "
        "population and the remaining ones are buffered or discarded "
        "depending on the value of the cache_niched_designs flag.  Note that "
        "like other niching operators, this one will not discard an extreme "
        "design."
        );
    return ret;
}

GeneticAlgorithmOperator*
MaxDesignsNichePressureApplicator::Create(
    GeneticAlgorithm& algorithm
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return new MaxDesignsNichePressureApplicator(algorithm);
}







/*
================================================================================
Subclass Visible Methods
================================================================================
*/


JEGA::DoubleVector
MaxDesignsNichePressureApplicator::ComputeCutoffDistances(
    const eddy::utilities::extremes<obj_val_t>& exts
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    typedef eddy::utilities::extremes<obj_val_t> extremes_t;

    // the cutoff distance is a percentage of the range of the objective
    // considering only the non-dominated designs.
    const std::size_t nof = this->GetDesignTarget().GetNOF();

    JEGAIFLOG_CF_II_F(nof != exts.size(), GetLogger(), this,
        ostream_entry(lfatal(), this->GetName() + ": Extremes contain "
            "record of ") << exts.size() << " objectives for an "
            << nof << " objective problem."
        )

    // Prepare a vector for return.
    JEGA::DoubleVector ret(nof);

    for(extremes_t::size_type i=0; i<nof; ++i)
        ret[i] = Math::Abs(
            this->GetDistancePercentage(i) * exts.get_range(i)
            );

    // return the square route of the sum of squares.
    return ret;
}

double
MaxDesignsNichePressureApplicator::ComputeObjectiveDistance(
    const Design& des1,
    const Design& des2,
    size_t of
    )
{
    EDDY_FUNC_DEBUGSCOPE
    return Math::Abs(des1.GetObjective(of) - des2.GetObjective(of));
}






/*
================================================================================
Subclass Overridable Methods
================================================================================
*/


string
MaxDesignsNichePressureApplicator::GetName(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return MaxDesignsNichePressureApplicator::Name();
}

string
MaxDesignsNichePressureApplicator::GetDescription(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return MaxDesignsNichePressureApplicator::Description();
}

GeneticAlgorithmOperator*
MaxDesignsNichePressureApplicator::Clone(
    GeneticAlgorithm& algorithm
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return new MaxDesignsNichePressureApplicator(*this, algorithm);
}

bool
MaxDesignsNichePressureApplicator::PollForParameters(
    const ParameterDatabase& db
    )
{
    EDDY_FUNC_DEBUGSCOPE

    bool success = ParameterExtractor::GetDoubleVectorFromDB(
        db, "method.jega.niche_vector", this->_distPcts
        );

    // If we did not find the distance percentages, warn about it and use the
    // default values.  Note that if !success, then _distPcts has not been
    // altered.
    JEGAIFLOG_CF_II(!success, this->GetLogger(), lverbose(), this,
        text_entry(lverbose(), this->GetName() + ": The distance percentages "
            "were not found in the parameter database.  Using the current "
            "values.")
        )

    this->SetDistancePercentages(this->_distPcts);

    success = ParameterExtractor::GetSizeTypeFromDB(
        db, "method.jega.max_designs", this->_maxDesigns
        );

    // If we did not find the maximum number of designs, warn about it and move
    // on to trying the population size.  Note that if !success, then
    // _maxDesigns has not been altered.
    JEGAIFLOG_CF_II(!success, this->GetLogger(), lverbose(), this,
        text_entry(lverbose(), this->GetName() + ": The maximum post niching "
            "design count was not found in the parameter database.  Attempting "
            "to find the population size to use it.")
        )

    if(!success) success = ParameterExtractor::GetSizeTypeFromDB(
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
MaxDesignsNichePressureApplicator::PreSelection(
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
MaxDesignsNichePressureApplicator::ApplyNichePressure(
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

    // Now continue by extracting the populations objective function extremes
    eddy::utilities::extremes<obj_val_t> popExtremes(
        DesignStatistician::GetObjectiveFunctionExtremes(popByOf)
        );

    NicheCountMap ncts(this->ComputeNicheCounts(popByOf, popExtremes));

    // We want to make sure not to discard any "Pareto" extremes.  So find them
    // all and keep them.  Then go into the loop below.
    DesignOFSortSet pareto(GetBest(popByOf, fitnesses));

    // Now continue by extracting the Pareto extremes
    this->TagTabooNicheDesigns(pareto);

    double minFit = DBL_MAX;
    double maxFit = DBL_MIN;

    for(DesignOFSortSet::const_iterator it(popByOf.begin());
        it!=popByOf.end(); ++it)
    {
        const double fit = fitnesses.GetFitness(**it);
        if(fit < minFit) minFit = fit;
        if(fit > maxFit) maxFit = fit;
    }

    const double fitRng = maxFit - minFit;

    const size_t maxNCT = ncts.GetMaxValue();
    const size_t nctRng = maxNCT - ncts.GetMinValue();

    FitnessRecord nicheFits(popByOf.size());
    vector<const Design*> allDesVec;
    allDesVec.reserve(popByOf.size());

    for(DesignOFSortSet::const_iterator curr(popByOf.begin());
        curr!=popByOf.end(); ++curr)
    {
        const Design& des = **curr;
        double normFit = (fitnesses.GetFitness(des) - minFit) / fitRng;
        double normNCV = double(maxNCT - ncts.GetValue(des)) / nctRng;
        nicheFits.AddFitness(&des, normFit + normNCV);
        allDesVec.push_back(&des);
    }

    sort(
        allDesVec.begin(), allDesVec.end(),
        GeneticAlgorithmSelector::FitnessPred(nicheFits)
        );

    // Start at the end of the list and remove the required number all the while
    // skipping any taboo.
    size_t n2Remove = allDesVec.size() - n2Keep;
    size_t index = allDesVec.size()-1;

    while(n2Remove > 0)
    {
        const Design* des = allDesVec[index--];
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
MaxDesignsNichePressureApplicator::NicheCountMap
MaxDesignsNichePressureApplicator::ComputeNicheCounts(
    const DesignOFSortSet& designs,
    const eddy::utilities::extremes<obj_val_t>& exts
    ) const
{
    NicheCountMap ncm(designs.size());
    ncm.SuspendStatistics();

    JEGA::DoubleVector dists(this->ComputeCutoffDistances(exts));

    const size_t nof = this->GetDesignTarget().GetNOF();

    for(
        DesignOFSortSet::const_iterator iit(designs.begin());
        iit!=designs.end(); ++iit
        )
    {
        size_t ncAfter = 1; // 1 to count self
        const Design& id = **iit;

        DesignOFSortSet::const_iterator jit(iit);
        for(++jit; jit!=designs.end(); ++jit)
        {
            const Design& jd = **jit;

            // If jit is too close to iit, then we increment ncAfter and
            // increment the count for jit stored in ncm.  That way, ncm will
            // always have record of those that are too close prior to iit at
            // the beginning of an iit loop.
            const double obj0Dist = this->ComputeObjectiveDistance(id, jd, 0);

            // If the distance at obj0 is large enough, we can get out of this
            // inner loop and move onto the next "iit".  This is b/c of the
            // hierarchical sorting by obj0.
            if(obj0Dist > dists[0]) break;

            // prepare to store whether or not jit is too close.
            bool tooClose = true;

            // We need to see if the distances are all too small.  If any are
            // larger than the cutoff, then jit is far enough away.
            for(size_t of=1; of<nof; ++of)
                if(this->ComputeObjectiveDistance(id, jd, of) > dists[of])
                { tooClose = false; break; }

            // If it is too close, we mark it as such.  Otherwise we move on.
            if(tooClose)
            {
                ++ncAfter;
                ncm.AddToValue(jd, 1);
            }
        }

        ncm.AddToValue(id, ncAfter);
    }

    ncm.ResumeStatistics(true);
    return ncm;
}






/*
================================================================================
Structors
================================================================================
*/



MaxDesignsNichePressureApplicator::MaxDesignsNichePressureApplicator(
    GeneticAlgorithm& algorithm
    ) :
        GeneticAlgorithmNichePressureApplicator(algorithm),
        _maxDesigns(DEFAULT_MAX_DESIGNS)
{
    EDDY_FUNC_DEBUGSCOPE
}

MaxDesignsNichePressureApplicator::MaxDesignsNichePressureApplicator(
    const MaxDesignsNichePressureApplicator& copy
    ) :
        GeneticAlgorithmNichePressureApplicator(copy),
        _maxDesigns(copy._maxDesigns)
{
    EDDY_FUNC_DEBUGSCOPE
}

MaxDesignsNichePressureApplicator::MaxDesignsNichePressureApplicator(
    const MaxDesignsNichePressureApplicator& copy,
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

