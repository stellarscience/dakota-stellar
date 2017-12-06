/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Inline methods of class MultiObjectiveStatistician.

    NOTES:

        See notes of MultiObjectiveStatistician.hpp.

    PROGRAMMERS:

        John Eddy (jpeddy@sandia.gov) (JE)

    ORGANIZATION:

        Sandia National Laboratories

    COPYRIGHT:

        See the LICENSE file in the top level JEGA directory.

    VERSION:

        1.0.0

    CHANGES:

        Thu May 22 08:02:26 2003 - Original Version (JE)

================================================================================
*/




/*
================================================================================
Document This File
================================================================================
*/
/** \file
 * \brief Contains the inline methods of the MultiObjectiveStatistician class.
 */




/*
================================================================================
Includes
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
File Scope Utility Class Implementations
================================================================================
*/
inline
bool
DominationPredicate::operator ()(
    const Design& des
    ) const
{
    return this->operator ()(&des);
}

inline
bool
DominationPredicate::operator ()(
    const Design* des
    ) const
{
    return MultiObjectiveStatistician::DominationCompare(
            *des, this->_des, this->_infos
            ) == this->_value;
}

inline
DominatesPredicate::DominatesPredicate(
    const Design& des
    ) :
        DominationPredicate(des, 1)
{
} // DominatesPredicate::DominatesPredicate

inline
DominatesPredicate::DominatesPredicate(
    const DominatesPredicate& copy
    ) :
        DominationPredicate(copy)
{
};

inline
DominatedPredicate::DominatedPredicate(
    const Design& des
    ) :
        DominationPredicate(des, -1)
{
} // DominatedPredicate::DominatedPredicate

inline
DominatedPredicate::DominatedPredicate(
    const DominatedPredicate& copy
    ) :
        DominationPredicate(copy)
{
} // DominatedPredicate::DominatedPredicate



/*
================================================================================
Inline Mutators
================================================================================
*/







/*
================================================================================
Inline Accessors
================================================================================
*/







/*
================================================================================
Inline Public Methods
================================================================================
*/
inline
int
MultiObjectiveStatistician::DominationCompare(
    const Design& des1,
    const Design& des2
    )
{
    return DominationCompare(
        des1, des2, des1.GetDesignTarget().GetObjectiveFunctionInfos()
        );

} // MultiObjectiveStatistician::DominationCompare

inline
bool
MultiObjectiveStatistician::IsDominatedByAtLeast1(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in,
    DesignOFSortSet::const_iterator where
    )
{
    return FindDominatingDesign(des, in, where) != in.end();

} // MultiObjectiveStatistician::IsDominatedByAtLeast1

inline
bool
MultiObjectiveStatistician::IsDominatedByAtLeast1(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in
    )
{
    return IsDominatedByAtLeast1(des, in, in.end());

} // MultiObjectiveStatistician::IsDominatedByAtLeast1

inline
DesignOFSortSet::const_iterator
MultiObjectiveStatistician::FindDominatingDesign(
    const Design& des,
    const DesignOFSortSet& in
    )
{
    return FindDominatingDesign(des, in, in.end());

} // MultiObjectiveStatistician::FindDominatingDesign

inline
bool
MultiObjectiveStatistician::DominatesAtLeast1(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in
    )
{
    return DominatesAtLeast1(des, in, in.begin());

} // MultiObjectiveStatistician::DominatesAtLeast1

inline
bool
MultiObjectiveStatistician::DominatesAtLeast1(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in,
    DesignOFSortSet::const_iterator where
    )
{
    return FindDominatedDesign(des, in, where) != in.end();

} // MultiObjectiveStatistician::DominatesAtLeast1

inline
DesignOFSortSet::const_iterator
MultiObjectiveStatistician::FindDominatedDesign(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in
    )
{
    return FindDominatedDesign(des, in, in.begin());

} // MultiObjectiveStatistician::FindDominatedDesign

inline
DesignOFSortSet::const_iterator
MultiObjectiveStatistician::FindDominatedDesign(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in,
    bool hasInfeasible
    )
{
    return FindDominatedDesign(des, in, in.begin(), hasInfeasible);

} // MultiObjectiveStatistician::FindDominatedDesign

inline
std::size_t
MultiObjectiveStatistician::CountNumDominatedBy(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in,
    int cutoff
    )
{
    return CountNumDominatedBy(des, in, in.begin(), cutoff);

} // MultiObjectiveStatistician::CountNumDominatedBy

inline
std::size_t
MultiObjectiveStatistician::CountNumDominatedBy(
    const JEGA::Utilities::Design& des,
    const DesignOFSortSet& in,
    bool hasInfeasible,
    int cutoff
    )
{
    return CountNumDominatedBy(des, in, in.begin(), hasInfeasible, cutoff);

} // MultiObjectiveStatistician::CountNumDominatedBy

inline
std::size_t
MultiObjectiveStatistician::CountNumDominatedBy(
    const Design& des,
    const DesignOFSortSet& in,
    DesignOFSortSet::const_iterator where,
    int cutoff
    )
{
    return CountNumDominatedBy(
        des, in, where,
        DesignStatistician::FindNonFeasible(in) != in.end(), cutoff
        );

} // MultiObjectiveStatistician::CountNumDominatedBy
//
//inline
//std::size_t
//MultiObjectiveStatistician::CountNumDominating(
//    const JEGA::Utilities::Design& des,
//    const DesignOFSortSet& in,
//    int cutoff
//    )
//{
//    return CountNumDominating(des, in, in.end(), cutoff);
//
//} // MultiObjectiveStatistician::CountNumDominating
//
//inline
//std::size_t
//MultiObjectiveStatistician::CountNumDominating(
//    const JEGA::Utilities::Design& des,
//    const DesignOFSortSet& in,
//    DesignOFSortSet::const_iterator where,
//    int cutoff
//    )
//{
//    return CountNumDominating(des, in, where, cutoff);
//
//} // MultiObjectiveStatistician::CountNumDominating

inline
DesignOFSortSet
MultiObjectiveStatistician::GetNonDominated(
    const DesignOFSortSet& from
    )
{
    // prepare a return container which will have the non-dominated designs
    DesignOFSortSet ret;

    for(DesignOFSortSet::const_iterator it(from.begin()); it!=from.end(); ++it)
        if(!IsDominatedByAtLeast1(*(*it), from, it)) ret.insert(*it);

    // return the non-dominated.
    return ret;
}

inline
DesignOFSortSet
MultiObjectiveStatistician::GetDominated(
    const DesignOFSortSet& from
    )
{
    // prepare a return container which will have the dominated designs
    DesignOFSortSet ret;

    for(DesignOFSortSet::const_iterator it(from.begin()); it!=from.end(); ++it)
        if(IsDominatedByAtLeast1(*(*it), from, it)) ret.insert(*it);

    // return the dominated.
    return ret;
}

inline
DesignOFSortSet
MultiObjectiveStatistician::SeparateNonDominated(
    DesignOFSortSet& from
    )
{
    // This MUST be carried out in reverse so that we don't remove designs
    // that might dominate other as yet untested designs.  This sucks because
    // reverse iterators are very poorly supported.  Instead, we will have
    // to take the brute force method of separating the dominated and then
    // swapping.
    DesignOFSortSet dom(SeparateDominated(from));

    // at this point, from contains the non dominated and nonDom contains the
    // dominated.  So we have to swap them.
    std::swap(from, dom);

    // now we can return the nonDominated and leave the dominated in from.
    return dom;
}

inline
DesignOFSortSet
MultiObjectiveStatistician::SeparateDominated(
    DesignOFSortSet& from
    )
{
    // prepare a return container which will have the dominated designs
    DesignOFSortSet ret;

    for(DesignOFSortSet::iterator it(from.begin()); it!=from.end();)
    {
        if(IsDominatedByAtLeast1(*(*it), from, it))
        {
            ret.insert(*it);
            from.erase(it++);
        }
        else ++it;
    }

    // return the dominated.
    return ret;
}

inline
bool
MultiObjectiveStatistician::IsDominatedByAtLeastN(
    const Design& des,
    const DesignOFSortSet& in,
    DesignOFSortSet::size_type n
    )
{
    return IsDominatedByAtLeastN(
            des,
            in,
            n,
            des.IsFeasible() ?
                in.lower_bound(const_cast<Design*>(&des)) :
                in.end()
            );
}

inline
bool
MultiObjectiveStatistician::IsDominatedByAtLeastN(
    const Design& des,
    const DesignOFSortSet& in,
    DesignOFSortSet::size_type n,
    DesignOFSortSet::const_iterator where
    )
{
    return DesignStatistician::CountDesigns(
        in.begin(), where, DominatedPredicate(des), static_cast<int>(n)
        ) >= n;
}

inline
std::size_t
MultiObjectiveStatistician::CountNumDominating(
    const Design& des,
    const DesignOFSortSet& in,
    int cutoff
    )
{
    return CountNumDominating(
        des,
        in,
        des.IsFeasible() ?
            in.lower_bound(const_cast<Design*>(&des)) :
            in.end(),
        cutoff
        );
}

template <typename DesContT>
std::size_t
MultiObjectiveStatistician::TagParetoExtremeDesigns(
    const DesContT& designs,
    const std::size_t& tag
    )
{
    EDDY_FUNC_DEBUGSCOPE

    if(designs.empty()) return 0;

    // The algorithm is as follows:
    //
    // For each objective dimension, find the set of all designs that share
    // the "best" objective function values.  From that subset, visit each
    // other objective function and find the one that has the lowest value in
    // that objective.  It is protected.  If more than one has that best lowest
    // value, then take whatever one, doesn't matter.
    const Design* des = *(designs.begin());
    const DesignTarget& target = des->GetDesignTarget();
    const std::size_t nof = target.GetNOF();

    const ObjectiveFunctionInfoVector& ofInfos =
        target.GetObjectiveFunctionInfos();

    std::size_t nMarked = 0;

    for(size_t oof=0; oof<nof; ++oof)
    {
        std::vector<Design*> bests(GetBestDesigns(designs, *ofInfos[oof]));

        for(size_t iof=0; iof<nof; ++iof)
        {
            if(iof == oof) continue;
            nMarked += MarkOneOfBestDesigns(bests, *ofInfos[iof], tag);
        }
    }

    return nMarked;
}

/*
================================================================================
Inline Subclass Visible Methods
================================================================================
*/








/*
================================================================================
Inline Private Methods
================================================================================
*/


template <typename DesContT>
std::size_t
MultiObjectiveStatistician::MarkOneOfBestDesigns(
    const DesContT& designs,
    const ObjectiveFunctionInfo& ofInfo,
    const std::size_t& tag
    )
{
    EDDY_FUNC_DEBUGSCOPE

    const typename DesContT::const_iterator e(designs.end());
    typename DesContT::const_iterator it(designs.begin());
    Design* best = *it;

    for(++it; it!=e; ++it)
        if(ofInfo.GetPreferredDesign(*best, **it) == *it) best = *it;

    best->ModifyAttribute(tag, true);

    return 1;
}

template <typename DesContT>
std::vector<JEGA::Utilities::Design*>
MultiObjectiveStatistician::GetBestDesigns(
    const DesContT& designs,
    const ObjectiveFunctionInfo& ofInfo
    )
{
    EDDY_FUNC_DEBUGSCOPE

    std::vector<JEGA::Utilities::Design*> ret;
    ret.reserve(4);

    const typename DesContT::const_iterator e(designs.end());
    typename DesContT::const_iterator it(designs.begin());
    ret.push_back(*it);

    for(++it; it!=e; ++it)
    {
        const Design* prefDes = ofInfo.GetPreferredDesign(*ret[0], **it);
        if(prefDes == *it)
        {
            ret.clear();
            ret.push_back(*it);
        }
        else if(prefDes == 0x0)
        {
            ret.push_back(*it);
        }
    }

    return ret;
}





/*
================================================================================
Inline Structors
================================================================================
*/








/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA
