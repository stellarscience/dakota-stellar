/*
================================================================================
    PROJECT:

        John Eddy's Genetic Algorithms (JEGA)

    CONTENTS:

        Implementation of class DiscreteDesignVariableNature.

    NOTES:

        See notes of DiscreteDesignVariableNature.hpp.

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
 * \brief Contains the implementation of the DiscreteDesignVariableNature
 *        class.
 */




/*
================================================================================
Includes
================================================================================
*/
// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>

#include <cfloat>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <utilities/include/Math.hpp>
#include <utilities/include/EDDY_DebugScope.hpp>
#include <utilities/include/RandomNumberGenerator.hpp>
#include <../Utilities/include/DiscreteDesignVariableNature.hpp>

/*
================================================================================
Namespace Using Directives
================================================================================
*/
using namespace std;
using namespace JEGA;
using namespace eddy::utilities;







/*
================================================================================
Begin Namespace
================================================================================
*/
namespace JEGA {
    namespace Utilities {


double
relative_difference(
    double arg_a,
    double arg_b
    )
{
    const double min_val = std::numeric_limits<double>::min();
    const double max_val = std::numeric_limits<double>::max();

    // Screen out NaN's first, if either value is a NaN then the distance is
    // "infinite":
    if(std::isnan(arg_a) || std::isnan(arg_b)) return max_val;

    // Screen out infinites:
    if(fabs(arg_b) > max_val)
    {
        if(fabs(arg_a) > max_val)
            return (arg_a < 0) == (arg_b < 0) ? 0 : max_val;
        else
            return max_val;
    }
    else if(fabs(arg_a) > max_val)
        return max_val;

    // If the values have different signs, treat as infinite difference:
    if(((arg_a < 0.0) != (arg_b < 0.0)) && (arg_a != 0.0) && (arg_b != 0.0))
        return max_val;

    arg_a = fabs(arg_a);
    arg_b = fabs(arg_b);

    // Now deal with zero's, if one value is zero (or denorm) then treat it the
    // same as min_val for the purposes of the calculation that follows:
    if(arg_a < min_val) arg_a = min_val;
    if(arg_b < min_val) arg_b = min_val;

    return std::max(
        fabs((arg_a - arg_b) / arg_a), fabs((arg_a - arg_b) / arg_b
        ));
}

/**
 * \brief A base class for the Min and Max predicates.
 */
template <typename Comp>
class CutoffPred :
    public unary_function<double, bool>
{
    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    protected:

        /// Cutoff value
        argument_type _value;

        /// The comparison operator to use on the passed in values.
        Comp _comp;

    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        /// Comparison operator
        /**
         * \param val The value to test against the stored value.
         * \return True if \a val is greater than or equal to the stored value.
         */
        result_type
        operator ()(
            argument_type val
            )
        {
            EDDY_FUNC_DEBUGSCOPE
            return this->_comp(val, this->_value);
        };

    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        /// Constructs a CutoffPred to compare passed in values to \a val.
        /**
         * \param val The value for this operator to compare passed in values
         *            to.
         */
        CutoffPred(
            argument_type val
            ) :
                _value(val),
                _comp()
        {
            EDDY_FUNC_DEBUGSCOPE
        };

}; // class DiscreteDesignVariableNature::CutoffPred

/**
 * \brief A class to test passed in values to see if they are less
 *        than or equal to a stored value.
 */
class MaxPred :
	public CutoffPred<greater_equal<double> >
{
    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        /// Constructs a MaxPred to compare passed in values to \a val.
        /**
         * \param val The value for this operator to compare passed in values
         *            to using greater than or equal to.
         */
        MaxPred(
            argument_type val
            ) :
                CutoffPred<greater_equal<double> >(val)
        {
            EDDY_FUNC_DEBUGSCOPE
        };

}; // class DiscreteDesignVariableNature::MaxPred


/**
 * \brief A class to test passed in values to see if they are greater
 *        than or equal to a stored value.
 */
class MinPred :
    public CutoffPred<less_equal<double> >
{
    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        /// Constructs a MinPred to compare passed in values to \a val.
        /**
         * \param val The value for this operator to compare passed in values
         *            to using less than or equal to.
         */
        MinPred(
            argument_type val
            ) :
                CutoffPred<less_equal<double> >(val)
        {
            EDDY_FUNC_DEBUGSCOPE
        };

}; // class DiscreteDesignVariableNature::MinPred


class EqComp :
    public unary_function<double, bool>
{
    /*
    ============================================================================
    Member Data Declarations
    ============================================================================
    */
    protected:

        double _value;

        double _threshold;

    /*
    ============================================================================
    Public Methods
    ============================================================================
    */
    public:

        result_type
        operator ()(
            argument_type val
            )
        {
            EDDY_FUNC_DEBUGSCOPE
            return relative_difference(val, this->_value) < this->_threshold;
        };

    /*
    ============================================================================
    Structors
    ============================================================================
    */
    public:

        /// Constructs an EqComp to compare passed in values for equality.
        /**
         * \param threshold The relative difference below which equality is
         *                  obtained in accordance with this operator.
         */
        EqComp(
            double value,
            double threshold
            ) :
                unary_function<double, bool>(),
                _value(value),
                _threshold(threshold)
        {
            EDDY_FUNC_DEBUGSCOPE
        };

}; // class DiscreteDesignVariableNature::MinPred

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
DiscreteDesignVariableNature::ToString(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return string("Discrete");
}

DesignVariableNatureBase*
DiscreteDesignVariableNature::Clone(
    DesignVariableTypeBase& forType
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return new DiscreteDesignVariableNature(*this, forType);
}

var_rep_t
DiscreteDesignVariableNature::GetMaxRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return static_cast<var_rep_t>(this->_disVals.size()-1);
}

var_rep_t
DiscreteDesignVariableNature::GetMinRep(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_disVals.empty() ? var_rep_t(-1) : var_rep_t(0);
}

var_rep_t
DiscreteDesignVariableNature::GetRandomRep(
    var_rep_t lb,
    var_rep_t ub
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return static_cast<var_rep_t>(
        RandomNumberGenerator::UniformInt<DoubleVector::size_type>(
            // if anything, shrink the range.
            static_cast<size_t>(Math::Ceil(static_cast<double>(lb))),
            static_cast<size_t>(Math::Floor(static_cast<double>(ub)))
            )
        );
}

var_rep_t
DiscreteDesignVariableNature::GetRepOf(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    DoubleVector::const_iterator it(find_if(
        this->_disVals.begin(), this->_disVals.end(), EqComp(value, 0.000000001)
        ));

    return (it==this->_disVals.end()) ?
        -numeric_limits<var_rep_t>::max() :
        static_cast<var_rep_t>(distance(this->_disVals.begin(), it));
}

double
DiscreteDesignVariableNature::GetRandomValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    DoubleVector::size_type elem =
        static_cast<DoubleVector::size_type>(
            this->DesignVariableNatureBase::GetRandomRep()
            );
    return this->_disVals[elem];
}

double
DiscreteDesignVariableNature::GetMaxValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_disVals.back();
}

double
DiscreteDesignVariableNature::GetMinValue(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->_disVals.front();
}

double
DiscreteDesignVariableNature::GetValueOf(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    EDDY_ASSERT(this->IsValidRep(rep));
    DoubleVector::size_type loc =
        static_cast<DoubleVector::size_type>(Math::Round(rep));
    return this->IsValidRep(rep) ?
		this->_disVals[loc] : -numeric_limits<double>::max();
}

double
DiscreteDesignVariableNature::GetNearestValidValue(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

	typedef DoubleVector::const_iterator dvci;

    // go through the values (which are in sorted order)
    // and find the closest one.

    // prepare to use the beginning and end of the vector over
    // and over again.
    dvci b(this->_disVals.begin());
    dvci e(this->_disVals.end());

    // bound the passed in value
    pair<dvci, dvci> p = equal_range(b, e, value);

    // if the value actually exists in the list, the lowerbound will be
    // the correct value.
    if(*(p.first) == value) return value;

    // otherwise, if the upperbound is the beginning of the vector, then
    // the value passed in is smaller than all others and we return
    // the lowerbound
    if(p.second == b) return this->_disVals.front();

    // otherwise, if the lowerbound is equal to the end of the vector,
    // then our value is larger than all _disVals and we return the
    // last (largest) value.
    if(p.first == e) return this->_disVals.back();

    // if we make it here, we have a value that is in-between
    // two legitimate values.  To handle this, we have to back up the
    // lowerbound and see which is furthest from the value (lb or ub).
    EDDY_ASSERT(p.first != b);
    --p.first;

    const double ubd = Math::Abs(*(p.second) - value);
    const double lbd = Math::Abs(*(p.first) - value);

    return (ubd > lbd) ? *(p.first) : *(p.second);
}

var_rep_t
DiscreteDesignVariableNature::GetNearestValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    if(rep == -numeric_limits<var_rep_t>::max()) return rep;
    const var_rep_t temp = Math::Round(rep);
    return Math::Max(Math::Min(this->GetMaxRep(), temp), this->GetMinRep());
}

var_rep_t
DiscreteDesignVariableNature::GetDistanceBetweenReps(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return var_rep_t(1);
}

bool
DiscreteDesignVariableNature::AddDiscreteValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE

    DoubleVector::const_iterator it(find_if(
        this->_disVals.begin(), this->_disVals.end(), EqComp(value, 0.000000001)
        ));

    EDDY_DEBUG(it!=this->_disVals.end(),
        "Attempt to add duplicate discrete value failed");

    if(it == this->_disVals.end())
    {
        this->_disVals.insert(lower_bound(
            this->_disVals.begin(), this->_disVals.end(), value
            ), value);
        return true;
    }
    return false;
}

void
DiscreteDesignVariableNature::ClearDiscreteValues(
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->_disVals.clear();
}

bool
DiscreteDesignVariableNature::RemoveDiscreteValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE

    DoubleVector::size_type osize = this->_disVals.size();
    this->_disVals.erase(remove(
        this->_disVals.begin(), this->_disVals.end(), value
        ), this->_disVals.end());
    return this->_disVals.size() != osize;
}

void
DiscreteDesignVariableNature::SetMaxValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
     this->_disVals.erase(
         remove_if(
            this->_disVals.begin(), this->_disVals.end(), MaxPred(value)
            ),
         this->_disVals.end()
         );
    this->AddDiscreteValue(value);
}

void
DiscreteDesignVariableNature::SetMinValue(
    double value
    )
{
    EDDY_FUNC_DEBUGSCOPE
    this->_disVals.erase(
        remove_if(this->_disVals.begin(), this->_disVals.end(), MinPred(value)),
        this->_disVals.end()
        );
    this->AddDiscreteValue(value);
}

bool
DiscreteDesignVariableNature::IsDiscreteValueLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return false;
}

bool
DiscreteDesignVariableNature::IsValueInBounds(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE

    return find_if(
        this->_disVals.begin(), this->_disVals.end(), EqComp(value, 0.000000001)
        ) != this->_disVals.end();
}

bool
DiscreteDesignVariableNature::IsRepInBounds(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return Math::IsWhole(rep) && (rep >= 0.0) && (rep < this->_disVals.size());
}

bool
DiscreteDesignVariableNature::IsOutOfBoundsDefined(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return false;
}

bool
DiscreteDesignVariableNature::IsPrecisionLocked(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return true;
}

bool
DiscreteDesignVariableNature::IsValidValue(
    double value
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->DesignVariableNatureBase::IsValidValue(value) &&
           this->IsValueInBounds(value);
}

bool
DiscreteDesignVariableNature::IsValidRep(
    var_rep_t rep
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return this->DesignVariableNatureBase::IsValidRep(rep) &&
           this->IsRepInBounds(rep);
}

bool
DiscreteDesignVariableNature::IsDiscrete(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return true;
}

bool
DiscreteDesignVariableNature::IsContinuum(
    ) const
{
    EDDY_FUNC_DEBUGSCOPE
    return false;
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
DiscreteDesignVariableNature::DiscreteDesignVariableNature(
    DesignVariableTypeBase& type
    ) :
        DesignVariableNatureBase(type)
{
    EDDY_FUNC_DEBUGSCOPE
}

DiscreteDesignVariableNature::DiscreteDesignVariableNature(
    const DiscreteDesignVariableNature& copy,
    DesignVariableTypeBase& type
    ) :
        DesignVariableNatureBase(copy, type),
        _disVals(copy._disVals)
{
    EDDY_FUNC_DEBUGSCOPE
}

DiscreteDesignVariableNature::~DiscreteDesignVariableNature(
    )
{
    EDDY_FUNC_DEBUGSCOPE
}







/*
================================================================================
End Namespace
================================================================================
*/
    } // namespace Utilities
} // namespace JEGA
