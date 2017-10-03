// $Id: HOPSPACK_SolveLinConstrProj.cpp 167 2010-03-24 00:07:00Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_SolveLinConstrProj.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2010 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This file is part of HOPSPACK.
//
// HOPSPACK is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see http://www.gnu.org/licenses/.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov)
//                 or Todd Plantenga (tplante@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/*!
  @file HOPSPACK_SolveLinConstrProj.cpp
  @brief Implement HOPSPACK::SolveLinConstrProj
*/

#include <math.h>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_SolveLinConstrProj.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//---- ENABLE THIS TO HELP DEBUG THE ACTIVE SET ALGORITHM.
static const bool  bDEBUG = false;


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
SolveLinConstrProj::SolveLinConstrProj (void)
{
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
SolveLinConstrProj::~SolveLinConstrProj (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method solve
//----------------------------------------------------------------------
bool SolveLinConstrProj::solve (const ProblemDef &  cProbDef,
                                const LinConstr  &  cLinConstr,
                                const Vector     &  vX,
                                      Vector     &  vProjection)
{
    if (cLinConstr.hasLinearConstraints() == false)
    {
        //---- NO GENERAL CONSTRAINTS, SO JUST PROJECT TO VARIABLE BOUNDS.
        vProjection = vX;
        if (cProbDef.isBndsFeasible (vProjection) == false)
            cProbDef.makeBndsFeasible (-1.0, vProjection);
        return( true );
    }
    _dActiveTol = cLinConstr.getActiveTol();

    Vector  vScaledX = vX;
    cLinConstr.scale (vScaledX);

    //---- cLinConstr PROVIDES THE LINEAR CONSTRAINTS AND METHODS THAT TEST
    //---- FEASIBILITY OF A POINT.  CONSTRAINTS ARE SCALED AND SHIFTED,
    //---- AND ALREADY TREAT VARIABLE BOUNDS AS THE FIRST INEQUALITIES.
    const Matrix &  mMatEqs = cLinConstr.getAtildeEq();
    const Vector &  vRhsEqs = cLinConstr.getBtildeEq();
    const Matrix &  mMatIneqs = cLinConstr.getAhat();
    const Vector &  vLoIneqs = cLinConstr.getBhatLower();
    const Vector &  vUpIneqs = cLinConstr.getBhatUpper();
    if (bDEBUG)
    {
        cout << "===== BEGIN Debugging for active set method =====" << endl;
        cout << "ActSet problem definition, shifted and scaled:";
        mMatEqs.formattedPrint ("  Eqs ", cout);  cout << endl;
        cout << "  Eqs RHS: ";  vRhsEqs.leftshift (cout);
        mMatIneqs.formattedPrint ("  Bnds+Ineqs ", cout);  cout << endl;
        cout << "  Ineq LBnds: ";  vLoIneqs.leftshift (cout);  cout << endl;
        cout << "  Ineq UBnds: ";  vUpIneqs.leftshift (cout);  cout << endl;
        cout << "  Target Point:    ";  vScaledX.leftshift (cout);  cout << endl;
    }

    //---- PROJECT TO VARIABLE BOUNDS.
    vProjection = vX;
    if (cProbDef.isBndsFeasible (vProjection) == false)
        cProbDef.makeBndsFeasible (-1.0, vProjection);
    Vector  vScaledProjection = vProjection;
    cLinConstr.scale (vScaledProjection);
    if (bDEBUG)
    {
        cout << "  Bnds feas Point: ";
        vScaledProjection.leftshift (cout);  cout << endl;
        cout << endl;
    }

    //---- FIND AN INITIAL POINT FEASIBLE WITH RESPECT TO INEQUALITIES.
    //---- THIS CAN FAIL IF INEQUALITIES CANNOT BE SATISFIED, OR CONSTRAINTS
    //---- ARE DEPENDENT.
    if (findFeasibleIneqPoint_ (cLinConstr,
                                mMatIneqs, vLoIneqs, vUpIneqs,
                                vScaledProjection) == false)
    {
        if (bDEBUG)
            cout << "===== END Debugging for active set method =====" << endl;
        return( false );
    }
    if (bDEBUG)
    {
        cout << "  Ineq feas Point: ";
        vScaledProjection.leftshift (cout);  cout << endl;
        cout << endl;
    }

    //---- FROM THE INITIAL FEASIBLE POINT, FIND THE CLOSEST FEASIBLE POINT.
    //---- THIS CAN FAIL IF EQUALITIES CANNOT BE SATISFIED, OR CONSTRAINTS
    //---- ARE DEPENDENT.
    if (findClosestPoint_ (mMatEqs, vRhsEqs,
                           mMatIneqs, vLoIneqs, vUpIneqs,
                           vScaledX, vScaledProjection) == false)
    {
        if (bDEBUG)
            cout << "===== END Debugging for active set method =====" << endl;
        return( false );
    }

    //---- TEST FEASIBILITY USING THE TOLERANCES DEFINED BY THE USER.
    //---- IF IT FAILS, THEN TOLERANCES SHOULD BE LOOSENED OR SCALING IMPROVED.
    vProjection = vScaledProjection;
    cLinConstr.unscale (vProjection);
    if (bDEBUG)
    {
        cout << "ActSet final unscaled point: " << endl;
        cout << "   ";  vProjection.leftshift (cout);  cout << endl;
    }
    if (cLinConstr.isFeasible (vProjection, false) == false)
    {
        cerr << "WARNING: Active set point declared infeasible by"
             << " LinConstr.isFeasible()" << endl;
        cLinConstr.isFeasible (vProjection, true);
        cerr << "         Check problem scaling and the Active Tolerance"
             << " parameter" << endl;
        if (bDEBUG)
            cout << "===== END Debugging for active set method =====" << endl;
        return( false );
    }

    if (bDEBUG)
        cout << "===== END Debugging for active set method =====" << endl;

    return( true );
}


//----------------------------------------------------------------------
//  Private method findFeasiblePoint_
//----------------------------------------------------------------------
bool SolveLinConstrProj::findFeasibleIneqPoint_
         (const LinConstr &  cLinConstr,
          const Matrix    &  mMatIneqs,
          const Vector    &  vLoIneqs,
          const Vector    &  vUpIneqs,
                Vector    &  vX) const
{
    //---- FORM A LINEAR SYSTEM WITH BOUNDS AND INEQUALITIES ONLY.
    LinConstr  cLinConstrNoEqs (cLinConstr, true);

    //---- CHECK IF THE INITIAL POINT IS ALREADY INEQUALITY FEASIBLE.
    //---- THE OBJECT WILL BE USED LATER TO VERIFY THE FINAL POINT.
    Vector  vUnscaledX = vX;
    cLinConstrNoEqs.unscale (vUnscaledX);
    if (cLinConstrNoEqs.isFeasible (vUnscaledX, false) == true)
        return( true );

    //---- THE FIRST N INEQUALITIES ARE DOUBLE-SIDED BOUND CONSTRAINTS,
    //---- WHICH ARE ASSUMED TO ALREADY BE SATISFIED.  ASSIGN A NEW
    //---- VIOLATION VARIABLE FOR EACH INEQUALITY AFTER THE BOUNDS.
    int  nNumVars = mMatIneqs.getNcols();
    int  nNumViolVars = 0;
    for (int  i = nNumVars; i < mMatIneqs.getNrows(); i++)
    {
        if (exists (vLoIneqs[i]))
            nNumViolVars++;
        if (exists (vUpIneqs[i]))
            nNumViolVars++;
    }

    //---- THE SUBPROBLEM INEQUALITY MATRIX HAS VARIABLE BOUNDS, THEN
    //---- VIOLATION VARIABLE BOUNDS, AND THEN INEQUALITIES.  ANY DOUBLE-SIDED
    //---- INEQUALITIES ARE CONVERTED TO A PAIR OF SUBPROBLEM INEQUALITIES.
    Matrix  mSubprobIneqs;
    Vector  vSubprobLoIneqs (nNumVars + (2 * nNumViolVars));
    Vector  vSubprobUpIneqs (nNumVars + (2 * nNumViolVars));
    for (int  i = 0; i < nNumVars; i++)
    {
        Vector  vTmp (mMatIneqs.getRow (i));
        vTmp.append (nNumViolVars, 0.0);
        mSubprobIneqs.addRow (vTmp);
        vSubprobLoIneqs[i] = vLoIneqs[i];
        vSubprobUpIneqs[i] = vUpIneqs[i];
    }
    for (int  i = 0; i < nNumViolVars; i++)
    {
        Vector  vTmp (nNumVars + nNumViolVars, 0.0);
        vTmp[nNumVars + i] = 1.0;
        mSubprobIneqs.addRow (vTmp);
        vSubprobLoIneqs[nNumVars + i] = 0.0;
        vSubprobUpIneqs[nNumVars + i] = dne();
    }

    //---- INITIALIZE VIOLATION VARIABLES AND DEFINE THEIR INEQUALITIES.
    Vector  vInitPt (nNumVars + nNumViolVars);
    for (int  i = 0; i < nNumVars; i++)
        vInitPt[i] = vX[i];
    int  k = nNumVars + nNumViolVars;
    for (int  i = nNumVars; i < mMatIneqs.getNrows(); i++)
    {
        //---- NOTE: IF computeActiveSolution_ WERE CHANGED TO ACCEPT
        //---- AN INITIAL GUESS OF THE ACTIVE SET, THEN THE INEQUALITY
        //---- SHOULD BE ACTIVE IF VIOLATION IS POSITIVE, AND THE
        //---- VIOLATION VARIABLE BOUND CONSTRAINT OTHERWISE.

        const Vector &  vIneqRow = mMatIneqs.getRow (i);
        double  dLhs = vIneqRow.dot (vX);
        if (exists (vLoIneqs[i]))
        {
            Vector  vTmp (mMatIneqs.getRow (i));
            vTmp.append (nNumViolVars, 0.0);
            vTmp[k - nNumViolVars] = 1.0;
            mSubprobIneqs.addRow (vTmp);
            vSubprobLoIneqs[k] = vLoIneqs[i];
            vSubprobUpIneqs[k] = dne();

            double  dViolationValue = vLoIneqs[i] - dLhs;
            if (dViolationValue > 0.0)
                vInitPt[k - nNumViolVars] = dViolationValue;
            else
                vInitPt[k - nNumViolVars] = 0.0;
            k++;
        }
        if (exists (vUpIneqs[i]))
        {
            Vector  vTmp (mMatIneqs.getRow (i));
            vTmp.append (nNumViolVars, 0.0);
            vTmp[k - nNumViolVars] = -1.0;
            mSubprobIneqs.addRow (vTmp);
            vSubprobLoIneqs[k] = dne();
            vSubprobUpIneqs[k] = vUpIneqs[i];

            double  dViolationValue = dLhs - vUpIneqs[i];
            if (dViolationValue > 0.0)
                vInitPt[k - nNumViolVars] = dViolationValue;
            else
                vInitPt[k - nNumViolVars] = 0.0;
            k++;
        }
    }

    if (bDEBUG)
    {
        cout << "ActSet feasible point subproblem, shifted and scaled:";
        mSubprobIneqs.formattedPrint ("  Bnds+Ineqs ", cout);  cout << endl;
        cout << "  Ineq LBnds: ";  vSubprobLoIneqs.leftshift (cout);
        cout << endl;
        cout << "  Ineq UBnds: ";  vSubprobUpIneqs.leftshift (cout);
        cout << endl;
        cout << "  init Point:    ";  vInitPt.leftshift (cout);  cout << endl;
    }

    //---- DEFINE THE OBJECTIVE FUNCTION TO BE ||v||_2 AND SOLVE.
    //---- IF THE PRIMAL VARIABLES x ARE NOT IN THE OBJECTIVE AND ARE NOT
    //---- CONSTRAINED, THEN LAPACK dgglse IS FREE TO SET x VALUES.
    //---- SOMETIMES (FOR INSTANCE, test_ineqs_3() ON MAC OSX) IT CHOOSES
    //---- GARBAGE VALUES.  THEREFORE, INCLUDE A SMALL PENALTY TERM ON x
    //---- TO KEEP IT NEAR THE STARTING POINT.
    const double  dSMALL_PENALTY = 10.0 * getMachineEpsilon();
    Vector  vZero (nNumVars + nNumViolVars, 0.0);
    Vector  vObj (nNumVars + nNumViolVars, 0.0);
    for (int  i = 0; i < nNumVars; i++)
    {
        vZero[i] = dSMALL_PENALTY * vInitPt[i];
        vObj[i]  = dSMALL_PENALTY;
    }
    for (int  i = nNumVars; i < nNumVars + nNumViolVars; i++)
    {
        vObj[i] = 1.0;
    }
    Matrix  mEmptyEqs;
    Vector  vEmptyRhsEqs;
    Vector  vFeasPoint (nNumVars + nNumViolVars);

    bool  bOK = computeActiveSetSolution_ (vZero,
                                           vObj,
                                           vInitPt,
                                           mEmptyEqs, vEmptyRhsEqs,
                                           mSubprobIneqs,
                                           vSubprobLoIneqs, vSubprobUpIneqs,
                                           vFeasPoint);
    if (bOK == false)
    {
        cerr << "ERROR: Could not find a feasible start point" << endl;
        cerr << "       Active set subproblem failed" << endl;
        return( false );
    }

    //---- CHECK THAT VIOLATION VARIABLES WERE REDUCED TO ZERO.
    double  dTmp = 0.0;
    for (int  i = 0; i < nNumViolVars; i++)
        dTmp = dTmp + vFeasPoint[nNumVars + i];
    if (dTmp > (2.0 * getMachineEpsilon()))
    {
        cerr << "ERROR: Cannot find a feasible start point" << endl;
        cerr << "       Linear constraints may be inconsistent" << endl;
        return( false );
    }

    for (int  i = 0; i < nNumVars; i++)
    {
        vX[i] = vFeasPoint[i];
        vUnscaledX[i] = vX[i];
    }
    cLinConstrNoEqs.unscale (vUnscaledX);
    if (cLinConstrNoEqs.isFeasible (vUnscaledX, false) == false)
    {
        cerr << "ERROR: Cannot find a feasible start point" << endl;
        cerr << "       findFeasiblePoint_ has a point, but failed the feasibility test" << endl;
        return( false );
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private method findClosestPoint_
//----------------------------------------------------------------------
bool SolveLinConstrProj::findClosestPoint_ (const Matrix &  mMatEqs,
                                            const Vector &  vRhsEqs,
                                            const Matrix &  mMatIneqs,
                                            const Vector &  vLoIneqs,
                                            const Vector &  vUpIneqs,
                                            const Vector &  vXtarget,
                                                  Vector &  vX) const
{
    Vector  vXinitial = vX;
    Vector  vIdentity (vX.size());
    for (int  i = 0; i < vIdentity.size(); i++)
        vIdentity[i] = 1.0;

    bool  bOK = computeActiveSetSolution_ (vXtarget,
                                           vIdentity,
                                           vXinitial,
                                           mMatEqs, vRhsEqs,
                                           mMatIneqs, vLoIneqs, vUpIneqs,
                                           vX);
    return( bOK );
}


//----------------------------------------------------------------------
//  Private method computeActiveSetSolution_
//----------------------------------------------------------------------
bool  SolveLinConstrProj::computeActiveSetSolution_ (const Vector &  vC,
                                                     const Vector &  vD,
                                                     const Vector &  vXinit,
                                                     const Matrix &  mMatEqs,
                                                     const Vector &  vRhsEqs,
                                                     const Matrix &  mMatIneqs,
                                                     const Vector &  vLoIneqs,
                                                     const Vector &  vUpIneqs,
                                                           Vector &  vXsol) const
{
    //---- HANDLE THE TRIVIAL CASE OF NO CONSTRAINTS.
    if ((mMatEqs.getNrows() + mMatIneqs.getNrows()) == 0)
    {
        calcUnconstrainedSolution_ (vC, vD, vXsol);
        return( true );
    }

    if (bDEBUG)
    {
        //---- TEST THE INPUT TO HELP DEBUG WHEN THE CALLER MESSES UP,
        //---- BUT THIS IS NOT AN EXACT TEST.
        if (isIneqFeasible_ (vXinit, mMatIneqs, vLoIneqs, vUpIneqs) == false)
        {
            cerr << "WARNING: Call to computeActiveSetSolution_"
                 << " has infeasible start point" << endl;
            cerr << "         Linear constraints may be dependent" << endl;
        }
    }

    int  nNumEqs = mMatEqs.getNrows();
    int  nNumIneqs = mMatIneqs.getNrows();

    bool *  pbaIsActive = new bool[nNumIneqs];
    bool *  pbaIsLowerBnd = new bool[nNumIneqs];

    //---- INITIAL ACTIVE SET IS JUST THE EQUALITY CONSTRAINTS.
    for (int  i = 0; i < nNumIneqs; i++)
        pbaIsActive[i] = false;

    //---- ITERATE UNTIL FINISHED.
    //---- THE ITERATION LIMIT IS JUST A SAFEGUARD AGAINST INFINITE LOOPS;
    //---- HOPEFULLY, THEY NEVER HAPPEN.
    Matrix  mActiveCon = mMatEqs;
    Vector  vActiveRHS = vRhsEqs;
    Vector  vXcurrent = vXinit;
    bool    bFoundSolution = false;
    int     nMaxIters = 3 * (vXinit.size() + nNumEqs + nNumIneqs);
    int     nCurrentIter = 0;
    while (nCurrentIter < nMaxIters)
    {
        nCurrentIter++;

        //---- SET UP THE CURRENT EQP AND CALL LAPACK TO SOLVE IT.
        int  nNumActiveIneqs = mActiveCon.getNrows();
        for (int  i = nNumEqs; i < nNumActiveIneqs; i++)
            mActiveCon.deleteRow (nNumEqs);
        for (int  i = vActiveRHS.size() - 1; i >= nNumEqs; i--)
            vActiveRHS.erase (i);
        for (int  i = 0; i < nNumIneqs; i++)
        {
            if (pbaIsActive[i] == true)
            {
                //---- ADD THE ACTIVE INEQUALITY.
                if (pbaIsLowerBnd[i] == true)
                {
                    mActiveCon.addRow (mMatIneqs.getRow (i));
                    vActiveRHS.push_back (vLoIneqs[i]);
                    if (bDEBUG)
                        cout << "ActSet subproblem adding inequality " << i
                             << " with lower bound" << endl;
                }
                else
                {
                    //---- FLIP THE SIGN ON UPPER BOUNDS SO THE SUBPROBLEM
                    //---- INEQUALITY IS IN THE FORM c(x) > b.  THEN ALL
                    //---- ACTIVE MULTIPLIERS WILL BE POSITIVE.
                    const Vector  vTmp = mMatIneqs.getRow (i);
                    Vector  vNew = vTmp;
                    vNew.scale (-1.0);
                    mActiveCon.addRow (vNew);
                    vActiveRHS.push_back (-1.0 * vUpIneqs[i]);
                    if (bDEBUG)
                        cout << "ActSet subproblem adding inequality " << i
                             << " with upper bound" << endl;
                }
            }
        }
        //---- ALWAYS CALL LSQR WITH THE UNCONSTRAINED SOLUTION IN CASE THE
        //---- CONSTRAINT SET IS EMPTY.
        calcUnconstrainedSolution_ (vC, vD, vXsol);
        if (bDEBUG)
        {
            cout << "ActSet calling LSQR (LAPACK.dgglse) for new point with";
            mActiveCon.formattedPrint ("  Cons ", cout);  cout << endl;
            cout << "  RHS: ";  vActiveRHS.leftshift (cout);  cout << endl;
            cout << "  c: ";  vC.leftshift (cout);  cout << endl;
            cout << "  d: ";  vD.leftshift (cout);  cout << endl;
        }
        if (mActiveCon.generalConstrainedLSQR (vC, vD,
                                               vActiveRHS, vXsol) == false)
        {
            cerr << "WARNING: Call to solve LSQR subproblem failed" << endl;
            cerr << "         Linear constraints may be dependent" << endl;
            bFoundSolution = false;
            break;
        }
        if (bDEBUG)
        {
            cout << "  LSQR returns: ";  vXsol.leftshift (cout);
            cout << endl << endl;
        }


        //---- FIND THE FIRST INEQUALITY (IF ANY) VIOLATED WHEN MOVING TO
        //---- THE EQP SOLUTION.
        double  dMaxStep = 1.0;
        int     nFirstViolatedIneq = -1;
        bool    bIsLowerViolated = false;
        for (int  i = 0; i < nNumIneqs; i++)
        {
            if (pbaIsActive[i] == false)
            {
                double  dLHS = (mMatIneqs.getRow(i)).dot (vXsol);
                if (exists (vLoIneqs[i]) && (dLHS < vLoIneqs[i]))
                {
                    //---- THE LOWER BOUND OF THIS INEQUALITY IS VIOLATED.
                    double  dInit = (mMatIneqs.getRow(i)).dot (vXcurrent);
                    double  dStep = (dInit - vLoIneqs[i]) / (dInit - dLHS);
                    if (dStep < dMaxStep)
                    {
                        dMaxStep = dStep;
                        nFirstViolatedIneq = i;
                        bIsLowerViolated = true;
                    }
                }
                else if (exists (vUpIneqs[i]) && (dLHS > vUpIneqs[i]))
                {
                    //---- THE UPPER BOUND OF THIS INEQUALITY IS VIOLATED.
                    double  dInit = (mMatIneqs.getRow(i)).dot (vXcurrent);
                    double  dStep = (vUpIneqs[i] - dInit) / (dLHS - dInit);
                    if (dStep < dMaxStep)
                    {
                        dMaxStep = dStep;
                        nFirstViolatedIneq = i;
                        bIsLowerViolated = false;
                    }
                }
            }
        }
        if (nFirstViolatedIneq >= 0)
        {
            if ( (dMaxStep < 0.0) && (dMaxStep > -_dActiveTol) )
                dMaxStep = 0.0;
            if (dMaxStep < 0.0)
            {
                //---- THIS CAN HAPPEN IF THE vXinit IS INFEASIBLE TO
                //---- BEGIN WITH, BUT SO BADLY ILL-CONDITIONED THAT
                //---- LinConstr.isFeasible() CANNOT DETECT IT; FOR INSTANCE,
                //---- IN test_ineqs_5().
                cerr << "ERROR: computeActiveSetSolution_ became infeasible"
                     << endl;
                bFoundSolution = false;
                break;
            }

            //---- ADD THE CONSTRAINT AND LOOP TO THE NEXT EQP.
            pbaIsActive[nFirstViolatedIneq] = true;
            pbaIsLowerBnd[nFirstViolatedIneq] = bIsLowerViolated;
            for (int  k = 0; k < vXcurrent.size(); k++)
                vXcurrent[k] += dMaxStep * (vXsol[k] - vXcurrent[k]);
            if (bDEBUG)
            {
                cout << "ActSet step to violated ineq = " << dMaxStep << endl;
                cout << "  New x: ";  vXcurrent.leftshift (cout);
                cout << endl << endl;
            }

            continue;
        }


        //---- COMPUTE MULTIPLIERS TO SEE IF ANY IMPROVEMENT ON THE POINT
        //---- CAN BE MADE.  IF YES, THEN DROP ONE CONSTRAINT AND CONTINUE.
        int  nIndexWorstIneqMultiplier = 0;
        if (computeMultipliers_ (vC, vD, mActiveCon, nNumEqs, vXsol,
                                 nIndexWorstIneqMultiplier) == false)
        {
            bFoundSolution = false;
            break;
        }
        if (nIndexWorstIneqMultiplier == -1)
        {
            //---- MULTIPLIERS ARE FINE, EXIT WITH SUCCESS.
            bFoundSolution = true;
            break;
        }
        int  nNextActiveIneq = 0;
        for (int  i = 0; i < nNumIneqs; i++)
        {
            if (pbaIsActive[i] == true)
            {
                if (nNextActiveIneq == nIndexWorstIneqMultiplier)
                {
                    pbaIsActive[i] = false;
                    if (bDEBUG)
                    {
                        cout << "ActSet subproblem dropping inequality " << i
                             << endl;
                    }
                    break;
                }
                else
                    nNextActiveIneq++;
            }
        }

    }    //-- END ACTIVE SET LOOP
    if (nCurrentIter >= nMaxIters)
    {
        if (bDEBUG)
            cout << "ActSet subproblem halting due to interation limit" << endl;
        bFoundSolution = false;
    }

    delete[] pbaIsActive;
    delete[] pbaIsLowerBnd;

    return( bFoundSolution );
}


//----------------------------------------------------------------------
//  Private method calcUnconstrainedSolution_
//----------------------------------------------------------------------
void  SolveLinConstrProj::calcUnconstrainedSolution_
          (const Vector &  vC,
           const Vector &  vD,
                 Vector &  vXsol) const
{
    const double  dNEAR_ZERO = getMachineEpsilon() * getMachineEpsilon();

    for (int  i = 0; i < vD.size(); i++)
    {
        if (fabs (vD[i]) < dNEAR_ZERO)
            vXsol[i] = vC[i];
        else
            vXsol[i] = vC[i] / vD[i];
    }
    return;
}


//----------------------------------------------------------------------
//  Private method isIneqFeasible_
//----------------------------------------------------------------------
bool  SolveLinConstrProj::isIneqFeasible_ (const Vector &  vX,
                                           const Matrix &  mMatIneqs,
                                           const Vector &  vLoIneqs,
                                           const Vector &  vUpIneqs) const
{
    //---- THIS TOLERANCE IS AN ARBITRARY, SOMEWHAT LOOSE NUMBER,
    //---- THAT ASSUMES PROBLEMS ARE SCALED.
    const double  dTOLERANCE = sqrt (getMachineEpsilon());

    for (int  i = 0; i < mMatIneqs.getNrows(); i++)
    {
        double  dLHS = (mMatIneqs.getRow(i)).dot (vX);
        if (exists (vLoIneqs[i]) && ((vLoIneqs[i] - dLHS) > dTOLERANCE))
            return( false );
        if (exists (vUpIneqs[i]) && ((dLHS - vUpIneqs[i]) > dTOLERANCE))
            return( false );
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private method computeMultipliers_
//----------------------------------------------------------------------
bool  SolveLinConstrProj::computeMultipliers_ (const Vector &  vC,
                                               const Vector &  vD,
                                               const Matrix &  mMatCons,
                                               const int       nNumEqs,
                                               const Vector &  vX,
                                                     int    &  nConIndex) const
{
    if (nNumEqs == mMatCons.getNrows())
    {
        //---- SKIP CALCULATION BECAUSE THERE ARE NO INEQUALITIES.
        nConIndex = -1;
        return( true );
    }

    Matrix  mConsTrans;
    mConsTrans.transpose (mMatCons);
    Vector  vMults (mConsTrans.getNcols());

    //---- GRADIENT OF THE OBJECTIVE = 2 D^T (D^T x - c).
    Vector  vGradObj (mConsTrans.getNrows());
    for (int  i = 0; i < vGradObj.size(); i++)
        vGradObj[i] = 2.0 * vD[i] * (vD[i] * vX[i] - vC[i]);

    if (bDEBUG)
    {
        cout << "ActSet calling LS (LAPACK.dgelss) for multipliers with";
        mConsTrans.formattedPrint ("  ConsTrans ", cout);  cout << endl;
        cout << "  RHS: ";  vGradObj.leftshift (cout);  cout << endl;
    }
    if (mConsTrans.generalLS (vGradObj, vMults) == false)
    {
        cerr << "WARNING: Call to solve LS subproblem failed" << endl;
        return( false );
    }
    if (bDEBUG)
    {
        cout << "  LS returns: ";  vMults.leftshift (cout);  cout << endl;
    }

    //---- FIND THE MOST NEGATIVE INEQUALITY MULTIPLIER.
    nConIndex = -1;
    double  dMostNeg = - _dActiveTol;
    for (int  i = nNumEqs; i < mMatCons.getNrows(); i++)
    {
        if (vMults[i] < dMostNeg)
        {
            dMostNeg = vMults[i];
            nConIndex = i - nNumEqs;
        }
    }

    return( true );
}


}     //-- namespace HOPSPACK
