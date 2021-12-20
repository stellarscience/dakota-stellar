#ifndef ODE_H_
#define ODE_H_

#include <boost/property_tree/ptree.hpp>

#include <cvodes/cvodes.h>               /* access to CVODE                 */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sunnonlinsol/sunnonlinsol_newton.h>     /* access to the newton SUNNonlinearSolver      */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> /* access to the fixed point SUNNonlinearSolver */


#include "MUQ/Modeling/ModPiece.h"

namespace muq {
namespace Modeling {


  /**
  @brief A class for integrating ODEs with CVODES
  @details

  <B>Configuration Parameters:</B>
  | Parameter Key | Type | Default Value | Description |
  | ------------- | ------------- | ------------- | ------------- |
  | "NonlinearSolver"  | ptree | - | A child tree of options for the nonlinear solver. See ODEPiece::NonlinearSolverOptions. |
  | "LinearSolver" | ptree | - | A child tree of options for the linear solver. See ODEPiece::LinearSolverOptions. |
  | "Integrator" | ptree | - | A child tree of options for the nonlinear solver See ODEPiece::IntegratorOptions. |
  */
  class ODE : public ModPiece {
  public:


    ODE(std::shared_ptr<ModPiece>   const& rhsIn,
        Eigen::VectorXd             const& evalTimesIn,
        boost::property_tree::ptree        opts = boost::property_tree::ptree()) : ODE(rhsIn, evalTimesIn, true, opts){};

    ODE(std::shared_ptr<ModPiece>   const& rhsIn,
        Eigen::VectorXd             const& evalTimesIn,
        bool                               isAuto,
        boost::property_tree::ptree        opts = boost::property_tree::ptree());


    virtual ~ODE() = default;

    /// Does the RHS explicitly depend on time?
    const bool isAutonomous;

    /// The time when this simulation starts -- either set in the opts ptree or set to the minimum of evalTimes
    const double startTime;

    protected:

      /**
        @brief Integrates the ODE forward in time
        @details
       */
      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void GradientImpl(unsigned int                       outWrt,
                                unsigned int                       inWrt,
                                ref_vector<Eigen::VectorXd> const& inputs,
                                Eigen::VectorXd const&             sens) override;

      virtual void JacobianImpl(unsigned int outWrt,
                                unsigned int inWrt,
                                ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual void ApplyJacobianImpl(unsigned int outWrt,
                                      unsigned int inWrt,
                                      ref_vector<Eigen::VectorXd> const& inputs,
                                      Eigen::VectorXd const& vec) override;

      /** @class ODE::NonlinearSolverOptions
        Contains options for the nonlinear solver used by CVODES.

        <B>Configuration Parameters:</B>
        | Parameter Key | Type | Default Value | Description |
        | ------------- | ------------- | ------------- | ------------- |
        | "Method"  | string | "Newton" | (Optional) Defines the Sundials nonlinear solver type.  Either "Newton" or "Iter". |
        | "MaximumIterations" | int | -1 | (Optional) Defines the maximum number of iterations allowed in the nonlienar solve.  Ignored if negative. |

        */
      struct NonlinearSolverOptions{
        NonlinearSolverOptions() = default;
        NonlinearSolverOptions(boost::property_tree::ptree const& opts);

        void SetOptions(boost::property_tree::ptree const& opts);

        std::string method; // e.g., Iter, Newton
        int maxIters;
      };

      /** @class ODE::LinearSolverOptions
      Contains options for the linear solver used in the nonlinear solver.

      <B>Configuration Parameters:</B>
      | Parameter Key | Type | Default Value | Description |
      | ------------- | ------------- | ------------- | ------------- |
      | "Method"  | string | "Dense" | (Optional) Defines the Sundials linear solver type. Options are "Dense", "SPGMR", "SPFGMR", "SPBCGS", "SPTFQMR", "PCG". |
      | "MaximumIterations" | int | -1 | (Optional) Maximum number of iterations allowed in iterative solvers. |
      */
      struct LinearSolverOptions{
        LinearSolverOptions() = default;
        LinearSolverOptions(boost::property_tree::ptree const& opts);
        void SetOptions(boost::property_tree::ptree const& opts);

        std::string method; // e.g., Dense, SPGMR, SPFGMR, SPBCGS, SPTFQMR, PCG
        int maxIters; // maximum number of iterations allowed for iterative solvers

      };

      /** @class ODE::IntegratorOptions
      Contains options for the time integrator.

      <B>Configuration Parameters:</B>
      | Parameter Key | Type | Default Value | Description |
      | ------------- | ------------- | ------------- | ------------- |
      | "Method"  | string | "BDF" | (Optional) Defines the CVODES time integration algorithm.  Either "Adams" or "BDF".  |
      | "RelativeTolerance"  | double | 1e-5 | (Optional) Relative tolerance in ODE integration.  |
      | "AbsoluteTolerance"  | double | 1e-5 | (Optional) Absolute tolerance in ODE integration.  |
      | "MaximumStepSize"  | double | -1 | (Optional) Maximum stepsize.  Ignored if negative.  |
      | "MaximumSteps"  | int | -1 | (Optional) Maximum number of steps.  Ignored if negative.  |
      | "CheckPointSteps" | int | 100 | (Optional) Number of steps between checkpoints. |
      */
      struct IntegratorOptions{
        IntegratorOptions() = default;
        IntegratorOptions(boost::property_tree::ptree const& opts);

        void SetOptions(boost::property_tree::ptree const& opts);

        std::string method; // e.g., Adams or BDF
        double reltol; // relative tolerance
        double abstol; // absolute tolerance
        double maxStepSize; // maximum allowed step inputSize
        int maxNumSteps; // maximum allowed steps
        int checkPtGap; // gap between check points
      };

      /** Holds static functions needed for CVODES to interface with ModPiece RHS.
      */
      struct Interface{

        static int RHS(realtype time, N_Vector state, N_Vector deriv, void *user_data);
        static int RHSJacobian(realtype time, N_Vector state, N_Vector rhs, SUNMatrix jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
        static int RHSJacobianAction(N_Vector v, N_Vector Jv, realtype time, N_Vector state, N_Vector rhs, void *user_data, N_Vector tmp);

        /** Evaluates the RHS of the sensitivity equations, which are given by
        $$
        \dot{s}_i = \frac{\partial f}{\partial y}s_i + \frac{\partial f}{\partial p_i}
        $$
        See page 5 of the CVODES TOMS paper (https://computing.llnl.gov/sites/default/files/public/toms_cvodes_with_covers.pdf).
        */
        static int SensitivityRHS(int Ns, realtype time, N_Vector y, N_Vector ydot, N_Vector *ys, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2);

        /** Evaluates the RHS of the adjoint equation, which is given by
        $$
        \dot{\lambda} = -\left(\frac{\partial f}{\partial y}\right)^T \lambda - \left(\frac{\partial g}{\partial y}\right)
        $$
        For details, see Equation (15) in Section 2.3 of https://computing.llnl.gov/sites/default/files/public/toms_cvodes_with_covers.pdf
        */
        static int AdjointRHS(double time, N_Vector state, N_Vector lambda, N_Vector deriv, void *user_data);
        static int AdjointJacobian(double time, N_Vector state, N_Vector lambda, N_Vector rhs, SUNMatrix jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

        static int AdjointJacobianAction(N_Vector target, N_Vector output, double time, N_Vector state, N_Vector lambda, N_Vector adjRhs, void *user_data, N_Vector tmp);

        /** Function for accumulating sensitivity wrt parameters during adjoint solve.
        See Equation (16) of https://computing.llnl.gov/sites/default/files/public/toms_cvodes_with_covers.pdf
        */
        static int AdjointQuad(realtype time, N_Vector state, N_Vector lambda, N_Vector quadRhs, void *user_data);
      };

      NonlinearSolverOptions nonlinOpts;
      LinearSolverOptions linOpts;
      IntegratorOptions intOpts;

      // Store the times when where we want to return the solution of the ODE
      const Eigen::VectorXd evalTimes;

      // Store a ModPiece that evaluates the RHS function
      std::shared_ptr<ModPiece>  rhs;

    private:

      /** Creates a linear solver and, if necessary, a dense matrix for storing
          the Jacobian.  If the nonlinear solver is fixed point iteration, then
          no linear solver is necessary and this function returns an empty
          linear solver.  If no Jacobian is needed, either because a fixed point
          nonlinear solver is specified or because an iterative linear solver is
          specified, then NULL is returned. The options read by the constructor
          and stored in `linOpts` are used to construct the solver.
      */
      std::pair<SUNLinearSolver, SUNMatrix> CreateLinearSolver(void *cvode_mem, N_Vector const& state);

      SUNNonlinearSolver CreateNonlinearSolver(void *cvode_mem, N_Vector const& state);

      std::pair<SUNLinearSolver, SUNMatrix> CreateAdjointLinearSolver(void *cvode_mem, N_Vector const& state, int indexB);

      SUNNonlinearSolver CreateAdjointNonlinearSolver(void *cvode_mem, N_Vector const& state, int indexB);


      static Eigen::VectorXi GetInputSizes(std::shared_ptr<ModPiece> const& rhsIn,
                                           bool                             isAuto);

      static Eigen::VectorXi GetOutputSizes(std::shared_ptr<ModPiece> const& rhsIn,
                                            Eigen::VectorXd           const& evalTimes);

      static double GetStartTime(Eigen::VectorXd             const& evalTimesIn,
                                 boost::property_tree::ptree const& opts);

      /// Checks the return value produced by CVODES.  Throws and exception if error occured.
      static void CheckFlag(void *returnvalue, const char *funcname, int opt);
    };

} // namespace Modeling
} // namespace muq

#endif // #ifndef ODE_H_
