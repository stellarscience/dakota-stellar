#ifndef DRKERNEL_H_
#define DRKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/Utilities/VectorSlice.h"

#include <set>

namespace muq {
  namespace SamplingAlgorithms {

    /**
      @ingroup MCMCKernels
      @class DRKernel
      @brief An implementation of the delayed rejection kernel.
      @details
      This class provides an implementation the delayed rejection kernel described
      in "On Metropolis-Hastings algorithms with delayed rejection" by Antonietta Mira.
      This algorithm uses a sequence of proposal stages to increase the probability of
      accepting a move.  In the first stage, a normal Metropolis-Hastings step
      is attempted with a proposal \f$q_0\f$.  However, if the proposed point is
      rejected, a second proposed move is generated with a proposal \f$q_1\f$.
      This proposed point is accepted or rejected with an acceptance probability
      that has been adjusted to maintain detailed balance.  The process continues
      for subsequent stages until either (1) a proposal is accepted, at which point
      the proposal is immediately returned as the next state in the Markov chain or (2)
      the maximum number of delayed rejection stages is met, at which point the
      current state is repeated in the Markov chain.

      The delayed rejection implementation in MUQ allows the proposal sequence to
      be defined in two ways: either by specifying the proposals for each stage
      individually, or by specifying a single proposal that is then scaled at
      each stage.  The latter approach is more restrictive, but can yield
      slightly more efficient implementations of standard methods like DRAM.

      More precisely, a scaling \f$a_i\f$ of some base proposal density \f$q(x)\f$
      defines a new proposal density \f$q_{a_i}(y) = a_i q\left(a_iy\right)\f$, where the
      new variable \f$y\f$ is equal in distribution to a scaling of \f$x\f$, i.e.,
      \f$y=\frac{x}{a_i}\f$.  We assume the scale \f$a_i\f$ for stage \f$i\f$ is
      given by one of two functions of the stage index \f$i\f$.  Either the **linear
      scaling**
      \f[
      a_i = \frac{b}{i+1},
      \f]
      or the **power scaling**
      \f[
      a_i = \frac{b}{2^i}.
      \f]
      Note that in both cases, \f$i\in\{0,1,...,N-1\}\f$, where \f$N\f$ is the number
      of stages.


      <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Proposal"   | string | - | Either the name of another block in the ptree that defines the base proposal \f$q\f$ that will be scaled for each stage, or a comma separated list of blocks that define the proposals for each stage.  |
      "NumStages"   | integer       | -   | Number of stages to use if scalings of a single proposal are to be used.  If a list of proposals is specified, then this option does nothing. |
      "ScaleFunction"   | string        | "Power" |  Either "Linear" or "Power" to dictate the type of scaling function to use. |
      "Scale"       | double        | 2.0 |  Value of \f$b\f$ in scaling equation.  Value must be positive and will generally be greater than 1 to ensure the proposal shrinks with increasing stage \f$i\f$. |

     */
    class DRKernel : public TransitionKernel {
    public:

      DRKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem);

      DRKernel(boost::property_tree::ptree         const& pt,
               std::shared_ptr<AbstractSamplingProblem>   problem,
               std::vector<std::shared_ptr<MCMCProposal>> proposalsIn,
               std::vector<double>                        scales);

      virtual ~DRKernel() = default;

      /**
      Return a  vector of the MCMC proposals used in each stage.
      */
      virtual inline std::vector<std::shared_ptr<MCMCProposal>> Proposals() {return proposals;};

      virtual void PostStep(unsigned int const t,
                            std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t,
                                                               std::shared_ptr<SamplingState> prevState) override;

      /**
      Print the status of this kernel to std::cout.
      @param[in] prefix A string that should be appended to the output.  Typically
                        used for handling indentation in nested kernels.
      */
      virtual void PrintStatus(std::string prefix) const override;

      /**
      For each stage, we define the acceptance rate as the ratio of the number
      of proposal calls to the number of times the proposal from that stage was
      accepted.  Note that later stages will be called less often than earlier
      stages.
      @return An Eigen::VectorXd containing acceptance rates for each stage.
      */
      virtual inline Eigen::VectorXd AcceptanceRates() const{return numProposalAccepts.cast<double>().array()/numProposalCalls.cast<double>().array();};


      /** Generates a sample of a stage proposal at the point x.  This function
          handles all necessary scaling.
      */
      std::shared_ptr<SamplingState> SampleProposal(unsigned int                          stage,
                                                    std::shared_ptr<SamplingState> const& state) const;

      /** Evaluates the log density of a stage proposal at point y of the proposal
          located at point x.  This function takes care of all necessary scaling.
      */
      double EvaluateProposal(unsigned int                          stage,
                              std::shared_ptr<SamplingState> const& x,
                              std::shared_ptr<SamplingState> const& y) const;


      /** Returns a vector with the scaling used for each proposal stage. */
      std::vector<double> GetScales() const{return propScales;};

    protected:

      /** Extracts information from the property tree and creates MCMC proposals.
      */
      static std::vector<std::shared_ptr<MCMCProposal>> CreateProposals(boost::property_tree::ptree const& pt,
                                                                        std::shared_ptr<AbstractSamplingProblem> const& problem);

      /** In some cases, the proposal at each stage is just a scaled version of
      some base proposal distribution.  In these cases, we can store and adapt a
      a single proposal if we explicitly store the scales.  This function reads
      a property tree and creates a vector of proposal scalings.  See the class
      level documentation for the necessary ptree format.
      */
      static std::vector<double> CreateScales(boost::property_tree::ptree const& pt);

      // A vector containing one proposal for each stage
      std::vector<std::shared_ptr<MCMCProposal>> proposals;

      /* A set containing unique proposals (the proposals vector might have
        duplicates but with different scales, e.g., for DRAM).  This unique set
        is used to make sure that adaptive proposals that used in multiple
        stages are updated correctly.
      */
      std::set<std::shared_ptr<MCMCProposal>> uniqueProps;

      // Scales for each stage
      std::vector<double> propScales;
      bool isScaled = false;

      // store how many times each proposal has been called
      Eigen::VectorXi numProposalCalls;

      // store how many times each proposal has been accepted
      Eigen::VectorXi numProposalAccepts;



      template<typename VecType1, typename VecType2>
      double Alpha(VecType1& likelies, VecType2& proposed_points) const
      {
        // create subcontainers of input variables to use in recursion

        int stage = likelies.size() - 1;

        double a1 =  1.0;
        double a2 = 1.0;

        for (int k = 1; k < stage; ++k) {
          auto slice1 = muq::Utilities::GetSlice(likelies, stage, stage - k-1,-1);
          auto slice2 = muq::Utilities::GetSlice(proposed_points, stage, stage - k-1, -1);
          a2 *= (1.0 - Alpha(slice1, slice2));

          slice1 = muq::Utilities::GetSlice(likelies, 0, k+1);
          slice2 = muq::Utilities::GetSlice(proposed_points, 0, k+1);
          a1 *= (1.0 - Alpha(slice1, slice2)); // forward

          // if a2==0, there is no chance of getting accepted, so just give up and try again
          if (a2 == 0)
            return 0.0;
        }

        double q = 0.0;

        //now put in qs

        for (int k = 1; k <= stage; ++k) {
          q +=  QFun(muq::Utilities::GetSlice(proposed_points, stage, stage - k-1, -1));
          q -= QFun(muq::Utilities::GetSlice(proposed_points, 0, k+1)); // forward probabilities
        }

        return std::min<double>(1.0, exp(likelies[stage] - likelies[0] + q) * a2 / a1);
      }

      template<typename VecType>
      double QFun(VecType const& proposed_points) const
      {
        //use the proposal that generated this distance move
        int stage = proposed_points.size() - 1;
        return EvaluateProposal(stage - 1, proposed_points[0], proposed_points[stage]);
      }


    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
