#ifndef REMOTEMIPROPOSAL_H_
#define REMOTEMIPROPOSAL_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif


#include <boost/property_tree/ptree.hpp>
#include "MUQ/SamplingAlgorithms/Phonebook.h"

#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

		/**
		 * @brief Proposal retrieving samples from other ranks.
		 * @details This is particularly relevant for parallel MIMCMC type methods.
		 * It allows drawing proposals from coarser chains being computed on other
		 * processes. Which process to draw a sample from is determined via phonebook.
		 */
		class RemoteMIProposal : public MCMCProposal {
		public:
			RemoteMIProposal (pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> prob, std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<MultiIndex> remoteIndex, std::shared_ptr<MultiIndex> sourceIndex, std::shared_ptr<PhonebookClient> phonebookClient)
				: MCMCProposal(pt,prob),
				//subsampling(pt.get<int>("Subsampling")),
				comm(comm),
		    remoteIndex(remoteIndex),
		    sourceIndex(sourceIndex),
				phonebookClient(phonebookClient)
			{
			}

			std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) {


				int remoteRank = phonebookClient->Query(remoteIndex, sourceIndex, true);

				comm->Send(ControlFlag::SAMPLE, remoteRank, ControlTag);
				Eigen::VectorXd remoteState = comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag);
				double remoteLogTarget = comm->Recv<double>(remoteRank, ControlTag);
				Eigen::VectorXd remoteQOI = comm->Recv<Eigen::VectorXd>(remoteRank, ControlTag);

				auto proposal = std::make_shared<SamplingState>(remoteState); // TODO: Support non-QOI samples!!
				proposal->meta["QOI"] = std::make_shared<SamplingState>(remoteQOI);
				proposal->meta["LogTarget"] = remoteLogTarget;

				return proposal;
			}

			double LogDensity(std::shared_ptr<SamplingState> const& currState,
                        std::shared_ptr<SamplingState> const& propState) {
				return 0;
			}

		private:
			int sampleID = 0;
			int sampleWeight = 0;
			//const int subsampling;
			std::shared_ptr<parcer::Communicator> comm;
		  std::shared_ptr<MultiIndex> remoteIndex;
		  std::shared_ptr<MultiIndex> sourceIndex;
			std::shared_ptr<PhonebookClient> phonebookClient;
		};
	}
}

#endif

#endif
