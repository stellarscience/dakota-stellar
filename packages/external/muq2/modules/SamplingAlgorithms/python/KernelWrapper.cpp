#include "AllClassWrappers.h"

#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/DRKernel.h"
#include "MUQ/SamplingAlgorithms/DILIKernel.h"
#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;
namespace py = pybind11;

void PythonBindings::KernelWrapper(py::module &m) {
  py::class_<TransitionKernel, std::shared_ptr<TransitionKernel>> transKern(m, "TransitionKernel");
  transKern
    .def_static("Construct", [](py::dict d, std::shared_ptr<AbstractSamplingProblem> problem)->std::shared_ptr<TransitionKernel>{return TransitionKernel::Construct(ConvertDictToPtree(d), problem);})
    .def("PreStep", &TransitionKernel::PreStep)
    .def("PostStep", &TransitionKernel::PostStep)
    .def("Step", &TransitionKernel::Step)
    .def_readonly("blockInd", &TransitionKernel::blockInd);


  py::class_<MHKernel, TransitionKernel, std::shared_ptr<MHKernel>> mhKern(m, "MHKernel");
  mhKern
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> problem) {return new MHKernel(ConvertDictToPtree(d), problem);}))
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn) {return new MHKernel(ConvertDictToPtree(d), problem, proposalIn);}))
    .def("Proposal", &MHKernel::Proposal)
    .def("PostStep", &MHKernel::PostStep)
    .def("Step", &MHKernel::Step)
    .def("AcceptanceRate", &MHKernel::AcceptanceRate);

  py::class_<DRKernel, TransitionKernel, std::shared_ptr<DRKernel>>(m, "DRKernel")
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> problem) {return new DRKernel(ConvertDictToPtree(d), problem);}))
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> problem, std::vector<std::shared_ptr<MCMCProposal>> proposals, std::vector<double> scales) {return new DRKernel(ConvertDictToPtree(d), problem, proposals, scales);}))
    .def("PostStep", &DRKernel::PostStep)
    .def("Step", &DRKernel::Step)
    .def("Proposals",&DRKernel::Proposals)
    .def("AcceptanceRates",&DRKernel::AcceptanceRates)
    .def("SampleProposal", &DRKernel::SampleProposal)
    .def("EvaluateProposal", &DRKernel::EvaluateProposal);

  py::class_<DILIKernel, TransitionKernel, std::shared_ptr<DILIKernel>>(m,"DILIKernel")
    .def(py::init( [](py::dict d,
                      std::shared_ptr<AbstractSamplingProblem> problem){return new DILIKernel(ConvertDictToPtree(d),problem);}))
    .def(py::init( [](py::dict d,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      std::shared_ptr<muq::Modeling::GaussianBase> const& prior,
                      std::shared_ptr<muq::Modeling::ModPiece> const& noiseMod,
                      std::shared_ptr<muq::Modeling::ModPiece> const& likelihood){return new DILIKernel(ConvertDictToPtree(d),problem,prior,noiseMod, likelihood);}))
    .def(py::init( [](py::dict d,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      std::shared_ptr<muq::Modeling::GaussianBase> const& prior,
                      std::shared_ptr<muq::Modeling::ModPiece> const& likelihood){return new DILIKernel(ConvertDictToPtree(d),problem,prior,likelihood);}))
    .def(py::init( [](py::dict d,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      Eigen::VectorXd const& vals,
                      Eigen::MatrixXd const& vecs){return new DILIKernel(ConvertDictToPtree(d),problem, vals, vecs);}))
    .def(py::init( [](py::dict d,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      std::shared_ptr<muq::Modeling::GaussianBase> const& prior,
                      std::shared_ptr<muq::Modeling::ModPiece> const& noiseMod,
                      std::shared_ptr<muq::Modeling::ModPiece> const& likelihood,
                      Eigen::VectorXd const& vals,
                      Eigen::MatrixXd const& vecs){return new DILIKernel(ConvertDictToPtree(d),problem,prior,noiseMod, likelihood, vals, vecs);}))
    .def(py::init( [](py::dict d,
                      std::shared_ptr<AbstractSamplingProblem> problem,
                      std::shared_ptr<muq::Modeling::GaussianBase> const& prior,
                      std::shared_ptr<muq::Modeling::ModPiece> const& likelihood,
                      Eigen::VectorXd const& vals,
                      Eigen::MatrixXd const& vecs){return new DILIKernel(ConvertDictToPtree(d),problem,prior,likelihood,vals,vecs);}))
    .def("LISKernel", &DILIKernel::LISKernel)
    .def("CSKernel", &DILIKernel::CSKernel)
    .def_static("ExtractLikelihood", &DILIKernel::ExtractLikelihood)
    .def_static("ExtractNoiseModel", &DILIKernel::ExtractNoiseModel)
    .def_static("ExtractForwardModel", &DILIKernel::ExtractForwardModel)
    .def_static("CreateLikelihood", &DILIKernel::CreateLikelihood)
    .def("LISVecs", &DILIKernel::LISVecs)
    .def("LISVals", &DILIKernel::LISVals);
}
