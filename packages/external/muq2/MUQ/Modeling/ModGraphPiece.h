#ifndef MODGRAPHPIECE_H_
#define MODGRAPHPIECE_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

#include <map>

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Modeling/NodeNameFinder.h"
#include "MUQ/Modeling/WorkGraph.h"

namespace muq {
  namespace Modeling {

    /// A filtered graph that only has nodes downstream of a specified input
    typedef boost::filtered_graph<Graph, DependentEdgePredicate, DependentPredicate> FilteredGraph;

    /// A muq::Modeling::ModPiece created from a muq::Modeling::WorkGraph
    class ModGraphPiece : public ModPiece {
    public:

      /// Construct a muq::Modeling::ModGraphPiece
      /**
      	 Typically muq::Modeling::ModGraphPiece's are constructed by calling muq::Modeling::WorkGraph::CreateModPiece
      	 @param[in] graph From inputs to the output-of-interest
      	 @param[in] constantPieces Pointers to the muq::Modeling::ConstantVector's that hold the graph's inputs
      	 @param[in] inputNames The names of each node in the graph corresponding to the constantPieces
      	 @param[in] outputNode The muq::Modeling::ModPiece that we ultimately want to evaluate
       */
      ModGraphPiece(std::shared_ptr<WorkGraph>                           graph,
                    std::vector<std::shared_ptr<ConstantVector> > const& constantPieces,
                    std::vector<std::string>                      const& inputNames,
                    std::shared_ptr<ModPiece>                            outputNode);

      /// Default destructor
      virtual ~ModGraphPiece() = default;

      std::shared_ptr<WorkGraph> GetGraph(){return wgraph;};

      std::vector<std::shared_ptr<ConstantVector> > GetConstantPieces(){return constantPieces;};

      /** @brief Returns another ModGraphPiece containing only part of the graph
                 making up this model.
          @details Consider the graph for a posterior density.  You might have 5
                   nodes in the graph: "Parameters", "Prior", "Likelihood", "Forward Model",
                   and "Posterior".   This function returns a model using only
                   a subset of these nodes.  For example, calling GetSubModel with
                   the node name "Likelihood" will return a new ModGraphPiece containing
                   three nodes: "Parameters", "Forward Model", and "Likelihood".
                   The resulting model will only evaluate the likelihood.
          @param[in] nodeName The name of a node in the graph that will serve as
                              the new "output" node.
          @return A ModGraphPiece with the "nodeName" node as output.
      */
      std::shared_ptr<ModGraphPiece> GetSubModel(std::string const& nodeName) const;

      /**
        @brief Matches inputs with another ModGraphPiece by looking at node names and edges.

        @details Assume this ModGraphPiece has three inputs call \f$x_1\f$, \f$x_2\f$, and \f$y\f$,
        where \f$x_i\f$ denotes the \f$i^{th}\f$ input of a node \f$x\f$.  Let there
        be another ModGraphPiece with inputs \f$x_2\f$, \f$y\f$, and \f$z\f$.  This function
        tries to match the inputs between the two ModGraphPieces by looking at the node
        names and input indices.  It returns the input indices of *this that correspond to
        the inputs of the other piece, with the value "-1" reserved for nonoverlapping inputs.
        For the example above, the result would be \f$[1,2,-1]\f$.

        @params[in] otherPiece The other ModGraphPiece that we want to match with the inputs to this piece.
        @returns A std::vector containing indices into the inputs of this ModGraphPiece that correspond with the inputs of the otherPiece.  The length of this vector is equal to the number of inputs of otherPiece.
      */
      std::vector<int> MatchInputs(std::shared_ptr<ModGraphPiece> otherPiece) const;

      /**
      Returns a ModGraphPiece that, when evaluated, returns the gradient of this
      ModPiece.  Note that the returned ModPiece will have an additional input
      for the sensitivity vector used in the Gradient call.
      */
      std::shared_ptr<ModGraphPiece> GradientGraph(unsigned int                const  outputDimWrt,
                                                   unsigned int                const  inputDimWrt);

      /**
        Returns a ModGraphPiece that, when evaluated, returns the action of the Jacobian of this
        ModPiece on a vector.  Note that the returned ModPiece will have an
        additional input for the vector that we want to apply the Jacobian to.
      */
      std::shared_ptr<ModGraphPiece> JacobianGraph(unsigned int                const  outputDimWrt,
                                                   unsigned int                const  inputDimWrt);

      /** Returns the name of the output node in the graph.*/
      std::string GetOutputName() const{return wgraph->GetName(outputPiece);};

      /** Returns the output ModPiece. */
      std::shared_ptr<ModPiece> GetOutputPiece() const{return outputPiece;};

    private:

      // Indices are [outWrt][inWrt]
      std::map<std::pair<unsigned int, unsigned int>, std::shared_ptr<ModGraphPiece>> gradientPieces;
      std::map<std::pair<unsigned int, unsigned int>, std::shared_ptr<ModGraphPiece>> jacobianPieces;
      std::map<std::tuple<unsigned int, unsigned int, unsigned int>, std::shared_ptr<ModGraphPiece>> hessianPieces;

      /** Returns the input index of this ModGraphPiece that corresponds to the
          the specified placeholder ModPiece.  If the ModPiece is not an input,
          than -1 is returned.
      */
      int GetInputIndex(std::shared_ptr<WorkPiece> const& piece) const;

      static Eigen::VectorXi ConstructInputSizes(std::vector<std::shared_ptr<ConstantVector> > const& constantPiecesIn);

      /// Evaluate each muq::Modeling::WorkPiece in the graph
      /**
	       @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
      */
      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Compute the action of the Jacobian transpose for this muq::Modeling::WorkGraphPiece using the chain rule
      /**
      @param[in] wrtIn We are taking the Jacobian with respect to this input
      @param[in] wrtOut We are taking the Jacobian of this output
      @param[in] vec We are applying the Jacobian transpose to this object
      @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
       */
      virtual void GradientImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd             const& sensitivity) override;

      /// Compute the Jacobian for this muq::Modeling::WorkGraphPiece using the chain rule
      /**
         @param[in] wrtIn We are taking the Jacobian with respect to this input
         @param[in] wrtOut We are taking the Jacobian of this output
         @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
      */
      virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input) override;

      /// Compute the action of the Jacobian for this muq::Modeling::WorkGraphPiece using the chain rule
      /**
         @param[in] wrtIn We are taking the Jacobian with respect to this input
         @param[in] wrtOut We are taking the Jacobian of this output
         @param[in] vec We are applying the Jacobian to this object
         @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
       */
      virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                     unsigned int                const  inputDimWrt,
                                     ref_vector<Eigen::VectorXd> const& input,
                                     Eigen::VectorXd             const& vec) override;


     virtual void ApplyHessianImpl(unsigned int                const  outWrt,
                                   unsigned int                const  inWrt1,
                                   unsigned int                const  inWrt2,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& sens,
                                   Eigen::VectorXd             const& vec) override;


      /// Get the required outputs for a node in one of the filtered graphs
      /**
	       @param[in] node We want the outputs of this node
	       @param[in] wrtIn The input whose downstream nodes we care about
	       @param[in] wrtOut The output we are ultimately trying to differentiate wrt
	       @return The output nodes --- tuple: the output WorkPiece ID, the output number, and the input number
      */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > RequiredOutputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn, unsigned int wrtOut) const;

      /// Get the required inputs for a node in one of the filtered graphs
      /**
	       @param[in] node We want the inputs of this node
	       @param[in] wrtIn The input whose downstream nodes we care about
	       @return The input nodes --- tuple: the input WorkPiece ID, the output number, and the input number
       */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > RequiredInputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn) const;

      /// Fill the map from each node's muq::Modeling::WorkPiece::ID to its outputs
      void FillOutputMap();

      /// Set the inputs
      /**
	       Set the inputs in each the muq::Modeling::ConstantPiece.
      */
      void SetInputs(ref_vector<Eigen::VectorXd> const& inputs);

      /// Get the inputs from muq::Modeling::WorkGraphPiece::valMap to a specified node in the graph
      /**
	       @param[in] id The ID of the node of interest
	       @return A reference vector of inputs to that node
       */
      ref_vector<Eigen::VectorXd> GetNodeInputs(boost::graph_traits<Graph>::vertex_descriptor node) const;

      /// Get a the input nodes for a node
      /**
	       @param[in] node We want the input nodes for this node
	       @return A map from the input node's ID to the input/output number
      */
      std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > InputNodes(boost::graph_traits<Graph>::vertex_descriptor const& node) const;

      /// Run order computed during construction (input->output order)
      std::deque<boost::graph_traits<Graph>::vertex_descriptor> runOrder;

      // Run order for computing the derivatives of this muq::Modeling::WorkGraphPiece
      /**
	       Like muq::Modeling::WorkGraphPiece::runOrder, but specific to which input node is used (also in output->input order)
      */
      std::vector<std::deque<boost::graph_traits<Graph>::vertex_descriptor> > adjointRunOrders;

      /// The WorkGraph associated with this WorkGraphPiece
      std::shared_ptr<WorkGraph> wgraph;

      std::vector<std::shared_ptr<FilteredGraph> > filtered_graphs;

      /// A the map from each node's muq::Modeling::WorkPiece::ID to its outputs
      std::unordered_map<unsigned int, ref_vector<Eigen::VectorXd> > valMap;

      /// The ID of the WorkPiece corresponding to the output node
      unsigned int outputID;

      std::shared_ptr<ModPiece> outputPiece;

      /// The muq::Modeling::ConstantVector's that store the inputs
      std::vector<std::shared_ptr<ConstantVector> > constantPieces;

    };
  } // namespace Modeling
} // namespace muq

#endif
