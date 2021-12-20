#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/GradientPiece.h"
#include "MUQ/Modeling/JacobianPiece.h"
#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Modeling/SumPiece.h"

#include <boost/range/adaptor/reversed.hpp>
#include <boost/graph/topological_sort.hpp>

#include <Eigen/Core>

using namespace muq::Modeling;
using namespace muq::Utilities;



ModGraphPiece::ModGraphPiece(std::shared_ptr<WorkGraph>                           graph,
                             std::vector<std::shared_ptr<ConstantVector> > const& constantPiecesIn,
                             std::vector<std::string>                      const& inputNames,
                             std::shared_ptr<ModPiece>                            outputPieceIn) : ModPiece(ConstructInputSizes(constantPiecesIn),
                                                                                                          outputPieceIn->outputSizes),
                                                                                                 wgraph(graph),
                                                                                                 outputID(outputPieceIn->ID()),
                                                                                                 outputPiece(outputPieceIn),
                                                                                                 constantPieces(constantPiecesIn)  {

  // build the run order
  assert(graph);
  boost::topological_sort(wgraph->graph, std::front_inserter(runOrder));

  // each input only needs to loop over its downstream nodes when computing derivatives
  adjointRunOrders.resize(inputSizes.size());
  filtered_graphs.resize(inputSizes.size());

  // compute a run order for each of the inputs so we only have to loop over their downstream nodes
  assert(numInputs==inputNames.size());

  for( unsigned int i=0; i<numInputs; ++i ) { // loop through the inputs

    // get iterators to the begining and end of the graph
    boost::graph_traits<Graph>::vertex_iterator v, v_end;
    boost::tie(v, v_end) = vertices(wgraph->graph);

    // determine the downstream nodes of this input
    DependentPredicate nFilt(*std::find_if(v, v_end, NodeNameFinder(inputNames[i], wgraph->graph)), wgraph->graph);
    DependentEdgePredicate eFilt(nFilt, wgraph->graph);

    // filter the graph, we only care about downstream nodes of this input
    filtered_graphs[i] = std::make_shared<boost::filtered_graph<Graph, DependentEdgePredicate, DependentPredicate> >(wgraph->graph, eFilt, nFilt);

    // build specialized run order for each input dimension
    boost::topological_sort(*filtered_graphs[i], std::back_inserter(adjointRunOrders[i]));
  }
}

int ModGraphPiece::GetInputIndex(std::shared_ptr<WorkPiece> const& piece) const
{
  for(int i = 0; i<constantPieces.size(); ++i){
    if(piece == constantPieces.at(i))
      return i;
  }

  return -1;
}

std::shared_ptr<ModGraphPiece> ModGraphPiece::GradientGraph(unsigned int                const  outputDimWrt,
                                                            unsigned int                const  inputDimWrt)
{
  assert(outputDimWrt<outputSizes.size());
  assert(inputDimWrt<inputSizes.size());

  WorkGraph gradGraph;
  auto& filtGraph = *filtered_graphs[inputDimWrt];

  std::vector<std::string> inputNames(constantPieces.size()+1);
  std::vector<boost::graph_traits<Graph>::vertex_descriptor> inputNodes(constantPieces.size());

  std::string outputName = filtGraph[adjointRunOrders[inputDimWrt][0]]->name;
  inputNames.at(constantPieces.size()) = outputName + "_Sensitivity";

  bool needsFinalSum = false;

  // Add the forward model
  for(auto node = runOrder.begin(); node != runOrder.end()-1; ++node){
    auto piece = wgraph->graph[*node]->piece;

    int placeholderInd = GetInputIndex(piece);
    if(placeholderInd>=0){
      auto pieceMod = std::dynamic_pointer_cast<ModPiece>(piece);
      assert(pieceMod);
      gradGraph.AddNode(std::make_shared<IdentityOperator>(pieceMod->outputSizes(0)), wgraph->graph[*node]->name);
      inputNames.at(placeholderInd) = wgraph->graph[*node]->name;
      inputNodes.at(placeholderInd) = *node;
    }else{
      gradGraph.AddNode(piece, wgraph->graph[*node]->name);

      // Add the input edges
      for(auto es = in_edges(*node, wgraph->graph); es.first!=es.second; ++es.first){
        auto source = boost::source(*es.first, wgraph->graph);
        gradGraph.AddEdge(wgraph->graph[source]->name, wgraph->graph[*es.first]->outputDim, wgraph->graph[*node]->name, wgraph->graph[*es.first]->inputDim);
      }
    }
  }

  // Add a node in the graph to hold the sensitivity input
  gradGraph.AddNode(std::make_shared<IdentityOperator>(outputSizes(outputDimWrt)),
                    outputName+"_Sensitivity");


  // Add the adjoint components
  for(auto node : adjointRunOrders[inputDimWrt]){
    std::string baseName = wgraph->graph[node]->name;

    auto piece = std::dynamic_pointer_cast<ModPiece>(filtGraph[node]->piece);
    assert(piece);

    // outEdges stores information about where information for each output of this piece goes
    // outEdges[node output index][edge index]
    std::vector<std::vector<std::pair<boost::graph_traits<boost::filtered_graph<Graph, DependentEdgePredicate, DependentPredicate>>::vertex_descriptor, unsigned int>>> outEdges(piece->outputSizes.size()); // Pairs are node name and input id.
    for(auto eout = out_edges(node, filtGraph); eout.first!=eout.second; ++eout.first){
      auto target = boost::target(*eout.first, filtGraph);
      outEdges.at( filtGraph[*eout.first]->outputDim ).push_back(std::make_pair(target, filtGraph[*eout.first]->inputDim ));
    }

    // Figure out if the output of the node is used by more than one downstream node, which means we'll have to add up the gradients wrt this nodes output
    std::vector<unsigned int> edgeCounts(outEdges.size());
    for(unsigned int outInd=0; outInd<outEdges.size(); ++outInd){

      // Get the total number of gradient pieces that might be flowing into this sum
      edgeCounts[outInd] = 0;
      for(unsigned int i = 0; i<outEdges.at(outInd).size(); ++i){
        auto mod = std::dynamic_pointer_cast<ModPiece>(wgraph->graph[outEdges.at(outInd).at(i).first]->piece);
        assert(mod);
        edgeCounts[outInd] += mod->outputSizes.size();
      }

      // If there is more than one downstream node, we need to add up their sensitivities to this node
      if(edgeCounts[outInd]>1){

        std::stringstream sumName;
        sumName << baseName << "[" << outInd << "]_SensitivitySum";
        gradGraph.AddNode(std::make_shared<SumPiece>(piece->outputSizes(outInd), edgeCounts[outInd]), sumName.str());

        if(baseName==inputNames.at(inputDimWrt)){
          needsFinalSum = true;
        }
      }
    } // end of ouput index loop for counting edges

    // for each combination of incominging and outgoing edges for the forward node in the filtered graph, add a Gradient piece.
    for( auto e =in_edges(node, filtGraph); e.first!=e.second; ++e.first ) {
      auto source = boost::source(*e.first,filtGraph);

      // If this is the first node in the adjoint run order, then it is the last node in the forward graph
      if(node==adjointRunOrders[inputDimWrt][0]){
        std::stringstream nodeName;

        unsigned int inDim = filtGraph[*e.first]->inputDim;
        auto gradPiece = std::make_shared<GradientPiece>(piece, outputDimWrt, inDim);
        nodeName << filtGraph[adjointRunOrders[inputDimWrt][0]]->name << "_Gradient[" << outputDimWrt << "," << inDim << "]";
        try{
          gradGraph.AddNode(gradPiece, nodeName.str());
        }catch(const std::logic_error&){}

        gradGraph.AddEdge(outputName +"_Sensitivity", 0 , nodeName.str(), gradPiece->inputSizes.size()-1);

        // Add edges for the inputs
        for( auto e2 =in_edges(node, wgraph->graph); e2.first!=e2.second; ++e2.first ) {
          auto source = boost::source(*e2.first,wgraph->graph);
          gradGraph.AddEdge(wgraph->graph[source]->name, wgraph->graph[*e2.first]->outputDim,
                            nodeName.str(),              wgraph->graph[*e2.first]->inputDim);
        }

      }else{

        // Loop through all the edges and add the adjoint nodes and edges
        for(unsigned int outInd=0; outInd<outEdges.size(); ++outInd){

          unsigned int inDim = filtGraph[*e.first]->inputDim;
          auto gradPiece = std::make_shared<GradientPiece>(piece, outInd, inDim);

          std::stringstream nodeName;
          nodeName << filtGraph[node]->name << "_Gradient[" << outInd << "," << inDim << "]";
          gradGraph.AddNode(gradPiece, nodeName.str());

          // Add edges for the inputs
          for( auto e2 =in_edges(node, wgraph->graph); e2.first!=e2.second; ++e2.first ) {
            auto source = boost::source(*e2.first,wgraph->graph);
            gradGraph.AddEdge(wgraph->graph[source]->name, wgraph->graph[*e2.first]->outputDim,
                              nodeName.str(),              wgraph->graph[*e2.first]->inputDim);
          }

          // Add edges for the sensitivities.
          /* If the output of any node is used in more than one other ModPiece, to
             compute the gradient, we will need to accumulate the sensitivities from
             all outputs using a SumPiece.  This loop creates the SumPieces
          */

          // If there is more than one downstream node, we need to add up their sensitivities to this node
          if(edgeCounts[outInd]>1){
            std::stringstream sumName;
            sumName << baseName << "[" << outInd << "]_SensitivitySum";

            // Add edges for the gradient pieces that flow into the sum
            unsigned int edgeInd = 0;
            for(unsigned int i = 0; i<outEdges.at(outInd).size(); ++i){
              auto nextNode = filtGraph[outEdges.at(outInd).at(i).first];
              auto nextMod = std::dynamic_pointer_cast<ModPiece>(nextNode->piece);
              int nextNumOutputs = nextMod->outputSizes.size();

              for(unsigned int nextOutInd=0; nextOutInd<nextNumOutputs; ++nextOutInd){
                std::stringstream otherNodeName;
                otherNodeName << nextNode->name << "_Gradient[" << nextOutInd << "," << outEdges.at(outInd).at(i).second << "]";
                gradGraph.AddEdge(otherNodeName.str(),0, sumName.str(), edgeInd);
                edgeInd++;
              }
            }

            // Connect the SensitivitySum to the gradient pice
            gradGraph.AddEdge(sumName.str(), 0, nodeName.str(), gradPiece->inputSizes.size()-1);


          }else if(outEdges.at(outInd).size()>0){ // There is only one edge flowing into the gradient
            auto es=out_edges(outEdges.at(outInd).at(0).first, filtGraph);
            if(outEdges.at(outInd).at(0).first == adjointRunOrders[inputDimWrt][0]){
              std::stringstream otherNodeName;
              otherNodeName << filtGraph[adjointRunOrders[inputDimWrt][0]]->name << "_Gradient[" << outputDimWrt << "," << outEdges.at(outInd).at(0).second << "]";

              gradGraph.AddEdge(otherNodeName.str(), 0, nodeName.str(), gradPiece->inputSizes.size()-1);

            }else if(es.first != es.second){
              auto nextNode = filtGraph[outEdges.at(outInd).at(0).first];

              std::stringstream otherNodeName;
              otherNodeName << nextNode->name << "_Gradient[" << filtGraph[*es.first]->outputDim << "," << outEdges.at(outInd).at(0).second << "]";

              gradGraph.AddEdge(otherNodeName.str(), 0, nodeName.str(), gradPiece->inputSizes.size()-1);
            }
          }

        } // end of outgoing edge loop
      }
    } // end of incoming edge

  } // end of adjoint run order loop

  // Add any necessary edges for the final gradient
  std::string outNodeName;
  if(needsFinalSum){

    outNodeName = inputNames.at(inputDimWrt) + "[0]_SensitivitySum";

    for( auto eout=out_edges(inputNodes.at(inputDimWrt), filtGraph); eout.first!=eout.second; ++eout.first )
    {
      auto target = boost::target(*eout.first,filtGraph);
      auto targetPiece = std::dynamic_pointer_cast<ModPiece>(wgraph->graph[target]->piece);
      if(targetPiece){
        // Loop over all the output edges
        for(unsigned int otherOutInd=0; otherOutInd<targetPiece->outputSizes.size(); ++otherOutInd){
          std::stringstream otherName;
          otherName << wgraph->graph[target]->name << "_Gradient[" << otherOutInd << "," << wgraph->graph[*eout.first]->inputDim << "]";
          gradGraph.AddEdge(otherName.str(), 0, outNodeName, otherOutInd);
        }
      }

    }

  }else{

    // At this point, there should only be one node with an unconnected output.
    for(auto vs = boost::vertices(gradGraph.graph); vs.first!=vs.second; ++vs.first){
      auto pieceMod =  std::dynamic_pointer_cast<ModPiece>(gradGraph.graph[*vs.first]->piece);
      if(pieceMod){
        if(out_degree(*vs.first, gradGraph.graph) <pieceMod->outputSizes.size()){
          outNodeName = gradGraph.graph[*vs.first]->name;
          break;
        }
      }
    }
  }

  return gradGraph.CreateModPiece(outNodeName, inputNames);
}

std::shared_ptr<ModGraphPiece> ModGraphPiece::JacobianGraph(unsigned int                const  outputDimWrt,
                                                            unsigned int                const  inputDimWrt)
{
  assert(outputDimWrt<outputSizes.size());
  assert(inputDimWrt<inputSizes.size());

  WorkGraph jacGraph;
  auto& filtGraph = *filtered_graphs[inputDimWrt];

  std::vector<std::string> inputNames(constantPieces.size()+1);

  // Add the forward model
  for(auto node = runOrder.begin(); node != runOrder.end()-1; ++node){
    auto piece = wgraph->graph[*node]->piece;

    int placeholderInd = GetInputIndex(piece);
    if(placeholderInd>=0){
      auto pieceMod = std::dynamic_pointer_cast<ModPiece>(piece);
      assert(pieceMod);
      jacGraph.AddNode(std::make_shared<IdentityOperator>(pieceMod->outputSizes(0)), wgraph->graph[*node]->name);
      inputNames.at(placeholderInd) = wgraph->graph[*node]->name;

    }else{
      jacGraph.AddNode(piece, wgraph->graph[*node]->name);

      // Add the input edges
      for(auto es = in_edges(*node, wgraph->graph); es.first!=es.second; ++es.first){
        auto source = boost::source(*es.first, wgraph->graph);
        jacGraph.AddEdge(wgraph->graph[source]->name, wgraph->graph[*es.first]->outputDim, wgraph->graph[*node]->name, wgraph->graph[*es.first]->inputDim);
      }
    }
  }

  // Add a node in the graph to hold the vector input
  jacGraph.AddNode(std::make_shared<IdentityOperator>(inputSizes(inputDimWrt)), inputNames.at(inputDimWrt) + "_Vec");

  // Add the Jacobian components
  std::unordered_map<std::string, std::vector<std::string>> cumNames;

  for(auto node = adjointRunOrders[inputDimWrt].rbegin(); node!=adjointRunOrders[inputDimWrt].rend(); ++node){
    std::string baseName = wgraph->graph[*node]->name;
    auto piece = std::dynamic_pointer_cast<ModPiece>(filtGraph[*node]->piece);

    int inDegree = boost::in_degree(*node, filtGraph);
    int outDegree = boost::out_degree(*node, filtGraph);

    auto modCast = std::dynamic_pointer_cast<ModPiece>(filtGraph[*node]->piece);
    assert(modCast);

    /** If the node has multiple inputs, than we have to sum the Jacobians for each of them.
        Here, we add a SumPiece if necessary, and store the name of the node returning the cumulative Jacobian action.
    */
    cumNames[baseName] = std::vector<std::string>(out_degree(*node,filtGraph));
    for(auto eout = boost::out_edges(*node, filtGraph); eout.first!=eout.second; ++eout.first){
      if(inDegree>1){
        std::stringstream jacName;
        jacName << baseName << "_CumulativeJacobian[" << filtGraph[*eout.first]->outputDim << "]";
        cumNames[baseName].at(filtGraph[*eout.first]->outputDim) = jacName.str();

        jacGraph.AddNode(std::make_shared<SumPiece>(modCast->outputSizes(filtGraph[*eout.first]->outputDim),inDegree), jacName.str());
      }else if(inDegree==1){
        std::stringstream jacName;
        jacName << baseName << "_Jacobian[" << filtGraph[*eout.first]->outputDim <<  ",0]";
        cumNames[baseName].at(filtGraph[*eout.first]->outputDim) = jacName.str();
      }
    }

    // For each output in the filtered graph, add a jacobian term for each input in the filtered graph
    for(auto ein=boost::in_edges(*node, filtGraph); ein.first!=ein.second; ++ein.first){
      // If we're the last node....
      if(*node == adjointRunOrders[inputDimWrt][0]){
        std::stringstream jacName;
        jacName << baseName << "_Jacobian[" << outputDimWrt << "," << filtGraph[*ein.first]->inputDim << "]";
        auto jacPiece = std::make_shared<JacobianPiece>(piece, outputDimWrt, filtGraph[*ein.first]->inputDim);

        jacGraph.AddNode(jacPiece, jacName.str());

        // Add input edges from the original forward graph
        for(auto es=boost::in_edges(*node, wgraph->graph); es.first!=es.second; ++es.first){
          auto source = boost::source(*es.first, wgraph->graph);
          jacGraph.AddEdge(wgraph->graph[source]->name, wgraph->graph[*es.first]->outputDim, jacName.str(), wgraph->graph[*es.first]->inputDim);
        }

        auto source = boost::source(*ein.first,filtGraph);
        if(filtGraph[source]->name == inputNames.at(inputDimWrt)){
          jacGraph.AddEdge(inputNames.at(inputDimWrt)+ "_Vec", 0, jacName.str(), jacPiece->inputSizes.size()-1);
        }else if(boost::in_degree(source,filtGraph)>0){
          jacGraph.AddEdge(cumNames[filtGraph[source]->name].at(filtGraph[*ein.first]->outputDim), 0, jacName.str(),  jacPiece->inputSizes.size()-1);
        }

        // If there are multiple inputs, we have to add them up

      }else{

        for(auto eout = boost::out_edges(*node, filtGraph); eout.first!=eout.second; ++eout.first){

          std::stringstream jacName;
          jacName << baseName << "_Jacobian[" << filtGraph[*eout.first]->outputDim << "," << filtGraph[*ein.first]->inputDim << "]";

          // Try statement needed because we might end up adding the same piece multiple times
          try{
            auto jacPiece = std::make_shared<JacobianPiece>(piece, filtGraph[*eout.first]->outputDim, filtGraph[*ein.first]->inputDim);
            jacGraph.AddNode(jacPiece, jacName.str());

            // Add input edges from the original forward graph
            for(auto es=boost::in_edges(*node, wgraph->graph); es.first!=es.second; ++es.first){
              auto source = boost::source(*es.first, wgraph->graph);
              jacGraph.AddEdge(wgraph->graph[source]->name, wgraph->graph[*es.first]->outputDim, jacName.str(), wgraph->graph[*es.first]->inputDim);
            }

            // Add an input edge for the upstream Jacobian
            auto source = boost::source(*ein.first,filtGraph);
            if(filtGraph[source]->name == inputNames.at(inputDimWrt)){
              jacGraph.AddEdge(inputNames.at(inputDimWrt)+ "_Vec", 0, jacName.str(), jacPiece->inputSizes.size()-1);
            }else if(boost::in_degree(source,filtGraph)>0){
              jacGraph.AddEdge(cumNames[filtGraph[source]->name].at(filtGraph[*ein.first]->outputDim), 0, jacName.str(),  jacPiece->inputSizes.size()-1);
            }
          }catch(const std::logic_error& e){
             //std::cout << "duplicate..." << std::endl;
          }
        } // Loop over output edges
      }
    } // Loop over input edges

    // Potentially add up the Jacobians from multiple inputs
    for(auto eout = boost::out_edges(*node, filtGraph); eout.first!=eout.second; ++eout.first){
      if(inDegree>1){
        std::stringstream cumName;
        cumName << baseName << "_CumulativeJacobian[" << filtGraph[*eout.first]->outputDim << "]";

        // For each input ...
        unsigned int tempInd = 0;
        for(auto ein = boost::in_edges(*node, filtGraph); ein.first!=ein.second; ++ein.first){
          std::stringstream jacName;
          jacName << baseName << "_Jacobian[" << filtGraph[*eout.first]->outputDim <<  "," << filtGraph[*ein.first]->inputDim << "]";
          jacGraph.AddEdge(jacName.str(),0, cumName.str(), tempInd);
          tempInd++;
        }
      }// if(inDegree>1)
    }// loop over outputs

    // If this is the last node and it has multiple inputs, then we'll need to sum them up
    if(*node == adjointRunOrders[inputDimWrt][0]){
      unsigned int tempInd = 0;

      if(inDegree>1){
        std::stringstream jacName;
        jacName << baseName << "_CumulativeJacobian[" << outputDimWrt << "]";
        jacGraph.AddNode(std::make_shared<SumPiece>(modCast->outputSizes(outputDimWrt), inDegree), jacName.str());

        for(auto ein = boost::in_edges(*node, filtGraph); ein.first!=ein.second; ++ein.first){
          std::stringstream otherName;
          otherName << baseName << "_Jacobian[" << outputDimWrt << "," << filtGraph[*ein.first]->inputDim << "]";
          jacGraph.AddEdge(otherName.str(), 0, jacName.str(), tempInd);
          tempInd++;
        }
      }
    }

  } // Loop over adjoint run order

  // At this point, there should only be one node with an unconnected output.
  std::string outNodeName;
  for(auto vs = boost::vertices(jacGraph.graph); vs.first!=vs.second; ++vs.first){
    auto pieceMod =  std::dynamic_pointer_cast<ModPiece>(jacGraph.graph[*vs.first]->piece);
    if(pieceMod){
      if(out_degree(*vs.first, jacGraph.graph) <pieceMod->outputSizes.size()){
        outNodeName = jacGraph.graph[*vs.first]->name;
        break;
      }
    }
  }

  inputNames.at(constantPieces.size()) = inputNames.at(inputDimWrt) + "_Vec";
  return jacGraph.CreateModPiece(outNodeName, inputNames);
}

Eigen::VectorXi ModGraphPiece::ConstructInputSizes(std::vector<std::shared_ptr<ConstantVector> > const& constantPiecesIn)
{
  Eigen::VectorXi sizes(constantPiecesIn.size());
  for(int i=0; i<constantPiecesIn.size(); ++i){
    sizes(i)  = constantPiecesIn.at(i)->outputSizes(0);
    assert(constantPiecesIn.at(i)->outputSizes.size()==1);
  }
  return sizes;
}

void ModGraphPiece::ApplyHessianImpl(unsigned int                const  outWrt,
                                     unsigned int                const  inWrt1,
                                     unsigned int                const  inWrt2,
                                     ref_vector<Eigen::VectorXd> const& input,
                                     Eigen::VectorXd             const& sens,
                                     Eigen::VectorXd             const& vec)
{
  ref_vector<Eigen::VectorXd> newInputs(input.begin(),input.end());
  newInputs.push_back( std::cref(sens));
  newInputs.push_back( std::cref(vec));

  std::shared_ptr<ModGraphPiece> hessGraph;
  std::tuple<unsigned int, unsigned int, unsigned int> wrts = std::make_tuple(outWrt,inWrt1,inWrt2);

  auto iter = hessianPieces.find(wrts);
  if(iter==hessianPieces.end())
    hessianPieces[wrts] = GradientGraph(outWrt,inWrt1)->JacobianGraph(0,inWrt2);

  hessAction = hessianPieces[wrts]->Evaluate(newInputs).at(0);
}

void ModGraphPiece::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                      unsigned int                const  inputDimWrt,
                                      ref_vector<Eigen::VectorXd> const& input,
                                      Eigen::VectorXd             const& vec)
{
  ref_vector<Eigen::VectorXd> newInputs(input.begin(),input.end());
  newInputs.push_back( std::cref(vec));

  // Check to see if this input-output pair has been computed before.  If so, we don't need to recomptue the gradient graph
  std::pair<unsigned int, unsigned int> indexPair = std::make_pair(outputDimWrt, inputDimWrt);
  auto iter = jacobianPieces.find(indexPair);
  if(iter==jacobianPieces.end())
    jacobianPieces[indexPair] = JacobianGraph(outputDimWrt,inputDimWrt);

  // Evaluate the gradient
  jacobianAction = jacobianPieces[indexPair]->Evaluate(newInputs).at(0);
}

void ModGraphPiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  // set the inputs
  SetInputs(inputs);

  // fill the map from the WorkPiece ID to its outputs
  FillOutputMap();

  // store the result in the output vector
  outputs.resize(valMap[outputID].size());

  for(int i=0; i<outputs.size(); ++i) {
    outputs.at(i) = valMap[outputID].at(i).get();
  }
}

void ModGraphPiece::GradientImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& sensitivity) {

  ref_vector<Eigen::VectorXd> newInputs(input.begin(),input.end());
  newInputs.push_back( std::cref(sensitivity));

  // Check to see if this input-output pair has been computed before.  If so, we don't need to recomptue the gradient graph
  std::pair<unsigned int, unsigned int> indexPair = std::make_pair(outputDimWrt, inputDimWrt);
  auto iter = gradientPieces.find(indexPair);
  if(iter==gradientPieces.end()){
    gradientPieces[indexPair] = GradientGraph(outputDimWrt,inputDimWrt);
  }

  // Evaluate the gradient
  gradient = gradientPieces[indexPair]->Evaluate(newInputs).at(0);
}

void ModGraphPiece::JacobianImpl(unsigned int                const  wrtOut,
                                 unsigned int                const  wrtIn,
                                 ref_vector<Eigen::VectorXd> const& inputs) {

  // set the inputs
  SetInputs(inputs);

  // fill the map from the WorkPiece ID to its outputs
  FillOutputMap();

  // a map from the WorkPiece ID to a vector holding the cumulative jacobians of that output wrt the specified input
  std::map<unsigned int, std::vector<Eigen::MatrixXd> > jacMap;

  // loop through each downstream node
  for( auto node : boost::adaptors::reverse(adjointRunOrders[wrtIn]) ) {

    std::shared_ptr<ModPiece> nodePiece = std::dynamic_pointer_cast<ModPiece>(filtered_graphs[wrtIn]->operator[](node)->piece);
    assert(nodePiece);

    // the ID of the current node
    const unsigned int nodeID = nodePiece->ID();

    // Initialize the jacobian map for this node
    jacMap[nodeID] = std::vector<Eigen::MatrixXd>(nodePiece->outputSizes.size());

    // get the outputs of this node that impact the specified output node
    const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > & requiredOutNodes = RequiredOutputs(node, wrtIn, wrtOut);

    // remove duplicates
    std::vector<unsigned int> requiredOuts;
    requiredOuts.reserve(requiredOutNodes.size());
    for( auto out : requiredOutNodes ) {
      auto it = std::find(requiredOuts.begin(), requiredOuts.end(), std::get<1>(out));
      if( it==requiredOuts.end() ) {
	       requiredOuts.push_back(std::get<1>(out));
      }
    }

    // Initialize the jacobians that will be stored
    for(int i=0; i<requiredOuts.size(); ++i)
      jacMap[nodeID].at(requiredOuts.at(i)) = Eigen::MatrixXd::Zero(nodePiece->outputSizes(requiredOuts.at(i)), inputSizes(wrtIn));

    // get the inputs for this node --- the input WorkPiece ID, the output number, and the input number
    const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredIns = RequiredInputs(node, wrtIn);

    // if there are no inputs, it is an input node, so set the Jacobian to the identity
    if(requiredIns.size()==0){
      jacMap[nodeID][0] = Eigen::MatrixXd::Identity(nodePiece->outputSizes(0), nodePiece->outputSizes(0));
      assert(jacMap[nodeID].size()==1);
    }else{

      // the inputs to this WorkPiece
      const ref_vector<Eigen::VectorXd>& ins = GetNodeInputs(node);

      // compute the jacobian of each required output wrt each input
      for( auto out : requiredOuts ) {
        // To compute the Jacobian of out, we need to add the add the combination from each input
  	    for( auto in : requiredIns )
  	      jacMap[nodeID][out] += nodePiece->Jacobian(out, std::get<2>(in), ins) * jacMap[std::get<0>(in)][std::get<1>(in)];
      }
    }

  } // loop over run order

  // set the Jacobian for this WorkPiece
  jacobian = jacMap[outputID][wrtOut];
}
//
// void WorkGraphPiece::JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) {
//   // set the inputs
//   SetInputs(inputs);
//
//   // fill the map from the WorkPiece ID to its outputs
//   FillOutputMap();
//
//   // a map from the WorkPiece ID a map from the output number to the action of the jacobian of that output wrt the specified input
//   std::map<unsigned int, std::map<unsigned int, boost::any> > jacActionMap;
//
//   // loop through each downstream node
//   for( auto node : boost::adaptors::reverse(derivRunOrders[wrtIn]) ) {
//     // the ID of the current node
//     const unsigned int nodeID = filtered_graphs[wrtIn]->operator[](node)->piece->ID();
//
//     // get the outputs of this node that depend on the specified input
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredOutNodes = RequiredOutputs(node, wrtIn, wrtOut);
//     // remove duplicates
//     std::vector<unsigned int> requiredOuts;
//     requiredOuts.reserve(requiredOutNodes.size());
//     for( auto out : requiredOutNodes ) {
//       auto it = std::find(requiredOuts.begin(), requiredOuts.end(), std::get<1>(out));
//       if( it==requiredOuts.end() ) {
// 	requiredOuts.push_back(std::get<1>(out));
//       }
//     }
//
//     // get the inputs for this node --- the input WorkPiece ID, the output number, and the input number
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredIns = RequiredInputs(node, wrtIn);
//
//     // the inputs to this WorkPiece
//     const ref_vector<boost::any>& ins = Inputs(node);
//
//     // compute the jacobian of each required output wrt each input
//     for( auto out : requiredOuts ) {
//       if( requiredIns.size()==0 ) {
// 	// if there are no inputs, it is the input so set the Jacobian to the identity
// 	jacActionMap[nodeID][out] = vec;
//       } else {
// 	// initize the jacobian to nothing
// 	jacActionMap[nodeID][out] = boost::none;
//
// 	for( auto in : requiredIns ) {
// 	  // compute the Jacobian with respect to each required input
// 	  graph->operator[](node)->piece->JacobianAction(std::get<2>(in), out, jacActionMap[std::get<0>(in)][std::get<1>(in)], ins);
//
// 	  // use chain rule to get the jacobian wrt to the required input
// 	  jacActionMap[nodeID][out] = algebra->Add(jacActionMap[nodeID][out], *(graph->operator[](node)->piece->jacobianAction));
// 	}
//       }
//     }
//   }
//
//   // set the action of the Jacobian for this WorkPiece
//   jacobianAction = jacActionMap[outputID][wrtOut];
// }
//
// void WorkGraphPiece::JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) {
//     // set the inputs
//   SetInputs(inputs);
//
//   // fill the map from the WorkPiece ID to its outputs
//   FillOutputMap();
//
//   // a map from the WorkPiece ID a map from the output number to the action of the jacobian of that output wrt the specified input
//   std::map<unsigned int, std::map<unsigned int, boost::any> > jacTransActionMap;
//
//   // loop through each downstream node
//   for( auto node : derivRunOrders[wrtIn] ) {
//     // the ID of the current node
//     const unsigned int nodeID = filtered_graphs[wrtIn]->operator[](node)->piece->ID();
//
//     // get the outputs of this node that depend on the specified input
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredOuts = RequiredOutputs(node, wrtIn, wrtOut);
//
//     // get the inputs for this node --- the input WorkPiece ID, the output number, and the input number
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredIns = RequiredInputs(node, wrtIn);
//
//     // the inputs to this WorkPiece
//     const ref_vector<boost::any>& ins = Inputs(node);
//
//     for( auto in : requiredIns ) {
//       if( nodeID==outputID ) {
// 	assert(requiredOuts.size()==1);
// 	assert(std::get<1>(requiredOuts[0])==wrtOut);
//
// 	// compute the Jacobian transpose action of the output node
// 	graph->operator[](node)->piece->JacobianTransposeAction(std::get<2>(in), wrtOut, vec, ins);
// 	jacTransActionMap[nodeID][std::get<2>(in)] = *(graph->operator[](node)->piece->jacobianTransposeAction);
//       } else {
// 	// initialize the jacobian transpose action to nothing
// 	jacTransActionMap[nodeID][std::get<2>(in)] = boost::none;
//
// 	// loop through the outputs
// 	for( auto out : requiredOuts ) {
// 	  // compute the jacobian transpose action for this output
// 	  graph->operator[](node)->piece->JacobianTransposeAction(std::get<2>(in), std::get<1>(out), jacTransActionMap[std::get<0>(out)][std::get<2>(out)], ins);
// 	  // add it (chain rule)
// 	  jacTransActionMap[nodeID][std::get<2>(in)] = algebra->Add(jacTransActionMap[nodeID][std::get<2>(in)], *(graph->operator[](node)->piece->jacobianTransposeAction));
// 	}
//       }
//     }
//
//     // if this is the input node
//     if( requiredIns.size()==0 ) {
//       // loop though the outputs
//       for( auto out : requiredOuts ) {
// 	if( jacobianTransposeAction ) { // if the jacobian transpose action has not be initilized ...
// 	  // it is equal to the action of the output
// 	  *jacobianTransposeAction = algebra->Add(*jacobianTransposeAction, jacTransActionMap[std::get<0>(out)][std::get<1>(out)]);
// 	} else {
// 	  // add it to the existing jacobian transpose action (chain rule)
// 	  jacobianTransposeAction = jacTransActionMap[std::get<0>(out)][std::get<1>(out)];
// 	}
//       }
//     }
//   }
// }

void ModGraphPiece::SetInputs(ref_vector<Eigen::VectorXd> const& inputs) {
  // get the inputs and set them to the ConstantPiece nodes
  assert(inputs.size()==constantPieces.size());
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    constantPieces[i]->SetValue(inputs.at(i));
  }
}

std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > ModGraphPiece::InputNodes(boost::graph_traits<Graph>::vertex_descriptor const& node) const {
  // the map of input nodes
  std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > inMap;

  // loop though the input nodes
  boost::graph_traits<Graph>::in_edge_iterator e, e_end;
  for( tie(e, e_end)=boost::in_edges(node, wgraph->graph); e!=e_end; ++e ) {
    // get the WorkPiece id number, the output that it supplies, and the input that receives it
    const unsigned int id = wgraph->graph[boost::source(*e, wgraph->graph)]->piece->ID();
    const unsigned int inNum = wgraph->graph[*e]->inputDim;
    const unsigned int outNum = wgraph->graph[*e]->outputDim;

    // try to find the WorkPiece in the other upstream nodes
    auto it = inMap.find(id);

    if( it==inMap.end() ) { // if we have not yet needed this WorkPiece ...
      // ... add it to the list and store the input/output pair
      inMap[id] = std::vector<std::pair<unsigned int, unsigned int> >(1, std::pair<unsigned int, unsigned int>(inNum, outNum));
    } else { // we have needed this WorkPiece
      // ... add the input/output pair
      inMap[id].push_back(std::pair<unsigned int, unsigned int>(inNum, outNum));
    }
  }

  return inMap;
}

void ModGraphPiece::FillOutputMap() {
  // clear the map
  valMap.clear();

  // loop over the run order
  for( auto it : runOrder ) {

    // the inputs to this WorkPiece
    const ref_vector<Eigen::VectorXd>& ins = GetNodeInputs(it);

    // evaluate the current map and store a vector of references to its output
    auto wPiece = wgraph->GetPiece(it);
    auto currPiece = std::dynamic_pointer_cast<ModPiece>(wPiece);

    if(!currPiece){

      // If it can't be cast to a ModPiece, check to see if the output can be cast to an Eigen vector
      ref_vector<Eigen::VectorXd> output;
      std::vector<boost::any> anyIns(ins.size());
      for(int i=0; i<ins.size(); ++i)
        anyIns.at(i) = boost::any(ins.at(i));

      std::vector<boost::any> const& anyOut = wPiece->Evaluate(anyIns);

      for(int i=0; i<anyOut.size(); ++i){
        Eigen::VectorXd const& temp = AnyConstCast(anyOut.at(i));
        output.push_back(std::cref(temp));
      }
      valMap[wgraph->GetPiece(it)->ID()] = output;

    }else{
      assert(currPiece);
      valMap[wgraph->GetPiece(it)->ID()] = ToRefVector(currPiece->Evaluate(ins));
    }
  }
}

ref_vector<Eigen::VectorXd> ModGraphPiece::GetNodeInputs(boost::graph_traits<Graph>::vertex_descriptor node) const {

  // how many inputs does this node require?
  const int numIns = wgraph->GetPiece(node)->numInputs;

  // get the inputs for this node
  const std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > >& inMap = InputNodes(node);

  Eigen::VectorXd empty;
  ref_vector<Eigen::VectorXd> ins(numIns, std::cref(empty));

  // loop through the edges again, now we know which outputs supply which inputs
  for( auto edge : inMap ) {
    // loop over the input/output pairs supplied by this input
    for( auto in_out : edge.second ) {
      ins[in_out.first] = valMap.at(edge.first)[in_out.second];
    }
  }

  return ins;
}

std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > ModGraphPiece::RequiredOutputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn, unsigned int const wrtOut) const {
  // the ID of the current node
  const unsigned int nodeID = filtered_graphs[wrtIn]->operator[](node)->piece->ID();

  // get the outputs of this node that depend on the specified input
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > requiredOuts;

  if( nodeID==outputID ) { // if it is the output node ...
    // ... the user specifies the output derivative
    requiredOuts.push_back(std::tuple<unsigned int, unsigned int, unsigned int>(nodeID, wrtOut, wrtIn));

    return requiredOuts;
  }

  // loop though the output nodes
  boost::graph_traits<FilteredGraph>::out_edge_iterator eout, eout_end;
  for( tie(eout, eout_end)=boost::out_edges(node, *filtered_graphs[wrtIn]); eout!=eout_end; ++eout ) {
    // get the output number
    const unsigned int id = wgraph->GetPiece(boost::target(*eout, *filtered_graphs[wrtIn]))->ID();
    const unsigned int outNum = filtered_graphs[wrtIn]->operator[](*eout)->outputDim;
    const unsigned int inNum = filtered_graphs[wrtIn]->operator[](*eout)->inputDim;

    // if we have not already required this output, save it
    auto it = std::find(requiredOuts.begin(), requiredOuts.end(), std::tuple<unsigned int, unsigned int, unsigned int>(id, outNum, inNum));
    if( it==requiredOuts.end() ) {
      requiredOuts.push_back(std::tuple<unsigned int, unsigned int, unsigned int>(id, outNum, inNum));
    }
  }

  return requiredOuts;
}

std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > ModGraphPiece::RequiredInputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn) const {
  // how many inputs does this node require?
  const int numIns = filtered_graphs[wrtIn]->operator[](node)->piece->numInputs;

  std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > requiredIns;
  requiredIns.reserve(numIns);

  // loop though the output nodes
  boost::graph_traits<FilteredGraph>::in_edge_iterator ein, ein_end;
  for( tie(ein, ein_end)=boost::in_edges(node, *filtered_graphs[wrtIn]); ein!=ein_end; ++ein ) {
    // get the WorkPiece id number, the output that it supplies, and the input that receives it
    const unsigned int id = wgraph->GetPiece(boost::source(*ein, *filtered_graphs[wrtIn]))->ID();
    const unsigned int outNum = filtered_graphs[wrtIn]->operator[](*ein)->outputDim;
    const unsigned int inNum = filtered_graphs[wrtIn]->operator[](*ein)->inputDim;

    // store the requried input
    requiredIns.push_back(std::tuple<unsigned int, unsigned int, unsigned int>(id, outNum, inNum));
  }

  return requiredIns;
}

std::shared_ptr<ModGraphPiece> ModGraphPiece::GetSubModel(std::string const& nodeName) const
{

  // Create a new workgraph wit the nodes necessary to evaluate nodeName
  auto subGraph = wgraph->Clone();
  for(auto piece : constantPieces)
    subGraph->RemoveNode(wgraph->GetName(piece));

  return subGraph->CreateModPiece(nodeName);
}

std::vector<int> ModGraphPiece::MatchInputs(std::shared_ptr<ModGraphPiece> otherPiece) const
{

  std::vector<std::shared_ptr<ConstantVector>> otherIns = otherPiece->GetConstantPieces();
  std::shared_ptr<WorkGraph> otherGraph = otherPiece->GetGraph();

  std::vector<int> outputs(otherIns.size());

  for(int i=0; i<otherIns.size(); ++i)
  {
    // get the downstream node and input index corresponding to this constant piece
    std::string constName = otherGraph->GetName( otherIns.at(i) );
    std::string sharedName = otherGraph->GetChildren( constName ).at(0);

    // Now try to find the same node and input in *this graph
    if(wgraph->HasNode(sharedName)){
      int inputIndex = otherGraph->GetEdges( constName, sharedName ).at(0).second;

      std::string upstreamName = wgraph->GetParent( sharedName, inputIndex);
      assert(upstreamName.size()>0);

      // Now, get the parent piece and check it against all of the constant pieces
      auto iter = std::find(constantPieces.begin(), constantPieces.end(), wgraph->GetPiece(upstreamName));
      if(iter != constantPieces.end()){
        outputs.at(i) = std::distance(constantPieces.begin(), iter);
      }else{
        outputs.at(i) = -1;
      }

    }else{
      outputs.at(i) = -1;
    }
  }

  return outputs;
}
