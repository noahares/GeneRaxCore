#include "ReconciliationBLEstimator.hpp"

#include <algorithm>
#include <memory>

#include <maths/ModelParameters.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <routines/Routines.hpp>

static void getAverageDepthRec(corax_rnode_t *speciesNode, double currentDepth,
                               double &sumDepths, unsigned int &count) {
  currentDepth += speciesNode->length;
  if (!speciesNode->left) {
    count++;
    sumDepths += currentDepth;
  } else {
    getAverageDepthRec(speciesNode->left, currentDepth, sumDepths, count);
    getAverageDepthRec(speciesNode->right, currentDepth, sumDepths, count);
  }
}

static double getAverageDepth(corax_rnode_t *speciesNode) {
  double sumDepths = 0.0;
  unsigned int count = 0;
  // we do not want to count this node's length
  getAverageDepthRec(speciesNode, -speciesNode->length, sumDepths, count);
  assert(count);
  return sumDepths / double(count);
}

static void balanceRoot(PLLRootedTree &speciesTree) {
  auto root = speciesTree.getRoot();
  auto left = root->left;
  auto right = root->right;
  auto initialLength = left->length + right->length;
  auto leftDepth = getAverageDepth(left);
  auto rightDepth = getAverageDepth(right);
  auto diff = leftDepth - rightDepth;
  left->length -= diff / 2.0;
  right->length += diff / 2.0;
  double epsilon = 0.0000001;
  left->length =
      std::min(std::max(epsilon, left->length), initialLength - epsilon);
  right->length =
      std::min(std::max(epsilon, right->length), initialLength - epsilon);
}

/**
 *  Fill branch lengths between consequent speciation or leaf events that
 *  happened at directly related nodes of the species tree
 *  Speciation-loss events are not used because they happen along a gene
 *  branch, and not on a gene node, and thus do not hold the time
 *  of speciation
 */
static void estimateBLRecursive(
    corax_rtree_t *speciesTree, corax_unode_t *node, unsigned int depth,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    unsigned int ancestralSpeciesId, double lengthToAncestralSpecies,
    double familyWeight, std::vector<double> &speciesSumBL,
    std::vector<double> &speciesSumWeight) {
  const auto &event = geneToEvents[node->node_index].back(); // terminal event
  bool isSpeciation = (event.type == ReconciliationEventType::EVENT_S ||
                       event.type == ReconciliationEventType::EVENT_None);
  // divide the root BL by two to place the root
  // at the middle of this branch
  lengthToAncestralSpecies += (depth == 1) ? node->length / 2.0 : node->length;
  if (isSpeciation) {
    auto currentSpeciesId = event.speciesNode;
    auto currentSpeciesNode = speciesTree->nodes[currentSpeciesId];
    bool isDirectSpeciation =
        (currentSpeciesNode->parent &&
         (currentSpeciesNode->parent->node_index == ancestralSpeciesId));
    if (isDirectSpeciation) {
      speciesSumBL[currentSpeciesId] += lengthToAncestralSpecies * familyWeight;
      speciesSumWeight[currentSpeciesId] += familyWeight;
    }
    ancestralSpeciesId = currentSpeciesId;
    lengthToAncestralSpecies = 0.0;
  }
  if (node->next) {
    auto left = (depth == 0) ? node->next : node->next->back;
    auto right = (depth == 0) ? node->next->back : node->next->next->back;
    estimateBLRecursive(speciesTree, left, depth + 1, geneToEvents,
                        ancestralSpeciesId, lengthToAncestralSpecies,
                        familyWeight, speciesSumBL, speciesSumWeight);
    estimateBLRecursive(speciesTree, right, depth + 1, geneToEvents,
                        ancestralSpeciesId, lengthToAncestralSpecies,
                        familyWeight, speciesSumBL, speciesSumWeight);
  }
}

static void estimateBLForFamily(const Scenario &scenario, double familyWeight,
                                std::vector<double> &speciesSumBL,
                                std::vector<double> &speciesSumWeight) {
  auto geneRoot = scenario.getGeneRoot();
  corax_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = scenario.getVirtualRootIndex();
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  auto ancestralSpeciesId = static_cast<unsigned int>(-1);
  double lengthToAncestralSpecies = 0.0;
  estimateBLRecursive(scenario.getSpeciesTree(), &virtualRoot, 0,
                      scenario.getGeneIdToEvents(), ancestralSpeciesId,
                      lengthToAncestralSpecies, familyWeight, speciesSumBL,
                      speciesSumWeight);
}

static double getFamilyWeight(const FamilyInfo &info) {
  if (info.alignmentFile == "" || info.libpllModel == "" ||
      info.libpllModel == "true") {
    return 1.0;
  }
  try {
    Model model(info.libpllModel);
  } catch (...) {
    return 1.0;
  }
  auto alignmentLength =
      LibpllParsers::getMSAEntropy(info.alignmentFile, info.libpllModel);
  if (alignmentLength == 0.0) {
    return 1.0;
  }
  return alignmentLength;
}

void ReconciliationBLEstimator::estimate(
    const std::string &speciesTreeFile, const Families &families,
    const ModelParameters &modelParameters) {
  PLLRootedTree speciesTree(speciesTreeFile, true);
  PerCoreGeneTrees geneTrees(families);
  const unsigned int samples = 0;
  const bool optimizeRates = false;
  std::vector<std::shared_ptr<Scenario>> scenarios;
  Logger::timed << std::endl;
  Logger::timed << "[Species BL estimation] Infering reconciliation scenarios"
                << std::endl;
  Routines::inferAndGetReconciliationScenarios(speciesTree, geneTrees,
                                               modelParameters, samples,
                                               optimizeRates, scenarios);
  double speciesNodesNumber = speciesTree.getNodeNumber();
  std::vector<double> speciesSumBL(speciesNodesNumber, 0.0);
  std::vector<double> speciesSumWeight(speciesNodesNumber, 0.0);

  Logger::timed
      << "[Species BL estimation] Infering branch lengths from gene trees"
      << std::endl;
  for (unsigned int i = 0; i < geneTrees.getTrees().size(); ++i) {
    double familyWeight =
        getFamilyWeight(families[geneTrees.getTrees()[i].familyIndex]);
    estimateBLForFamily(*scenarios[i], familyWeight, speciesSumBL,
                        speciesSumWeight);
  }
  ParallelContext::sumVectorDouble(speciesSumBL);
  ParallelContext::sumVectorDouble(speciesSumWeight);
  for (unsigned int e = 0; e < speciesSumBL.size(); ++e) {
    auto length = 0.0;
    if (speciesSumWeight[e] != 0.0) {
      length = speciesSumBL[e] / speciesSumWeight[e];
    }
    speciesTree.getNode(e)->length = length;
  }
  balanceRoot(speciesTree);
  if (ParallelContext::getRank() == 0) {
    speciesTree.save(speciesTreeFile);
  }
  Logger::timed << "[Species BL estimation] Done" << std::endl;
}
