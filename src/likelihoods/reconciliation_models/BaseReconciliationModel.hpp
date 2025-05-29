#pragma once

#include <cmath>
#include <memory>
#include <unordered_set>

#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <maths/Random.hpp>
#include <maths/ScaledValue.hpp>
#include <trees/PLLRootedTree.hpp>
#include <util/RecModelInfo.hpp>
#include <util/Scenario.hpp>
#include <util/enums.hpp>
#include <util/types.hpp>

#define IS_PROBA(x) ((REAL(0.0) <= REAL(x)) && (REAL(x) <= REAL(1.0)))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));

inline double solveSecondDegreePolynome(double a, double b, double c) {
  return (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
}

/**
 *  Common implementations for all reconciliation likelihood
 *  computation classes
 */
class BaseReconciliationModel {
public:
  BaseReconciliationModel(PLLRootedTree &speciesTree,
                          const GeneSpeciesMapping &geneSpeciesMapping,
                          const RecModelInfo &recModelInfo);

  virtual ~BaseReconciliationModel() {}

  /**
   *  Set the per-species branch rates
   */
  virtual void setRates(const RatesVector &rates) = 0;

  /**
   *  Should be called after changing speciation order on a fixed
   *  species tree topology
   */
  virtual void onSpeciesDatesChange() { invalidateAllSpeciesNodes(); }

  /**
   *  Should be called after each change in the species tree topology
   */
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);

  /**
   *  Make CLV components to be recomputed for all species nodes upon a CLV
   * update
   */
  virtual void invalidateAllSpeciesNodes() { _allSpeciesNodesInvalid = true; }

  /**
   *  (Incrementally) compute and return the reconciliation likelihood
   */
  virtual double computeLogLikelihood() = 0;

  /**
   *  Fill scenario with the maximum likelihood set of
   *  events that would lead to the current tree.
   *  Return true in case of success
   */
  virtual bool inferMLScenario(Scenario &scenario) = 0;

  /**
   *  Sample scenarios and add them to the scenarios vector.
   *  Return true in case of success
   */
  virtual bool
  sampleReconciliations(unsigned int samples,
                        std::vector<std::shared_ptr<Scenario>> &scenarios) = 0;

protected:
  /**
   *  Given the gene clades, fill _geneToSpecies, _speciesCoverage
   *  and _numberOfCoveredSpecies and build the pruned species tree
   *  representation based on the species coverage
   */
  virtual void mapGenesToSpecies() = 0;

  /**
   *  Accessors to the current species tree state
   */
  PLLRootedTree &getSpeciesTree() { return _speciesTree; }
  unsigned int getAllSpeciesNodeNumber() const {
    return _allSpeciesNodes.size();
  }
  unsigned int getPrunedSpeciesNodeNumber() const {
    return _prunedSpeciesNodes.size();
  }
  std::vector<corax_rnode_t *> &getAllSpeciesNodes() {
    return _allSpeciesNodes;
  }
  std::vector<corax_rnode_t *> &getPrunedSpeciesNodes() {
    return _prunedSpeciesNodes;
  }
  bool prunedMode() const { return _info.pruneSpeciesTree; }
  size_t getSpeciesTreeHash() const;

  /**
   *  Accessors to the internal species tree representation
   */
  corax_rnode_t *getSpeciesLeft(corax_rnode_t *node) {
    return _speciesLeft[node->node_index];
  }
  corax_rnode_t *getSpeciesRight(corax_rnode_t *node) {
    return _speciesRight[node->node_index];
  }
  corax_rnode_t *getSpeciesParent(corax_rnode_t *node) {
    return _speciesParent[node->node_index];
  }
  corax_rnode_t *getPrunedRoot() { return _prunedRoot; }

  /**
   *  Callback to be always called at the start of recomputing CLVs
   */
  void beforeComputeCLVs() {
    if (_allSpeciesNodesInvalid || _invalidatedSpeciesNodes.size()) {
      recomputeSpeciesProbabilities();
    }
  }

private:
  /**
   *  Init all structures describing the species tree, in particular the
   *  structure representing the pruned species tree in the pruned mode.
   *  This function should be called once at the start
   */
  void initSpeciesTree();

  /**
   *  Fill the nodes vector with all the children of the given node based
   *  on the pointers of the node.
   *  Nodes are filled in the postorder fashion
   */
  void fillNodesPostOrder(corax_rnode_t *node,
                          std::vector<corax_rnode_t *> &nodes);

  /**
   *  Fill the nodes vector with all the children of the given node based
   *  on the model's species tree representation. Use to fill
   * _prunedSpeciesNodes. Nodes are filled in the postorder fashion
   */
  void fillPrunedNodesPostOrder(corax_rnode_t *node,
                                std::vector<corax_rnode_t *> &nodes);

  /**
   *  Set the path to a file containing information about the
   *  proportions of missing genes for each extant species
   */
  void setFractionMissingGenes(const std::string &fractionMissingFile);

  /**
   *  Recompute probability values depending on the species tree only
   *  (not on the gene tree), such as the exctinction probability or
   *  the per-species event probabilities.
   *  This function should be typically called after changing the DTL
   *  rates or after updating the species tree
   */
  virtual void recomputeSpeciesProbabilities() = 0;

  size_t getTreeHashRec(const corax_rnode_t *node, size_t i) const;

protected:
  // description of the recmodel
  RecModelInfo _info;
  // reference to the species tree
  PLLRootedTree &_speciesTree;
  // list of all species nodes in postorder used for the likelihood computation
  std::vector<corax_rnode_t *> _allSpeciesNodes;
  std::vector<corax_rnode_t *> _prunedSpeciesNodes;
  // map species leaf names to species leaf indices. Species leaf indices run
  // from 0 to (_speciesTree.getLeafNumber() - 1)
  std::map<std::string, unsigned int> _speciesNameToId;
  // map gene leaf names to species leaf names
  std::map<std::string, std::string> _geneNameToSpeciesName;
  // map gene leaf indices to species leaf indices. Not computed by this class
  std::map<unsigned int, unsigned int> _geneToSpecies;
  // number of gene copies covering each species leaf. Not computed by this
  // class
  std::vector<unsigned int> _speciesCoverage;
  // number of species leaves covered by this gene family. Not computed by this
  // class
  unsigned int _numberOfCoveredSpecies;
  // fraction of missing genes, indexed by species leaf indices
  std::vector<double> _fm;
  // if true, updating a CLV will recompute its values for all species nodes
  bool _allSpeciesNodesInvalid;
  // species nodes for which values of a CLV will be recomputed on its update
  std::unordered_set<corax_rnode_t *> _invalidatedSpeciesNodes;
  // internal representation of the current species tree. Always use these
  // pointers to be compliant with the pruned species tree mode
  std::vector<corax_rnode_t *> _speciesLeft;
  std::vector<corax_rnode_t *> _speciesRight;
  std::vector<corax_rnode_t *> _speciesParent;
  corax_rnode_t *_prunedRoot;
};
