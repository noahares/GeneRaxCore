#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <maths/Parameters.hpp>
#include <search/Moves.hpp>
#include <sstream>
#include <stack>
#include <trees/PLLRootedTree.hpp>
#include <util/GeneRaxCheckpoint.hpp>
#include <util/enums.hpp>

struct RecModelInfo;
void printLibpllNode(corax_unode_s *node, std::ostream &os, bool isRoot);
void printLibpllTreeRooted(corax_unode_t *root, std::ostream &os);

class JointTree {
public:
  JointTree(const std::string &newickString, const std::string &alignment_file,
            const std::string &speciestree_file,
            const std::string &geneSpeciesMapfile,
            const std::string &substitutionModel,
            const RecModelInfo &recModelInfo, RecOpt reconciliationOpt,
            bool madRooting, double supportThreshold, double recWeight,
            bool safeMode, bool optimizeDTLRates, const Parameters &ratesVector,
            const std::string &checkpointPath);
  JointTree(const JointTree &) = delete;
  JointTree &operator=(const JointTree &) = delete;
  JointTree(JointTree &&) = delete;
  JointTree &operator=(JointTree &&) = delete;

  virtual ~JointTree();
  void printLibpllTree() const;
  void optimizeParameters(bool felsenstein = true, bool reconciliation = true);
  double computeLibpllLoglk(bool incremental = false);
  double computeReconciliationLoglk();
  double computeJointLoglk();
  void printLoglk(bool libpll = true, bool rec = true, bool joint = true,
                  Logger &os = Logger::info);
  corax_unode_t *getNode(unsigned int index);
  void applyMove(SPRMove &move);
  void optimizeMove(SPRMove &move);
  void reOptimizeMove(SPRMove &move);

  void invalidateCLV(corax_unode_s *node);
  void printAllNodes(std::ostream &os);
  void printInfo();
  void rollbackLastMove();
  void save(const std::string &fileName, bool append);
  corax_treeinfo_t *getTreeInfo();
  void setRates(const Parameters &ratesVector);
  PLLRootedTree &getSpeciesTree() { return _speciesTree; }
  size_t getUnrootedTreeHash();
  ReconciliationEvaluation &getReconciliationEvaluation() {
    return *reconciliationEvaluation_;
  }
  std::shared_ptr<ReconciliationEvaluation> getReconciliationEvaluationPtr() {
    return reconciliationEvaluation_;
  }

  corax_unode_t *getRoot() const {
    return reconciliationEvaluation_->getRoot();
  }
  void setRoot(corax_unode_t *root) {
    reconciliationEvaluation_->setRoot(root);
  }
  const Parameters &getRatesVector() const { return _ratesVector; }
  void inferMLScenario(Scenario &scenario) {
    reconciliationEvaluation_->inferMLScenario(scenario);
  }
  bool isSafeMode() { return _safeMode; }
  void enableReconciliation(bool enable) { _enableReconciliation = enable; }
  void enableLibpll(bool enable) { _enableLibpll = enable; }
  unsigned int getGeneTaxaNumber() { return getTreeInfo()->tip_count; }
  PLLUnrootedTree &getGeneTree() { return _libpllEvaluation.getGeneTree(); }
  Model &getModel() { return _libpllEvaluation.getModel(); }
  const GeneSpeciesMapping &getMappings() const { return _geneSpeciesMap; }
  double getSupportThreshold() const { return _supportThreshold; }
  /**
   *  Are we allowed to apply a move that would move a prune node that is
   *  on one side of branch to the other side of branch?
   *  This might be constrained by the enforced gene root for instance
   */
  bool canSPRCrossBranch(const corax_unode_t *branch) const;
  void saveCheckpoint();
  const GeneRaxCheckpoint &getCheckpoint() const { return _checkpoint; }

private:
  GeneRaxCheckpoint _checkpoint;
  LibpllEvaluation _libpllEvaluation;
  std::shared_ptr<ReconciliationEvaluation> reconciliationEvaluation_;
  PLLRootedTree _speciesTree;
  GeneSpeciesMapping _geneSpeciesMap;
  Parameters _ratesVector;
  std::stack<std::shared_ptr<SPRRollback>> _rollbacks;
  bool _optimizeDTLRates;
  bool _safeMode;
  bool _enableReconciliation;
  bool _enableLibpll;
  RecOpt _recOpt;
  double _recWeight;
  double _supportThreshold;
  bool _madRooting;
  std::string _enforcedRootedGeneTree;
};
