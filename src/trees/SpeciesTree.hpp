#pragma once

#include <memory>

#include "DatedTree.hpp"
#include "PLLRootedTree.hpp"
#include <IO/Families.hpp>
#include <util/types.hpp>

class SpeciesTree {
public:
  SpeciesTree(const std::string &str, bool isFile, bool useBLs);
  SpeciesTree(const std::unordered_set<std::string> &labels);
  SpeciesTree(const Families &families);

  // forbid copy
  SpeciesTree(const SpeciesTree &) = delete;
  SpeciesTree &operator=(const SpeciesTree &) = delete;
  SpeciesTree(SpeciesTree &&) = delete;
  SpeciesTree &operator=(SpeciesTree &&) = delete;

  std::unique_ptr<SpeciesTree> buildRandomTree() const;

  void saveToFile(const std::string &fileName, bool masterRankOnly) const;
  std::string toString() const { return _speciesTree.getNewickString(); }

  corax_rnode_t *getRoot() const { return _speciesTree.getRoot(); }
  corax_rnode_t *getNode(unsigned int nodeIndex) const {
    return _speciesTree.getNode(nodeIndex);
  }

  const PLLRootedTree &getTree() const { return _speciesTree; }
  const DatedTree &getDatedTree() const { return _datedTree; }
  PLLRootedTree &getTree() { return _speciesTree; }
  DatedTree &getDatedTree() { return _datedTree; }

  void getLabelToId(StringToUint &labelToId) const;

  size_t getHash() const;
  size_t getNodeIndexHash() const;

  class Listener {
  public:
    virtual ~Listener() {}
    virtual void onSpeciesDatesChange() = 0;
    virtual void onSpeciesTreeChange(
        const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) = 0;
  };
  void addListener(Listener *listener);
  void removeListener(Listener *listener);

  // should be called every time after changing the tree node dates
  void onSpeciesDatesChange();
  // should be called every time after changing the tree topology
  void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);

private:
  PLLRootedTree _speciesTree;
  DatedTree _datedTree;
  std::vector<Listener *> _listeners;
};

class SpeciesTreeOperator {
public:
  static void restoreDates(SpeciesTree &speciesTree, const DatedBackup &backup);
  static bool canChangeRoot(const SpeciesTree &speciesTree,
                            unsigned int direction);
  /**
   * Change the root to the neighboring branch described by direction where
   * direction is in [0:4[
   */
  static void changeRoot(SpeciesTree &speciesTree, unsigned int direction);
  static void revertChangeRoot(SpeciesTree &speciesTree,
                               unsigned int direction);
  static bool canApplySPRMove(SpeciesTree &speciesTree, unsigned int prune,
                              unsigned int regraft);
  /**
   *  Add to affectedBranches all branches that would be affected (whose
   * bipartition would change) if we prune prune and regraft it to regraft on
   * speciesTree
   */
  static void getAffectedBranches(SpeciesTree &speciesTree, unsigned int prune,
                                  unsigned int regraft,
                                  std::vector<unsigned int> &affectedBranches);
  static unsigned int applySPRMove(SpeciesTree &speciesTree, unsigned int prune,
                                   unsigned int regraft);
  static void reverseSPRMove(SpeciesTree &speciesTree, unsigned int prune,
                             unsigned int applySPRMoveReturnValue);
  static void getPossiblePrunes(SpeciesTree &speciesTree,
                                std::vector<unsigned int> &prunes,
                                std::vector<double> support, double maxSupport);
  static void getPossibleRegrafts(SpeciesTree &speciesTree, unsigned int prune,
                                  unsigned int radius,
                                  std::vector<unsigned int> &regrafts);
};
