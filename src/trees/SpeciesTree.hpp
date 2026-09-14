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

  ~SpeciesTree() = default;

  // forbid copy and move
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

  size_t getHash() const { return _speciesTree.getTreeHash(); }
  size_t getNodeIndexHash() const { return _speciesTree.getTreeHash(false); }

  const PLLRootedTree &getTree() const { return _speciesTree; }
  const DatedTree &getDatedTree() const { return _datedTree; }
  PLLRootedTree &getTree() { return _speciesTree; }
  DatedTree &getDatedTree() { return _datedTree; }

  void getLabelToId(StringToUint &labelToId) const;

  class Listener {
  public:
    virtual ~Listener() = default;
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
  SpeciesTreeOperator() = delete;

  /**
   *  restoreDates operation:
   *  speciesTree branch lenghts, node ranks and speciation order are changed
   *  in accordance with the ranks in backup
   */
  static void restoreDates(SpeciesTree &speciesTree, const DatedBackup &backup);

  /**
   *  changeRoot operation:
   *  speciesTree is changed in such a way that one of root's children becomes
   *  the new root, one of the new root's children is swapped with root and
   *  itself takes the place of root's original child
   *
   *  direction shows on which branch to place the new root relative to root:
   *  0 == right-right, 1 == left-right, 2 == right-left, 3 == left-left
   */
  static bool canChangeRoot(const SpeciesTree &speciesTree,
                            unsigned int direction);
  static void changeRoot(SpeciesTree &speciesTree, unsigned int direction);
  static void revertChangeRoot(SpeciesTree &speciesTree,
                               unsigned int direction);

  /**
   *  applySPRMove operation:
   *  speciesTree is changed in such a way that prune's brother takes the place
   *  of prune's father, prune's father takes the place of regraft, and regraft
   *  becomes prune's new brother
   */
  static bool canApplySPRMove(SpeciesTree &speciesTree, unsigned int prune,
                              unsigned int regraft);
  static unsigned int applySPRMove(SpeciesTree &speciesTree, unsigned int prune,
                                   unsigned int regraft);
  static void reverseSPRMove(SpeciesTree &speciesTree, unsigned int prune,
                             unsigned int applySPRMoveReturnValue);

  /**
   *  Get nodes that can be moved in speciesTree based on precomputed branch
   *  quartet supports
   */
  static void getPossiblePrunes(const SpeciesTree &speciesTree,
                                const std::vector<double> &supportValues,
                                double maxSupport,
                                std::vector<unsigned int> &prunes);

  /**
   *  Get possible destination nodes to move prune to in speciesTree based on
   *  the SPR radius and basic topological constraints
   */
  static void getPossibleRegrafts(const SpeciesTree &speciesTree,
                                  unsigned int prune, unsigned int radius,
                                  std::vector<unsigned int> &regrafts);

  /**
   *  Add to affectedBranches all branches whose bipartition would change if
   *  we pruned prune and regrafted it to regraft in speciesTree
   */
  static void getAffectedBranches(SpeciesTree &speciesTree, unsigned int prune,
                                  unsigned int regraft,
                                  std::vector<unsigned int> &affectedBranches);
};
