#pragma once

#include "PLLRootedTree.hpp"

/**
 *  Wrapper around PLLRootedTree handling the order of speciations
 */
class DatedTree {
public:
  DatedTree(PLLRootedTree &rootedTree, bool useBLs);

  const PLLRootedTree &getRootedTree() const { return _rootedTree; }
  PLLRootedTree &getRootedTree() { return _rootedTree; }

  bool isDated() const { return _fromBL; }

  const std::vector<corax_rnode_t *> &getOrderedSpeciations() const {
    return _orderedSpeciations;
  }
  const std::vector<unsigned int> &getOrderedSpeciesRanks() const {
    return _ranks;
  }
  unsigned int getRank(corax_rnode_t *node) const {
    return _ranks[node->node_index];
  }
  DatedBackup getBackup() const { return _ranks; }

  void updateSpeciationOrderAndRanks();
  void rescaleBranchLengths();

  bool moveUp(unsigned int rank);
  bool moveDown(unsigned int rank);

  void restore(const DatedBackup &backup);

  bool canTransferUnderRelDated(unsigned int e, unsigned int d) const;

  void randomize();

  /**
   *  Hash value characterizing the current order of speciations
   */
  size_t getOrderingHash(size_t startingHash = 42) const;

private:
  void checkRanks() const;

private:
  const bool _fromBL;
  PLLRootedTree &_rootedTree;
  // all nodes, from the root to the most recent speciation followed by leaves
  std::vector<corax_rnode_t *> _orderedSpeciations;
  // ranks for all nodes, parents always have a lower rank than their children
  std::vector<unsigned int> _ranks;
};
