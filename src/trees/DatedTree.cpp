#include "DatedTree.hpp"

#include <cassert>
#include <functional>

#include <maths/Random.hpp>

DatedTree::DatedTree(PLLRootedTree &rootedTree, bool useBLs)
    : _fromBL(useBLs), _rootedTree(rootedTree) {
  // get _orderedSpeciations and _ranks either from tree topology
  // or from branch lengths
  _ranks.resize(_rootedTree.getNodeNumber(), 0);
  updateSpeciationOrderAndRanks();
  // standardize branch lengths either by default or from _ranks
  rescaleBranchLengths();
}

void DatedTree::updateSpeciationOrderAndRanks() {
  // update _orderedSpeciations
  _orderedSpeciations.clear();
  if (!_fromBL) { // order the speciations in reverse postorder
    auto nodes = _rootedTree.getPostOrderNodes();
    for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
      auto node = *it;
      _orderedSpeciations.push_back(node);
    }
  } else { // order the speciations according to branch lengths
    _orderedSpeciations = _rootedTree.getOrderedSpeciations();
  }
  // update _ranks
  unsigned int rank = 0;
  for (auto node : _orderedSpeciations) {
    _ranks[node->node_index] = rank++;
  }
}

void DatedTree::rescaleBranchLengths() {
  // check _ranks
  checkRanks();
  // update branch lengths
  if (!_fromBL) { // set all lengths to a standard value
    _rootedTree.equalizeBranchLengths();
  } else { // set lengths according to _ranks
    double treeHeight = 0.0;
    for (auto node : _orderedSpeciations) {
      if (!node->parent || !node->left) { // the root or a leaf
        node->length = 1.0;               // not included in treeHeight
        continue;
      }
      auto e = node->node_index;
      auto p = node->parent->node_index;
      node->length = double(_ranks[e]) - double(_ranks[p]);
      treeHeight = double(_ranks[e]);
    }
    treeHeight += 1.0; // account for the leaves
    for (auto leaf : _rootedTree.getLeaves()) {
      auto p = leaf->parent->node_index;
      leaf->length = treeHeight - double(_ranks[p]);
    }
  }
}

bool DatedTree::moveUp(unsigned int rank) {
  assert(_fromBL);
  if (rank == 0) {
    return false;
  }
  // move the previous node one rank away from the root
  return moveDown(rank - 1);
}

bool DatedTree::moveDown(unsigned int rank) {
  assert(_fromBL);
  if (rank > _orderedSpeciations.size() - 2) {
    return false;
  }
  auto n1 = _orderedSpeciations[rank];
  auto n2 = _orderedSpeciations[rank + 1];
  // n1 has a lower rank than n2. We want to swap them
  if (!n1->left || !n2->left || n2->parent == n1) {
    return false;
  }
  _orderedSpeciations[rank + 1] = n1;
  _orderedSpeciations[rank] = n2;
  _ranks[n1->node_index]++;
  _ranks[n2->node_index]--;
  return true;
}

void DatedTree::checkRanks() const {
  // check that _ranks are consistent with _orderedSpeciations
  for (unsigned int i = 0; i < _orderedSpeciations.size() - 1; ++i) {
    auto n1 = _orderedSpeciations[i];
    auto n2 = _orderedSpeciations[i + 1];
    assert(_ranks[n1->node_index] + 1 == _ranks[n2->node_index]);
  }
  // check that _ranks are consistent with the tree topology
  for (auto node : _rootedTree.getNodes()) {
    if (node->parent) {
      auto e = node->node_index;
      auto p = node->parent->node_index;
      assert(_ranks[p] < _ranks[e]);
    }
  }
}

void DatedTree::restore(const DatedBackup &backup) {
  _ranks = backup;
  auto speciations = _orderedSpeciations;
  for (auto node : speciations) {
    _orderedSpeciations[_ranks[node->node_index]] = node;
  }
}

// TODO: Too high hash collision rate! Fix somehow before using!
// taken from https://stackoverflow.com/a/27952689
static size_t hash_combine(size_t lhs, size_t rhs) {
  lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  return lhs;
}

size_t DatedTree::getOrderingHash(size_t startingHash) const {
  assert(_fromBL);
  std::hash<size_t> hash_fn;
  size_t hash = hash_fn(startingHash);
  for (auto rank : _ranks) {
    hash = hash_combine(rank, hash);
  }
  return hash;
}

bool DatedTree::canTransferUnderRelDated(unsigned int e, unsigned int d) const {
  // the destination species (d) should be younger than
  // the parent of the source species (e)
  if (d == e) {
    return false;
  }
  auto srcSpeciesNode = _rootedTree.getNode(e);
  if (!srcSpeciesNode->parent) {
    return true;
  }
  auto p = srcSpeciesNode->parent->node_index;
  return _ranks[d] > _ranks[p];
}

void DatedTree::randomize() {
  assert(_fromBL);
  std::vector<corax_rnode_t *> toAdd;
  toAdd.push_back(_rootedTree.getRoot());
  unsigned int currentRank = 0;
  while (toAdd.size()) {
    unsigned int i = Random::getInt() % toAdd.size();
    auto node = toAdd[i];
    if (node->left) {
      _orderedSpeciations[currentRank] = node;
      _ranks[node->node_index] = currentRank;
      currentRank++;
      toAdd[i] = node->left;
      toAdd.push_back(node->right);
    } else {
      toAdd[i] = toAdd.back();
      toAdd.pop_back();
    }
  }
}
