#include "SpeciesTree.hpp"

#include <algorithm>
#include <cassert>

#include <IO/GeneSpeciesMapping.hpp>
#include <parallelization/ParallelContext.hpp>

static std::unordered_set<std::string>
getLabelsFromFamilies(const Families &families) {
  std::unordered_set<std::string> leafLabels;
  for (const auto &family : families) {
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    for (const auto &p : mapping.getMap()) {
      leafLabels.insert(p.second);
    }
  }
  return leafLabels;
}

SpeciesTree::SpeciesTree(const std::string &str, bool isFile, bool useBLs)
    : _speciesTree(str, isFile), _datedTree(_speciesTree, useBLs) {}

SpeciesTree::SpeciesTree(const std::unordered_set<std::string> &labels)
    : _speciesTree(labels), _datedTree(_speciesTree, false) {}

SpeciesTree::SpeciesTree(const Families &families)
    : _speciesTree(getLabelsFromFamilies(families)),
      _datedTree(_speciesTree, false) {}

std::unique_ptr<SpeciesTree> SpeciesTree::buildRandomTree() const {
  return std::make_unique<SpeciesTree>(_speciesTree.getLabels(true));
}

void SpeciesTree::saveToFile(const std::string &fileName,
                             bool masterRankOnly) const {
  if (masterRankOnly && ParallelContext::getRank()) {
    return;
  }
  _speciesTree.save(fileName);
}

void SpeciesTree::getLabelToId(StringToUint &labelToId) const {
  labelToId.clear();
  for (auto node : getTree().getNodes()) {
    labelToId[std::string(node->label)] = node->node_index;
  }
}

void SpeciesTree::addListener(Listener *listener) {
  _listeners.push_back(listener);
}

void SpeciesTree::removeListener(Listener *listener) {
  _listeners.erase(std::remove(_listeners.begin(), _listeners.end(), listener),
                   _listeners.end());
}

void SpeciesTree::onSpeciesDatesChange() {
  _datedTree.rescaleBranchLengths(); // update branch lengths
  for (auto listener : _listeners) {
    listener->onSpeciesDatesChange();
  }
}

void SpeciesTree::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  _datedTree.rescaleBranchLengths();                   // update branch lengths
  _speciesTree.onSpeciesTreeChange(nodesToInvalidate); // update labels and lcas
  for (auto listener : _listeners) {
    listener->onSpeciesTreeChange(nodesToInvalidate);
  }
}

void SpeciesTreeOperator::restoreDates(SpeciesTree &speciesTree,
                                       const DatedBackup &backup) {
  speciesTree.getDatedTree().restore(backup);
  speciesTree.onSpeciesDatesChange();
}

static void setRootAux(SpeciesTree &speciesTree, corax_rnode_t *root) {
  speciesTree.getTree().getRawPtr()->root = root;
  root->parent = nullptr;
}

// direction: 0 == RR, 1 == LR, 2 == RL, 3 == LL
bool SpeciesTreeOperator::canChangeRoot(const SpeciesTree &speciesTree,
                                        unsigned int direction) {
  assert(direction < 4);
  auto root = speciesTree.getRoot();
  auto newRoot = (direction % 2) ? root->left : root->right;
  return newRoot->left && newRoot->right;
}

void SpeciesTreeOperator::changeRoot(SpeciesTree &speciesTree,
                                     unsigned int direction) {
  assert(canChangeRoot(speciesTree, direction));
  std::unordered_set<corax_rnode_t *> nodesToInvalidate;
  // reroot
  auto root = speciesTree.getRoot();
  auto rootLeft = root->left;
  auto rootRight = root->right;
  auto A = rootLeft->left;
  auto B = rootLeft->right;
  auto C = rootRight->left;
  auto D = rootRight->right;
  auto newRoot = (direction % 2) ? rootLeft : rootRight;
  nodesToInvalidate.insert(root);
  setRootAux(speciesTree, newRoot);
  switch (direction) {
  case 0: // RR
    PLLRootedTree::setSon(rootRight, root, true);
    PLLRootedTree::setSon(rootRight, D, false);
    PLLRootedTree::setSon(root, rootLeft, true);
    PLLRootedTree::setSon(root, C, false);
    break;
  case 1: // LR
    PLLRootedTree::setSon(rootLeft, B, true);
    PLLRootedTree::setSon(rootLeft, root, false);
    PLLRootedTree::setSon(root, rootRight, true);
    PLLRootedTree::setSon(root, A, false);
    break;
  case 2: // RL
    PLLRootedTree::setSon(rootRight, root, true);
    PLLRootedTree::setSon(rootRight, C, false);
    PLLRootedTree::setSon(root, D, true);
    PLLRootedTree::setSon(root, rootLeft, false);
    break;
  case 3: // LL
    PLLRootedTree::setSon(rootLeft, A, true);
    PLLRootedTree::setSon(rootLeft, root, false);
    PLLRootedTree::setSon(root, B, true);
    PLLRootedTree::setSon(root, rootRight, false);
    break;
  default:
    assert(false);
  }
  // update ranks, branches, and labels
  auto &datedTree = speciesTree.getDatedTree();
  if (datedTree.isDated()) {
    while (datedTree.moveUp(datedTree.getRank(newRoot))); // move newRoot rank
  } else {
    datedTree.updateSpeciationOrderAndRanks(); // get ranks from topology
  }
  speciesTree.onSpeciesTreeChange(&nodesToInvalidate);
}

void SpeciesTreeOperator::revertChangeRoot(SpeciesTree &speciesTree,
                                           unsigned int direction) {
  changeRoot(speciesTree, 3 - direction);
}

static corax_rnode_t *getBrother(corax_rnode_t *node) {
  auto father = node->parent;
  if (father) {
    return (father->left == node) ? father->right : father->left;
  }
  return nullptr;
}

bool SpeciesTreeOperator::canApplySPRMove(SpeciesTree &speciesTree,
                                          unsigned int prune,
                                          unsigned int regraft) {
  auto pruneNode = speciesTree.getNode(prune);
  auto pruneFatherNode = pruneNode->parent;
  auto pruneBrotherNode = getBrother(pruneNode);
  auto regraftNode = speciesTree.getNode(regraft);
  // check prune is not the root
  if (!pruneFatherNode) {
    return false;
  }
  // check regraft is not adjacent to prune
  if (regraftNode == pruneFatherNode || regraftNode == pruneBrotherNode) {
    return false;
  }
  // check regraft is not a prune's descendant (or prune itself)
  // this should be fast, as the ancestors are cached
  if (speciesTree.getTree().isAncestorOf(prune, regraft)) {
    return false;
  }
  return true;
}

unsigned int SpeciesTreeOperator::applySPRMove(SpeciesTree &speciesTree,
                                               unsigned int prune,
                                               unsigned int regraft) {
  assert(canApplySPRMove(speciesTree, prune, regraft));
  std::unordered_set<corax_rnode_t *> nodesToInvalidate;
  // prune
  auto pruneNode = speciesTree.getNode(prune);
  auto pruneFatherNode = pruneNode->parent;
  assert(pruneFatherNode);
  auto pruneGrandFatherNode = pruneFatherNode->parent;
  auto pruneBrotherNode = getBrother(pruneNode);
  nodesToInvalidate.insert(pruneFatherNode);
  if (!pruneGrandFatherNode) { // prune is at the root
    setRootAux(speciesTree, pruneBrotherNode);
  } else {
    nodesToInvalidate.insert(pruneGrandFatherNode);
    PLLRootedTree::setSon(pruneGrandFatherNode, pruneBrotherNode,
                          pruneGrandFatherNode->left == pruneFatherNode);
  }
  // regraft
  auto regraftNode = speciesTree.getNode(regraft);
  auto regraftParentNode = regraftNode->parent;
  if (!regraftParentNode) { // regraft is the root
    setRootAux(speciesTree, pruneFatherNode);
    PLLRootedTree::setSon(pruneFatherNode, regraftNode,
                          pruneFatherNode->left != pruneNode);
  } else {
    nodesToInvalidate.insert(regraftParentNode);
    PLLRootedTree::setSon(regraftParentNode, pruneFatherNode,
                          regraftParentNode->left == regraftNode);
    PLLRootedTree::setSon(pruneFatherNode, regraftNode,
                          pruneFatherNode->left != pruneNode);
  }
  // update ranks, branches, and labels
  auto &datedTree = speciesTree.getDatedTree();
  assert(!datedTree.isDated());
  datedTree.updateSpeciationOrderAndRanks(); // get ranks from topology
  speciesTree.onSpeciesTreeChange(&nodesToInvalidate);
  // return info for a rollback
  return pruneBrotherNode->node_index;
}

void SpeciesTreeOperator::reverseSPRMove(SpeciesTree &speciesTree,
                                         unsigned int prune,
                                         unsigned int applySPRMoveReturnValue) {
  applySPRMove(speciesTree, prune, applySPRMoveReturnValue);
}

void SpeciesTreeOperator::getPossiblePrunes(
    const SpeciesTree &speciesTree, const std::vector<double> &supportValues,
    double maxSupport, std::vector<unsigned int> &prunes) {
  // we can take potentially all nodes (even the root), as the root
  // may be changed during a SPR search round
  for (auto node : speciesTree.getTree().getNodes()) {
    auto parentNode = node->parent;
    if (supportValues.size() && parentNode) {
      double parentSupport = supportValues[parentNode->node_index];
      if (parentSupport > maxSupport) {
        continue; // but we do not want to break well-supported clades
      }
    }
    prunes.push_back(node->node_index);
  }
}

// direction: 0 == from parent, 1 == from left, 2 == from right
static void recursiveGetNodes(corax_rnode_t *node, unsigned int direction,
                              unsigned int radius,
                              std::vector<unsigned int> &nodes, bool addNode) {
  if (!node) {
    return;
  }
  if (addNode) {
    nodes.push_back(node->node_index);
  }
  if (radius) {
    radius--;
    switch (direction) {
    case 0:
      recursiveGetNodes(node->left, 0, radius, nodes, true);
      recursiveGetNodes(node->right, 0, radius, nodes, true);
      break;
    case 1:
    case 2:
      recursiveGetNodes((direction == 1) ? node->right : node->left, 0, radius,
                        nodes, true);
      if (node->parent) {
        recursiveGetNodes(node->parent, (node->parent->left == node) ? 1 : 2,
                          radius, nodes, true);
      }
      break;
    default:
      assert(false);
    }
  }
}

// Warning:
// This function performs a modified regraft-sampling algorithm, optimized
// for the exhaustive SPR search that uses all possible prunes; compare:
// - standard radius-1 SPR:
// Takes as possible regrafts prune's
// GrandFather, Uncle, and Brother's children (4 nodes in total).
// - modified radius-1 SPR:
// Takes as possible regrafts prune's
// GreatGrandFather, Uncle, and Brother's grandchildren (6 nodes in total).
// Of the standard radius-1 regrafts, only prune's Uncle is taken, while the
// other regrafts are omitted to ensure that moves are not repeated, as the
// same topologies can be obtained from different prunes, namely:
// prune -> Brother's child1 = Brother's child1 -> prune (which is his uncle)
// prune -> Brother's child2 = Brother's child2 -> prune (which is his uncle)
// prune -> GrandFather = Brother -> Uncle (which is his uncle too)
void SpeciesTreeOperator::getPossibleRegrafts(
    const SpeciesTree &speciesTree, unsigned int prune, unsigned int radius,
    std::vector<unsigned int> &regrafts) {
  auto pruneNode = speciesTree.getNode(prune);
  auto pruneFatherNode = pruneNode->parent;
  if (!pruneFatherNode) {
    return;
  }
  auto pruneGrandFatherNode = pruneFatherNode->parent;
  if (pruneGrandFatherNode) {
    recursiveGetNodes(pruneGrandFatherNode,
                      (pruneGrandFatherNode->left == pruneFatherNode) ? 1 : 2,
                      radius, regrafts, false);
  }
  auto pruneBrotherNode = getBrother(pruneNode);
  recursiveGetNodes(pruneBrotherNode->left, 0, radius, regrafts, false);
  recursiveGetNodes(pruneBrotherNode->right, 0, radius, regrafts, false);
}

void SpeciesTreeOperator::getAffectedBranches(
    SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft,
    std::vector<unsigned int> &affectedBranches) {
  // we return the nodes between p and r (excluding p, r, and their LCA)
  // this should be fast, as the LCAs are cached
  auto lca = speciesTree.getTree().getLCA(prune, regraft);
  auto p = speciesTree.getNode(prune);
  auto r = speciesTree.getNode(regraft);
  if (p != lca) {
    p = p->parent; // skip the prune node
  }
  if (r != lca) {
    r = r->parent; // skip the regraft node
  }
  while (p != lca) {
    affectedBranches.push_back(p->node_index);
    p = p->parent;
  }
  while (r != lca) {
    affectedBranches.push_back(r->node_index);
    r = r->parent;
  }
}
