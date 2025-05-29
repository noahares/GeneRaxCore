#include "BaseReconciliationModel.hpp"

#include <functional>

BaseReconciliationModel::BaseReconciliationModel(
    PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMapping,
    const RecModelInfo &recModelInfo)
    : _info(recModelInfo), _speciesTree(speciesTree),
      _geneNameToSpeciesName(geneSpeciesMapping.getMap()),
      _allSpeciesNodesInvalid(true) {
  initSpeciesTree();
  setFractionMissingGenes(_info.fractionMissingFile);
}

void BaseReconciliationModel::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) {
  if (!nodesToInvalidate) {
    _allSpeciesNodesInvalid = true;
  } else {
    assert(nodesToInvalidate->size());
    for (auto speciesNode : *nodesToInvalidate) {
      while (speciesNode) {
        _invalidatedSpeciesNodes.insert(speciesNode);
        speciesNode = speciesNode->parent;
      }
    }
  }
  _allSpeciesNodes.clear();
  fillNodesPostOrder(_speciesTree.getRoot(), _allSpeciesNodes);
  for (auto speciesNode : getAllSpeciesNodes()) {
    auto e = speciesNode->node_index;
    _speciesLeft[e] = speciesNode->left;
    _speciesRight[e] = speciesNode->right;
    _speciesParent[e] = speciesNode->parent;
  }
  _prunedRoot = _speciesTree.getRoot();
  if (prunedMode() && _speciesCoverage.size()) {
    std::vector<corax_rnode_t *> pruned(getAllSpeciesNodeNumber(), nullptr);
    for (auto speciesNode : getAllSpeciesNodes()) {
      auto e = speciesNode->node_index;
      if (!speciesNode->left) {
        pruned[e] = (_speciesCoverage[e] > 0) ? speciesNode : nullptr;
      } else {
        auto left = _speciesLeft[e];
        auto right = _speciesRight[e];
        auto prunedLeft = pruned[left->node_index];
        auto prunedRight = pruned[right->node_index];
        if (prunedLeft && prunedRight) {
          _speciesLeft[e] = prunedLeft;
          _speciesRight[e] = prunedRight;
          pruned[e] = speciesNode;
          _speciesParent[prunedLeft->node_index] = speciesNode;
          _speciesParent[prunedRight->node_index] = speciesNode;
          _prunedRoot = speciesNode;
        } else if (!prunedLeft && prunedRight) {
          pruned[e] = prunedRight;
        } else if (prunedLeft && !prunedRight) {
          pruned[e] = prunedLeft;
        }
      }
    }
  }
  _prunedSpeciesNodes.clear();
  fillPrunedNodesPostOrder(getPrunedRoot(), _prunedSpeciesNodes);
  assert(getAllSpeciesNodeNumber());
  assert(getPrunedSpeciesNodeNumber());
}

void BaseReconciliationModel::initSpeciesTree() {
  // fill the list of the species nodes
  _allSpeciesNodes.clear();
  fillNodesPostOrder(_speciesTree.getRoot(), _allSpeciesNodes);
  assert(getAllSpeciesNodeNumber());
  // build the species tree structure representation
  _speciesLeft =
      std::vector<corax_rnode_t *>(getAllSpeciesNodeNumber(), nullptr);
  _speciesRight =
      std::vector<corax_rnode_t *>(getAllSpeciesNodeNumber(), nullptr);
  _speciesParent =
      std::vector<corax_rnode_t *>(getAllSpeciesNodeNumber(), nullptr);
  onSpeciesTreeChange(nullptr);
  // fill _speciesNameToId
  _speciesNameToId.clear();
  for (auto speciesNode : getAllSpeciesNodes()) {
    if (!speciesNode->left) {
      _speciesNameToId[speciesNode->label] = speciesNode->node_index;
    }
  }
}

void BaseReconciliationModel::fillNodesPostOrder(
    corax_rnode_t *node, std::vector<corax_rnode_t *> &nodes) {
  if (node->left) {
    assert(node->right);
    fillNodesPostOrder(node->left, nodes);
    fillNodesPostOrder(node->right, nodes);
  }
  nodes.push_back(node);
}

void BaseReconciliationModel::fillPrunedNodesPostOrder(
    corax_rnode_t *node, std::vector<corax_rnode_t *> &nodes) {
  if (getSpeciesLeft(node)) {
    assert(getSpeciesRight(node));
    fillPrunedNodesPostOrder(getSpeciesLeft(node), nodes);
    fillPrunedNodesPostOrder(getSpeciesRight(node), nodes);
  }
  nodes.push_back(node);
}

void BaseReconciliationModel::setFractionMissingGenes(
    const std::string &fractionMissingFile) {
  if (!fractionMissingFile.size()) {
    _fm = std::vector<double>(_speciesTree.getLeafNumber(), 0.0);
    return;
  }
  _fm = std::vector<double>(_speciesTree.getLeafNumber(), -1.0);
  std::ifstream is(fractionMissingFile);
  std::string species;
  double fm;
  while (is >> species >> fm) {
    if (_speciesNameToId.find(species) == _speciesNameToId.end()) {
      Logger::error << "Error: species " << species
                    << " from the fraction missing file " << fractionMissingFile
                    << " is not in the species tree" << std::endl;
      assert(false);
    }
    _fm[_speciesNameToId[species]] = fm;
  }
  for (unsigned int i = 0; i < _speciesTree.getLeafNumber(); ++i) {
    if (_fm[i] == -1.0) {
      Logger::error << "Error: the fraction missing file "
                    << fractionMissingFile << " does not cover species "
                    << _speciesTree.getNode(i)->label << std::endl;
      assert(false);
    }
  }
}

size_t BaseReconciliationModel::getSpeciesTreeHash() const {
  if (!_prunedRoot) {
    return 0;
  }
  return getTreeHashRec(_prunedRoot, 0);
}

size_t BaseReconciliationModel::getTreeHashRec(const corax_rnode_t *node,
                                               size_t i) const {
  assert(node);
  std::hash<size_t> hash_fn;
  if (i == 0) {
    i = 1;
  }
  if (!node->left) {
    return hash_fn(node->node_index);
  }
  auto hash1 = getTreeHashRec(_speciesLeft[node->node_index], i + 1);
  auto hash2 = getTreeHashRec(_speciesRight[node->node_index], i + 1);
  auto m = std::min(hash1, hash2);
  auto M = std::max(hash1, hash2);
  auto res = hash_fn(m * i + M);
  res = hash_fn(res * i + node->node_index);
  return res;
}
