#include "util/Scenario.hpp"
#include "util/enums.hpp"
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <IO/ReconciliationWriter.hpp>
#include <map>
#include <sstream>
#include <string>
#include <vector>

const char *Scenario::eventNames[] = {"S",  "SL", "D",    "DL",     "T",
                                      "TL", "L",  "Leaf", "Invalid"};

const unsigned int Scenario::EVENT_TYPE_NUMBER = 9;

struct TransferPair {
  double count;
  unsigned int src;
  unsigned int dest;
  TransferPair(double count, unsigned int src, unsigned int dest)
      : count(count), src(src), dest(dest) {}
  bool operator<(const TransferPair &other) const {
    if (count != other.count) {
      return count < other.count;
    } else if (src != other.src) {
      return src < other.src;
    } else {
      return dest < other.dest;
    }
  }
};

void PerSpeciesEvents::parallelSum() {
  for (auto &speciesEvents : events) {
    ParallelContext::sumUInt(speciesEvents.LeafCount);
    ParallelContext::sumUInt(speciesEvents.DCount);
    ParallelContext::sumUInt(speciesEvents.DLCount);
    ParallelContext::sumUInt(speciesEvents.SCount);
    ParallelContext::sumUInt(speciesEvents.SLCount);
    ParallelContext::sumUInt(speciesEvents.TCount);
    ParallelContext::sumUInt(speciesEvents.TLCount);
  }
}

void Scenario::addEvent(ReconciliationEventType type, unsigned int geneNode,
                        unsigned int speciesNode,
                        unsigned int destSpeciesNode) {

  addTransfer(type, geneNode, speciesNode, destSpeciesNode);
}

void Scenario::addTransfer(ReconciliationEventType type, unsigned int geneNode,
                           unsigned int speciesNode,
                           unsigned int destSpeciesNode) {
  Event event;
  event.type = type;
  event.geneNode = geneNode;
  event.speciesNode = speciesNode;
  event.destSpeciesNode = destSpeciesNode;
  addEvent(event);
}

void Scenario::addEvent(const Event &event) {
  _events.push_back(event);
  assert(static_cast<int>(event.type) >= 0);
  _eventsCount[static_cast<unsigned int>(event.type)]++;
  if (_geneIdToEvents.size() <= static_cast<size_t>(event.geneNode)) {
    _geneIdToEvents.resize(event.geneNode + 1);
  }
  _geneIdToEvents[event.geneNode].push_back(event);
}

void Scenario::writeEventCountsHeader(ParallelOfstream &os, bool useSample) {
  if (useSample) {
    os << "sample\t";
  }
  for (unsigned int i = 0;
       i < static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid);
       ++i) {
    if (i != 0)
      os << "\t";
    os << eventNames[i];
  }
  os << std::endl;
}

void Scenario::saveEventCounts(ParallelOfstream &os, int sample) const {
  if (sample > -1) {
    os << sample << "\t";
  }
  for (unsigned int i = 0;
       i < static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid);
       ++i) {
    if (i != 0)
      os << "\t";
    os << _eventsCount[i];
  }
  os << std::endl;
}

void Scenario::writeTransfersHeader(ParallelOfstream &os, bool useSample) {
  if (useSample) {
    os << "sample\t";
  }
  os << "src\tdest\ttransfers" << std::endl;
}

void Scenario::saveTransfers(ParallelOfstream &os, int sample) const {
  auto N = _speciesTree->tip_count + _speciesTree->inner_count;
  MatrixUint transferCounts(N, VectorUint(N, 0));
  for (const auto &event : _events) {
    if (event.type == ReconciliationEventType::EVENT_T ||
        event.type == ReconciliationEventType::EVENT_TL) {
      transferCounts[event.speciesNode][event.destSpeciesNode]++;
    }
  }
  std::vector<TransferPair> transferPairs;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      double transferCount = transferCounts[i][j];
      if (transferCount > 0.0) {
        transferPairs.push_back(TransferPair(transferCount, i, j));
      }
    }
  }
  std::sort(transferPairs.rbegin(), transferPairs.rend());
  for (const auto &tp : transferPairs) {
    if (sample > -1) {
      os << sample << "\t";
    }
    os << _speciesTree->nodes[tp.src]->label << "\t"
       << _speciesTree->nodes[tp.dest]->label << "\t" << tp.count << std::endl;
  }
}

std::vector<unsigned int> Scenario::getPerSpeciesCopies() const {
  std::vector<unsigned int> copies(_rootedTree->getNodeNumber(), 0);
  for (const auto &event : _events) {
    if (event.type == ReconciliationEventType::EVENT_S ||
        event.type == ReconciliationEventType::EVENT_SL ||
        event.type == ReconciliationEventType::EVENT_None) {
      copies[event.speciesNode]++;
    }
  }
  return copies;
}

void Scenario::gatherReconciliationStatistics(
    PerSpeciesEvents &perSpeciesEvents) const {
  for (auto &event : _events) {
    auto e = event.speciesNode;
    auto &speciesEvents = perSpeciesEvents.events[e];
    switch (event.type) {
    case ReconciliationEventType::EVENT_S:
      speciesEvents.SCount++;
      break;
    case ReconciliationEventType::EVENT_None:
      speciesEvents.LeafCount++;
      break;
    case ReconciliationEventType::EVENT_SL:
      speciesEvents.SLCount++;
      break;
    case ReconciliationEventType::EVENT_D:
      speciesEvents.DCount++;
      break;
    case ReconciliationEventType::EVENT_DL:
      speciesEvents.DLCount++;
      break;
    case ReconciliationEventType::EVENT_T:
      speciesEvents.TCount++;
      break;
    case ReconciliationEventType::EVENT_TL:
      speciesEvents.TLCount++;
      break;
    default:
      assert(false);
    }
  }
}

static void
dumpSpeciesEventCountVector(ParallelOfstream &os,
                            const std::vector<double> &eventCounts,
                            const std::vector<std::string> &indexToLabel) {
  auto N = indexToLabel.size();
  assert(Scenario::EVENT_TYPE_NUMBER * N == eventCounts.size());
  for (unsigned int i = 0; i < N; ++i) {
    bool print = false;
    for (unsigned int e = 0; e < Scenario::EVENT_TYPE_NUMBER; ++e) {
      print |= (eventCounts[i * Scenario::EVENT_TYPE_NUMBER + e] > 0.0);
    }
    if (print) {
      os << indexToLabel[i];
      for (unsigned int e = 0; e < Scenario::EVENT_TYPE_NUMBER; ++e) {
        os << "\t" << eventCounts[i * Scenario::EVENT_TYPE_NUMBER + e];
      }
      os << std::endl;
    }
  }
}

void Scenario::writePerSpeciesEventCountsHeader(ParallelOfstream &os,
                                                bool useSample) {
  if (useSample) {
    os << "sample\t";
  }
  os << "species\tspeciations\tduplications\tlosses\ttransfers\t"
     << "presence\torigination\tcopies\tsingletons\ttransfers_to" << std::endl;
}

void Scenario::savePerSpeciesEventCounts(ParallelOfstream &os,
                                         int sample) const {
  std::unordered_map<std::string, std::vector<double>> speciesToEventCount;
  std::vector<double> defaultCount(EVENT_TYPE_NUMBER, 0);
  auto N = _speciesTree->tip_count + _speciesTree->inner_count;
  for (unsigned int e = 0; e < N; ++e) {
    assert(_speciesTree->nodes[e]->label);
    speciesToEventCount.insert(
        {std::string(_speciesTree->nodes[e]->label), defaultCount});
  }
  // 0: speciations (number of gene copies)
  // 1: duplications
  // 2: losses
  // 3: transfers
  // 4: presence (1 or 0)
  // 5: origination (1 for the origination species, 0 otherwise)
  // 6: gene copies
  // 7: singletons
  // 8: transfers to
  for (const auto &event : _events) {
    auto &eventCount =
        speciesToEventCount.at(_speciesTree->nodes[event.speciesNode]->label);
    switch (event.type) {
    case ReconciliationEventType::EVENT_SL:
      speciesToEventCount.at(event.pllLostSpeciesNode->label)[2]++;
      /* fallthrough */
    case ReconciliationEventType::EVENT_S:
    case ReconciliationEventType::EVENT_None:
      eventCount[0]++;
      eventCount[4] = 1;
      eventCount[6]++;
      if (event.previousType == ReconciliationEventType::EVENT_S ||
          event.previousType == ReconciliationEventType::EVENT_SL) {
        eventCount[7]++;
      }
      break;
    case ReconciliationEventType::EVENT_D:
      eventCount[1]++;
      break;
    case ReconciliationEventType::EVENT_TL:
      eventCount[2]++;
      /* fallthrough */
    case ReconciliationEventType::EVENT_T:
      eventCount[3]++;
      speciesToEventCount.at(
          _speciesTree->nodes[event.destSpeciesNode]->label)[8]++;
      break;
    default:
      // no other event type should ever occur in a Scenario
      assert(false);
    }
  }
  auto originationNode = getOriginationSpecies();
  speciesToEventCount[originationNode->label][5] = 1.0;
  // write to os
  for (const auto &p : speciesToEventCount) {
    if (defaultCount == p.second) {
      // do not write species without any event
      continue;
    }
    assert(p.second.size() == Scenario::EVENT_TYPE_NUMBER);
    if (sample > -1) {
      os << sample << "\t";
    }
    os << p.first;
    for (auto v : p.second) {
      os << "\t" << v;
    }
    os << std::endl;
  }
}

void Scenario::mergeTransfers(const PLLRootedTree &speciesTree,
                              const std::string &fileName,
                              const std::vector<std::string> &fileNames,
                              bool parallel, bool normalize,
                              unsigned int samples) {
  ParallelOfstream os(fileName, parallel);
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int N = idToLabel.size();
  MatrixDouble transferCounts(N, VectorDouble(N, 0.0));
  for (const auto &subfile : fileNames) {
    std::ifstream is(subfile);
    std::string line;
    std::getline(is, line); // skip header
    while (std::getline(is, line)) {
      std::istringstream iss(line);
      if (samples > 0) { // 0 just indicates no sample field
        std::string sample;
        iss >> sample;
      }
      std::string sp1;
      std::string sp2;
      double count = 0.0;
      iss >> sp1 >> sp2 >> count;
      transferCounts[labelToId.at(sp1)][labelToId.at(sp2)] += count;
    }
  }
  std::vector<TransferPair> transferPairs;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      double transferCount = transferCounts[i][j];
      if (parallel) {
        ParallelContext::sumDouble(transferCount);
      }
      if (transferCount > 0.0) {
        transferPairs.push_back(TransferPair(transferCount, i, j));
      }
    }
  }
  if (normalize) {
    unsigned int subfileCount = fileNames.size();
    if (parallel) {
      ParallelContext::sumUInt(subfileCount);
    }
    if (samples == 0) {
      samples = 1;
    }
    for (auto &tp : transferPairs) {
      tp.count /= static_cast<double>(subfileCount * samples);
    }
  }
  std::sort(transferPairs.rbegin(), transferPairs.rend());
  writeTransfersHeader(os);
  for (const auto &tp : transferPairs) {
    os << idToLabel[tp.src] << "\t" << idToLabel[tp.dest] << "\t" << tp.count
       << std::endl;
  }
}

void Scenario::mergePerSpeciesEventCounts(
    const PLLRootedTree &speciesTree, const std::string &fileName,
    const std::vector<std::string> &fileNames, bool parallel, bool normalize,
    unsigned int samples) {
  ParallelOfstream os(fileName, parallel);
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int N = idToLabel.size();
  VectorDouble speciesEventCounts(EVENT_TYPE_NUMBER * N, 0.0);
  for (const auto &subfile : fileNames) {
    std::ifstream is(subfile);
    std::string line;
    std::getline(is, line); // skip header
    while (std::getline(is, line)) {
      std::istringstream iss(line);
      if (samples > 0) { // 0 just indicates no sample field
        std::string sample;
        iss >> sample;
      }
      std::string sp;
      iss >> sp;
      auto i = labelToId.at(sp);
      auto begin = i * EVENT_TYPE_NUMBER;
      for (unsigned int j = begin; j < begin + EVENT_TYPE_NUMBER; ++j) {
        double count = 0.0;
        iss >> count;
        speciesEventCounts[j] += count;
      }
    }
  }
  if (parallel) {
    ParallelContext::sumVectorDouble(speciesEventCounts);
  }
  if (normalize) {
    unsigned int subfileCount = fileNames.size();
    if (parallel) {
      ParallelContext::sumUInt(subfileCount);
    }
    if (samples == 0) {
      samples = 1;
    }
    for (auto &count : speciesEventCounts) {
      count /= static_cast<double>(subfileCount * samples);
    }
  }
  writePerSpeciesEventCountsHeader(os);
  dumpSpeciesEventCountVector(os, speciesEventCounts, idToLabel);
}

void Scenario::saveReconciliation(const std::string &fileName,
                                  ReconciliationFormat format,
                                  bool masterRankOnly) const {
  ParallelOfstream os(fileName, masterRankOnly);
  saveReconciliation(os, format);
}

void Scenario::saveReconciliation(ParallelOfstream &os,
                                  ReconciliationFormat format) const {
  switch (format) {
  case ReconciliationFormat::NHX:
    ReconciliationWriter::saveReconciliationNHX(
        _speciesTree, _geneRoot, _virtualRootIndex, _geneIdToEvents, os);
    break;
  case ReconciliationFormat::ALE:
    ReconciliationWriter::saveReconciliationALE(
        _speciesTree, _geneRoot, _virtualRootIndex, _geneIdToEvents, os);
    break;
  case ReconciliationFormat::RecPhyloXML:
    ReconciliationWriter::saveReconciliationRecPhyloXML(
        _speciesTree, _virtualRootIndex, _geneIdToEvents, os);
    break;
  case ReconciliationFormat::NewickEvents:
    ReconciliationWriter::saveReconciliationNewickEvents(
        _geneRoot, _virtualRootIndex, _geneIdToEvents, os);
    break;
  }
}

void Scenario::countTransfers(const StringToUint &labelToId,
                              MatrixUint &transferCounts) const {
  char *labelFrom = nullptr;
  char *labelTo = nullptr;
  for (const auto &event : _events) {
    switch (event.type) {
    case ReconciliationEventType::EVENT_T:
    case ReconciliationEventType::EVENT_TL:
      labelFrom = _speciesTree->nodes[event.speciesNode]->label;
      labelTo = _speciesTree->nodes[event.destSpeciesNode]->label;
      transferCounts[labelToId.at(labelFrom)][labelToId.at(labelTo)]++;
      break;
    default:
      break;
    }
  }
}

void Scenario::countOrigins(const StringToUint &labelToId, VectorUint &fromO,
                            VectorUint &fromS, VectorUint &fromSButL,
                            MatrixUint &transferCounts) const {
  // origins come from O, S, SL, T, and TL events
  char *origLabel = getOriginationSpecies()->label;
  fromO[labelToId.at(origLabel)]++;
  char *destLabel1 = nullptr;
  char *destLabel2 = nullptr;
  for (const auto &event : _events) {
    switch (event.type) {
    case ReconciliationEventType::EVENT_S:
      destLabel1 = _speciesTree->nodes[event.speciesNode]->left->label;
      destLabel2 = _speciesTree->nodes[event.speciesNode]->right->label;
      fromS[labelToId.at(destLabel1)]++;
      fromS[labelToId.at(destLabel2)]++;
      break;
    case ReconciliationEventType::EVENT_SL:
      destLabel1 = event.pllDestSpeciesNode->label;
      destLabel2 = event.pllLostSpeciesNode->label;
      fromS[labelToId.at(destLabel1)]++;
      fromSButL[labelToId.at(destLabel2)]++;
      break;
    default:
      break;
    }
  }
  countTransfers(labelToId, transferCounts);
}

unsigned int Scenario::countTransfers(const std::string &from,
                                      const std::string &to) const {
  unsigned int transfers = 0;
  for (const auto &event : _events) {
    if (event.type == ReconciliationEventType::EVENT_T ||
        event.type == ReconciliationEventType::EVENT_TL) {
      if (_speciesTree->nodes[event.speciesNode]->label == from &&
          _speciesTree->nodes[event.destSpeciesNode]->label == to) {
        transfers++;
      }
    }
  }
  return transfers;
}

void Scenario::saveLargestOrthoGroup(const std::string &fileName,
                                     bool masterRankOnly) const {
  ParallelOfstream os(fileName, masterRankOnly);
  corax_unode_t virtualRoot;
  virtualRoot.next = _geneRoot;
  virtualRoot.node_index = _virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  std::unique_ptr<OrthoGroup> orthoGroup(
      getLargestOrthoGroupRec(&virtualRoot, true));
  for (auto value : *orthoGroup) {
    os << value << std::endl;
  }
  os << "-" << std::endl;
}

void Scenario::saveAllOrthoGroups(const std::string &fileName,
                                  bool masterRankOnly) const {
  ParallelOfstream os(fileName, masterRankOnly);
  corax_unode_t virtualRoot;
  virtualRoot.next = _geneRoot;
  virtualRoot.node_index = _virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  OrthoGroupPtr currentOrthoGroup = std::make_shared<OrthoGroup>();
  OrthoGroups orthoGroups;
  getAllOrthoGroupRec(&virtualRoot, orthoGroups, currentOrthoGroup, true);
  orthoGroups.push_back(currentOrthoGroup);
  for (auto orthoGroup : orthoGroups) {
    for (auto value : *orthoGroup) {
      os << value << std::endl;
    }
    os << "-" << std::endl;
  }
}

void Scenario::saveTransfersGlobal(const PLLRootedTree &speciesTree,
                                   const MatrixUint &transferCounts,
                                   const std::string &fileName,
                                   unsigned int samples) {
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int N = idToLabel.size();
  std::vector<TransferPair> transferPairs;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      double transferCount = transferCounts[i][j];
      ParallelContext::sumDouble(transferCount);
      if (transferCount > 0.0) {
        if (samples > 1) {
          transferCount /= static_cast<double>(samples);
        }
        transferPairs.push_back(TransferPair(transferCount, i, j));
      }
    }
  }
  unsigned int parentTransfers = 0;
  for (auto speciesNode : speciesTree.getNodes()) {
    auto from = labelToId.at(speciesNode->label);
    auto parentNode = speciesNode;
    while (parentNode) {
      auto to = labelToId.at(parentNode->label);
      parentTransfers += transferCounts[from][to];
      parentNode = parentNode->parent;
    }
  }
  ParallelContext::sumUInt(parentTransfers);
  std::sort(transferPairs.rbegin(), transferPairs.rend());
  ParallelOfstream os(fileName);
  writeTransfersHeader(os);
  for (const auto &tp : transferPairs) {
    os << idToLabel[tp.src] << "\t" << idToLabel[tp.dest] << "\t" << tp.count
       << std::endl;
  }
}

void Scenario::saveOriginsGlobal(const PLLRootedTree &speciesTree,
                                 const VectorUint &fromO,
                                 const VectorUint &fromS,
                                 const VectorUint &fromSButL,
                                 const MatrixUint &transferCounts,
                                 const std::string &outputDir,
                                 unsigned int samples) {
  const auto labelToNode = speciesTree.getLabelToNode(false);
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int N = idToLabel.size();
  // iterate over all dest species and compute their origins
  for (unsigned int j = 0; j < N; ++j) {
    // origins from O events
    auto noLabel = std::string("-");
    double originationCount = fromO[j];
    ParallelContext::sumDouble(originationCount);
    if (samples > 1) {
      originationCount /= static_cast<double>(samples);
    }
    // origins from S and SL events
    auto speciesLabel = idToLabel[j];
    auto parentNode = labelToNode.at(speciesLabel)->parent;
    auto parentLabel = (parentNode) ? std::string(parentNode->label) : noLabel;
    double verticalCount = fromS[j] + fromSButL[j];
    ParallelContext::sumDouble(verticalCount);
    if (samples > 1) {
      verticalCount /= static_cast<double>(samples);
    }
    // origins from T and TL events
    std::vector<TransferPair> transferPairs;
    for (unsigned int i = 0; i < N; ++i) { // all transfer sources
      double transferCount = transferCounts[i][j];
      ParallelContext::sumDouble(transferCount);
      if (transferCount > 0.0) {
        if (samples > 1) {
          transferCount /= static_cast<double>(samples);
        }
        transferPairs.push_back(TransferPair(transferCount, i, j));
      }
    }
    std::sort(transferPairs.rbegin(), transferPairs.rend());
    // save origins for the current species
    auto fileName = FileSystem::joinPaths(outputDir, speciesLabel + ".tsv");
    ParallelOfstream os(fileName);
    os << "src\ttype\torigins" << std::endl;
    os << noLabel << "\torigination\t" << originationCount << std::endl;
    os << parentLabel << "\tvertical\t" << verticalCount << std::endl;
    for (const auto &tp : transferPairs) {
      os << idToLabel[tp.src] << "\thorizontal\t" << tp.count << std::endl;
    }
  }
}

OrthoGroup *Scenario::getLargestOrthoGroupRec(corax_unode_t *geneNode,
                                              bool isVirtualRoot) const {
  auto &events = _geneIdToEvents[geneNode->node_index];
  for (auto &event : events) {
    if (event.type == ReconciliationEventType::EVENT_TL) {
      return new OrthoGroup();
    }
  }
  if (geneNode->next == nullptr) {
    auto res = new OrthoGroup();
    res->insert(std::string(geneNode->label));
    return res;
  } else {
    auto left = geneNode->next->back;
    auto right = geneNode->next->next->back;
    if (isVirtualRoot) {
      left = geneNode->next;
      right = geneNode->next->back;
    }
    auto leftOrthoGroup = getLargestOrthoGroupRec(left, false);
    auto rightOrthoGroup = getLargestOrthoGroupRec(right, false);
    auto &event = events.back();
    switch (event.type) {
    case ReconciliationEventType::EVENT_S:
      // keep both groups
      for (auto &v : *leftOrthoGroup) {
        rightOrthoGroup->insert(v);
      }
      delete leftOrthoGroup;
      return rightOrthoGroup;
    case ReconciliationEventType::EVENT_T:
      // only keep the non transfered gene
      if (event.rightGeneIndex == left->node_index) {
        delete leftOrthoGroup;
        return rightOrthoGroup;
      } else if (event.rightGeneIndex == right->node_index) {
        delete rightOrthoGroup;
        return leftOrthoGroup;
      } else {
        assert(false);
      }
    case ReconciliationEventType::EVENT_TL:
    case ReconciliationEventType::EVENT_DL:
    case ReconciliationEventType::EVENT_None:
    case ReconciliationEventType::EVENT_SL:
    case ReconciliationEventType::EVENT_L:
      // we already handled these cases
      // or they should not be the last event
      // attached to this gene node
      assert(false);
      break;
    case ReconciliationEventType::EVENT_Invalid:
      delete leftOrthoGroup;
      delete rightOrthoGroup;
      break;
    case ReconciliationEventType::EVENT_D:
      // keep the bigger group!
      if (leftOrthoGroup->size() < rightOrthoGroup->size()) {
        std::swap(leftOrthoGroup, rightOrthoGroup);
      }
      delete rightOrthoGroup;
      return leftOrthoGroup;
    }
    return new OrthoGroup();
  }
}

static void appendOrtho(OrthoGroupPtr &orthoGroups,
                        const OrthoGroupPtr &orthoGroupsToAppend) {
  orthoGroups->insert(orthoGroupsToAppend->begin(), orthoGroupsToAppend->end());
}

void Scenario::getAllOrthoGroupRec(corax_unode_t *geneNode,
                                   OrthoGroups &orthoGroups,
                                   OrthoGroupPtr &currentOrthoGroup,
                                   bool isVirtualRoot) const {
  auto &events = _geneIdToEvents[geneNode->node_index];
  bool underTL = false;
  for (auto &event : events) {
    if (event.type == ReconciliationEventType::EVENT_TL) {
      underTL = true;
    }
  }
  if (geneNode->next == nullptr) {
    currentOrthoGroup->insert(std::string(geneNode->label));
  } else {
    auto left = geneNode->next->back;
    auto right = geneNode->next->next->back;
    if (isVirtualRoot) {
      left = geneNode->next;
      right = geneNode->next->back;
    }
    OrthoGroupPtr leftOrthoGroup = std::make_shared<OrthoGroup>();
    OrthoGroupPtr rightOrthoGroup = std::make_shared<OrthoGroup>();
    getAllOrthoGroupRec(left, orthoGroups, leftOrthoGroup, false);
    getAllOrthoGroupRec(right, orthoGroups, rightOrthoGroup, false);
    auto &event = events.back();
    switch (event.type) {
    case ReconciliationEventType::EVENT_S:
      // merge both groups
      appendOrtho(currentOrthoGroup, leftOrthoGroup);
      appendOrtho(currentOrthoGroup, rightOrthoGroup);
      break;
    case ReconciliationEventType::EVENT_T:
      // save the orthoGroup from the transfered gene,
      // and keep filling the orthoGroup on the non transfered gene
      if (event.rightGeneIndex == left->node_index) {
        appendOrtho(currentOrthoGroup, rightOrthoGroup);
        orthoGroups.push_back(leftOrthoGroup);
      } else if (event.rightGeneIndex == right->node_index) {
        appendOrtho(currentOrthoGroup, leftOrthoGroup);
        orthoGroups.push_back(rightOrthoGroup);
      } else {
        assert(false);
      }
      break;
    case ReconciliationEventType::EVENT_TL:
    case ReconciliationEventType::EVENT_DL:
    case ReconciliationEventType::EVENT_None:
    case ReconciliationEventType::EVENT_SL:
    case ReconciliationEventType::EVENT_L:
      // we already handled these cases
      // or they should not be the last event
      // attached to this gene node
      assert(false);
      break;
    case ReconciliationEventType::EVENT_Invalid:
      orthoGroups.push_back(leftOrthoGroup);
      orthoGroups.push_back(rightOrthoGroup);
      break;
    case ReconciliationEventType::EVENT_D:
      // save the biggest group and continue working on the other
      if (leftOrthoGroup->size() < rightOrthoGroup->size()) {
        std::swap(leftOrthoGroup, rightOrthoGroup);
      }
      orthoGroups.push_back(rightOrthoGroup);
      appendOrtho(currentOrthoGroup, leftOrthoGroup);
      break;
    }
  }
  if (underTL) {
    // we got transfered: separate the current orthogroup
    orthoGroups.push_back(currentOrthoGroup);
    currentOrthoGroup->clear();
  }
}

void Scenario::initBlackList(unsigned int genesNumber,
                             unsigned int speciesNumber) {
  std::vector<int> blackListAux(speciesNumber, false);
  _blacklist = std::make_unique<ScenarioBlackList>(genesNumber, blackListAux);
}

void Scenario::blackList(unsigned int geneNode, unsigned int speciesNode) {
  if (_blacklist &&
      _blacklist->size() > geneNode) { // not true for virtual nodes
    (*_blacklist)[geneNode][speciesNode] = true;
  }
}
bool Scenario::isBlacklisted(unsigned int geneNode, unsigned int speciesNode) {
  if (_blacklist &&
      _blacklist->size() > geneNode) { // not true for virtual nodes
    return (*_blacklist)[geneNode][speciesNode];
  }
  return false;
}

void Scenario::resetBlackList() {
  for (auto &aux : *_blacklist) {
    for (unsigned int i = 0; i < aux.size(); ++i) {
      aux[i] = false;
    }
  }
}

static corax_unode_t *createNode(unsigned int index) {
  auto node = new corax_unode_t();
  node->next = nullptr;
  node->back = nullptr;
  node->data = nullptr;
  node->label = nullptr;
  node->node_index = index;
  return node;
}

static void backLink(corax_unode_t *n1, corax_unode_t *n2) {
  n1->back = n2;
  n2->back = n1;
}

static void nextLink(corax_unode_t *n1, corax_unode_t *n2, corax_unode_t *n3) {
  n1->next = n2;
  n2->next = n3;
  n3->next = n1;
}

corax_unode_t *Scenario::generateGeneRoot() {
  auto node1 = createNode(0);
  _geneNodeBuffer.push_back(node1);
  return node1;
}

void Scenario::generateGeneChildren(corax_unode_t *geneNode,
                                    corax_unode_t *&leftGeneNode,
                                    corax_unode_t *&rightGeneNode) {
  if (_geneNodeBuffer.size() == 1) {
    // we are at the root
    leftGeneNode = _geneNodeBuffer.back();
    rightGeneNode = createNode(_geneNodeBuffer.size());
    backLink(leftGeneNode, rightGeneNode);
    _geneNodeBuffer.push_back(rightGeneNode);
    return;
  }
  auto next1 = createNode(geneNode->node_index);
  auto next2 = createNode(geneNode->node_index);
  nextLink(geneNode, next1, next2);
  leftGeneNode = createNode(_geneNodeBuffer.size());
  backLink(next1, leftGeneNode);
  _geneNodeBuffer.push_back(leftGeneNode);
  rightGeneNode = createNode(_geneNodeBuffer.size());
  backLink(next2, rightGeneNode);
  _geneNodeBuffer.push_back(rightGeneNode);
}

corax_rnode_t *Scenario::getOriginationSpecies() const {
  const auto &rootEvents = _geneIdToEvents[_virtualRootIndex];
  assert(rootEvents.size());
  return _speciesTree->nodes[rootEvents[0].speciesNode];
}
