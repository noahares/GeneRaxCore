#pragma once

#include <corax/corax.h>
#include <memory>
#include <string>
#include <trees/PLLRootedTree.hpp>
#include <unordered_set>
#include <util/enums.hpp>
#include <util/types.hpp>
#include <vector>

/*
typedef struct corax_utree_s corax_utree_t;
typedef struct corax_unode_s corax_unode_t;
typedef struct corax_rtree_s corax_rtree_t;
typedef struct corax_rnode_s corax_rnode_t;
*/
class ParallelOfstream;
typedef std::unordered_set<std::string> OrthoGroup;
typedef std::shared_ptr<OrthoGroup> OrthoGroupPtr;
typedef std::vector<OrthoGroupPtr> OrthoGroups;

struct SpeciesEvents {
  unsigned int LeafCount;
  unsigned int DCount;
  unsigned int DLCount;
  unsigned int SCount;
  unsigned int SLCount;
  unsigned int TCount;
  unsigned int TLCount;
  SpeciesEvents()
      : LeafCount(0), DCount(0), DLCount(0), SCount(0), SLCount(0), TCount(0),
        TLCount(0) {}

  double speciesFrequency() const { return SCount + SLCount + LeafCount; }
};

struct PerSpeciesEvents {
  std::vector<SpeciesEvents> events;
  PerSpeciesEvents(unsigned int speciesNodeCount = 0)
      : events(speciesNodeCount) {}
  void parallelSum();
};

/*
 * Store the set of events that reconciles a gene tree with a species tree
 */
class Scenario {
public:
  // value representing an invalid node ID (0 is a valid node ID)
  static const unsigned int INVALID_NODE_ID = static_cast<unsigned int>(-1);

  /**
   *  Reconciliation event description: type of event, and
   *  involved gene and species nodes
   */
  struct Event {
    Event()
        : type(ReconciliationEventType::EVENT_Invalid),
          previousType(ReconciliationEventType::EVENT_Invalid),
          geneNode(INVALID_NODE_ID), speciesNode(INVALID_NODE_ID),
          destSpeciesNode(INVALID_NODE_ID), lostSpeciesNode(INVALID_NODE_ID),
          leftGeneIndex(INVALID_NODE_ID), rightGeneIndex(INVALID_NODE_ID),
          pllTransferedGeneNode(nullptr), pllDestSpeciesNode(nullptr),
          pllLostSpeciesNode(nullptr) {}
    ReconciliationEventType type;
    ReconciliationEventType previousType;
    unsigned int geneNode;
    unsigned int speciesNode;
    unsigned int destSpeciesNode; // for T and TL events only
    unsigned int lostSpeciesNode; // for SL events only

    // left and gene indices:
    // - for speciation: left gene goes to left species
    // - for transfer: left gene goes to source species and right gene
    // goes to recieving species
    // - for duplication, the order doesn't matter
    // - for other events, those members are irrelevant
    unsigned int leftGeneIndex;
    unsigned int rightGeneIndex;

    // temporary variables for event inference
    corax_unode_t *pllTransferedGeneNode;
    corax_rnode_t *pllDestSpeciesNode;
    corax_rnode_t *pllLostSpeciesNode;

    std::string label;

    /**
     *  Is this event valid? (if not, something went wrong)
     */
    bool isValid() const { return speciesNode != INVALID_NODE_ID; }
    bool isLeaf() const { return type == ReconciliationEventType::EVENT_None; }
  };

  /**
   * Default constructor
   */
  Scenario()
      : _lastEventType(ReconciliationEventType::EVENT_S),
        _eventsCount(
            static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid),
            0),
        _geneRoot(nullptr), _rootedTree(nullptr), _speciesTree(nullptr),
        _virtualRootIndex(INVALID_NODE_ID) {}

  // forbid copy
  Scenario(const Scenario &) = delete;
  Scenario &operator=(const Scenario &) = delete;
  Scenario(Scenario &&) = delete;
  Scenario &operator=(Scenario &&) = delete;

  void setGeneRoot(corax_unode_t *geneRoot) { _geneRoot = geneRoot; }
  void setSpeciesTree(PLLRootedTree *speciesTree) {
    _rootedTree = speciesTree;
    _speciesTree = _rootedTree->getRawPtr();
  }
  void setVirtualRootIndex(unsigned int virtualRootIndex) {
    _virtualRootIndex = virtualRootIndex;
  }

  /**
   * Various methods to add an event in the Scnenario
   */
  void addEvent(const Event &event);
  void addEvent(ReconciliationEventType type, unsigned int geneNode,
                unsigned int speciesNode,
                unsigned int destSpeciesNode = INVALID_NODE_ID);
  void addTransfer(ReconciliationEventType type, unsigned int geneNode,
                   unsigned int speciesNode, unsigned int destSpeciesNode);

  /**
   * Various methods to save information from the Scenario
   */
  void saveEventCounts(ParallelOfstream &os, int sample = -1) const;
  static void writeEventCountsHeader(ParallelOfstream &os,
                                     bool useSample = false);
  void saveTransfers(ParallelOfstream &os, int sample = -1) const;
  static void writeTransfersHeader(ParallelOfstream &os,
                                   bool useSample = false);
  void savePerSpeciesEventCounts(ParallelOfstream &os, int sample = -1) const;
  static void writePerSpeciesEventCountsHeader(ParallelOfstream &os,
                                               bool useSample = false);

  static void mergeTransfers(const PLLRootedTree &speciesTree,
                             const std::string &fileName,
                             const std::vector<std::string> &fileNames,
                             bool parallel, bool normalize,
                             unsigned int samples = 0);
  static void mergePerSpeciesEventCounts(
      const PLLRootedTree &speciesTree, const std::string &fileName,
      const std::vector<std::string> &fileNames, bool parallel, bool normalize,
      unsigned int samples = 0);

  void saveReconciliation(const std::string &fileName,
                          ReconciliationFormat format,
                          bool masterRankOnly = true) const;
  void saveReconciliation(ParallelOfstream &os,
                          ReconciliationFormat format) const;

  static void saveTransfersGlobal(const PLLRootedTree &speciesTree,
                                  const MatrixUint &transferCounts,
                                  const std::string &fileName,
                                  unsigned int samples = 0);
  static void saveOriginsGlobal(const PLLRootedTree &speciesTree,
                                const VectorUint &fromS,
                                const VectorUint &fromSButL,
                                const MatrixUint &transferCounts,
                                const std::string &outputDir,
                                unsigned int samples = 0);

  void saveLargestOrthoGroup(const std::string &fileName,
                             bool masterRankOnly = true) const;
  void saveAllOrthoGroups(const std::string &fileName,
                          bool masterRankOnly = true) const;

  void gatherReconciliationStatistics(PerSpeciesEvents &perSpeciesEvents) const;

  void countTransfers(const StringToUint &labelToId,
                      MatrixUint &transferCounts) const;
  void countOrigins(const StringToUint &labelToId, VectorUint &fromS,
                    VectorUint &fromSButL, MatrixUint &transferCounts) const;

  unsigned int countTransfers(const std::string &from,
                              const std::string &to) const;

  std::vector<unsigned int> getPerSpeciesCopies() const;

  unsigned int getEventCount(ReconciliationEventType type) const {
    return _eventsCount[static_cast<unsigned int>(type)];
  }
  /**
   * The following methods help blacklisting pairs of gene and species nodes.
   * The motivation is that both ML and sampling reconciliation inference
   * algorithm could get stuck in an infinite loop of transfers, because of some
   * approximations in the likelihood computation. The blacklist can be used by
   * these algorithm to detect this infinite loop.
   */
  void initBlackList(unsigned int genesNumber, unsigned int speciesNumber);
  void blackList(unsigned int geneNode, unsigned int speciesNode);
  bool isBlacklisted(unsigned int geneNode, unsigned int speciesNode);
  void resetBlackList();

  std::vector<corax_unode_t *> &getGeneNodeBuffer() { return _geneNodeBuffer; }
  corax_unode_t *generateGeneRoot();
  void generateGeneChildren(corax_unode_t *geneNode,
                            corax_unode_t *&leftGeneNode,
                            corax_unode_t *&rightGeneNode);

  corax_unode_t *getGeneRoot() const { return _geneRoot; }
  unsigned int getVirtualRootIndex() const { return _virtualRootIndex; }
  corax_rtree_t *getSpeciesTree() const { return _speciesTree; }
  const std::vector<std::vector<Event>> &getGeneIdToEvents() const {
    return _geneIdToEvents;
  }
  corax_rnode_t *getOriginationSpecies() const;
  static const unsigned int EVENT_TYPE_NUMBER;

  void setLastEventType(ReconciliationEventType type) { _lastEventType = type; }
  ReconciliationEventType getLastEventType() const { return _lastEventType; }

private:
  static const char *eventNames[];
  ReconciliationEventType _lastEventType;
  std::vector<Event> _events;
  std::vector<unsigned int> _eventsCount;
  std::vector<std::vector<Event>> _geneIdToEvents;
  corax_unode_t *_geneRoot;
  PLLRootedTree *_rootedTree;
  corax_rtree_s *_speciesTree;
  unsigned int _virtualRootIndex;
  std::vector<corax_unode_t *> _geneNodeBuffer;
  typedef std::vector<std::vector<int>> ScenarioBlackList;
  std::unique_ptr<ScenarioBlackList> _blacklist;
  OrthoGroup *getLargestOrthoGroupRec(corax_unode_t *geneNode,
                                      bool isVirtualRoot) const;
  void getAllOrthoGroupRec(corax_unode_t *geneNode, OrthoGroups &orthogroups,
                           OrthoGroupPtr &currentOrthoGroup,
                           bool isVirtualRoot) const;
};
