#pragma once

#include <string>
#include <vector>

#include <util/Scenario.hpp>
#include <util/enums.hpp>

class ParallelOfstream;

class ReconciliationWriter {
public:
  ReconciliationWriter() = delete;

  /**
   *  Write a reconciliation into a stream using the NHX format
   */
  static void
  saveReconciliationNHX(corax_rtree_t *speciesTree, corax_unode_t *geneRoot,
                        unsigned int virtualRootIndex,
                        std::vector<std::vector<Scenario::Event>> &geneToEvent,
                        ParallelOfstream &os);

  /**
   *  Write a reconciliation into a stream using the AleRec format
   */
  static void
  saveReconciliationALE(corax_rtree_t *speciesTree, corax_unode_t *geneRoot,
                        unsigned int virtualRootIndex,
                        std::vector<std::vector<Scenario::Event>> &geneToEvents,
                        ParallelOfstream &os);

  /**
   *  Write a reconciliation into a stream using the RecPhyloXML format
   */
  static void saveReconciliationRecPhyloXML(
      corax_rtree_t *speciesTree, unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event>> &geneToEvent,
      ParallelOfstream &os);

  /**
   *  Write a reconciliation into a stream using the NewickEvents format
   */
  static void saveReconciliationNewickEvents(
      corax_unode_t *geneRoot, unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event>> &geneToEvent,
      ParallelOfstream &os);
};
