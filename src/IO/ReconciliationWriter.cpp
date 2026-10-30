#include "ReconciliationWriter.hpp"

#include <IO/ParallelOfstream.hpp>

/**
 *
 *  NHX format
 *
 */
static void printEventNHX(const Scenario::Event &event,
                          corax_rtree_t *speciesTree, double nodeBL,
                          ParallelOfstream &os) {
  assert(event.isValid());
  corax_rnode_t *species = speciesTree->nodes[event.speciesNode];
  corax_rnode_t *speciesDest = nullptr;
  assert(species->label);
  os << "[&&NHX";
  os << ":S=" << species->label;
  os << ":D=" << ((event.type == ReconciliationEventType::EVENT_D) ? "Y" : "N");
  os << ":H="
     << ((event.type == ReconciliationEventType::EVENT_T ||
          event.type == ReconciliationEventType::EVENT_TL)
             ? "Y"
             : "N");
  if (event.type == ReconciliationEventType::EVENT_T ||
      event.type == ReconciliationEventType::EVENT_TL) {
    speciesDest = speciesTree->nodes[event.destSpeciesNode];
    assert(speciesDest->label);
    os << "@" << species->label;
    os << "@" << speciesDest->label;
  }
  os << ":B=" << nodeBL;
  os << "]";
}

static void recursivelySaveReconciliationsNHX(
    corax_rtree_t *speciesTree, corax_unode_t *node, unsigned int depth,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  if (node->next) {
    auto left = (depth == 0) ? node->next : node->next->back;
    auto right = (depth == 0) ? node->next->back : node->next->next->back;
    os << "(";
    recursivelySaveReconciliationsNHX(speciesTree, left, depth + 1,
                                      geneToEvents, os);
    os << ",";
    recursivelySaveReconciliationsNHX(speciesTree, right, depth + 1,
                                      geneToEvents, os);
    os << ")";
  }
  if (node->label) {
    os << node->label;
  } else {
    os << "n" << node->node_index;
  }
  // divide the root BL by two to place the root
  // at the middle of this branch
  auto nodeBL = (depth == 1) ? node->length / 2.0 : node->length;
  if (depth > 0) {
    os << ":" << nodeBL;
  }
  printEventNHX(geneToEvents[node->node_index].back(), speciesTree, nodeBL, os);
}

void ReconciliationWriter::saveReconciliationNHX(
    corax_rtree_t *speciesTree, corax_unode_t *geneRoot,
    unsigned int virtualRootIndex,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  corax_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  recursivelySaveReconciliationsNHX(speciesTree, &virtualRoot, 0, geneToEvents,
                                    os);
  os << ";";
  os << std::endl;
}

/**
 *
 *  AleRec format
 *
 */
static void printEventALE(const Scenario::Event &event,
                          corax_rtree_t *speciesTree, ParallelOfstream &os) {
  assert(event.isValid());
  corax_rnode_t *species = speciesTree->nodes[event.speciesNode];
  corax_rnode_t *speciesDest = nullptr;
  assert(species->label);
  os << ".";
  switch (event.type) {
  case ReconciliationEventType::EVENT_S:
  case ReconciliationEventType::EVENT_SL:
    os << "S@" << species->label;
    break;
  case ReconciliationEventType::EVENT_D:
    os << "D@" << species->label;
    break;
  case ReconciliationEventType::EVENT_T:
  case ReconciliationEventType::EVENT_TL:
    speciesDest = speciesTree->nodes[event.destSpeciesNode];
    assert(speciesDest->label);
    os << "T@" << species->label << "->" << speciesDest->label;
    break;
  case ReconciliationEventType::EVENT_None:
    os << "Leaf@" << species->label << "...";
    break;
  default:
    break;
  }
}

static void recursivelySaveReconciliationsALE(
    corax_rtree_t *speciesTree, corax_unode_t *node, unsigned int depth,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  if (node->next) {
    auto left = (depth == 0) ? node->next : node->next->back;
    auto right = (depth == 0) ? node->next->back : node->next->next->back;
    os << "(";
    recursivelySaveReconciliationsALE(speciesTree, left, depth + 1,
                                      geneToEvents, os);
    os << ",";
    recursivelySaveReconciliationsALE(speciesTree, right, depth + 1,
                                      geneToEvents, os);
    os << ")";
  }
  for (const auto &event : geneToEvents[node->node_index]) {
    printEventALE(event, speciesTree, os);
  }
  if (!node->next) {
    os << (node->label ? node->label : "null");
  }
  // divide the root BL by two to place the root
  // at the middle of this branch
  auto nodeBL = (depth == 1) ? node->length / 2.0 : node->length;
  if (depth > 0) {
    os << ":" << nodeBL;
  }
}

void ReconciliationWriter::saveReconciliationALE(
    corax_rtree_t *speciesTree, corax_unode_t *geneRoot,
    unsigned int virtualRootIndex,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  corax_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  recursivelySaveReconciliationsALE(speciesTree, &virtualRoot, 0, geneToEvents,
                                    os);
  os << ";";
  os << std::endl;
}

/**
 *
 *  RecPhyloXML format
 *
 */
static void recursivelySaveSpeciesTreeRecPhyloXML(corax_rnode_t *species,
                                                  std::string &indent,
                                                  ParallelOfstream &os) {
  if (!species) {
    return;
  }
  os << indent << "<clade>" << std::endl;
  indent += "\t";
  assert(species->label);
  os << indent << "<name>" << species->label << "</name>" << std::endl;
  recursivelySaveSpeciesTreeRecPhyloXML(species->left, indent, os);
  recursivelySaveSpeciesTreeRecPhyloXML(species->right, indent, os);
  indent.pop_back();
  os << indent << "</clade>" << std::endl;
}

static void saveSpeciesTreeRecPhyloXML(corax_rtree_t *speciesTree,
                                       ParallelOfstream &os) {
  os << "<spTree>" << std::endl;
  os << "<phylogeny>" << std::endl;
  std::string indent = "";
  recursivelySaveSpeciesTreeRecPhyloXML(speciesTree->root, indent, os);
  os << "</phylogeny>" << std::endl;
  os << "</spTree>" << std::endl;
}

static void writeEventRecPhyloXML(corax_rtree_t *speciesTree,
                                  unsigned int geneIndex,
                                  const Scenario::Event &event,
                                  const Scenario::Event *previousEvent,
                                  std::string &indent, ParallelOfstream &os) {
  corax_rnode_t *species = speciesTree->nodes[event.speciesNode];
  assert(species->label);
  os << indent << "<eventsRec>" << std::endl;
  if (event.type != ReconciliationEventType::EVENT_L) {
    if ((previousEvent->type == ReconciliationEventType::EVENT_T &&
         geneIndex == previousEvent->rightGeneIndex) ||
        previousEvent->type == ReconciliationEventType::EVENT_TL) {
      corax_rnode_t *previousEventSpeciesDest =
          speciesTree->nodes[previousEvent->destSpeciesNode];
      assert(previousEventSpeciesDest->label);
      os << indent << "\t<transferBack destinationSpecies=\""
         << previousEventSpeciesDest->label << "\"/>" << std::endl;
    }
  }
  switch (event.type) {
  case ReconciliationEventType::EVENT_None:
    os << indent << "\t<leaf speciesLocation=\"" << species->label << "\"/>"
       << std::endl;
    break;
  case ReconciliationEventType::EVENT_S:
  case ReconciliationEventType::EVENT_SL:
    os << indent << "\t<speciation speciesLocation=\"" << species->label
       << "\"/>" << std::endl;
    break;
  case ReconciliationEventType::EVENT_D:
    os << indent << "\t<duplication speciesLocation=\"" << species->label
       << "\"/>" << std::endl;
    break;
  case ReconciliationEventType::EVENT_T:
  case ReconciliationEventType::EVENT_TL:
    os << indent << "\t<branchingOut speciesLocation=\"" << species->label
       << "\"/>" << std::endl;
    break;
  case ReconciliationEventType::EVENT_L:
    os << indent << "\t<loss speciesLocation=\"" << species->label << "\"/>"
       << std::endl;
    break;
  default:
    break;
  }
  os << indent << "</eventsRec>" << std::endl;
}

static void recursivelySaveGeneTreeRecPhyloXML(
    corax_rtree_t *speciesTree, unsigned int geneIndex,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    const Scenario::Event *previousEvent, std::string &indent,
    ParallelOfstream &os) {
  const auto &events = geneToEvents[geneIndex];
  // open new clades for loss events of the given geneIndex
  for (unsigned int i = 0; i < events.size() - 1; ++i) {
    os << indent << "<clade>" << std::endl;
    indent += "\t";
    const auto &event = events[i];
    assert(event.type == ReconciliationEventType::EVENT_SL ||
           event.type == ReconciliationEventType::EVENT_TL);
    os << indent << "<name>" << "NULL" << "</name>" << std::endl;
    writeEventRecPhyloXML(speciesTree, geneIndex, event, previousEvent, indent,
                          os);
    previousEvent = &event;
    // now we are one level further from the root:
    // the two child clades are lossEvent and the next event from events
    os << indent << "<clade>" << std::endl;
    indent += "\t";
    Scenario::Event lossEvent;
    lossEvent.type = ReconciliationEventType::EVENT_L;
    lossEvent.speciesNode = (event.type == ReconciliationEventType::EVENT_SL)
                                ? event.lostSpeciesNode
                                : event.speciesNode;
    os << indent << "<name>loss</name>" << std::endl;
    writeEventRecPhyloXML(speciesTree, geneIndex, lossEvent, previousEvent,
                          indent, os);
    indent.pop_back();
    os << indent << "</clade>" << std::endl;
  }
  // handle the last event of the given geneIndex
  os << indent << "<clade>" << std::endl;
  indent += "\t";
  const auto &event = events.back();
  auto label = event.label.size() ? event.label : "NULL";
  os << indent << "<name>" << label << "</name>" << std::endl;
  writeEventRecPhyloXML(speciesTree, geneIndex, event, previousEvent, indent,
                        os);
  if (!event.isLeaf()) {
    recursivelySaveGeneTreeRecPhyloXML(speciesTree, event.leftGeneIndex,
                                       geneToEvents, &event, indent, os);
    recursivelySaveGeneTreeRecPhyloXML(speciesTree, event.rightGeneIndex,
                                       geneToEvents, &event, indent, os);
  }
  indent.pop_back();
  os << indent << "</clade>" << std::endl;
  // close the clades for loss events
  for (unsigned int i = 0; i < events.size() - 1; ++i) {
    indent.pop_back();
    os << indent << "</clade>" << std::endl;
  }
}

static void saveGeneTreeRecPhyloXML(
    corax_rtree_t *speciesTree, unsigned int virtualRootIndex,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  os << "<recGeneTree>" << std::endl;
  os << "<phylogeny rooted=\"true\">" << std::endl;
  std::string indent;
  Scenario::Event previousEvent;
  previousEvent.type = ReconciliationEventType::EVENT_Invalid;
  recursivelySaveGeneTreeRecPhyloXML(speciesTree, virtualRootIndex,
                                     geneToEvents, &previousEvent, indent, os);
  os << "</phylogeny>" << std::endl;
  os << "</recGeneTree>" << std::endl;
}

void ReconciliationWriter::saveReconciliationRecPhyloXML(
    corax_rtree_t *speciesTree, unsigned int virtualRootIndex,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  os << "<recPhylo " << std::endl;
  os << "\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
     << std::endl;
  os << "\txsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\""
     << std::endl;
  os << "\txmlns=\"http://www.recg.org\">" << std::endl;
  saveSpeciesTreeRecPhyloXML(speciesTree, os);
  saveGeneTreeRecPhyloXML(speciesTree, virtualRootIndex, geneToEvents, os);
  os << "</recPhylo>" << std::endl;
}

/**
 *
 *  NewickEvents format
 *
 */
static void recursivelySaveReconciliationsNewickEvents(
    corax_unode_t *node, unsigned int depth,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  if (node->next) {
    auto left = (depth == 0) ? node->next : node->next->back;
    auto right = (depth == 0) ? node->next->back : node->next->next->back;
    os << "(";
    recursivelySaveReconciliationsNewickEvents(left, depth + 1, geneToEvents,
                                               os);
    os << ",";
    recursivelySaveReconciliationsNewickEvents(right, depth + 1, geneToEvents,
                                               os);
    os << ")";
  }
  if (!node->next) {
    os << (node->label ? node->label : "null");
  } else {
    os << Enums::getEventName(geneToEvents[node->node_index].back().type);
  }
  // divide the root BL by two to place the root
  // at the middle of this branch
  auto nodeBL = (depth == 1) ? node->length / 2.0 : node->length;
  if (depth > 0) {
    os << ":" << nodeBL;
  }
}

void ReconciliationWriter::saveReconciliationNewickEvents(
    corax_unode_t *geneRoot, unsigned int virtualRootIndex,
    const std::vector<std::vector<Scenario::Event>> &geneToEvents,
    ParallelOfstream &os) {
  corax_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  recursivelySaveReconciliationsNewickEvents(&virtualRoot, 0, geneToEvents, os);
  os << ";";
  os << std::endl;
}
