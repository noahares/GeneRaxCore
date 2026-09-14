#include "SpeciesTransferSearch.hpp"

#include <algorithm>
#include <cassert>

#include "SpeciesSearchCommon.hpp"
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <trees/SpeciesTree.hpp>
#include <util/Scenario.hpp>

void PerCorePotentialTransfers::addScenario(const Scenario &scenario) {
  auto famCopies = scenario.getPerSpeciesCopies();
  if (!copies.size()) {
    copies.resize(famCopies.size());
  }
  assert(copies.size() == famCopies.size());
  for (unsigned int species = 0; species < famCopies.size(); ++species) {
    copies[species].push_back(famCopies[species]);
  }
}

unsigned int
PerCorePotentialTransfers::countPotentialTransfers(unsigned int src,
                                                   unsigned int dest) const {
  unsigned int res = 0;
  if (copies.size()) {
    assert(copies[src].size() == copies[dest].size());
    for (unsigned int fam = 0; fam < copies[src].size(); ++fam) {
      if (copies[dest][fam]) {
        res += copies[src][fam];
      }
    }
  }
  ParallelContext::sumUInt(res);
  return res;
}

void SpeciesTransferSearch::getSortedTransferList(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    unsigned int minTransfers, const MovesBlackList &blacklist,
    std::vector<TransferMove> &transferMoves) {
  TransferFrequencies frequencies;
  PerSpeciesEvents perSpeciesEvents;
  PerCorePotentialTransfers potentialTransfers;
  evaluator.getTransferInformation(speciesTree, frequencies, perSpeciesEvents,
                                   potentialTransfers);
  unsigned int transfers = 0;
  StringToUint labelToId;
  speciesTree.getLabelToId(labelToId);
  auto N = frequencies.count.size();
  for (unsigned int from = 0; from < N; ++from) {
    for (unsigned int to = 0; to < N; ++to) {
      auto regraft = labelToId[frequencies.idToLabel[from]];
      auto prune = labelToId[frequencies.idToLabel[to]];
      auto count = frequencies.count[from][to];
      transfers += count;
      if (count < minTransfers) {
        continue;
      }
      TransferMove move(prune, regraft, count);
      if (true || !blacklist.isBlackListed(move)) { // disable the check for now
        // support: the proportion of regraft copies surviving in prune
        move.support = static_cast<double>(count);
        move.support /= static_cast<double>(
            potentialTransfers.countPotentialTransfers(regraft, prune));
        transferMoves.push_back(move);
      }
    }
  }
  std::sort(transferMoves.rbegin(), transferMoves.rend());
  Logger::timed << "Inferred transfers: " << transfers << std::endl;
  Logger::timed << "Selected transfer directions: " << transferMoves.size()
                << std::endl;
}

static bool transferRound(SpeciesTree &speciesTree,
                          SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                          SpeciesSearchState &searchState,
                          MovesBlackList &blacklist) {
  Logger::timed << "[Species search] Start new transfer-guided round"
                << std::endl;
  // update sprBoots on the current topology
  PerFamLL perFamLL;
  evaluator.computeLikelihood(&perFamLL);
  std::vector<unsigned int> affectedBranches;
  for (unsigned int i = 0; i < speciesTree.getTree().getNodeNumber(); ++i) {
    affectedBranches.push_back(i);
  }
  for (auto &bs : searchState.sprBoots) {
    bs.test(perFamLL, affectedBranches, true);
  }
  // do transfer-guided SPR search
  auto hash1 = speciesTree.getNodeIndexHash();
  const unsigned int minTransfers = 1;
  std::vector<TransferMove> transferMoves;
  SpeciesTransferSearch::getSortedTransferList(
      speciesTree, evaluator, minTransfers, blacklist, transferMoves);
  const unsigned int speciesNumber = speciesTree.getTree().getNodeNumber();
  const unsigned int minTrials = std::max(50u, speciesNumber / 2);
  const unsigned int stopAfterFailures = 50u;
  const unsigned int stopAfterImprovements = std::max(15u, speciesNumber / 4);
  unsigned int trials = 0;
  unsigned int failures = 0;
  unsigned int improvements = 0;
  std::unordered_set<unsigned int> alreadyPruned;
  for (const auto &move : transferMoves) {
    // check if this prune has already been regrafted with
    // a move that has a higher transfer support
    if (alreadyPruned.find(move.prune) != alreadyPruned.end()) {
      continue;
    }
    // as the tree topology can be changed during the round,
    // the precomputed transfer moves should be checked for validity
    if (!SpeciesTreeOperator::canApplySPRMove(speciesTree, move.prune,
                                              move.regraft)) {
      continue;
    }
    blacklist.blackList(move);
    trials++;
    /*
    Logger::info << "Test "
                 << speciesTree.getNode(move.prune)->label << " -> "
                 << speciesTree.getNode(move.regraft)->label << std::endl;
    */
    if (SpeciesSearchCommon::testSPR(speciesTree, evaluator, searchState,
                                     move.prune, move.regraft)) {
      failures = 0;
      improvements++;
      alreadyPruned.insert(move.prune);
      Logger::timed << "\tbetter tree "
                    << "(support: " << move.support << ", trial: " << trials
                    << ", LL=" << searchState.bestLL
                    << ", hash=" << speciesTree.getHash() << ") "
                    << speciesTree.getNode(move.prune)->label << " -> "
                    << speciesTree.getNode(move.regraft)->label << std::endl;
      hash1 = speciesTree.getNodeIndexHash();
      assert(ParallelContext::isIntEqual(hash1));
      if (!searchState.farFromPlausible) {
        SpeciesSearchCommon::veryLocalSearch(speciesTree, evaluator,
                                             searchState, move.prune);
      }
    } else {
      failures++;
    }
    // check if we've had enough failures to stop or enough improvements
    // to recompute the new transfers
    bool stop = (trials > minTrials) && (failures > stopAfterFailures);
    bool maxImprovementsReached = (improvements > stopAfterImprovements);
    stop |= maxImprovementsReached;
    if (stop) { // end this round
      if (searchState.farFromPlausible && !maxImprovementsReached) {
        // the inspected transfers do not point to a better tree, so the tree
        // is already near-optimal, we switch to a more thorough search mode,
        // this will affect the ENTIRE downstream analysis
        Logger::timed << "[Species search] Switch to hardToFindBetter mode"
                      << std::endl;
        searchState.farFromPlausible = false;
      }
      return improvements > 0;
    }
  }
  return improvements > 0;
}

bool SpeciesTransferSearch::transferSearch(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    SpeciesSearchState &searchState) {
  Logger::timed
      << "[Species search] Starting species tree transfer-guided search "
      << "(bestLL=" << searchState.bestLL << ", hash=" << speciesTree.getHash()
      << ")" << std::endl;
  assert(evaluator.computeLikelihood() == searchState.bestLL);
  bool better = false;
  bool tryAgain = false;
  MovesBlackList blacklist;
  do {
    // optimize model parameters
    searchState.bestLL = evaluator.optimizeModelRates();
    PerFamLL perFamLL;
    evaluator.computeLikelihood(&perFamLL);
    searchState.betterLikelihoodCallback(searchState.bestLL, perFamLL);
    // run a transfer-guided search round
    tryAgain = transferRound(speciesTree, evaluator, searchState, blacklist);
    better |= tryAgain;
  } while (tryAgain);
  // update bestLL due to possible precision changes
  PerFamLL perFamLL;
  double ll = evaluator.computeLikelihood(&perFamLL);
  searchState.betterLikelihoodCallback(ll, perFamLL);
  Logger::timed << "[Species search] After transfer-guided search: LL="
                << searchState.bestLL << std::endl;
  return better;
}
