#include "DatedSpeciesTreeSearch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "SpeciesSearchCommon.hpp"
#include "SpeciesTransferSearch.hpp"
#include <IO/Logger.hpp>
#include <maths/Random.hpp>
#include <parallelization/ParallelContext.hpp>
#include <trees/SpeciesTree.hpp>

// Search for the speciation order (dating) optimizing the score returned by
// the evaluator. If searchState is provided and the score gets higher than
// searchState.bestLL, save the new best tree and update searchState.bestLL
static double
optimizeDatesLocal(SpeciesTree &speciesTree,
                   SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                   SpeciesSearchState *searchState = nullptr) {
  bool verbose = evaluator.isVerbose();
  auto &datedTree = speciesTree.getDatedTree();
  auto bestLL = evaluator.computeLikelihood();
  if (verbose) {
    Logger::timed << "Starting new naive dating search from ll=" << bestLL
                  << std::endl;
  }
  bool tryAgain = false;
  auto maxRank = datedTree.getRootedTree().getInnerNodeNumber();
  do {
    auto initialItLL = bestLL;
    for (unsigned int rank = 0; rank < maxRank; ++rank) {
      if (!datedTree.moveUp(rank)) { // the node with rank gets rank-1
        continue;
      }
      speciesTree.onSpeciesDatesChange();
      PerFamLL perFamLL;
      auto ll = evaluator.computeLikelihood(&perFamLL);
      if (searchState && ll > searchState->bestLL) {
        // the tree is better than the last saved tree
        // update searchState to save the tree
        searchState->betterTreeCallback(ll, perFamLL);
      }
      if (ll > bestLL) {
        // the best tree over all performed iterations
        bestLL = ll;
        rank -= std::min(2u, rank);
      } else {
        datedTree.moveUp(rank); // reversal: the node with rank-1 gets rank
      }
    }
    // we'll run another iteration only if the improvement is above 1.0
    tryAgain = (bestLL - initialItLL > 1.0);
    if (verbose) {
      Logger::timed << " end of round, ll=" << bestLL << std::endl;
    }
  } while (tryAgain);
  speciesTree.onSpeciesDatesChange();
  if (verbose) {
    Logger::timed << "End of naive dating search, ll=" << bestLL << std::endl;
  }
  return bestLL;
}

// Randomly perturb the order of speciation events in the species tree.
// The number of perturbations is proportional to perturbation, which is
// typically between 0 and 1 (but can be greater)
static void perturbateDates(SpeciesTree &speciesTree,
                            double perturbation = 1.0) {
  assert(perturbation > 0.0);
  auto &datedTree = speciesTree.getDatedTree();
  auto N = datedTree.getRootedTree().getInnerNodeNumber();
  auto perturbations =
      static_cast<unsigned int>(double(N) * 2.0 * perturbation);
  auto maxDisplacement =
      static_cast<unsigned int>(std::sqrt(double(N)) * 2.0 * perturbation);
  maxDisplacement = std::max(2u, maxDisplacement);
  for (unsigned int i = 0; i < perturbations; ++i) {
    bool isUp = Random::getBool();
    unsigned int rank = Random::getInt() % N;
    unsigned int displacement = 1 + (Random::getInt() % maxDisplacement);
    unsigned int nodesToMove = 1 + (Random::getInt() % 10);
    bool ok;
    for (unsigned int k = 0; k < nodesToMove; ++k) {
      for (unsigned int j = 0; j < displacement; ++j) {
        if (isUp) {
          ok = datedTree.moveUp(rank + k - j);
        } else {
          ok = datedTree.moveDown(rank - k + j);
        }
        if (!ok) { // break both cycles
          k = nodesToMove, j = displacement;
        }
      }
    }
  }
  speciesTree.onSpeciesDatesChange();
}

double DatedSpeciesTreeSearch::optimizeDates(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    SpeciesSearchState &searchState, bool thorough) {
  // initial tree LL (it may differ from searchState.bestLL)
  PerFamLL perFamLL;
  auto initialLL = evaluator.computeLikelihood(&perFamLL);
  if (initialLL > searchState.bestLL) {
    searchState.betterTreeCallback(initialLL, perFamLL);
  }
  if (!evaluator.isDated()) {
    return initialLL;
  }
  Logger::timed << "[Species search] Optimizing dates, ll=" << initialLL
                << std::endl;
  // initial optimization
  auto bestLL = optimizeDatesLocal(speciesTree, evaluator, &searchState);
  // perturbation-optimization cycles
  const double perturbation = 0.1;
  const unsigned int maxTrials = 2;
  unsigned int unsuccessfulTrials = 0;
  while (thorough && unsuccessfulTrials < maxTrials) {
    auto backup = speciesTree.getDatedTree().getBackup();
    perturbateDates(speciesTree, perturbation);
    auto ll = optimizeDatesLocal(speciesTree, evaluator, &searchState);
    if (ll > bestLL) {
      bestLL = ll;
      unsuccessfulTrials = 0;
      Logger::timed << "[Species search]   better ll=" << bestLL << std::endl;
    } else {
      SpeciesTreeOperator::restoreDates(speciesTree, backup);
      unsuccessfulTrials++;
    }
  }
  Logger::timed << "[Species search]   After date opt, ll=" << bestLL
                << std::endl;
  return bestLL;
}

// Evaluate the current tree dating based on the share of precomputed
// undated transfer events that are supported by the dating. Better
// datings permit more precomputed transfers and get higher scores
static unsigned int getTransferScore(SpeciesTree &speciesTree,
                                     const TransferFrequencies &frequencies) {
  unsigned int score = 0;
  StringToUint labelToId;
  speciesTree.getLabelToId(labelToId);
  // parallelize across species for less computational redundancy
  auto N = frequencies.count.size();
  auto begin = ParallelContext::getBegin(N);
  auto end = ParallelContext::getEnd(N);
  for (unsigned int from = begin; from < end; ++from) {
    for (unsigned int to = 0; to < N; ++to) {
      auto count = frequencies.count[from][to];
      if (count) {
        // check if the current dating permits the precomputed transfer
        auto src = labelToId[frequencies.idToLabel[from]];
        auto dest = labelToId[frequencies.idToLabel[to]];
        if (speciesTree.getDatedTree().canTransferUnderRelDated(src, dest)) {
          score += count;
        }
      }
    }
  }
  ParallelContext::sumUInt(score);
  return score;
}

class TransferScoreEvaluator : public SpeciesTreeLikelihoodEvaluatorInterface {
private:
  SpeciesTree &_speciesTree;
  TransferFrequencies &_frequencies;

public:
  TransferScoreEvaluator(SpeciesTree &speciesTree,
                         TransferFrequencies &frequencies)
      : _speciesTree(speciesTree), _frequencies(frequencies) {}
  virtual double computeLikelihood(PerFamLL *perFamLL = nullptr) {
    (void)(perFamLL);
    return computeLikelihoodFast();
  }
  virtual double computeLikelihoodFast() {
    return getTransferScore(_speciesTree, _frequencies);
  }
  virtual bool providesFastLikelihoodImpl() const { assert(false); }
  virtual bool isDated() const { assert(false); }
  virtual double optimizeModelRates(bool) { assert(false); }
  virtual void pushRollback() { assert(false); }
  virtual void popAndApplyRollback() { assert(false); }
  virtual void getTransferInformation(SpeciesTree &, TransferFrequencies &,
                                      PerSpeciesEvents &,
                                      PerCorePotentialTransfers &) {
    assert(false);
  }
  virtual bool pruneSpeciesTree() const { assert(false); }
  virtual void onSpeciesDatesChange() { assert(false); }
  virtual void
  onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *) {
    assert(false);
  }
  virtual bool isVerbose() const { return false; }
};

ScoredBackups DatedSpeciesTreeSearch::getBestDatingsFromReconciliation(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator, unsigned int toTest,
    unsigned int toTake) {
  assert(toTake <= toTest);
  bool verbose = evaluator.isVerbose();
  auto &datedTree = speciesTree.getDatedTree();
  auto reconciliationDatingBackup = datedTree.getBackup();
  ScoredBackups scoredBackups;
  // get the transfers from reconciliations
  TransferFrequencies frequencies;
  PerSpeciesEvents perSpeciesEvents;
  PerCorePotentialTransfers potentialTransfers;
  evaluator.getTransferInformation(speciesTree, frequencies, perSpeciesEvents,
                                   potentialTransfers);
  TransferScoreEvaluator fakeEvaluator(speciesTree, frequencies);
  // start multiple searches from random datings
  for (unsigned int i = 0; i < toTest; ++i) {
    // we should replace this with anything that would produce
    // a random dating more efficiently
    datedTree.randomize();
    // first local search to get to a good starting tree
    auto bestScore = optimizeDatesLocal(speciesTree, fakeEvaluator);
    // Thorough round: at each step, randomly perturb the tree and
    // perform a local search. If no better tree is found, start
    // again with a greater perturbation, until maxTrials trials
    // without improvement. If there is an improvement, restart
    // the algorithm from the new best tree
    const unsigned int maxTrials = 20;
    unsigned int unsuccessfulTrials = 0;
    while (unsuccessfulTrials < maxTrials) {
      auto backup = datedTree.getBackup();
      // the perturbation parameter increases with the number of failures
      double perturbation = double(unsuccessfulTrials + 1) / double(maxTrials);
      perturbateDates(speciesTree, perturbation);
      auto score = optimizeDatesLocal(speciesTree, fakeEvaluator);
      if (score > bestScore) {
        // better tree found, reset the algorithm
        bestScore = score;
        unsuccessfulTrials = 0;
        // Logger::timed << " better score=" << bestScore << ", ll="
        //               << evaluator.computeLikelihood() << std::endl;
      } else {
        // this tree is worse than the best one, we rollback
        SpeciesTreeOperator::restoreDates(speciesTree, backup);
        unsuccessfulTrials++;
      }
    }
    scoredBackups.push_back(ScoredBackup(datedTree.getBackup(), bestScore));
    if (verbose) {
      Logger::timed << "End of iteration " << i << ", score=" << bestScore
                    << std::endl;
    }
  }
  // sort the datings by transfer score and take the best ones
  std::sort(scoredBackups.rbegin(), scoredBackups.rend());
  scoredBackups.resize(toTake);
  // for each dating compute the real likelihood (not the transfer score)
  // and set it as the dating score
  for (auto &sb : scoredBackups) {
    SpeciesTreeOperator::restoreDates(speciesTree, sb.backup);
    auto ll = evaluator.computeLikelihood();
    if (verbose) {
      Logger::info << "score=" << sb.score << ", ll=" << ll << std::endl;
    }
    sb.score = ll;
  }
  std::sort(scoredBackups.rbegin(), scoredBackups.rend());
  // reset the tree to its initial dating
  SpeciesTreeOperator::restoreDates(speciesTree, reconciliationDatingBackup);
  return scoredBackups;
}
