#include "DatedSpeciesTreeSearch.hpp"

#include "SpeciesSearchCommon.hpp"
#include "SpeciesTransferSearch.hpp"
#include <IO/Logger.hpp>
#include <trees/SpeciesTree.hpp>

static void saveDatedTree(SpeciesTree &speciesTree,
                          const SpeciesSearchState &state) {
  if (state.pathToBestSpeciesTree.size() == 0) {
    return;
  }
  speciesTree.getDatedTree().rescaleBranchLengths();
  speciesTree.saveToFile(state.pathToBestSpeciesTree, true);
  ParallelContext::barrier();
}

static double
optimizeDatesLocal(SpeciesTree &speciesTree, SpeciesSearchState &state,
                   SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                   bool verbose) {
  auto &tree = speciesTree.getDatedTree();
  unsigned int max = tree.getOrderedSpeciations().size();
  auto bestLL = evaluator.computeLikelihood();
  bool tryAgain = true;
  if (verbose) {
    Logger::timed << "Starting new naive dating search from ll=" << bestLL
                  << std::endl;
  }
  while (tryAgain) {
    tryAgain = false;
    const double initialBestLL = bestLL;
    for (unsigned int rank = 0; rank < max; ++rank) {
      if (!tree.moveUp(rank)) {
        continue;
      }
      evaluator.onSpeciesDatesChange();
      auto ll = evaluator.computeLikelihood();
      if (ll > state.bestLL) {
        // best tree over all search iterations, we save it!
        state.bestLL = ll;
        saveDatedTree(speciesTree, state);
      }
      if (ll > bestLL) {
        // best tree for this search iteration
        // (we might have found a better likelihood
        // in another iteration)
        // we only restart if the improvement is above 1.0
        if (ll - initialBestLL > 1.0) {
          tryAgain = true;
        }
        bestLL = ll;
        rank -= std::min((unsigned int)2, rank);
      } else {
        tree.moveUp(rank);
      }
      assert(tree.isConsistent());
    }
    if (verbose) {
      Logger::timed << " end of round, ll=" << bestLL << std::endl;
    }
  }
  evaluator.onSpeciesDatesChange();
  if (verbose) {
    Logger::timed << "End of naive dating search, ll=" << bestLL << std::endl;
  }
  return bestLL;
}

// Randomly perturbate the order of speciation events in the dated
// species tree. The level of perturbation is proportional with
// perturbation, which is typicall between 0 and 1 (but can be greater)
static void perturbateDates(SpeciesTree &speciesTree,
                            double perturbation = 1.0) {
  auto &tree = speciesTree.getDatedTree();
  size_t N = tree.getOrderedSpeciations().size();
  auto perturbations = size_t(double((2 * N) * perturbation));
  auto maxDisplacement = size_t(sqrt(double(N)) * 2.0 * perturbation);
  if (maxDisplacement < 2) {
    maxDisplacement = 2;
  }
  for (size_t i = 0; i < perturbations; ++i) {
    auto rank = Random::getInt() % N;
    bool upOrDown = Random::getBool();
    auto displacement = 1 + (Random::getInt() % maxDisplacement);
    auto nodesToMove = static_cast<size_t>((Random::getInt() % 10) + 1);
    bool ok = true;
    for (size_t k = 0; k < nodesToMove; ++k) {
      for (size_t j = 0; j < displacement; ++j) {
        if (upOrDown) {
          ok &= tree.moveUp(rank + k - j);
        } else {
          ok &= tree.moveDown(rank - k + j);
        }
        if (!ok) {
          break;
        }
      }
      if (!ok) {
        break;
      }
    }
  }
}

double DatedSpeciesTreeSearch::optimizeDates(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    SpeciesSearchState &state, bool thorough) {
  bool verbose = evaluator.isVerbose();
  auto bestLL = evaluator.computeLikelihood();
  if (!evaluator.isDated()) {
    return bestLL;
  }
  Logger::timed << "[Species search] Optimizing dates, ll=" << bestLL
                << std::endl;
  optimizeDatesLocal(speciesTree, state, evaluator, verbose);
  bestLL = evaluator.computeLikelihood();
  unsigned int unsuccessfulTrials = 0;
  const unsigned int maxTrials = 2;
  auto &tree = speciesTree.getDatedTree();
  while (thorough && unsuccessfulTrials < maxTrials) {
    auto backup = tree.getBackup();
    double perturbation =
        0.1; // double(unsuccessfulTrials + 1) / double(maxTrials);
    perturbateDates(speciesTree, perturbation);
    evaluator.onSpeciesDatesChange();
    auto ll = evaluator.computeLikelihood();
    optimizeDatesLocal(speciesTree, state, evaluator, verbose);
    ll = evaluator.computeLikelihood();
    if (ll <= bestLL) {
      tree.restore(backup);
      unsuccessfulTrials++;
    } else {
      unsuccessfulTrials = 0;
      bestLL = ll;
      Logger::timed << "[Species search]   better ll=" << bestLL << std::endl;
    }
  }
  Logger::timed << "[Species search]   After date opt, ll=" << bestLL
                << std::endl;
  return bestLL;
}

unsigned int getTransferScore(SpeciesTree &speciesTree,
                              const TransferFrequencies &frequencies) {
  std::unordered_map<std::string, unsigned int> labelsToIds;
  speciesTree.getLabelsToId(labelsToIds);
  unsigned int score = 0;
  auto N = frequencies.count.size();
  auto begin = ParallelContext::getBegin(N);
  auto end = ParallelContext::getEnd(N);
  for (unsigned int from = begin; from < end; ++from) {
    for (unsigned int to = 0; to < frequencies.count.size(); ++to) {
      auto count = frequencies.count[from][to];
      if (!count) {
        continue;
      }
      auto src = labelsToIds[frequencies.idToLabel[from]];
      auto dest = labelsToIds[frequencies.idToLabel[to]];
      if (speciesTree.getDatedTree().canTransferUnderRelDated(src, dest)) {
        score += count;
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
    assert(!perFamLL);
    return computeLikelihoodFast();
  }
  virtual double computeLikelihoodFast() {
    return getTransferScore(_speciesTree, _frequencies);
  }
  virtual bool providesFastLikelihoodImpl() const { return false; }
  virtual bool isDated() const { return true; }
  virtual double optimizeModelRates(bool) { return computeLikelihood(); }
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void getTransferInformation(SpeciesTree &, TransferFrequencies &,
                                      PerSpeciesEvents &,
                                      PerCorePotentialTransfers &) {
    assert(false);
  }
  virtual bool pruneSpeciesTree() const { return false; }
};

ScoredBackups DatedSpeciesTreeSearch::optimizeDatesFromReconciliation(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator, unsigned int searches,
    unsigned int toEvaluate) {
  bool verbose = evaluator.isVerbose();
  auto &tree = speciesTree.getDatedTree();
  auto reconciliationDatingBackup = tree.getBackup();
  ScoredBackups scoredBackups;
  // get the transfers from reconciliations
  TransferFrequencies frequencies;
  PerSpeciesEvents perSpeciesEvents;
  PerCorePotentialTransfers potentialTransfers;
  evaluator.getTransferInformation(speciesTree, frequencies, perSpeciesEvents,
                                   potentialTransfers);
  TransferScoreEvaluator fakeEvaluator(speciesTree, frequencies);
  // start multiple searches from random datings
  for (unsigned int i = 0; i < searches; ++i) {
    // we should replace this with anything that would produce
    // a random dating more efficiently
    tree.randomize();
    SpeciesSearchState fakeState(speciesTree, "", 0);
    fakeState.bestLL = fakeEvaluator.computeLikelihood();
    // first local search to get to a good starting tree
    auto bestScore =
        optimizeDatesLocal(speciesTree, fakeState, fakeEvaluator, false);
    unsigned int unsuccessfulTrials = 0;
    const unsigned int maxTrials = 20;
    // Thorough round: at each step, randomly perturbate the tree
    // and perform a local search. If no better tree is found, start
    // again with a greater perturbation, until maxTrials trials
    // without improvement. If there is an improvement, restart
    // the algorithm from the new best tree.
    while (unsuccessfulTrials < maxTrials) {
      auto backup = tree.getBackup();
      // the random perturbation increases with the number of failures
      double perturbation = double(unsuccessfulTrials + 1) / double(maxTrials);
      perturbateDates(speciesTree, perturbation);
      optimizeDatesLocal(speciesTree, fakeState, fakeEvaluator, false);
      auto score = fakeEvaluator.computeLikelihood();
      if (score <= bestScore) {
        // this tree is worse than the best one, we rollback
        tree.restore(backup);
        unsuccessfulTrials++;
      } else {
        // better tree found, reset the algorithm
        unsuccessfulTrials = 0;
        bestScore = score;
        // Logger::timed << " better score=" << bestScore  << std::endl;//" ll="
        // << evaluator.computeLikelihood()<< std::endl;
      }
    }
    scoredBackups.push_back(ScoredBackup(tree, bestScore));
    if (verbose) {
      Logger::timed << "End of iteration " << i << ", score=" << bestScore
                    << std::endl;
      //            << ", ll=" << ll << std::endl;
    }
  }
  std::sort(scoredBackups.rbegin(), scoredBackups.rend());
  scoredBackups.resize(std::min<size_t>(scoredBackups.size(), toEvaluate));
  for (auto &sb : scoredBackups) {
    // now we compute the "real likelihood" (not the transfer score)
    // and map it to this dating
    speciesTree.getDatedTree().restore(sb.backup);
    evaluator.onSpeciesDatesChange();
    auto ll = evaluator.computeLikelihood();
    if (verbose) {
      Logger::info << "score=" << sb.score << ", ll=" << ll << std::endl;
    }
    sb.score = ll;
  }
  std::sort(scoredBackups.rbegin(), scoredBackups.rend());
  // reset the tree to its initial dating
  tree.restore(reconciliationDatingBackup);
  return scoredBackups;
}
