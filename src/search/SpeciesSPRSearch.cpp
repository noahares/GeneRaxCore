#include "SpeciesSPRSearch.hpp"

#include <cassert>

#include "SpeciesSearchCommon.hpp"
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <trees/SpeciesTree.hpp>

static bool SPRRound(SpeciesTree &speciesTree,
                     SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                     SpeciesSearchState &searchState, unsigned int radius) {
  Logger::timed << "[Species search] Start new SPR round, radius=" << radius
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
  // do local SPR search
  bool better = false;
  auto hash1 = speciesTree.getNodeIndexHash();
  double maxSupport = 0.2;                    // ignored for now
  auto supportValues = std::vector<double>(); //_getSupport();
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(speciesTree, supportValues, maxSupport,
                                         prunes);
  for (auto prune : prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(speciesTree, prune, radius,
                                             regrafts);
    for (auto regraft : regrafts) {
      // as the tree topology can be changed while working with a prune,
      // the precomputed regrafts should be checked for validity
      if (!SpeciesTreeOperator::canApplySPRMove(speciesTree, prune, regraft)) {
        continue;
      }
      if (SpeciesSearchCommon::testSPR(speciesTree, evaluator, searchState,
                                       prune, regraft)) {
        better = true;
        Logger::timed << "\tbetter tree "
                      << "(LL=" << searchState.bestLL
                      << ", hash=" << speciesTree.getHash() << ") "
                      << speciesTree.getNode(prune)->label << " -> "
                      << speciesTree.getNode(regraft)->label << std::endl;
        hash1 = speciesTree.getNodeIndexHash();
        assert(ParallelContext::isIntEqual(hash1));
        SpeciesSearchCommon::veryLocalSearch(speciesTree, evaluator,
                                             searchState, prune);
      }
    }
  }
  return better;
}

bool SpeciesSPRSearch::SPRSearch(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    SpeciesSearchState &searchState, unsigned int radius) {
  Logger::timed << "[Species search] Starting species tree local SPR search, "
                << "radius=" << radius << " (bestLL=" << searchState.bestLL
                << ", hash=" << speciesTree.getHash() << ")" << std::endl;
  assert(evaluator.computeLikelihood() == searchState.bestLL);
  bool better = false;
  bool tryAgain = false;
  do {
    // run a local SPR round
    tryAgain = SPRRound(speciesTree, evaluator, searchState, radius);
    better |= tryAgain;
  } while (tryAgain);
  // update bestLL due to possible precision changes
  PerFamLL perFamLL;
  double ll = evaluator.computeLikelihood(&perFamLL);
  searchState.betterLikelihoodCallback(ll, perFamLL);
  Logger::timed << "[Species search] After local SPR search: LL="
                << searchState.bestLL << std::endl;
  return better;
}
