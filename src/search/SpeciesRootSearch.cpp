#include "SpeciesRootSearch.hpp"

#include <cassert>

#include "DatedSpeciesTreeSearch.hpp"
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <trees/SpeciesTree.hpp>

static void rootSearchAux(SpeciesTree &speciesTree,
                          SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                          SpeciesSearchState &searchState,
                          std::vector<unsigned int> &movesHistory,
                          std::vector<unsigned int> &bestMovesHistory,
                          DatedBackup &bestDatedBackup, double bestStackLL,
                          unsigned int maxDepth,
                          RootLikelihoods *rootLikelihoods,
                          TreePerFamLLVec *treePerFamLLVec) {
  // update rootLikelihoods and treePerFamLLVec with the current rooting
  PerFamLL perFamLL;
  double ll = evaluator.computeLikelihood(&perFamLL);
  if (rootLikelihoods) {
    auto root = speciesTree.getRoot();
    rootLikelihoods->saveRootLikelihood(root, ll);
    rootLikelihoods->savePerFamilyLikelihoods(root, perFamLL);
  }
  if (treePerFamLLVec) {
    PerFamLL globalPerFamLL;
    ParallelContext::concatenateHeterogeneousDoubleVectors(perFamLL,
                                                           globalPerFamLL);
    auto newick = speciesTree.toString();
    treePerFamLLVec->push_back({newick, globalPerFamLL});
  }
  // generate new moves to apply or end the search
  if (movesHistory.size() >= maxDepth) {
    return;
  }
  std::vector<unsigned int> moves;
  for (unsigned int direction = 0; direction < 4; ++direction) {
    if (movesHistory.empty() || (movesHistory.back() % 2 == direction % 2)) {
      moves.push_back(direction);
    }
  }
  // apply the new moves independently to the current rooting
  for (auto direction : moves) {
    if (!SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      continue;
    }
    auto bestLL = searchState.bestLL; // best LL before applying the move
    movesHistory.push_back(direction);
    evaluator.pushRollback();
    auto backup = speciesTree.getDatedTree().getBackup();
    SpeciesTreeOperator::changeRoot(speciesTree, direction);
    double ll = DatedSpeciesTreeSearch::optimizeDates(
        speciesTree, evaluator, searchState, !searchState.farFromPlausible);
    if (ll > bestLL) {
      // found the best tree over all performed moves
      // already saved it to searchState in the optimizeDates function
      bestMovesHistory = movesHistory;
      bestDatedBackup = speciesTree.getDatedTree().getBackup();
      Logger::timed << "\tbetter root: LL=" << ll << std::endl;
    }
    auto newBestStackLL = bestStackLL;
    auto newMaxDepth = maxDepth;
    if (ll > bestStackLL) {
      // found the best tree in this history so far
      // switch to local search: find a better tree within a 2-move radius
      newBestStackLL = ll;
      newMaxDepth = movesHistory.size() + 2;
    }
    rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                  bestMovesHistory, bestDatedBackup, newBestStackLL,
                  newMaxDepth, rootLikelihoods, treePerFamLLVec);
    // rollback the changes before the next move or leaving the recursion
    SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
    SpeciesTreeOperator::restoreDates(speciesTree, backup);
    evaluator.popAndApplyRollback();
    movesHistory.pop_back();
  }
}

bool SpeciesRootSearch::rootSearch(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    SpeciesSearchState &searchState, unsigned int maxDepth,
    RootLikelihoods *rootLikelihoods, TreePerFamLLVec *treePerFamLLVec) {
  Logger::timed << "[Species search] Root search with depth=" << maxDepth
                << std::endl;
  assert(evaluator.computeLikelihood() == searchState.bestLL);
  // optimize initial tree node dates
  double initialLL = DatedSpeciesTreeSearch::optimizeDates(
      speciesTree, evaluator, searchState, !searchState.farFromPlausible);
  // run recursive root search
  if (treePerFamLLVec) {
    treePerFamLLVec->clear();
  }
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  auto bestDatedBackup = speciesTree.getDatedTree().getBackup();
  rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                bestMovesHistory, bestDatedBackup, initialLL, maxDepth,
                rootLikelihoods, treePerFamLLVec);
  // rearrange the initial tree (apply the best moves and the best dating)
  for (auto direction : bestMovesHistory) {
    SpeciesTreeOperator::changeRoot(speciesTree, direction);
  }
  SpeciesTreeOperator::restoreDates(speciesTree, bestDatedBackup);
  bool better = !bestMovesHistory.empty();
  // update bestLL due to possible precision changes
  PerFamLL perFamLL;
  double ll = evaluator.computeLikelihood(&perFamLL);
  searchState.betterLikelihoodCallback(ll, perFamLL);
  Logger::timed << "[Species search] After root search: LL="
                << searchState.bestLL << std::endl;
  return better;
}
