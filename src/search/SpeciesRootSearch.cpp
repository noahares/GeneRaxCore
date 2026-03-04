#include "SpeciesRootSearch.hpp"

#include "DatedSpeciesTreeSearch.hpp"
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <trees/SpeciesTree.hpp>

static void rootSearchAux(SpeciesTree &speciesTree,
                          SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                          SpeciesSearchState &searchState,
                          std::vector<unsigned int> &movesHistory,
                          std::vector<unsigned int> &bestMovesHistory,
                          DatedBackup &bestDatedBackup, double &bestLL,
                          double bestLLStack, unsigned int maxDepth,
                          RootLikelihoods *rootLikelihoods,
                          TreePerFamLLVec *treePerFamLLVec) {
  if (movesHistory.size() > maxDepth) {
    return;
  }
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  for (auto direction : moves) {
    if (!SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      continue;
    }
    movesHistory.push_back(direction);
    evaluator.pushRollback();
    auto backup = speciesTree.getDatedTree().getBackup();
    SpeciesTreeOperator::changeRoot(speciesTree, direction);
    double ll = DatedSpeciesTreeSearch::optimizeDates(
        speciesTree, evaluator, searchState, !searchState.farFromPlausible);
    PerFamLL perFamLL;
    ll = evaluator.computeLikelihood(&perFamLL);
    if (treePerFamLLVec) {
      PerFamLL globalPerFamLL;
      ParallelContext::concatenateHeterogeneousDoubleVectors(perFamLL,
                                                             globalPerFamLL);
      auto newick = speciesTree.toString();
      treePerFamLLVec->push_back({newick, globalPerFamLL});
    }
    if (rootLikelihoods) {
      auto root = speciesTree.getRoot();
      rootLikelihoods->saveRootLikelihood(root, ll);
      rootLikelihoods->savePerFamilyLikelihoods(root, perFamLL);
    }
    auto newMaxDepth = maxDepth;
    if (ll > bestLLStack) {
      bestLLStack = ll;
      newMaxDepth = movesHistory.size() + 2;
    }
    if (ll > bestLL) {
      bestLL = ll;
      bestMovesHistory = movesHistory;
      bestDatedBackup = speciesTree.getDatedTree().getBackup();
      Logger::timed << "\tbetter root: LL=" << ll << std::endl;
    }
    rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                  bestMovesHistory, bestDatedBackup, bestLL, bestLLStack,
                  newMaxDepth, rootLikelihoods, treePerFamLLVec);
    SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
    SpeciesTreeOperator::restoreDates(speciesTree, backup);
    evaluator.popAndApplyRollback();
    movesHistory.pop_back();
  }
}

double SpeciesRootSearch::rootSearch(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    SpeciesSearchState &searchState, unsigned int maxDepth,
    RootLikelihoods *rootLikelihoods, TreePerFamLLVec *treePerFamLLVec) {
  Logger::timed << "[Species search] Root search with depth=" << maxDepth
                << std::endl;
  PerFamLL perFamLL;
  double initialLL = evaluator.computeLikelihood(&perFamLL);
  if (treePerFamLLVec) {
    treePerFamLLVec->clear();
    PerFamLL globalPerFamLL;
    ParallelContext::concatenateHeterogeneousDoubleVectors(perFamLL,
                                                           globalPerFamLL);
    auto newick = speciesTree.toString();
    treePerFamLLVec->push_back({newick, globalPerFamLL});
  }
  if (rootLikelihoods) {
    auto root = speciesTree.getRoot();
    rootLikelihoods->saveRootLikelihood(root, initialLL);
    rootLikelihoods->savePerFamilyLikelihoods(root, perFamLL);
  }
  double bestLL = initialLL;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  auto bestDatedBackup = speciesTree.getDatedTree().getBackup();
  movesHistory.push_back(1);
  rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                bestMovesHistory, bestDatedBackup, bestLL, initialLL, maxDepth,
                rootLikelihoods, treePerFamLLVec);
  movesHistory[0] = 0;
  rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                bestMovesHistory, bestDatedBackup, bestLL, initialLL, maxDepth,
                rootLikelihoods, treePerFamLLVec);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(speciesTree, bestMovesHistory[i]);
  }
  SpeciesTreeOperator::restoreDates(speciesTree, bestDatedBackup);
  Logger::timed << "[Species search] After root search: LL=" << bestLL
                << std::endl;
  return bestLL;
}
