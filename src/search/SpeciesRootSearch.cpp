#include "SpeciesRootSearch.hpp"

#include "DatedSpeciesTreeSearch.hpp"
#include <numeric>
#include <trees/PLLRootedTree.hpp>
#include <trees/SpeciesTree.hpp>

static void rootSearchAux(SpeciesTree &speciesTree,
                          SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                          SpeciesSearchState &searchState,
                          std::vector<unsigned int> &movesHistory,
                          std::vector<unsigned int> &bestMovesHistory,
                          DatedTree::Backup &bestDatedBackup, double &bestLL,
                          double bestLLStack, unsigned int &visits,
                          unsigned int maxDepth,
                          RootLikelihoods *rootLikelihoods,
                          TreePerFamLLVec *treePerFamLLVec) {
  if (movesHistory.size() > maxDepth) {
    return;
  }
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  auto backup = speciesTree.getDatedTree().getBackup();
  for (auto direction : moves) {
    if (!SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      continue;
    }
    movesHistory.push_back(direction);
    evaluator.pushRollback();
    SpeciesTreeOperator::changeRoot(speciesTree, direction);
    double ll = DatedSpeciesTreeSearch::optimizeDates(
        speciesTree, evaluator, searchState, searchState.farFromPlausible);
    PerFamLL perFamLL;
    ll = evaluator.computeLikelihood(&perFamLL);
    if (treePerFamLLVec) {
      PerFamLL globalPerFamLL;
      ParallelContext::concatenateHetherogeneousDoubleVectors(perFamLL,
                                                              globalPerFamLL);
      auto newick = speciesTree.getTree().getNewickString();
      treePerFamLLVec->push_back({newick, globalPerFamLL});
    }
    auto root = speciesTree.getRoot();
    if (rootLikelihoods) {
      rootLikelihoods->saveRootLikelihood(root, ll);
      rootLikelihoods->savePerFamilyLikelihoods(root, perFamLL);
    }
    visits++;
    unsigned int additionalDepth = 0;
    if (ll > bestLLStack) {
      additionalDepth = 2;
      bestLLStack = ll;
    }
    if (ll > bestLL) {
      bestLL = ll;
      bestMovesHistory = movesHistory;
      bestDatedBackup = speciesTree.getDatedTree().getBackup();
      Logger::timed << "\tbetter root: LL=" << ll << std::endl;
      searchState.betterTreeCallback(ll, perFamLL);
    }
    auto newMaxDepth = maxDepth;
    if (additionalDepth) {
      newMaxDepth = movesHistory.size() + additionalDepth;
    }
    rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                  bestMovesHistory, bestDatedBackup, bestLL, bestLLStack,
                  visits, newMaxDepth, rootLikelihoods, treePerFamLLVec);
    SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
    evaluator.popAndApplyRollback();
    movesHistory.pop_back();
    speciesTree.getDatedTree().restore(backup);
  }
}

double SpeciesRootSearch::rootSearch(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
    SpeciesSearchState &searchState, unsigned int maxDepth,
    RootLikelihoods *rootLikelihoods, TreePerFamLLVec *treePerFamLLVec) {
  Logger::timed << "[Species search] Root search with depth=" << maxDepth
                << std::endl;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  PerFamLL perFamLL;
  double bestLL = evaluator.computeLikelihood(&perFamLL);
  if (treePerFamLLVec) {
    PerFamLL globalPerFamLL;
    ParallelContext::concatenateHetherogeneousDoubleVectors(perFamLL,
                                                            globalPerFamLL);
    treePerFamLLVec->clear();
    auto newick = speciesTree.getTree().getNewickString();
    treePerFamLLVec->push_back({newick, globalPerFamLL});
  }
  auto root = speciesTree.getRoot();
  if (rootLikelihoods) {
    rootLikelihoods->saveRootLikelihood(root, bestLL);
    rootLikelihoods->savePerFamilyLikelihoods(root, perFamLL);
  }
  unsigned int visits = 1;
  movesHistory.push_back(1);
  double temp1 = bestLL;
  double temp2 = bestLL;
  auto bestDatedBackup = speciesTree.getDatedTree().getBackup();
  rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                bestMovesHistory, bestDatedBackup, bestLL, temp1, visits,
                maxDepth, rootLikelihoods, treePerFamLLVec);
  movesHistory[0] = 0;
  rootSearchAux(speciesTree, evaluator, searchState, movesHistory,
                bestMovesHistory, bestDatedBackup, bestLL, temp2, visits,
                maxDepth, rootLikelihoods, treePerFamLLVec);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(speciesTree, bestMovesHistory[i]);
  }
  speciesTree.getDatedTree().restore(bestDatedBackup);
  if (rootLikelihoods) {
    auto newick = speciesTree.getTree().getNewickString();
    PLLRootedTree tree(newick, false);
    rootLikelihoods->fillTree(tree);
  }
  Logger::timed << "[Species search] After root search: LL=" << bestLL << " "
                << std::endl;
  return bestLL;
}
