#pragma once

#include <util/types.hpp>

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class SpeciesSearchState;

struct ScoredBackup {
  DatedBackup backup;
  double score;
  ScoredBackup() : score(0.0) {}
  ScoredBackup(const DatedBackup &backup, double score)
      : backup(backup), score(score) {}
  bool operator<(const ScoredBackup &other) const {
    return score < other.score;
  }
};
using ScoredBackups = std::vector<ScoredBackup>;

class DatedSpeciesTreeSearch {
public:
  /**
   *  Optimize the speciation order (dating) of the current species
   *  tree. Save the tree and update searchState.bestLL on finding
   *  a dating with likelihood higher than searchState.bestLL
   *
   *  If the current species tree is not the current best tree, then
   *  the optimization may result in a tree having likelihood lower
   *  than searchState.bestLL (desired in SpeciesRootSearch class)
   *
   *  If thorough is not set, only apply one naive round. Otherwise,
   *  additionally conduct search with random dating perturbations
   */
  static double
  optimizeDates(SpeciesTree &speciesTree,
                SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                SpeciesSearchState &searchState, bool thorough);

  /**
   *  Generate and test random datings based on their transfer scores
   *  and return the best of them with computed LLs
   *  @param toTest The number of random datings to test
   *  @param toTake The number of best datings to return
   *
   *  The input species tree remains unchanged
   */
  static ScoredBackups getBestDatingsFromReconciliation(
      SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluatorInterface &evaluator, unsigned int toTest,
      unsigned int toTake);
};
