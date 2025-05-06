#pragma once

#include <trees/DatedTree.hpp>

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class SpeciesSearchState;

struct ScoredBackup {
  DatedTree::Backup backup;
  double score;
  ScoredBackup() : score(0.0) {}
  ScoredBackup(DatedTree &datedTree, double score)
      : backup(datedTree.getBackup()), score(score) {}
  bool operator<(const ScoredBackup &other) const {
    return score < other.score;
  }
};

using ScoredBackups = std::vector<ScoredBackup>;

class DatedSpeciesTreeSearch {
public:
  /**
   *  Search of the speciation order (dating) that optimizes
   *  the score returned by the evaluator. If the score gets higher
   *  than searchState.bestLL and if searchState.pathToBestSpeciesTree
   *  is set, save the new best tree and update bestLL.
   *
   *  If thorough is not set, we only apply one naive round.
   *  Otherwise, we conduct a more thorough search
   */
  static double
  optimizeDates(SpeciesTree &speciesTree,
                SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                SpeciesSearchState &searchState, bool thorough);

  /**
   *
   */
  static ScoredBackups optimizeDatesFromReconciliation(
      SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluatorInterface &evaluator, unsigned int searches,
      unsigned int toEvaluate = 10);
};
