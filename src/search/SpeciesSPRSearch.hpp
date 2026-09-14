#pragma once

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class SpeciesSearchState;

class SpeciesSPRSearch {
public:
  SpeciesSPRSearch() = delete;

  /**
   *  Search for the ML topology of the current species tree
   *  with local SPR moves
   */
  static bool SPRSearch(SpeciesTree &speciesTree,
                        SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                        SpeciesSearchState &searchState, unsigned int radius);
};
