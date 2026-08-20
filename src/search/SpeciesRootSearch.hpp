#pragma once

#include "SpeciesSearchCommon.hpp"

class SpeciesRootSearch {
public:
  SpeciesRootSearch() = delete;

  /**
   *  Search for the ML root of the current species tree
   */
  static bool rootSearch(SpeciesTree &speciesTree,
                         SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                         SpeciesSearchState &searchState, unsigned int maxDepth,
                         RootLikelihoods *rootLikelihoods = nullptr,
                         TreePerFamLLVec *treePerFamLLVec = nullptr);
};
