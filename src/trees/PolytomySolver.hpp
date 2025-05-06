#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <vector>

class PLLRootedTree;
class PLLUnrootedTree;

class PolytomySolver {
public:
  static void solveSimpleInterface(
      PLLRootedTree &speciesTree,
      std::map<std::string, unsigned int> &speciesLabelsToSolve);
};
