#pragma once

#include <IO/Families.hpp>

class ModelParameters;

class ReconciliationBLEstimator {
public:
  ReconciliationBLEstimator() = delete;

  /**
   *  Estimate the branch lengths of a species tree by
   *  reconciling the gene trees (from the families) using
   *  the model of reconciliation described in modelParameters
   *  and then averaging over the lengths of gene tree branches
   *  that represent a branch in the species tree.
   *  The exact procedure is formally described in the SpeciesRax
   *  paper
   */
  static void estimate(const std::string &speciesTreeFile,
                       const Families &families,
                       const ModelParameters &modelParameters);
};
