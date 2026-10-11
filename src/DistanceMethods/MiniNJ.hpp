#pragma once

#include "NeighborJoining.hpp"

class PLLUnrootedTree;
class GeneSpeciesMapping;

/*
 * Naive NJ implementation.
 *
 * Both distance matrix and NJ tree computations could be
 * implemented more efficiently if needed.
 */
class MiniNJ {
public:
  MiniNJ() = delete;

  /**
   *  Run the original NJst algorithm
   */
  static std::unique_ptr<PLLRootedTree> runNJst(const Families &families);

  /**
   *  Infer a NJ tree, using the gene tree internode distances to compute
   *  the distance matrix. The distance matrix is very similar to the one
   *  built in NJst
   */
  static std::unique_ptr<PLLRootedTree> runMiniNJ(const Families &families);
  static std::unique_ptr<PLLRootedTree> runWMiniNJ(const Families &families);
  static std::unique_ptr<PLLRootedTree> runUstar(const Families &families);

  static std::unique_ptr<PLLRootedTree>
  applyNJ(DistanceMatrix &distanceMatrix,
          std::vector<std::string> &speciesIdToSpeciesString,
          StringToUint &speciesStringToSpeciesId);

  static void
  computeDistanceMatrix(const Families &families, bool minMode, bool reweight,
                        bool ustar, double contractBranchUnder,
                        DistanceMatrix &distanceMatrix,
                        std::vector<std::string> &speciesIdToSpeciesString,
                        StringToUint &speciesStringToSpeciesId);

  static void geneDistancesFromGeneTree(PLLUnrootedTree &geneTree,
                                        GeneSpeciesMapping &mapping,
                                        StringToUint &speciesStringToSpeciesId,
                                        DistanceMatrix &distances,
                                        DistanceMatrix &distancesDenominator,
                                        bool minMode, bool reweight, bool ustar,
                                        double contractBranchUnder = 0.0000011);

private:
  static std::unique_ptr<PLLRootedTree>
  geneTreeNJ(const Families &families, bool minAlgo, bool ustarAlgo = false,
             bool reweight = false, double contractBranchUnder = 0.0000011);
};
