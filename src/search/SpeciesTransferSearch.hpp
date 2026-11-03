#pragma once

#include <functional>

#include <util/types.hpp>

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class SpeciesSearchState;
class Scenario;

struct TransferMove {
  unsigned int prune;
  unsigned int regraft;
  unsigned int transfers;
  double support;
  TransferMove() : prune(0), regraft(0), transfers(0), support(0.0) {}
  TransferMove(unsigned int p, unsigned int r, unsigned int t)
      : prune(p), regraft(r), transfers(t), support(0.0) {}
  // used for sorting
  bool operator<(const TransferMove &other) const {
    if (support != other.support) {
      return support < other.support;
    } else if (regraft != other.regraft) {
      return regraft < other.regraft;
    } else {
      return prune < other.prune;
    }
  }
  // used for blacklisting
  bool operator==(const TransferMove &other) const {
    return (prune == other.prune) && (regraft == other.regraft) &&
           (transfers == other.transfers);
  }
};

namespace std {
template <> struct hash<TransferMove> {
  size_t operator()(const TransferMove &move) const {
    auto hashints = [](unsigned int a, unsigned int b) {
      return (a + b) * (a + b + 1) / 2 + b;
    };
    return hash<int>()(static_cast<int>(
        hashints(hashints(move.prune, move.regraft), move.transfers)));
  }
};
} // namespace std

struct MovesBlackList {
  std::unordered_set<TransferMove> _blacklist;
  bool isBlackListed(const TransferMove &move) const {
    return _blacklist.find(move) != _blacklist.end();
  }
  void blackList(const TransferMove &move) { _blacklist.insert(move); }
};

struct PerCorePotentialTransfers {
  MatrixUint copies; // copies[species][family]
  void addScenario(const Scenario &scenario);
  unsigned int countPotentialTransfers(unsigned int src,
                                       unsigned int dest) const;
};

class SpeciesTransferSearch {
public:
  SpeciesTransferSearch() = delete;

  /**
   *  Return HGT directions sorted in the descending order by
   *  the number of transfer events across all gene families
   *  @param minTransfers The min number of transfers to include a move
   *  @param blacklist The moves that should not be included in any case
   *
   *  Note: in TransferMove, regraft is the HGT src and prune is the HGT dest
   */
  static void
  getSortedTransferList(SpeciesTree &speciesTree,
                        SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                        unsigned int minTransfers,
                        const MovesBlackList &blacklist,
                        std::vector<TransferMove> &transferMoves);

  /**
   *  Search for the ML topology of the current species tree
   *  with SPR moves along the most frequent HGT directions
   */
  static bool transferSearch(SpeciesTree &speciesTree,
                             SpeciesTreeLikelihoodEvaluatorInterface &evaluator,
                             SpeciesSearchState &searchState);
};
