#pragma once

#include <iostream>

#include <util/types.hpp>

class GeneSpeciesMapping;
class PLLRootedTree;
class PLLUnrootedTree;

using CladeSet = std::set<size_t>;

class Clade {
public:
  Clade();

  void mergeWith(const Clade &clade);
  void addId(unsigned int id);

  size_t getHash() const { return _hash; }
  friend std::ostream &operator<<(std::ostream &os, const Clade &c) {
    os << "[(";
    for (auto id : c._ids) {
      os << id << ",";
    }
    os << ")" << c._hash << "]";
    return os;
  }

  Clade getComplement(Clade &allTaxa);

  static Clade getMaximumClade(PLLRootedTree &tree);

  static CladeSet buildCladeSet(PLLRootedTree &tree);

  static CladeSet buildCladeSet(PLLUnrootedTree &tree,
                                const GeneSpeciesMapping &geneSpeciesMapping,
                                const StringToUint &speciesLabelToInt);

private:
  void _recomputeHash();

  std::set<unsigned int> _ids;
  size_t _hash;
};
