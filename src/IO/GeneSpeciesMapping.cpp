#include "GeneSpeciesMapping.hpp"
#include <IO/FileSystem.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <corax/corax.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <trees/PLLUnrootedTree.hpp>

void GeneSpeciesMapping::fill(const std::string &mappingFile,
                              const std::string &geneTreeStrOrFile) {
  if (mappingFile.size()) {
    buildFromMappingFile(mappingFile);
  } else {
    buildFromTrees(geneTreeStrOrFile);
  }
}

void GeneSpeciesMapping::fill(const GeneSpeciesMapping &mapping) {
  auto &m = mapping.getMap();
  for (auto &pair : m) {
    _map[pair.first] = pair.first;
  }
}

bool GeneSpeciesMapping::check(corax_utree_t *geneTree,
                               corax_rtree_t *speciesTree) {

  assert(false);
  std::unordered_set<std::string> geneLeaves;
  std::unordered_set<std::string> speciesLeaves;
  LibpllParsers::fillLeavesFromUtree(geneTree, geneLeaves);
  LibpllParsers::fillLeavesFromRtree(speciesTree, speciesLeaves);
  return check(geneLeaves, speciesLeaves);
}

bool GeneSpeciesMapping::check(
    const std::unordered_set<std::string> &geneLeaves,
    const std::unordered_set<std::string> &speciesLeaves) {
  bool ok = true;
  for (auto pair : getMap()) {
    auto &gene = pair.first;
    auto &species = pair.second;
    if (geneLeaves.find(gene) == geneLeaves.end()) {
      std::cerr << "[Error] Invalid mapping '" << gene << "' <-> '" << species
                << "': can't find the gene '" << gene << "' in the gene tree"
                << std::endl;
      ok = false;
    }
    if (speciesLeaves.find(species) == speciesLeaves.end()) {
      std::cerr << "[Error] Invalid mapping '" << gene << "' <-> '" << species
                << "': can't find the species '" << species
                << "' in the species tree" << std::endl;
      ok = false;
    }
  }
  for (auto &gene : geneLeaves) {
    auto speciesIt = getMap().find(gene);
    if (speciesIt == getMap().end()) {
      std::cerr << "[Error] Gene tree leaf '" << gene
                << "' is not mapped to any species" << std::endl;
      ok = false;
    }
  }
  return ok;
}

void GeneSpeciesMapping::buildFromMappingFile(const std::string &mappingFile) {
  std::ifstream f(mappingFile);
  std::string line;
  getline(f, line);
  f.seekg(0, std::ios::beg);
  if (line.find(":") != std::string::npos) {
    buildFromPhyldogMapping(f);
  } else {
    buildFromTreerecsMapping(f);
  }
}

static bool isBlanck(const std::string &s) {
  return s.empty() || std::all_of(s.begin(), s.end(),
                                  [](char c) { return std::isspace(c); });
}

void GeneSpeciesMapping::buildFromPhyldogMapping(std::ifstream &f) {
  /*
    species1:gene1;gene2;gene3
    species2:gene4;gene5
  */
  std::string line;
  while (getline(f, line)) {
    if (isBlanck(line)) {
      continue;
    }
    std::stringstream ss(line);
    std::string species;
    std::string gene;
    getline(ss, species, ':');
    while (getline(ss, gene, ';')) {
      _map[gene] = species;
    }
  }
}

void GeneSpeciesMapping::buildFromTreerecsMapping(std::ifstream &f) {
  /*
    gene1 species1
    gene2 species2
    gene3 species1
  */
  std::string line;
  while (getline(f, line)) {
    if (isBlanck(line)) {
      continue;
    }
    std::stringstream ss(line);
    std::string species;
    std::string gene;
    ss >> gene;
    ss >> species;
    _map[gene] = species;
  }
}

void GeneSpeciesMapping::buildFromTrees(const std::string &geneTreeStrOrFile) {
  std::string str = geneTreeStrOrFile;
  std::ifstream ifs(geneTreeStrOrFile);
  if (ifs.good()) {
    if (ifs.peek() == '#') {
      std::string temp;
      std::getline(ifs, temp);
      std::getline(ifs, str);
    } else {
      str.assign((std::istreambuf_iterator<char>(ifs)),
                 (std::istreambuf_iterator<char>()));
    }
    str = str.substr(0, str.find(";") + 1);
  }
  // now str is a newick string
  PLLUnrootedTree tree(str, false);
  auto labels = tree.getLabels();
  fillFromGeneLabels(labels);
}

void GeneSpeciesMapping::fillFromGeneLabels(
    const std::unordered_set<std::string> &labels) {
  for (const auto &label : labels) {
    std::string species;
    std::string gene;
    auto pos = label.find_first_of('_');
    species = label.substr(0, pos);
    gene = label; // label.substr(pos + 1);
    _map[gene] = species;
  }
}

std::unordered_set<std::string> GeneSpeciesMapping::getCoveredSpecies() const {
  std::unordered_set<std::string> res;
  for (auto &geneSpecies : _map) {
    res.insert(geneSpecies.second);
  }
  return res;
}
