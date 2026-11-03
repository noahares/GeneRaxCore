#pragma once

#include <string>
#include <unordered_set>
#include <vector>

struct FamilyInfo {
  std::string name;
  std::string startingGeneTree;
  std::string ccpFile;
  std::string alignmentFile;
  std::string mappingFile;
  std::string libpllModel;
  std::string statsFile;
  std::string likelihoodFile;
  unsigned int color;
  FamilyInfo() { reset(); }
  void reset() {
    name = "";
    startingGeneTree = "__random__";
    ccpFile = "";
    alignmentFile = "";
    mappingFile = "";
    libpllModel = "GTR";
    statsFile = "";
    color = 0;
    likelihoodFile = "";
  }
};
using Families = std::vector<FamilyInfo>;

class Family {
public:
  Family() = delete;
  static void filterFamilies(Families &families,
                             const std::string &speciesTreeFile,
                             bool checkAlignments, bool checkSpeciesTree);
  static void printStats(const Families &families,
                         const std::string &speciesTreeFile,
                         const std::string &coverageFile,
                         const std::string &fractionMissingFile);
  static std::unordered_set<std::string>
  getLabelsFromFamilies(const Families &families);
};
