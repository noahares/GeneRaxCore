#pragma once

#include <IO/FamiliesFileParser.hpp>
#include <string>
#include <util/enums.hpp>
#include <vector>

class Parameters;
struct RecModelInfo;

class GeneRaxMaster {
public:
  GeneRaxMaster() = delete;

  static void optimizeGeneTrees(
      Families &families, const RecModelInfo &recModelInfo, Parameters &rates,
      const std::string &output, const std::string &resultName,
      const std::string &execPath, const std::string &speciesTreePath,
      RecOpt reconciliationOpt, bool madRooting, double supportThreshold,
      double recWeight, bool enableRec, bool enableLibpll,
      unsigned int sprRadius, unsigned int iteration, bool schedulerSplitImplem,
      long &elapsed, bool inPlace = false);
};
