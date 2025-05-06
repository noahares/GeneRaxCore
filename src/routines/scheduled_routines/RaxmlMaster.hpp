#pragma once

#include <IO/FamiliesFileParser.hpp>
#include <string>
#include <vector>

class RaxmlMaster {
public:
  RaxmlMaster() = delete;
  static void runRaxmlOptimization(Families &families,
                                   const std::string &output,
                                   const std::string &execPath,
                                   unsigned int iteration, bool splitImplem,
                                   long &sumElapsedSec);
};
