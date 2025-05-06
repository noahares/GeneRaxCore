#pragma once

#include <corax/corax.h>
#include <maths/Parameters.hpp>
#include <string>
class JointTree;

struct GeneRaxCheckpoint {
  std::string path;
  bool checkpointExists;
  bool modelParamDone;
  Parameters ratesVector;
  std::string substModelStr;
  std::string geneTreeNewickStr;

  GeneRaxCheckpoint(const std::string checkpointPath);

  void save(bool masterRankOnly);
};
