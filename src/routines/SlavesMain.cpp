#include "SlavesMain.hpp"

#include <cassert>
#include <string>

#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <routines/scheduled_routines/GeneRaxSlave.hpp>
#include <routines/scheduled_routines/RaxmlSlave.hpp>

bool SlavesMain::isSlave(int argc, char **argv) {
  if (argc < 2) {
    return false;
  }
  std::string key(argv[1]);
  return key == "optimizeGeneTrees" || key == "raxmlLight";
}

extern "C" int static_scheduled_main(int argc, char **argv, void *commPtr) {
  int res = 1;
  assert(commPtr);
  ParallelContext::init(commPtr);
  bool isNewProc = (*static_cast<int *>(commPtr) == -1);
  if (isNewProc) { // new process, init Logger
    Logger::init();
  } else { // original process, Logger is already inited
    Logger::enableLogFile(false);
  }
  std::string key(argv[1]);
  if (key == "optimizeGeneTrees") {
    res = GeneRaxSlave::optimizeGeneTreesMain(argc, argv);
  } else if (key == "raxmlLight") {
    res = RaxmlSlave::runRaxmlOptimization(argc, argv);
  } else {
    assert(false);
  }
  if (isNewProc) {
    Logger::close();
  } else {
    Logger::enableLogFile(true);
  }
  ParallelContext::finalize();
  return res;
}
