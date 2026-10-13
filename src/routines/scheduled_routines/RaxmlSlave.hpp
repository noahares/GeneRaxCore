#pragma once

class RaxmlSlave {
public:
  RaxmlSlave() = delete;
  /**
   *  Runs a job scheduled by RaxmlMaster with MPIScheduler
   *  See RaxmlMaster code to see the list of input parmeters
   *  in argv.
   *  @param argc Number of parameters
   *  @param argv Parameters
   */
  static int runRaxmlOptimization(int argc, char **argv);
};
