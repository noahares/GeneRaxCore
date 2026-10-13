#pragma once

/**
 *  Entry point for slave jobs scheduled with the class Scheduler
 *  (through the external dependency MPIScheduler)
 */
extern "C" int static_scheduled_main(int argc, char **argv, void *commPtr);

class SlavesMain {
public:
  SlavesMain() = delete;

  /**
   *  Indicates whether the generax main is called by the user
   *  (the main run) or by the scheduler (a slave run).
   *  This is done through scanning the arguments in argv
   */
  static bool isSlave(int argc, char **argv);
};
